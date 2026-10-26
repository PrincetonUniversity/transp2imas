#!/usr/bin/python
"""Check curation.detectors on SYNTHETIC traces only (seeded random numbers; no IMAS, no real data).

These tests check the ALGORITHM's documented behaviour (including its false positives and negatives), not the
physical validity of any candidate.

    python examples/verify_curation_detectors.py
"""

import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from curation import detectors as D  # noqa: E402

DT = 1.0e-3
EXPECTED = {'EV', 'KEV'}


def trace(n=400, seed=1, noise=0.01, level=2.0):
    """Synthetic noisy flat trace in keV on a uniform 1 kHz clock."""
    rng = np.random.default_rng(seed)
    return np.arange(n) * DT, level + noise * rng.standard_normal(n)


def run(t, y, **kw):
    args = dict(units='keV', diagnostic='SYNTH_TE0', mode='causal', expected_units=EXPECTED, max_sample_interval=0.01)
    args.update(kw)
    return D.abrupt_drop(t, y, **args)


def drop_at(y, k, frac):
    """Multiply samples k.. by (1 - frac): a step drop between samples k-1 and k."""
    y = y.copy()
    y[k:] *= 1.0 - frac
    return y


class Detector(unittest.TestCase):
    def test_injected_drop_found_with_bounded_timing(self):
        t, y = trace()
        r = run(t, drop_at(y, 250, 0.3))
        self.assertEqual(r.status, 'ok')
        self.assertEqual(len(r.candidates), 1)
        c = r.candidates[0]
        self.assertTrue(c['start'] <= t[250] - DT / 2 <= c['end'])
        self.assertAlmostEqual(c['start'], t[249])
        self.assertAlmostEqual(c['end'], t[250])
        self.assertAlmostEqual(c['timing_uncertainty'], DT)
        self.assertEqual((c['type'], c['origin'], c['score_is_probability']), ('unclassified_transient',
                                                                              'automatic_candidate', False))
        self.assertIn('score', c)
        self.assertNotIn('probability', c)
        self.assertEqual((c['sample_before_index'], c['sample_after_index']), (249, 250))
        self.assertGreaterEqual(c['relative_drop'], 0.25)

    def test_padding_widens_interval_not_uncertainty(self):
        t, y = trace()
        c = run(t, drop_at(y, 250, 0.3), timing_padding=0.002).candidates[0]
        self.assertAlmostEqual(c['start'], t[249] - 0.002)
        self.assertAlmostEqual(c['end'], t[250] + 0.002)
        self.assertAlmostEqual(c['timing_uncertainty'], DT)

    def test_no_drop_ok_with_not_evidence_note(self):
        t, y = trace()
        r = run(t, y)
        self.assertEqual((r.status, r.candidates), ('ok', []))
        self.assertIn(D.NOT_EVIDENCE, r.as_dict()['notes'])

    def test_noise_only_false_positive_rate(self):
        flagged = steps = 0
        for seed in range(200):
            t, y = trace(seed=seed)
            r = run(t, y)
            flagged += len(r.candidates)
            steps += r.coverage['assessed']
        rate = flagged / steps
        # Observed on seeds 0..199 (78200 assessed steps): 0 candidates (rate 0). Documented bound: Gaussian noise
        # only, default threshold 5, expected well below 1e-3 candidates per assessed step.
        self.assertLess(rate, 1.0e-3, rate)

    def test_small_drop_is_documented_false_negative(self):
        t, y = trace(noise=0.001)
        r = run(t, drop_at(y, 250, 0.03))      # 3 % < min_relative_drop 5 %, though score is huge
        self.assertEqual(r.candidates, [])
        self.assertEqual(len(run(t, drop_at(y, 250, 0.03), min_relative_drop=0.02).candidates), 1)

    def test_drop_inside_nan_gap_not_bridged(self):
        t, y = trace()
        y = drop_at(y, 250, 0.3)
        y[248:252] = np.nan                    # the drop happens inside the outage
        r = run(t, y)
        self.assertEqual(r.candidates, [])
        self.assertTrue(any(lo <= t[248] and hi >= t[251] for lo, hi in r.coverage['unassessed_intervals']))
        self.assertEqual(r.coverage['unassessed_by_reason']['missing'], 5)

    def test_time_gap_splits_segments(self):
        t, y = trace()
        y = drop_at(y, 250, 0.3)
        t2, y2 = np.delete(t, range(245, 255)), np.delete(y, range(245, 255))
        r = run(t2, y2)
        self.assertEqual(r.candidates, [])
        self.assertEqual(r.coverage['unassessed_by_reason']['gap'], 1)

    def test_excluded_region_not_in_baseline_and_not_bridged(self):
        t, y = trace(noise=0.01)
        rng = np.random.default_rng(7)
        y[100:140] += 0.5 * rng.standard_normal(40)           # very noisy stretch, to be excluded by policy
        y = drop_at(y, 250, 0.3)
        excl = np.zeros(t.size, dtype=bool)
        excl[100:140] = True
        plain, masked = run(t, y), run(t, y, excluded=excl)
        self.assertEqual(masked.coverage['unassessed_by_reason']['excluded'], 41)
        self.assertEqual(len(masked.candidates), 1)
        # the noisy stretch inflates the baseline in the unexcluded run just after it
        near = lambda r: [c for c in r.candidates if c['start'] < 0.2]  # noqa: E731
        self.assertEqual(near(masked), [])
        # a drop placed inside the excluded region is not bridged
        y2 = y.copy()
        y2[120:] *= 0.5
        r = run(t, y2, excluded=excl)
        self.assertEqual([c for c in r.candidates if c['end'] < 0.15], [])
        # excluded values (e.g. NaN-free garbage) never feed the baseline: changing them changes nothing
        y3 = y.copy()
        y3[100:140] = 1.0e3
        a, b = run(t, y3, excluded=excl), masked
        self.assertEqual(a.as_dict()['candidates'], b.as_dict()['candidates'])
        self.assertEqual(plain.coverage['unassessed_by_reason']['excluded'], 0)

    def test_causal_equals_noncausal_when_stationary(self):
        t, y = trace(noise=0.01)
        y = drop_at(y, 250, 0.4)
        c, n = run(t, y), run(t, y, mode='noncausal')
        self.assertEqual([(x['start'], x['end']) for x in c.candidates], [(x['start'], x['end'])
                                                                         for x in n.candidates])
        self.assertEqual(len(c.candidates), 1)

    def test_variance_change_separates_modes(self):
        # noisy before the step, quiet after: the causal baseline is the noisy past (no flag); the noncausal
        # centered window also sees the quiet future and flags the step (checked for several seeds)
        t = np.arange(300) * DT
        for seed in (3, 4, 5, 6):
            rng = np.random.default_rng(seed)
            y = np.r_[2.0 + 0.1 * rng.standard_normal(150), 1.7 + 0.005 * rng.standard_normal(150)]
            causal, non = run(t, y), run(t, y, mode='noncausal')
            self.assertEqual(causal.candidates, [], seed)
            self.assertTrue(any(c['start'] <= t[149] and c['end'] >= t[150] for c in non.candidates), seed)
        self.assertIn('forecasting', ' '.join(non.notes))
        self.assertNotIn('forecasting', ' '.join(causal.notes))

    def test_causal_never_reads_future(self):
        t, y = trace()
        y = drop_at(y, 250, 0.3)
        k = 260
        full = run(t, y)
        cut = run(t[:k + 1], y[:k + 1])
        f = [c for c in full.candidates if c['sample_after_index'] <= k]
        self.assertEqual(f, cut.candidates)
        y2 = y.copy()
        y2[k + 1:] = np.random.default_rng(9).standard_normal(y.size - k - 1) * 50   # scramble the future
        c2 = [c for c in run(t, y2).candidates if c['sample_after_index'] <= k]
        self.assertEqual(c2, f)

    def test_low_cadence(self):
        t, y = trace()
        r = run(t, drop_at(y, 250, 0.3), max_sample_interval=DT / 2)
        self.assertEqual((r.status, r.candidates), ('insufficient_cadence', []))
        self.assertIn(D.NOT_EVIDENCE, r.notes)

    def test_status_missing_and_no_data_and_baseline(self):
        t, y = trace()
        self.assertEqual(run(None, None).status, 'missing_diagnostic')
        self.assertEqual(run(t, None).status, 'missing_diagnostic')
        self.assertEqual(run(t, np.full(t.size, np.nan)).status, 'no_data')
        self.assertEqual(run(t[:5], y[:5]).status, 'insufficient_baseline')
        self.assertEqual(run(t, np.full(t.size, 2.0)).status, 'insufficient_baseline')   # zero spread
        for s in (run(None, None), run(t[:5], y[:5])):
            self.assertIn(D.NOT_EVIDENCE, s.notes)

    def test_value_errors(self):
        t, y = trace()
        for bad in (dict(units=None), dict(units=''), dict(units='m'), dict(diagnostic=None), dict(diagnostic=' '),
                    dict(mode='both'), dict(mode=None), dict(time_units='ms'), dict(expected_units=set()),
                    dict(threshold=0.0), dict(window=3)):
            with self.assertRaises(ValueError, msg=str(bad)):
                run(t, y, **bad)
        with self.assertRaises(ValueError):
            run(t[::-1], y)
        with self.assertRaises(ValueError):
            run(np.r_[t[:10], t[9], t[11:]], y)
        with self.assertRaises(ValueError):
            run(t, y[:-1])
        with self.assertRaises(TypeError):
            D.abrupt_drop(t, y, 'keV', 'SYNTH', expected_units=EXPECTED, max_sample_interval=0.01)   # mode required
        self.assertEqual(run(t, y, units=' KeV ').status, 'ok')

    def test_adjacent_steps_merge(self):
        t, y = trace(noise=0.005)
        y = y.copy()
        y[250:] *= 0.8
        y[251:] *= 0.8                          # two consecutive large drops
        r = run(t, y)
        self.assertEqual(len(r.candidates), 1)
        c = r.candidates[0]
        self.assertAlmostEqual(c['start'], t[249])
        self.assertAlmostEqual(c['end'], t[251])
        self.assertAlmostEqual(c['timing_uncertainty'], 2 * DT)
        self.assertEqual((c['first_step_index'], c['last_step_index']), (250, 251))

    def test_masked_array_input(self):
        t, y = trace()
        y = drop_at(y, 250, 0.3)
        m = np.ma.masked_array(y, mask=np.zeros(y.size, dtype=bool))
        m.mask[100:110] = True
        m.data[100:110] = 1.0e6                 # hidden values must be ignored
        ref = y.copy()
        ref[100:110] = np.nan
        self.assertEqual(run(t, m).as_dict(), run(t, ref).as_dict())
        self.assertEqual(run(t, m).coverage['unassessed_by_reason']['missing'], 11)

    def test_identity_and_json(self):
        t, y = trace()
        y = drop_at(y, 250, 0.3)
        y[50:53] = np.nan
        excl = np.zeros(t.size, dtype=bool)
        excl[200:210] = True
        for mode in D.MODES:
            r = run(t, y, excluded=excl, mode=mode)
            cov = r.coverage
            self.assertEqual(cov['assessed'] + cov['unassessed'], cov['total_steps'], mode)
            self.assertEqual(cov['total_steps'], t.size - 1)
            self.assertEqual(sum(cov['unassessed_by_reason'].values()), cov['unassessed'])
            text = json.dumps(r.as_dict(), allow_nan=False)
            self.assertEqual(json.loads(text)['method'], 'abrupt_drop')
            self.assertEqual(json.loads(text)['method_version'], D.METHOD_VERSION)

    def test_inputs_not_mutated_and_deterministic(self):
        t, y = trace()
        y = drop_at(y, 250, 0.3)
        excl = np.zeros(t.size, dtype=bool)
        t.setflags(write=False)
        y.setflags(write=False)
        excl.setflags(write=False)
        a = run(t, y, excluded=excl).as_dict()
        b = run(t, y, excluded=excl).as_dict()
        self.assertEqual(json.dumps(a, sort_keys=True), json.dumps(b, sort_keys=True))


if __name__ == '__main__':
    unittest.main(verbosity=1)
