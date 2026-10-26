#!/usr/bin/python
"""Tests of the event catalog: origins, padded/support-aware membership, exact accounting, and the optional detector.

Synthetic data only (no IMAS, no experimental fixtures).  What is checked:

* membership uses the closed padded event window and the declared sample averaging interval, or labels the
  point-sample approximation; TIME and TIME3 stay distinct clocks;
* only events that may act as exclusions do (an unreviewed automatic candidate or a rejected event never does), and
  every use reports the event origin; commanded and modeled events are never independent experimental evidence;
* assessed + excluded + unassessed counts are exact, and an exclusion never suppresses an integrity failure;
* the detector runs on explicitly configured, unit-checked diagnostics, records candidates that exclude nothing,
  and no candidates is reported as not evidence of an event-free plasma.

    python examples/verify_curation_events.py
"""

import copy
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from curation import MANIFEST_SCHEMA
from curation import events as event_rules
from curation.manifest import validate_manifest
import verify_curation as base
from verify_curation import Variable, run_screen, variables

V3 = MANIFEST_SCHEMA


def event(identifier='e1', rules=('PROFILE.TE.*',), start=0.9, end=1.1, type='pellet', **extra):
    return {'id': identifier, 'type': type, 'start': start, 'end': end, 'rules': list(rules), **extra}


def canonical(event_dict):
    return validate_manifest({'schema_version': V3, 'events': [event_dict]})['events'][0]


class MembershipTests(unittest.TestCase):
    times = np.array([0.0, 1.0, 2.0, 3.0])

    def test_endpoints_are_inclusive_and_simultaneous_events_union(self):
        events = [canonical(event('a', start=1.0, end=1.0)), canonical(event('b', start=3.0, end=3.0))]
        mask, each, basis = event_rules.membership(events, 'PROFILE.TE.NEGATIVE', self.times)
        self.assertEqual(mask.tolist(), [False, True, False, True])
        self.assertEqual({k: v.tolist() for k, v in each.items()}, {'a': [False, True, False, False], 'b': [False, False, False, True]})
        self.assertEqual(basis, event_rules.POINT)

    def test_timing_uncertainty_and_padding_widen_the_window(self):
        for extra in ({'timing_uncertainty': 0.5}, {'padding': 0.5}, {'timing_uncertainty': 0.25, 'padding': 0.25}):
            mask, _each, _b = event_rules.membership([canonical(event(start=1.5, end=1.5, **extra))], 'PROFILE.TE.X', self.times)
            self.assertEqual(mask.tolist(), [False, True, True, False], extra)       # [1.0, 2.0] inclusive
        mask, _each, _b = event_rules.membership([canonical(event(start=1.5, end=1.5))], 'PROFILE.TE.X', self.times)
        self.assertFalse(mask.any())                                                  # point samples miss a point event

    def test_a_declared_averaging_interval_makes_a_sample_a_window_not_a_point(self):
        ev = [canonical(event(start=1.4, end=1.4))]
        centered = {'averaging_interval': 1.0, 'alignment': 'center'}                 # sample t covers [t-0.5, t+0.5]
        mask, _e, basis = event_rules.membership(ev, 'PROFILE.TE.X', self.times, centered)
        self.assertEqual((mask.tolist(), basis), ([False, True, False, False], event_rules.DECLARED))
        trailing = {'averaging_interval': 1.0, 'alignment': 'trailing'}               # sample t covers [t-1, t]
        mask, _e, _b = event_rules.membership(ev, 'PROFILE.TE.X', self.times, trailing)
        self.assertEqual(mask.tolist(), [False, False, True, False])

    def test_events_that_do_not_apply_to_the_rule_have_no_membership(self):
        mask, each, _b = event_rules.membership([canonical(event(rules=['NUMERICS.*']))], 'PROFILE.TE.X', self.times)
        self.assertFalse(mask.any())
        self.assertEqual(each, {})


class OriginAndEligibilityTests(unittest.TestCase):
    def test_only_events_that_may_act_as_exclusions_do(self):
        cases = (({'origin': 'analyst_declaration'}, True), ({'origin': 'actuator_command'}, True),
                 ({'origin': 'modeled_transp'}, True), ({'origin': 'experimental_observation'}, True),
                 ({'origin': 'automatic_candidate'}, False),
                 ({'origin': 'automatic_candidate', 'review_state': 'unreviewed'}, False),
                 ({'origin': 'automatic_candidate', 'review_state': 'reviewed_accepted'}, True),
                 ({'origin': 'analyst_declaration', 'review_state': 'reviewed_rejected'}, False))
        for extra, expected in cases:
            self.assertEqual(event_rules.acts_as_exclusion(canonical(event(**extra))), expected, extra)

    def negative_te(self):
        data = variables()
        te = np.array(data['TE'].values, dtype=float)
        te[1, 1] = -5.0
        data['TE'] = Variable(te, ('TIME3', 'X'), 'EV')
        return data

    def screen(self, **extra):
        return run_screen(self.negative_te(), {'schema_version': V3, 'events': [event(**extra)]})

    def test_an_unreviewed_candidate_excludes_nothing_but_is_inventoried(self):
        report, found = self.screen(origin='automatic_candidate', review_state='unreviewed', type='unclassified_transient')
        self.assertEqual(report.policy_exclusions, {})
        row = found['EVENTS.INVENTORY'].evidence['events'][0]
        self.assertEqual((row['acts_as_exclusion'], row['unused']), (False, True))
        self.assertEqual(found['EVENTS.INVENTORY'].status, 'warning')               # it matched nothing it could act on

    def test_an_accepted_candidate_excludes_and_the_origin_is_reported(self):
        report, _found = self.screen(origin='automatic_candidate', review_state='reviewed_accepted', type='unclassified_transient')
        record = report.policy_exclusions['PROFILE.TE.NEGATIVE']
        self.assertEqual(record['origins'], {'e1': 'automatic_candidate'})

    def test_a_commanded_or_modeled_event_is_never_independent_experimental_evidence(self):
        for origin, independent in (('actuator_command', False), ('modeled_transp', False), ('experimental_observation', True),
                                    ('analyst_declaration', False)):
            _report, found = self.screen(origin=origin)
            row = found['EVENTS.INVENTORY'].evidence['events'][0]
            self.assertEqual((row['origin'], row['independent_experimental_evidence']), (origin, independent), origin)

    def test_a_rejected_event_excludes_nothing(self):
        report, _found = self.screen(review_state='reviewed_rejected')
        self.assertEqual(report.policy_exclusions, {})

    def test_an_exclusion_never_suppresses_the_integrity_failure_it_covers(self):
        report, found = self.screen(rules=['*'], start=-5.0, end=5.0)
        self.assertEqual(found['PROFILE.TE.NEGATIVE'].status, 'fail')
        self.assertEqual(found['GRID.TIME3'].status, 'pass')
        for identifier in report.policy_exclusions:
            self.assertFalse(identifier.startswith('GRID.'), identifier)


class AccountingTests(unittest.TestCase):
    def test_assessed_excluded_and_unassessed_counts_are_exact(self):
        manifest = {'schema_version': V3, 'assessment_window': {'time_min': 0.0, 'time_max': 0.5},
                    'events': [event('inside', rules=['PROFILE.TE.FINITE'], start=0.0, end=0.0),
                               event('outside', rules=['PROFILE.TE.FINITE'], start=1.0, end=1.0)]}
        report, found = run_screen(manifest=manifest)
        record = report.policy_exclusions['PROFILE.TE.FINITE']
        scope = found['PROFILE.TE.FINITE'].scope
        self.assertEqual(scope['evaluated_count'], 3)                                # TIME3 = 0 only, three radii
        self.assertEqual(record['total_count'], 6)
        self.assertEqual((record['unassessed_count'], record['excluded_count'], record['eligible_count']), (3, 3, 0))
        self.assertEqual(record['excluded_count'] + record['eligible_count'], scope['evaluated_count'])
        self.assertEqual(record['total_count'], scope['evaluated_count'] + record['unassessed_count'])
        self.assertEqual(record['excluded_by_event'], {'inside': 3})                 # the other event hit no evaluated sample
        self.assertEqual(record['policy_status'], 'not_assessed')

    def test_time_and_time3_membership_are_distinct_even_with_equal_lengths(self):
        data = variables()
        data['TIME'] = Variable([0.0, 1.0], ('TIME',), 'SECONDS')                    # same length as TIME3, different clock
        data['TIME3'] = Variable([0.5, 1.5], ('TIME3',), 'SECONDS')
        data['NPTRITER'] = Variable([0, 3], ('TIME',))
        manifest = {'schema_version': V3, 'events': [event(rules=['PROFILE.TE.FINITE', 'NUMERICS.NPTRITER.FINITE'],
                                                            start=1.0, end=1.0)]}
        report, _found = run_screen(data, manifest)
        numerics = report.masks['NUMERICS_NPTRITER_FINITE__excluded']
        self.assertNotIn('PROFILE_TE_FINITE__excluded', report.masks)                 # TIME3 has no sample at 1.0
        self.assertEqual(numerics.tolist(), [False, True])                            # TIME does

    def test_the_support_basis_is_recorded_per_clock(self):
        manifest = {'schema_version': V3, 'events': [event(rules=['PROFILE.TE.FINITE'], start=0.9, end=1.1)]}
        report, _f = run_screen(manifest=manifest)
        self.assertEqual(report.policy_exclusions['PROFILE.TE.FINITE']['support_basis'], {'TIME3': event_rules.POINT})
        manifest['sample_support'] = {'TIME3': {'averaging_interval': 0.5}}
        report, _f = run_screen(manifest=manifest)
        self.assertEqual(report.policy_exclusions['PROFILE.TE.FINITE']['support_basis'], {'TIME3': event_rules.DECLARED})

    def test_padding_can_bring_an_otherwise_outside_event_into_the_source_time(self):
        far = event(start=1.2, end=1.3)                                              # the source ends at 1.0
        _r, found = run_screen(manifest={'schema_version': V3, 'events': [far]})
        self.assertTrue(found['EVENTS.INVENTORY'].evidence['events'][0]['outside_source_time'])
        _r, found = run_screen(manifest={'schema_version': V3, 'events': [{**far, 'timing_uncertainty': 0.3}]})
        self.assertFalse(found['EVENTS.INVENTORY'].evidence['events'][0]['outside_source_time'])


def detector_data(n=200, drop_at=100, drop=200.0, noise=2.0, seed=7, units='EV', include=True):
    rng = np.random.default_rng(seed)
    data = variables()
    time = np.arange(n) * 0.01
    data['TIME'] = Variable(time, ('TIME',), 'SECONDS')
    data['NPTRITER'] = Variable(np.zeros(n), ('TIME',))
    if include:
        level = 1000.0 + rng.normal(0.0, noise, n)
        level[drop_at:] -= drop
        data['TE0'] = Variable(level, ('TIME',), units)
    return data


DETECTOR = {'id': 'te0_drop', 'method': 'abrupt_drop', 'diagnostic': 'TE0', 'clock': 'TIME', 'mode': 'causal',
            'expected_units': ['EV'], 'max_sample_interval': 0.05}


class DetectorWiringTests(unittest.TestCase):
    def screen(self, data=None, detector=None, events=()):
        manifest = {'schema_version': V3, 'event_detectors': [{**DETECTOR, **(detector or {})}], 'events': list(events)}
        return run_screen(data if data is not None else detector_data(), manifest)

    def test_a_candidate_is_recorded_as_an_unreviewed_automatic_candidate(self):
        report, found = self.screen()
        finding = found['EVENTS.DETECT.te0_drop']
        self.assertEqual(finding.status, 'warning')
        self.assertEqual(len(report.detected_events), 1)
        candidate = report.detected_events[0]
        self.assertEqual((candidate['type'], candidate['origin'], candidate['review_state']),
                         ('unclassified_transient', 'automatic_candidate', 'unreviewed'))
        self.assertLessEqual(candidate['start'], 1.0)
        self.assertGreaterEqual(candidate['end'], 0.99)
        self.assertFalse(candidate['score_is_probability'])
        self.assertEqual((candidate['detector'], candidate['diagnostic'], candidate['clock']), ('te0_drop', 'TE0', 'TIME'))

    def test_a_candidate_excludes_nothing_on_its_own(self):
        report, _found = self.screen()
        self.assertEqual(report.policy_exclusions, {})

    def test_no_candidates_is_reported_as_not_evidence_of_an_event_free_plasma(self):
        report, found = self.screen(detector_data(drop=0.0))
        finding = found['EVENTS.DETECT.te0_drop']
        self.assertEqual(finding.status, 'pass')
        self.assertIn('not evidence of an event-free plasma', finding.message)
        self.assertEqual(report.detected_events, [])

    def test_a_missing_diagnostic_is_unassessed_never_a_pass(self):
        _r, found = self.screen(detector_data(include=False))
        finding = found['EVENTS.DETECT.te0_drop']
        self.assertEqual(finding.status, 'not_assessed')
        self.assertIn('missing_diagnostic', finding.message)

    def test_a_units_mismatch_is_an_error(self):
        _r, found = self.screen(detector_data(units='KEV'))
        self.assertEqual(found['EVENTS.DETECT.te0_drop'].status, 'error')
        self.assertIn('units', found['EVENTS.DETECT.te0_drop'].message.lower())

    def test_low_cadence_is_not_assessed(self):
        _r, found = self.screen(detector=None, data=detector_data())
        _r2, found2 = self.screen(detector={'max_sample_interval': 0.001}, data=detector_data())
        self.assertEqual(found['EVENTS.DETECT.te0_drop'].status, 'warning')
        self.assertEqual(found2['EVENTS.DETECT.te0_drop'].status, 'not_assessed')
        self.assertIn('insufficient_cadence', found2['EVENTS.DETECT.te0_drop'].message)

    def test_a_profile_needs_a_radial_index_and_no_axis_is_borrowed(self):
        data = detector_data()
        data['TE0'] = Variable(np.tile(data['TE0'].values[:, None], (1, 3)), ('TIME', 'X'), 'EV')
        _r, found = self.screen(data)
        self.assertEqual(found['EVENTS.DETECT.te0_drop'].status, 'error')            # a 2-D diagnostic without radial_index
        _r, found = self.screen(data, {'radial_index': 1})
        self.assertEqual(found['EVENTS.DETECT.te0_drop'].status, 'warning')
        _r, found = self.screen(detector_data(), {'clock': 'TIME3'})
        self.assertEqual(found['EVENTS.DETECT.te0_drop'].status, 'error')            # TE0 lives on TIME, not TIME3

    def test_a_declared_event_excludes_samples_from_the_detector_only_when_it_names_the_rule(self):
        window = event('outage', rules=['EVENTS.DETECT.te0_drop'], start=0.90, end=1.10, type='diagnostic_outage')
        _r, found = self.screen(events=[window])
        finding = found['EVENTS.DETECT.te0_drop']
        self.assertGreater(finding.evidence['excluded_by_declared_events'], 0)
        self.assertEqual(finding.evidence['detection']['coverage']['unassessed_by_reason']['excluded'] > 0, True)
        other = event('other', rules=['PROFILE.*'], start=0.90, end=1.10)
        _r, found = self.screen(events=[other])
        self.assertEqual(found['EVENTS.DETECT.te0_drop'].evidence['excluded_by_declared_events'], 0)

    def test_the_source_arrays_are_not_modified_by_detection(self):
        data = detector_data()
        before = {k: v.values.copy() for k, v in data.items()}
        self.screen(data)
        for name, values in before.items():
            self.assertTrue(np.array_equal(data[name].values, values), name)

    def test_detection_is_deterministic(self):
        a, _ = self.screen()
        b, _ = self.screen()
        self.assertEqual(a.detected_events, b.detected_events)

    def test_the_report_validates_and_records_detected_events(self):
        report, _found = self.screen()
        document = report.as_dict()
        from curation.schema import validate_structure
        validate_structure(document, 'report.v3')
        self.assertEqual(document['detected_events'][0]['origin'], 'automatic_candidate')


if __name__ == '__main__':
    unittest.main(verbosity=1)
