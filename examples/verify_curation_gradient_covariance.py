#!/usr/bin/python
"""Tests of the supplied-covariance gradient rule inside the engine, on explicitly synthetic profiles and covariances.

The numerical machinery (operator, covariance validation, sampling) is tested in ``verify_curation_uncertainty.py``.
Here: the rule reads a hash-checked resource, requires the native TIME3 axis exactly, scales units, reports its scope
and what it does not include, and never invents an uncertainty.  Nothing here is diagnostic data or evidence that a
real reconstruction's covariance has been calibrated.

    python examples/verify_curation_gradient_covariance.py
"""

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from curation import MANIFEST_SCHEMA, uncertainty
from curation.engine import CurationEngine
from curation.manifest import load_manifest
from curation.reporting import write_report
from curation.source import NativeSource
import verify_curation as base
from verify_curation import run_screen, variables

X = np.array([0.1, 0.5, 0.9])
H = 0.4                                                  # the synthetic grid spacing; TE = 1000 - 500 X, so dTE/dX = -500


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def resource(path, covariance, units='EV', time=(0.0, 1.0), **extra):
    np.savez(path, time=np.array(time), covariance=covariance)
    return {'id': 'cov_te', 'quantity': 'TE', 'kind': 'covariance', 'path': Path(path).name, 'sha256': sha(path),
            'units': units, 'coordinate': 'X', 'ordering': 'ascending_coordinate',
            'estimator': 'numpy_gradient_edge_order_1', 'provenance': 'synthetic fixture, not diagnostic data',
            'included_sources': ['synthetic statistical error only'], **extra}


def diagonal(sigma):
    return np.stack([np.diag(np.full(3, sigma ** 2))] * 2)


class CovarianceRuleTests(unittest.TestCase):
    def screen(self, tmp, covariance, data=None, **kwargs):
        spec = resource(Path(tmp) / 'cov.npz', covariance, **kwargs)
        manifest = load_manifest(self.write_sidecar(tmp, spec))
        return run_screen(data, manifest)

    def write_sidecar(self, tmp, spec):
        path = Path(tmp) / 'e.json'
        path.write_text(json.dumps({'schema_version': MANIFEST_SCHEMA, 'uncertainty': {'covariances': [spec]}}))
        return path

    def test_the_propagated_sigma_equals_the_analytic_value_and_the_scope_is_labelled(self):
        with tempfile.TemporaryDirectory() as tmp:
            report, found = self.screen(tmp, diagonal(10.0))
        finding = found['FIT.TE.GRADIENT_COVARIANCE']
        self.assertEqual(finding.status, 'pass')                    # |g| = 500 is far outside 2 sigma = 71 eV per unit X
        analytic = np.array([np.sqrt(2.0) * 10.0 / H, np.sqrt(2.0) * 10.0 / (2 * H), np.sqrt(2.0) * 10.0 / H])
        # g0 = (y1-y0)/h, g1 = (y2-y0)/(2h), g2 = (y2-y1)/h for independent unit-variance errors
        analytic[1] = 10.0 / H / np.sqrt(2.0)
        sigma = report.coordinates['uncertainty_cov_te__std_gradient']
        np.testing.assert_allclose(sigma, np.tile(analytic, (2, 1)), rtol=1e-12)
        self.assertEqual(finding.evidence['scope'], uncertainty.SCOPE_GRADIENT)
        self.assertEqual(finding.metrics['cross_time'], 'unmodeled')
        self.assertIn('radial-coordinate uncertainty', finding.evidence['not_included'])
        self.assertEqual(finding.evidence['included_sources'], ['synthetic statistical error only'])
        summary = report.uncertainty['cov_te']
        self.assertEqual((summary['slices_evaluated'], summary['slices_total']), (2, 2))
        np.testing.assert_allclose(summary['std_gradient_range'], [analytic.min(), analytic.max()], rtol=1e-12)

    def test_a_gradient_consistent_with_zero_is_flagged(self):
        with tempfile.TemporaryDirectory() as tmp:
            report, found = self.screen(tmp, diagonal(100.0))
        finding = found['FIT.TE.GRADIENT_COVARIANCE']
        self.assertEqual(finding.status, 'warning')
        # sigma_g = (354, 177, 354) per unit X, so |g| = 500 <= 2 sigma_g for the outer points only: 2 per slice, 2 slices
        self.assertEqual(finding.scope['flagged_count'], 4)
        self.assertEqual(finding.scope['flagged_index_examples'], [[0, 0], [0, 2], [1, 0], [1, 2]])
        self.assertIn('conditional on the supplied reconstruction', finding.message)

    def test_a_shared_additive_offset_contributes_no_gradient_variance(self):
        offset = np.stack([np.full((3, 3), 50.0 ** 2)] * 2)                          # C = sigma^2 * 1 1^T
        with tempfile.TemporaryDirectory() as tmp:
            report, found = self.screen(tmp, offset)
        sigma = report.coordinates['uncertainty_cov_te__std_gradient']
        self.assertLess(np.max(np.abs(sigma)), 1e-12 * 50.0 / H)
        self.assertEqual(found['FIT.TE.GRADIENT_COVARIANCE'].status, 'pass')

    def test_covariance_units_are_converted_by_the_square_of_the_factor(self):
        with tempfile.TemporaryDirectory() as tmp:
            report_ev, _ = self.screen(tmp, diagonal(10.0))
        with tempfile.TemporaryDirectory() as tmp:
            report_kev, _ = self.screen(tmp, diagonal(0.01), units='KEV')            # 0.01 keV = 10 eV
        np.testing.assert_allclose(report_kev.coordinates['uncertainty_cov_te__std_gradient'],
                                   report_ev.coordinates['uncertainty_cov_te__std_gradient'], rtol=1e-12)

    def test_unknown_units_are_an_error_not_a_guess(self):
        with tempfile.TemporaryDirectory() as tmp:
            _r, found = self.screen(tmp, diagonal(10.0), units='JOULE')
        self.assertEqual(found['FIT.TE.GRADIENT_COVARIANCE'].status, 'error')
        self.assertIn('JOULE', found['FIT.TE.GRADIENT_COVARIANCE'].message)

    def test_the_time_axis_must_be_the_native_time3_exactly(self):
        with tempfile.TemporaryDirectory() as tmp:
            _r, found = self.screen(tmp, diagonal(10.0), time=(0.0, 0.9999))
        self.assertEqual(found['FIT.TE.GRADIENT_COVARIANCE'].status, 'error')
        self.assertIn('no interpolation', found['FIT.TE.GRADIENT_COVARIANCE'].message)

    def test_a_wrong_shape_a_stacked_time_covariance_and_a_non_psd_matrix_are_errors(self):
        bad_psd = diagonal(10.0)
        bad_psd[0, 0, 0] = -1.0
        asymmetric = diagonal(10.0)
        asymmetric[1, 0, 1] = 50.0
        for label, covariance in (('shape', np.zeros((2, 4, 4))), ('stacked', np.zeros((1, 6, 6))), ('not psd', bad_psd),
                                  ('asymmetric', asymmetric), ('nan', np.full((2, 3, 3), np.nan))):
            with self.subTest(label), tempfile.TemporaryDirectory() as tmp:
                _r, found = self.screen(tmp, covariance)
                self.assertEqual(found['FIT.TE.GRADIENT_COVARIANCE'].status, 'error', label)

    def test_no_covariance_means_unknown_not_zero(self):
        _r, found = run_screen()
        self.assertNotIn('FIT.TE.GRADIENT_COVARIANCE', found)
        self.assertEqual(found['UNCERTAINTY.electron_temperature'].status, 'not_assessed')

    def test_a_missing_profile_is_an_error_finding(self):
        data = variables()
        del data['TE']
        with tempfile.TemporaryDirectory() as tmp:
            _r, found = self.screen(tmp, diagonal(10.0), data)
        self.assertEqual(found['FIT.TE.GRADIENT_COVARIANCE'].status, 'error')

    def test_the_resource_and_the_source_are_not_modified(self):
        data = variables()
        before = {k: v.values.copy() for k, v in data.items()}
        with tempfile.TemporaryDirectory() as tmp:
            digest = None
            spec = resource(Path(tmp) / 'cov.npz', diagonal(10.0))
            digest = sha(Path(tmp) / 'cov.npz')
            run_screen(data, load_manifest(self.write_sidecar(tmp, spec)))
            self.assertEqual(sha(Path(tmp) / 'cov.npz'), digest)
        for name, values in before.items():
            self.assertTrue(np.array_equal(data[name].values, values), name)

    def test_the_existing_supplied_sigma_path_is_unchanged(self):
        data = variables()
        data['DTE_DX_SIGMA'] = base.Variable(np.full((2, 3), 600.0), ('TIME3', 'X'), 'EV')
        manifest = {'schema_version': MANIFEST_SCHEMA, 'gradient_uncertainties': {
            'TE': {'sigma_variable': 'DTE_DX_SIGMA', 'coordinate': 'X', 'method': 'numpy_gradient_edge_order_1', 'source': 's'}}}
        _r, found = run_screen(data, manifest)
        self.assertEqual(found['FIT.TE.GRADIENT_UNCERTAINTY'].status, 'warning')
        self.assertNotIn('FIT.TE.GRADIENT_COVARIANCE', found)

    def test_a_policy_can_disable_the_rule_and_the_ledger_says_so(self):
        from curation.policy import Policy
        with tempfile.TemporaryDirectory() as tmp:
            spec = resource(Path(tmp) / 'cov.npz', diagonal(10.0))
            manifest = load_manifest(self.write_sidecar(tmp, spec))
            source = NativeSource(base.types.SimpleNamespace(variables=variables()))
            report = CurationEngine(source, manifest, None, Policy(disabled={'uncertainty.gradient_covariance': 'not needed'})).run()
        states = {e['rule']: e['state'] for e in report.ledger}
        self.assertEqual(states['uncertainty.gradient_covariance'], 'policy_disabled')
        self.assertNotIn('FIT.TE.GRADIENT_COVARIANCE', {f.id for f in report.findings})


class CovarianceCliTests(unittest.TestCase):
    def test_the_cli_runs_a_covariance_sidecar_resolved_beside_the_file_and_keeps_inputs_untouched(self):
        with tempfile.TemporaryDirectory() as tmp, tempfile.TemporaryDirectory() as elsewhere:
            tmp = Path(tmp)
            base.ReportingAndCliTests().write_cdf(tmp / 'run.CDF')
            spec = resource(tmp / 'cov.npz', diagonal(10.0))
            (tmp / 'e.json').write_text(json.dumps({'schema_version': MANIFEST_SCHEMA, 'uncertainty': {'covariances': [spec]}}))
            before = {n: sha(tmp / n) for n in ('run.CDF', 'cov.npz', 'e.json')}
            done = subprocess.run([sys.executable, str(ROOT / 'curate.py'), '-i', str(tmp / 'run.CDF'),
                                   '--evidence', str(tmp / 'e.json'), '--output-dir', str(tmp / 'out')],
                                  capture_output=True, text=True, cwd=elsewhere)
            self.assertEqual(done.returncode, 0, done.stderr + done.stdout)
            document = json.loads((tmp / 'out' / 'curation.json').read_text())
            found = {f['id']: f for f in document['findings']}
            self.assertEqual(found['FIT.TE.GRADIENT_COVARIANCE']['status'], 'pass')
            self.assertEqual(document['uncertainty']['cov_te']['scope'], uncertainty.SCOPE_GRADIENT)
            with np.load(tmp / 'out' / 'curation.npz', allow_pickle=False) as arrays:
                self.assertIn('uncertainty_cov_te__std_gradient', arrays.files)
            self.assertEqual(before, {n: sha(tmp / n) for n in before})
            (tmp / 'cov.npz').write_bytes(b'substituted')                           # a substituted resource fails clearly
            bad = subprocess.run([sys.executable, str(ROOT / 'curate.py'), '-i', str(tmp / 'run.CDF'),
                                  '--evidence', str(tmp / 'e.json'), '--output-dir', str(tmp / 'out2')],
                                 capture_output=True, text=True, cwd=elsewhere)
            self.assertEqual(bad.returncode, 2)
            self.assertIn('does not match the declared', bad.stderr)


if __name__ == '__main__':
    unittest.main(verbosity=1)
