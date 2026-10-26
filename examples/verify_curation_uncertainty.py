#!/usr/bin/python
"""Check curation.uncertainty on SYNTHETIC data only (made-up numbers, no calibration claim).

    python examples/verify_curation_uncertainty.py
"""

import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from curation import uncertainty as U  # noqa: E402

X = np.array([0.0, 0.1, 0.35, 0.5, 0.9])  # synthetic nonuniform grid


def spec(cov, x=X, sources=(), units='keV', estimator=U.ESTIMATOR):
    return U.CovarianceSpec('synthetic_te', units, x, cov, U.ORDERING, estimator, 'synthetic test', sources)


def corr_cov(x=X, sigma=0.3, length=0.4):
    return sigma ** 2 * np.exp(-np.abs(x[:, None] - x[None, :]) / length)


class OperatorTests(unittest.TestCase):
    def test_operator_equals_numpy_gradient_exactly(self):
        rng = np.random.default_rng(1)
        for _ in range(20):
            x = np.cumsum(rng.uniform(0.01, 0.5, 12))
            y = rng.normal(size=12)
            d = U.derivative_operator(x)
            np.testing.assert_allclose(d @ y, np.gradient(y, x, edge_order=1), rtol=1e-13, atol=1e-13)

    def test_edge_order_two_refused(self):
        with self.assertRaises(ValueError):
            U.derivative_operator(X, edge_order=2)

    def test_bad_grid_refused_not_sorted(self):
        for x in ([0.0, 0.5, 0.4], [0.0, 0.0, 1.0], [0.0, np.nan, 1.0], [1.0]):
            with self.assertRaises(ValueError):
                U.derivative_operator(x)

    def test_constant_gives_zero_row_sums(self):
        np.testing.assert_allclose(U.derivative_operator(X) @ np.ones(5), 0.0, atol=1e-13)


class PropagationTests(unittest.TestCase):
    def propagate(self, cov, y=None, x=X, scale=1.0):
        y = 2.0 - X if y is None else y
        return U.propagate_gradient([spec(cov, x)], y[None], x[None], unit_scale=scale)

    def test_linear_profile_exact_gradient_and_zero_covariance(self):
        r = self.propagate(np.zeros((5, 5)), y=3.0 - 2.0 * X)
        np.testing.assert_allclose(r.mean_gradient[0], -2.0, rtol=1e-13)
        np.testing.assert_array_equal(r.covariance[0], np.zeros((5, 5)))
        np.testing.assert_array_equal(r.std_gradient[0], 0.0)

    def test_diagonal_covariance_matches_hand_formula(self):
        var = np.array([1.0, 2.0, 3.0, 4.0, 5.0]) * 1e-2
        r = self.propagate(np.diag(var))
        # interior point i: sigma_g^2 = (h_l^2 var_{i+1} + (h_r^2-h_l^2)^2 var_i + h_r^2 var_{i-1})/(h_l+h_r)^2 / (h_l h_r)^2
        i = 2
        hl, hr = X[i] - X[i - 1], X[i + 1] - X[i]
        a = -hr / (hl * (hl + hr))
        b = (hr - hl) / (hl * hr)
        c = hl / (hr * (hl + hr))
        expected = a * a * var[i - 1] + b * b * var[i] + c * c * var[i + 1]
        self.assertAlmostEqual(r.covariance[0][i, i], expected, delta=1e-13 * expected)
        # edge point 0 with edge_order=1: (y1-y0)/h
        h = X[1] - X[0]
        self.assertAlmostEqual(r.covariance[0][0, 0], (var[0] + var[1]) / h ** 2, delta=1e-12)

    def test_correlated_covariance_matches_hand_computed_dcdt_small(self):
        x = np.array([0.0, 1.0, 3.0, 4.0])
        d = np.array([[-1.0, 1.0, 0.0, 0.0], [0.0] * 4, [0.0] * 4, [0.0, 0.0, -1.0, 1.0]])
        # interior rows by the 3-point nonuniform formula, computed by hand:
        # row1: h_l=1,h_r=2 -> a=-2/3, b=(2-1)/2=1/2, c=1/6 ; row2: h_l=2,h_r=1 -> a=-1/6, b=-1/2, c=2/3
        d[1] = [-2.0 / 3.0, 1.0 / 2.0, 1.0 / 6.0, 0.0]
        d[2] = [0.0, -1.0 / 6.0, -1.0 / 2.0, 2.0 / 3.0]
        c = np.array([[1.0, 0.5, 0.2, 0.0],
                      [0.5, 2.0, 0.3, 0.1],
                      [0.2, 0.3, 1.5, 0.4],
                      [0.0, 0.1, 0.4, 1.0]])
        r = U.propagate_gradient([spec(c, x)], np.array([[1.0, 2.0, 0.5, 4.0]]), x[None])
        np.testing.assert_allclose(U.derivative_operator(x), d, atol=1e-14)
        np.testing.assert_allclose(r.covariance[0], d @ c @ d.T, rtol=1e-13, atol=1e-14)
        np.testing.assert_allclose(r.mean_gradient[0], d @ np.array([1.0, 2.0, 0.5, 4.0]), atol=1e-13)

    def test_unit_scale_squares_on_covariance(self):
        c = corr_cov()
        a = self.propagate(c)
        b = self.propagate(c, scale=1e-3)
        np.testing.assert_allclose(b.covariance[0], 1e-6 * a.covariance[0], rtol=1e-13)
        np.testing.assert_allclose(b.mean_gradient, 1e-3 * a.mean_gradient, rtol=1e-13)

    def test_shared_additive_offset_gives_zero_gradient_variance(self):
        sigma = 0.7
        c = sigma ** 2 * np.ones((5, 5))
        r = self.propagate(c)
        self.assertLessEqual(np.max(np.abs(r.covariance[0])), 1e-12 * sigma ** 2 / np.min(np.diff(X)) ** 2)
        independent = c + 0.01 * np.eye(5)
        self.assertGreater(self.propagate(independent).std_gradient[0].min(), 0.05)

    def test_singular_but_valid_covariance_accepted(self):
        v = np.array([1.0, 0.5, -0.3, 0.2, 0.1])
        c = np.outer(v, v)  # rank 1, four zero eigenvalues
        self.assertEqual(self.propagate(c).evaluated.all(), True)

    def test_invalid_covariances_rejected(self):
        c = corr_cov()
        asym = c.copy()
        asym[0, 1] += 1e-3
        negative = np.diag([1.0, 1.0, 1.0, 1.0, -1e-3])
        nan = c.copy()
        nan[2, 2] = np.nan
        for bad in (asym, negative, nan, np.ones((4, 4)), np.ones(5)):
            with self.assertRaises(ValueError):
                U.validate_covariance(bad, 5)

    def test_psd_within_tolerance_accepted_without_modification(self):
        c = np.diag([1.0, 1.0, 1.0, 1.0, -1e-12])  # -1e-12 > -1e-10 * 1
        out, report = U.validate_covariance(c, 5)
        np.testing.assert_array_equal(out, c)
        self.assertEqual(report['regularization_absolute'], 0.0)

    def test_regularization_only_on_request_and_reported(self):
        c = np.outer(np.ones(5), np.ones(5))
        out, report = U.validate_covariance(c, 5, regularization=1e-8)
        self.assertEqual(report['regularization_relative'], 1e-8)
        np.testing.assert_allclose(out, c + 1e-8 * np.eye(5), atol=1e-15)
        with self.assertRaises(ValueError):
            U.validate_covariance(c, 5, regularization=1e-2)
        # regularization never rescues a non-PSD matrix
        with self.assertRaises(ValueError):
            U.validate_covariance(np.diag([1, 1, 1, 1, -0.5]), 5, regularization=1e-6)

    def test_missing_covariance_is_unknown_not_zero(self):
        with self.assertRaises(ValueError):
            U.validate_covariance(None, 5)
        r = U.propagate_gradient([None, spec(corr_cov())], np.tile(2 - X, (2, 1)), np.tile(X, (2, 1)))
        self.assertIsNone(r.covariance[0])
        self.assertTrue(np.all(np.isnan(r.std_gradient[0])))
        self.assertFalse(r.evaluated[0].any())
        self.assertTrue(r.evaluated[1].all())
        self.assertEqual(r.as_dict()['unknown_slices'], [0])
        self.assertIsNone(r.as_dict()['std_gradient'][0][0])

    def test_cross_time_not_implemented(self):
        y, x = np.tile(2 - X, (2, 1)), np.tile(X, (2, 1))
        with self.assertRaises(NotImplementedError):
            U.propagate_gradient([spec(corr_cov())] * 2, y, x, cross_time='correlated')
        with self.assertRaises(NotImplementedError):
            U.propagate_gradient([spec(np.eye(10))] * 2, y, x)

    def test_estimator_coordinate_and_count_mismatch_refused(self):
        y, x = (2 - X)[None], X[None]
        with self.assertRaises(ValueError):
            U.propagate_gradient([spec(corr_cov(), estimator='central')], y, x)
        with self.assertRaises(ValueError):
            U.propagate_gradient([spec(corr_cov(), x=X + 0.01)], y, x)
        with self.assertRaises(ValueError):
            U.propagate_gradient([spec(corr_cov())] * 2, y, x)

    def test_per_time_grids_use_their_own_operator(self):
        x2 = np.array([0.0, 0.3, 0.4, 0.8, 1.0])
        r = U.propagate_gradient([spec(corr_cov(X)), spec(corr_cov(x2), x=x2)],
                                 np.array([2 - X, 2 - x2]), np.array([X, x2]))
        np.testing.assert_allclose(r.covariance[1], U.derivative_operator(x2) @ corr_cov(x2) @ U.derivative_operator(x2).T)
        self.assertGreater(np.max(np.abs(r.covariance[0] - r.covariance[1])), 1e-3)

    def test_dense_size_limit(self):
        n = U.MAX_DENSE_POINTS + 1
        with self.assertRaises(ValueError):
            U.validate_covariance(np.zeros((n, n)), n)
        with self.assertRaises(ValueError):
            U.derivative_operator(np.arange(n, dtype=float))

    def test_scope_and_metadata(self):
        r = self.propagate(corr_cov())
        self.assertEqual(r.scope, 'profile_gradient_conditional_on_supplied_reconstruction')
        self.assertEqual(r.as_dict()['cross_time'], 'unmodeled')
        self.assertEqual(r.as_dict()['scope'], r.scope)


class SamplingTests(unittest.TestCase):
    def test_distribution_must_be_declared_and_gaussian(self):
        for dist in (None, 'student_t', 'Gaussian'):
            with self.assertRaises(ValueError):
                U.sample_from_covariance(np.zeros(5), corr_cov(), 10, 1, distribution=dist)

    def test_seed_required_and_validated(self):
        with self.assertRaises(TypeError):
            U.sample_from_covariance(np.zeros(5), corr_cov(), 10, distribution='gaussian')
        for seed in (-1, 1.5, True):
            with self.assertRaises(ValueError):
                U.sample_from_covariance(np.zeros(5), corr_cov(), 10, seed, distribution='gaussian')

    def test_same_seed_bit_identical_and_batch_independent(self):
        kw = dict(distribution='gaussian')
        a = U.sample_from_covariance(np.zeros(5), corr_cov(), 1000, 7, batch_size=1000, **kw)
        b = U.sample_from_covariance(np.zeros(5), corr_cov(), 1000, 7, batch_size=64, **kw)
        c = U.sample_from_covariance(np.zeros(5), corr_cov(), 1000, 8, **kw)
        self.assertTrue(np.array_equal(a.samples, b.samples))
        self.assertFalse(np.array_equal(a.samples, c.samples))
        self.assertEqual(a.record['seed'], 7)
        self.assertEqual(a.record['algorithm'], 'numpy.random.Generator(PCG64)')
        self.assertEqual(a.record['numpy_version'], np.__version__)
        self.assertEqual(a.record['failed_members'], [])
        json.dumps(a.as_dict())

    def test_limits_enforced(self):
        with self.assertRaises(ValueError):
            U.sample_from_covariance(np.zeros(5), corr_cov(), 11, 1, distribution='gaussian', max_draws=10)
        with self.assertRaises(ValueError):
            U.sample_from_covariance(np.zeros(5), corr_cov(), 10, 1, distribution='gaussian', max_bytes=399)
        U.sample_from_covariance(np.zeros(5), corr_cov(), 10, 1, distribution='gaussian', max_bytes=400)

    def test_singular_covariance_samples_stay_in_range_space(self):
        v = np.array([1.0, 0.5, -0.3, 0.2, 0.1])
        s = U.sample_from_covariance(np.zeros(5), np.outer(v, v), 200, 3, distribution='gaussian').samples
        # every member is a multiple of v: orthogonal complement component vanishes
        proj = s - np.outer(s @ v / (v @ v), v)
        self.assertLess(np.max(np.abs(proj)), 1e-12)

    def test_monte_carlo_agrees_with_analytic_within_standard_error(self):
        n, seed = 40000, 20240601
        c = corr_cov()
        mean = 2.0 - X
        s = U.sample_from_covariance(mean, c, n, seed, distribution='gaussian').samples
        g = U.gradient_from_samples(s, X)
        analytic = U.propagate_gradient([spec(c)], mean[None], X[None])
        d = U.derivative_operator(X)
        cg = analytic.covariance[0]
        # Gaussian: var(S_ij) ~ (C_ij^2 + C_ii C_jj)/(n-1); tolerance is 5 standard errors
        se = np.sqrt((cg ** 2 + np.outer(np.diag(cg), np.diag(cg))) / (n - 1))
        self.assertTrue(np.all(np.abs(g.covariance - cg) <= 5 * se), np.max(np.abs(g.covariance - cg) / se))
        se_mean = np.sqrt(np.diag(cg) / n)
        self.assertTrue(np.all(np.abs(g.mean_gradient - d @ mean) <= 5 * se_mean))
        self.assertEqual(g.scope, U.SCOPE_GRADIENT)

    def test_gradient_from_samples_records_rejections(self):
        s = np.tile(2 - X, (4, 1))
        s[1, 2] = np.nan
        s[3] += 0.1
        xs = np.tile(X, (4, 1))
        xs[2, 1], xs[2, 2] = xs[2, 2], xs[2, 1]  # displaced/unsorted coordinate: rejected, not sorted
        before = xs.copy()
        g = U.gradient_from_samples(s, xs, member_ids=[10, 11, 12, 13])
        self.assertEqual(g.member_ids_used.tolist(), [10, 13])
        self.assertEqual(g.rejected, [{'member_id': 11, 'reason': 'non-finite profile value'},
                                      {'member_id': 12, 'reason': 'coordinate not strictly increasing'}])
        np.testing.assert_array_equal(xs, before)
        np.testing.assert_allclose(g.mean_gradient, -1.0, atol=1e-13)
        self.assertEqual(g.as_dict()['rejected'], g.rejected)

    def test_gradient_from_samples_unknown_covariance_and_all_rejected(self):
        s = np.tile(2 - X, (2, 1))
        s[0, 0] = np.inf
        g = U.gradient_from_samples(s, X)
        self.assertIsNone(g.covariance)
        self.assertFalse(g.as_dict()['covariance_known'])
        with self.assertRaises(ValueError):
            U.gradient_from_samples(np.full((2, 5), np.nan), X)

    def test_temperatures_not_clipped(self):
        s = np.tile(X * 0 - 1.0, (3, 1))  # negative synthetic values
        g = U.gradient_from_samples(s, X)
        np.testing.assert_allclose(g.gradients, 0.0, atol=1e-13)
        self.assertEqual(g.gradients.shape, (3, 5))

    def test_covariance_from_samples_matches_numpy(self):
        s = np.random.default_rng(0).normal(size=(50, 5))
        np.testing.assert_allclose(U.covariance_from_samples(s), np.cov(s.T), rtol=1e-13)
        with self.assertRaises(ValueError):
            U.covariance_from_samples(s[:1])


class ScenarioTests(unittest.TestCase):
    def test_scenarios_and_double_counting(self):
        y = 2 - X
        sc = U.calibration_offset_scenarios(y, {'zero_shift': {'kind': 'additive', 'value': 0.1},
                                                'gain': {'kind': 'multiplicative', 'value': 1.05}})
        np.testing.assert_allclose(sc.scenarios['zero_shift'], y + 0.1)
        np.testing.assert_allclose(sc.scenarios['gain'], 1.05 * y)
        self.assertEqual(sc.scope, U.SCOPE_GRADIENT)
        with self.assertRaises(ValueError):
            U.calibration_offset_scenarios(y, {'zero_shift': {'kind': 'additive', 'value': 0.1}},
                                           included_sources=['zero_shift'])
        with self.assertRaises(ValueError):
            U.calibration_offset_scenarios(y, {'a': {'kind': 'additive', 'value': 0.1, 'source': 'cal'}},
                                           included_sources=('cal',))
        with self.assertRaises(ValueError):
            U.calibration_offset_scenarios(y, {'a': {'kind': 'other', 'value': 0.1}})

    def test_included_sources_carried_to_result(self):
        r = U.propagate_gradient([spec(corr_cov(), sources=('thomson_calibration',))],
                                 (2 - X)[None], X[None])
        self.assertEqual(r.included_sources, ['thomson_calibration'])


class RatioTests(unittest.TestCase):
    def ratio(self, num, den, **kw):
        return U.ratio_statistics(num, den, 'W/m^2', 'keV/m', 'synthetic heat flux over gradient', 1e-3, **kw)

    def test_basic_quantiles_no_clipping(self):
        den = np.array([1.0, 2.0, 4.0, 1e-4, 1.0])
        num = np.array([2.0, 2.0, 2.0, 2.0, 2.0])
        r = self.ratio(num, den)
        self.assertEqual(r.scope, 'frozen_physics_sensitivity_of_defined_quantity')
        self.assertEqual(r.n_defined, 5)
        self.assertAlmostEqual(r.near_zero_denominator_fraction, 0.2)
        np.testing.assert_allclose(r.quantiles, np.percentile(num / den, [5, 50, 95]))
        self.assertEqual(r.quantiles[2], np.percentile(num / den, 95))  # 2e4 member retained
        self.assertGreater(r.quantiles[2], 1e3)
        self.assertFalse(r.sign_ambiguous)
        self.assertEqual(r.units, '(W/m^2)/(keV/m)')

    def test_undefined_accounting(self):
        num = np.array([1.0, np.nan, 1.0, 1.0, 1.0, 1.0])
        den = np.array([1.0, 1.0, 0.0, np.nan, np.inf, 2.0])
        r = self.ratio(num, den, numerator_ids=np.arange(6), denominator_ids=np.arange(6))
        # undefined: nan numerator, den == 0, nan den, inf den
        self.assertEqual(r.n_defined, 2)
        self.assertEqual(r.undefined_ids, [1, 2, 3, 4])
        self.assertAlmostEqual(r.undefined_fraction, 4 / 6)
        d = r.as_dict()
        json.dumps(d)
        self.assertEqual(d['n_members'], 6)

    def test_all_undefined_gives_none_quantiles(self):
        r = self.ratio(np.ones(3), np.zeros(3))
        self.assertIsNone(r.quantiles)
        self.assertEqual(r.undefined_fraction, 1.0)

    def test_denominator_straddling_zero_is_sign_ambiguous(self):
        den = np.array([-1.0, 1.0, 2.0, 3.0])
        r = self.ratio(np.ones(4), den)
        self.assertTrue(r.sign_ambiguous)
        self.assertAlmostEqual(r.sign_ambiguity_fraction, 0.25)

    def test_pairing_and_definition_required(self):
        with self.assertRaises(ValueError):
            self.ratio(np.ones(3), np.ones(3), numerator_ids=[0, 1, 2], denominator_ids=[0, 2, 1])
        with self.assertRaises(ValueError):
            self.ratio(np.ones(3), np.ones(4))
        with self.assertRaises(ValueError):
            U.ratio_statistics(np.ones(3), np.ones(3), 'a', 'b', '  ', 0.1)
        with self.assertRaises(ValueError):
            U.ratio_statistics(np.ones(3), np.ones(3), '', 'b', 'defined', 0.1)


class IntegrityTests(unittest.TestCase):
    def test_inputs_not_mutated(self):
        y = 2.0 - X
        c = corr_cov()
        s = np.tile(y, (4, 1)) + 0.01 * np.arange(4)[:, None]
        for a in (y, c, s, X):
            a.setflags(write=False)
        before = [a.copy() for a in (y, c, s, X)]
        sp = spec(c)
        U.propagate_gradient([sp], y[None], X[None])
        U.gradient_from_samples(s, X)
        U.sample_from_covariance(y, c, 20, 1, distribution='gaussian')
        U.validate_covariance(c, 5, regularization=1e-9)
        U.calibration_offset_scenarios(y, {'o': {'kind': 'additive', 'value': 1.0}})
        for a, b in zip((y, c, s, X), before):
            np.testing.assert_array_equal(a, b)
        np.testing.assert_array_equal(sp.covariance, before[1])

    def test_npz_roundtrip_and_checksums(self):
        r = U.propagate_gradient([None, spec(corr_cov())], np.tile(2 - X, (2, 1)), np.tile(X, (2, 1)))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'grad.npz'
            sums = r.save_npz(path)
            with np.load(path, allow_pickle=False) as z:
                stored = json.loads(str(z['checksums_json']))
                self.assertEqual(stored, sums)
                for key, digest in stored.items():
                    self.assertEqual(U._sha256(z[key]), digest, key)
                np.testing.assert_array_equal(z['coordinate'], r.coordinate)
                np.testing.assert_array_equal(z['evaluated'], r.evaluated)
                np.testing.assert_array_equal(z['covariance_t0001'], r.covariance[1])
                self.assertNotIn('covariance_t0000', z.files)
                self.assertEqual(str(z['units']), 'keV')
                self.assertEqual(json.loads(str(z['meta_json']))['scope'], U.SCOPE_GRADIENT)

    def test_sample_results_save_npz(self):
        s = U.sample_from_covariance(np.zeros(5), corr_cov(), 30, 5, distribution='gaussian')
        with tempfile.TemporaryDirectory() as tmp:
            sums = s.save_npz(Path(tmp) / 's.npz')
            g = U.gradient_from_samples(s.samples, X)
            sums2 = g.save_npz(Path(tmp) / 'g.npz')
            self.assertEqual(U._sha256(s.samples), sums['samples'])
            self.assertIn('covariance', sums2)

    def test_synthetic_example_labelled(self):
        ex = U.synthetic_example()
        self.assertEqual(ex['label'], 'synthetic, not diagnostic data')
        self.assertEqual(ex['result'].scope, U.SCOPE_GRADIENT)
        json.dumps(ex['result'].as_dict())

    def test_no_forbidden_scope_labels(self):
        self.assertNotIn('ensemble', U.SCOPE_GRADIENT + U.SCOPE_RATIO)


if __name__ == '__main__':
    unittest.main(verbosity=1)
