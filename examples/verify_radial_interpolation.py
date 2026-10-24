#!/usr/bin/python
"""Verify the cached radial interpolation against the original per-slice ``np.interp`` loop.

1. unit cases: random monotone grids, targets inside, on the nodes and outside the source grid,
   the axis point (prepended only when the first node is positive), float32 input, per-time and
   scalar axis values, NaN input, non-monotone grids (fallback), shape guards, cache keyed by
   content (same shape but different grid must not share an operator), bounded cache,
2. (not with ``--quick``) every ``(TIME3, X)`` and ``(TIME3, XB)`` variable of every bundled CDF,
   mapped X -> XB and XB -> X, with and without an axis value, equals the original loop.

    python examples/verify_radial_interpolation.py [--quick]
"""

import argparse
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import radial_interpolation as ri

RUN = {'quick': False}
RTOL, ATOL = 1e-12, 1e-12


def reference(values, x_from, x_to, axis_value=None):
    """The implementation that was in nc2imas before the operator cache."""
    v, xf, xt = np.asarray(values), np.asarray(x_from), np.asarray(x_to)
    if v.size == 0 or xf.size == 0 or xt.size == 0 or v.ndim != 2:
        return v
    if xf.ndim != 2 or xt.ndim != 2 or v.shape[0] != xf.shape[0] or v.shape[0] != xt.shape[0]:
        return v
    out = np.empty((v.shape[0], xt.shape[1]), dtype=np.float64)
    av = None if axis_value is None else np.asarray(axis_value)
    for it in range(v.shape[0]):
        xx, yy = xf[it], v[it]
        if axis_value is not None and len(xx) and xx[0] > 0.0:
            a0 = float(av[it] if av.ndim else av)
            xx = np.concatenate(([0.0], xx))
            yy = np.concatenate(([a0], yy))
        out[it] = np.interp(xt[it], xx, yy)
    return out


def grids(rng, nt=6, n=40, m=41, first=0.0125):
    """Increasing zone-centre (X) and boundary (XB) grids that differ from slice to slice."""
    xb = np.cumsum(rng.uniform(0.5, 1.5, (nt, m)), axis=1)
    xb = xb / xb[:, -1:]
    x = 0.5 * (xb[:, :-1] + xb[:, 1:])
    return x[:, :n], xb


class UnitTests(unittest.TestCase):
    def setUp(self):
        ri.clear_cache()
        self.rng = np.random.default_rng(7)

    def same(self, values, x_from, x_to, axis=None):
        np.testing.assert_allclose(ri.interpolate_radial_profiles(values, x_from, x_to, axis),
                                   reference(values, x_from, x_to, axis), rtol=RTOL, atol=ATOL)

    def test_matches_np_interp_x_to_xb_and_back_with_and_without_axis(self):
        x, xb = grids(self.rng, n=40, m=41)
        values = self.rng.normal(size=x.shape)
        self.same(values, x, xb)                                 # target XB beyond the X range: clamped
        self.same(values, x, xb, axis=values[:, 0])
        self.same(values, x, xb, axis=1.7)                       # scalar axis value
        vb = self.rng.normal(size=xb.shape)
        self.same(vb, xb, x)
        self.same(vb, xb, x, axis=0.0)
        self.assertGreater(ri.STATS['misses'], 0)

    def test_targets_on_the_nodes_and_far_outside(self):
        x, _ = grids(self.rng)
        values = self.rng.normal(size=x.shape)
        targets = np.concatenate([x, x[:, :1] - 1.0, x[:, -1:] + 1.0], axis=1)       # nodes and beyond both ends
        self.same(values, x, targets)
        self.same(values, x, targets, axis=values[:, 0])

    def test_axis_point_only_where_the_first_node_is_positive(self):
        x, xb = grids(self.rng)
        x[1, 0] = 0.0                       # this slice already starts at the axis: no extra point
        x[1, 1:] = np.maximum(x[1, 1:], 1e-3)
        x[3, 0] = -0.2                      # and one below it
        values = self.rng.normal(size=x.shape)
        self.same(values, x, xb, axis=values[:, 0] + 3.0)

    def test_float32_and_integer_values(self):
        x, xb = grids(self.rng)
        self.same(self.rng.normal(size=x.shape).astype(np.float32), x.astype(np.float32), xb.astype(np.float32))
        self.same(self.rng.integers(0, 9, size=x.shape), x, xb)

    def test_nan_and_inf_inputs_keep_the_original_behaviour(self):
        x, xb = grids(self.rng)
        values = self.rng.normal(size=x.shape)
        values[2, 5] = np.nan
        values[4, 0] = np.inf
        got = ri.interpolate_radial_profiles(values, x, xb)
        want = reference(values, x, xb)
        np.testing.assert_array_equal(np.isnan(got), np.isnan(want))
        finite = np.isfinite(want)
        np.testing.assert_allclose(got[finite], want[finite], rtol=RTOL, atol=ATOL)
        xb_bad = xb.copy()
        xb_bad[1, 3] = np.nan
        np.testing.assert_array_equal(np.isnan(ri.interpolate_radial_profiles(values, x, xb_bad)),
                                      np.isnan(reference(values, x, xb_bad)))

    def test_non_monotone_grids_fall_back(self):
        x, xb = grids(self.rng)
        x_bad = x.copy()
        x_bad[2, 10], x_bad[2, 11] = x_bad[2, 11], x_bad[2, 10]
        values = self.rng.normal(size=x.shape)
        self.same(values, x_bad, xb)
        self.assertGreater(ri.STATS['fallbacks'], 0)
        x_dup = x.copy()
        x_dup[0, 5] = x_dup[0, 6]
        self.same(values, x_dup, xb)

    def test_guards_return_the_input_unchanged(self):
        x, xb = grids(self.rng)
        values = self.rng.normal(size=x.shape)
        for args in ((np.array([]), x, xb), (values[0], x, xb), (values, x[:3], xb), (values, x[0], xb)):
            self.assertIs(ri.interpolate_radial_profiles(*args), args[0]) if args[0].size == 0 or args[0].ndim != 2 else None
            np.testing.assert_array_equal(ri.interpolate_radial_profiles(*args), reference(*args))

    def test_cache_is_keyed_by_content_not_shape_and_is_bounded(self):
        x1, xb1 = grids(self.rng)
        x2, xb2 = grids(np.random.default_rng(99))
        values = self.rng.normal(size=x1.shape)
        a = ri.interpolate_radial_profiles(values, x1, xb1)
        b = ri.interpolate_radial_profiles(values, x2, xb2)           # same shapes, different grids
        self.assertEqual(ri.STATS['misses'], 2)
        np.testing.assert_allclose(a, reference(values, x1, xb1), rtol=RTOL, atol=ATOL)
        np.testing.assert_allclose(b, reference(values, x2, xb2), rtol=RTOL, atol=ATOL)
        self.assertFalse(np.allclose(a, b))
        ri.interpolate_radial_profiles(self.rng.normal(size=x1.shape), x1, xb1)   # same grids, other profile
        self.assertEqual(ri.STATS['hits'], 1)
        x1_mod = x1.copy()
        x1_mod[0, 3] += 1e-9                                          # a modified array is a different key
        ri.interpolate_radial_profiles(values, x1_mod, xb1)
        self.assertEqual(ri.STATS['misses'], 3)
        for seed in range(20):
            xs, xbs = grids(np.random.default_rng(1000 + seed))
            ri.interpolate_radial_profiles(self.rng.normal(size=xs.shape), xs, xbs)
        self.assertLessEqual(len(ri._CACHE), ri._CACHE_SIZE)

    def test_result_is_a_fresh_array(self):
        x, xb = grids(self.rng)
        values = self.rng.normal(size=x.shape)
        first = ri.interpolate_radial_profiles(values, x, xb)
        first[:] = 0.0
        second = ri.interpolate_radial_profiles(values, x, xb)
        np.testing.assert_allclose(second, reference(values, x, xb), rtol=RTOL, atol=ATOL)


@unittest.skipIf(False, '')
class RealDataTests(unittest.TestCase):
    def test_every_radial_variable_of_every_bundled_cdf(self):
        if RUN['quick']:
            self.skipTest('--quick')
        cdfs = sorted((ROOT / 'examples' / 'data').rglob('*.CDF'))
        self.assertTrue(cdfs)
        compared = 0
        for path in cdfs:
            with netcdf_file(str(path), 'r', mmap=False) as nc:
                if 'X' not in nc.variables or 'XB' not in nc.variables:
                    continue
                x, xb = np.array(nc.variables['X'].data), np.array(nc.variables['XB'].data)
                if x.ndim != 2 or xb.ndim != 2 or x.shape != xb.shape:
                    continue
                for name, var in nc.variables.items():
                    dims = tuple(var.dimensions)
                    if dims not in (('TIME3', 'X'), ('TIME3', 'XB')) or var.data.dtype.kind not in 'fiu':
                        continue
                    values = np.array(var.data)
                    src, dst = (x, xb) if dims[1] == 'X' else (xb, x)
                    for axis in (None, values[:, 0], 0.0):
                        got = ri.interpolate_radial_profiles(values, src, dst, axis)
                        want = reference(values, src, dst, axis)
                        np.testing.assert_allclose(got, want, rtol=RTOL, atol=ATOL, equal_nan=True,
                                                   err_msg=f'{path.name} {name} axis={axis is not None}')
                        compared += 1
        self.assertGreater(compared, 1000)
        print(f'\ncompared {compared} real interpolations over {len(cdfs)} CDFs; '
              f'cache {ri.STATS}')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='unit cases only')
    args, rest = parser.parse_known_args(argv)
    RUN['quick'] = args.quick
    result = unittest.main(argv=[sys.argv[0]] + rest, exit=False, verbosity=1).result
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
