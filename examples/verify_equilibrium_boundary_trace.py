#!/usr/bin/python
"""Verify the LCFS outline trace of the equilibrium writer on synthetic flux maps (no IMAS, no CDF).

``imas_equilibrium.trace_boundary`` contours the TRANSP rectangular flux map ``PSIRZ`` at the
enclosed-flux level ``PLFLXA``.  Some CDFs (DIII-D 147634T16) clamp ``PSIRZ`` just inside that flux
outside the plasma, so no contour exists at exactly that level, and which curve (if any) a contouring
library returned depended on the sign convention (DD3 versus DD4).  Pinned here:

* a level inside the data range gives the circle of the synthetic plasma and records no adjustment,
* a clamped map gives the same closed curve around the axis in both sign conventions, and records an
  ``adjusted`` entry that names the level and the range,
* a map with no usable flux range produces no outline rather than an exception.

    python examples/verify_equilibrium_boundary_trace.py
"""

import ast
import logging
from pathlib import Path
import sys
import types
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from contour_lines import line_contours
from writer_base import IdsWriter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mapping_status import ADJUSTED, MappingLedger

R0, A = 1.7, 0.5
EDGE = 1.2           # enclosed flux (Wb/rad)


# Registry entries exercised here: eq_boundary_outline
def trace_boundary_function():
    """The real method from the source, run without importing IMAS."""
    path = ROOT / 'imas_equilibrium.py'
    tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'imas_equilibrium')
    func = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == 'trace_boundary')
    namespace = {'np': np, 'plt': plt}
    exec(compile(ast.fix_missing_locations(ast.Module(body=[func], type_ignores=[])), str(path), 'exec'), namespace)
    return namespace['trace_boundary']


class Stub:
    pi = np.pi

    def __init__(self, psi_sign):
        self.psi_sign = psi_sign
        self.mapping = MappingLedger()
        self.mapping.writer = 'equilibrium'
        self.logger = logging.getLogger('verify-equilibrium-boundary-trace')

    @staticmethod
    def is_valid(*values, ind=0, name=''):
        return all(v is not None and np.ndim(v) > 0 and np.size(v) > 0 and len(v) > ind for v in values)

    def get_ncvar_on_time(self, name, target='TIME3', scale=1.0):
        values = {'PLFLXA': [EDGE / (2 * np.pi)], 'RAXIS': [R0 * 100.0], 'YAXIS': [0.0]}.get(name)
        return np.array([]) if values is None else np.asarray(values) * scale


def flux_map(plateau):
    """psi = EDGE * rho^2 inside the plasma (rho = distance / A), clamped at ``plateau`` outside."""
    r = np.linspace(1.0, 2.4, 141)
    z = np.linspace(-0.9, 0.9, 129)
    rr, zz = np.meshgrid(r, z)
    psi = np.minimum(EDGE * (((rr - R0) ** 2 + zz ** 2) / A ** 2), plateau)
    return r, z, psi


def trace(psi_sign, plateau):
    r, z, psi = flux_map(plateau)
    stub = Stub(psi_sign)
    writer = IdsWriter(stub)          # the views (cdf, signs, grid, fill, ctx) over the stub, as in production
    outline = trace_boundary_function()(writer, (psi_sign * psi).ravel(), r[None, :], z[None, :], 0, len(r), len(z))
    return outline, stub.mapping


class TraceTests(unittest.TestCase):
    def check_circle(self, outline, tolerance):
        r, z = (np.asarray(v) for v in outline)
        self.assertGreater(len(r), 50)
        self.assertAlmostEqual(r[0], r[-1])
        distance = np.hypot(r - R0, z)
        self.assertLess(np.max(np.abs(distance - A)) / A, tolerance)

    def test_a_level_inside_the_range_traces_the_plasma_and_records_nothing(self):
        for sign in (+1, -1):
            with self.subTest(psi_sign=sign):
                outline, ledger = trace(sign, plateau=1.5 * EDGE)
                self.check_circle(outline, 0.01)
                self.assertEqual([e for e in ledger.entries() if e['kind'] == ADJUSTED], [])

    def test_a_clamped_map_gives_the_same_curve_in_both_sign_conventions(self):
        plateau = (1.0 - 0.002) * EDGE                  # 0.2% below the edge flux, as in DIII-D 147634T16
        curves = {}
        for sign in (+1, -1):
            outline, ledger = trace(sign, plateau)
            # the clamp edge is at A*sqrt(0.998); the curve is only known to about one grid cell (1 cm = 2% of A)
            self.check_circle(outline, 0.03)
            adjusted = [e for e in ledger.entries('equilibrium') if e['kind'] == ADJUSTED]
            self.assertEqual(len(adjusted), 1)
            self.assertEqual((adjusted[0]['source'], adjusted[0]['target']), ('PSIRZ', 'time_slice/boundary/outline'))
            self.assertIn('outside the PSIRZ range', adjusted[0]['reason'])
            curves[sign] = np.column_stack(outline)
        # the same polyline, possibly traversed the other way round from another vertex: compare as point sets
        from scipy.spatial import cKDTree
        self.assertEqual(curves[+1].shape, curves[-1].shape)
        for a, b in ((+1, -1), (-1, +1)):
            nearest, _ = cKDTree(curves[b]).query(curves[a])
            self.assertLess(float(nearest.max()), 1e-9)

    def test_the_curve_is_closed_and_encloses_the_axis(self):
        from matplotlib.path import Path as MplPath
        for sign in (+1, -1):
            outline, _ = trace(sign, (1.0 - 0.002) * EDGE)
            points = np.column_stack(outline)
            self.assertTrue(np.allclose(points[0], points[-1]))
            self.assertTrue(MplPath(points).contains_point((R0, 0.0)))

    def test_an_unusable_flux_map_gives_no_outline_instead_of_an_error(self):
        stub = Stub(+1)
        writer = IdsWriter(stub)          # the views (cdf, signs, grid, fill, ctx) over the stub, as in production
        r, z, _ = flux_map(EDGE)
        flat = np.zeros(len(r) * len(z))               # no flux range at all
        outline = trace_boundary_function()(writer, flat, r[None, :], z[None, :], 0, len(r), len(z))
        self.assertEqual(len(outline[0]), 0)


class ContourEquivalenceTests(unittest.TestCase):
    """``line_contours`` returns the vertices ``Axes.contour`` does, without building a Figure."""

    @staticmethod
    def reference(r, z, psi, level):
        fig, ax = plt.subplots()
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')          # "No contour levels were found within the data range"
                return [np.array(segment) for segment in ax.contour(r, z, psi, levels=[level]).allsegs[0]]
        finally:
            plt.close(fig)

    def assert_same(self, r, z, psi, level, label=''):
        mine, theirs = line_contours(r, z, psi, level), self.reference(r, z, psi, level)
        self.assertEqual(len(mine), len(theirs), f'{label}: number of segments')
        for k, (a, b) in enumerate(zip(mine, theirs)):
            np.testing.assert_array_equal(a, b, err_msg=f'{label}: segment {k}')

    def test_synthetic_maps_inside_and_outside_the_data_range(self):
        r, z, psi = flux_map(plateau=EDGE * 0.999)
        for level in (0.2 * EDGE, 0.5 * EDGE, 0.9 * EDGE, EDGE * 0.9985, EDGE * 2.0, -1.0):
            self.assert_same(r, z, psi, level, f'level {level}')

    def test_non_finite_values_are_masked_like_matplotlib(self):
        r, z, psi = flux_map(plateau=EDGE * 2)
        psi = psi.copy()
        psi[10:14, 20:25] = np.nan
        self.assert_same(r, z, psi, 0.5 * EDGE, 'nan block')

    def test_a_flat_field_gives_the_same_empty_result(self):
        r, z, _ = flux_map(plateau=EDGE)
        flat = np.full((len(z), len(r)), 1.5)
        self.assert_same(r, z, flat, 1.5, 'flat field')

    def test_real_flux_maps_of_the_bundled_cases(self):
        from scipy.io import netcdf_file
        checked = 0
        for path in (ROOT / 'examples' / 'data' / 'iter' / '38530P94.CDF',
                     ROOT / 'examples' / 'data' / 'nstx' / '141711Y05' / '141711Y05.CDF',
                     ROOT / 'examples' / 'data' / 'd3d' / '147634T16' / '147634T16.CDF'):
            if not path.is_file():
                continue
            with netcdf_file(str(path), 'r', mmap=False) as nc:
                if 'PSIRZ' not in nc.variables:
                    continue
                rgrid = np.array(nc.variables['RGRID'].data, dtype=float) * 0.01
                zgrid = np.array(nc.variables['ZGRID'].data, dtype=float) * 0.01
                psirz = np.array(nc.variables['PSIRZ'].data, dtype=float)
            ntime = psirz.shape[0]
            for it in sorted({0, ntime // 2, ntime - 1}):
                r = rgrid[it] if rgrid.ndim > 1 else rgrid
                z = zgrid[it] if zgrid.ndim > 1 else zgrid
                psi = psirz[it].reshape((len(z), len(r)))
                lo, hi = float(psi.min()), float(psi.max())
                for fraction in (0.3, 0.7, 0.98):
                    self.assert_same(r, z, psi, lo + fraction * (hi - lo), f'{path.name} t{it} f{fraction}')
                    checked += 1
        self.assertGreater(checked, 0, 'no bundled CDF with PSIRZ was found')


if __name__ == '__main__':
    unittest.main(verbosity=1)
