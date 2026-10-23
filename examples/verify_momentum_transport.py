#!/usr/bin/python
"""Verify the predictive angular-momentum flux (AMTR_MOD) in the core_transport solver model.

``AMTR_MOD`` is a divergence of angular-momentum flux (N m cm^-3) on X, like the observed
``AMTR_OBS``.  The IMAS ``momentum_*/flux`` node needs the surface-averaged flux on XB, obtained
as ``cumsum(div * DVOL) / SURF`` (never by copying the divergence), with the IMAS toroidal
sign.  Only an AMTR_MOD with nonzero data counts: an all-zero array just says the momentum model
was inactive.

1. source checks (SciPy only): the writer's formula equals the neoclassical writer's verified
   ``divergence_to_surface_flux``; ITER is the bundled case with nonzero AMTR_MOD,
2. (not with ``--quick``) real conversions read back from HDF5 in DD 4.1.1 and 3.39.0: the ITER
   ``transport_solver`` model carries the flux and it equals an independent evaluation at
   every time and radius; NSTX 141711Y05 (all-zero AMTR_MOD) gets no momentum flux.

    python examples/verify_momentum_transport.py [--quick]
"""

import argparse
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
DATA = ROOT / 'examples' / 'data'
ITER = DATA / 'iter' / '38530P94.CDF'
NSTX = DATA / 'nstx' / '141711Y05' / '141711Y05.CDF'

from theory_transport_mapping import NEOCLASSICAL_CHANNELS, divergence_to_surface_flux

try:
    import imas
except ImportError:      # pragma: no cover - environment dependent
    imas = None

RUN = {'quick': False}


def read(path, names):
    with netcdf_file(path, mmap=False) as nc:
        return {n: np.array(nc.variables[n][:], dtype=float) for n in names}


def independent_flux(path):
    """Surface-averaged flux in SI from AMTR_MOD, written out with explicit unit conversions."""
    v = read(path, ('AMTR_MOD', 'DVOL', 'SURF'))
    divergence = v['AMTR_MOD'] * 1.0e6          # N m cm^-3 -> N m m^-3
    volume = v['DVOL'] * 1.0e-6                 # cm^3 -> m^3
    area = v['SURF'] * 1.0e-4                   # cm^2 -> m^2
    flux = np.zeros_like(divergence)
    for it in range(divergence.shape[0]):
        enclosed = 0.0
        for ix in range(divergence.shape[1]):
            enclosed += divergence[it, ix] * volume[it, ix]
            flux[it, ix] = enclosed / area[it, ix]
    return flux


class SourceTests(unittest.TestCase):
    def test_iter_is_the_case_with_predictive_momentum_flux(self):
        if not ITER.is_file():
            self.skipTest('ITER example not present')
        self.assertTrue(np.any(read(ITER, ('AMTR_MOD',))['AMTR_MOD']))
        if NSTX.is_file():
            self.assertFalse(np.any(read(NSTX, ('AMTR_MOD',))['AMTR_MOD']))

    def test_writer_formula_equals_the_verified_neoclassical_conversion(self):
        if not ITER.is_file():
            self.skipTest('ITER example not present')
        v = read(ITER, ('AMTR_MOD', 'DVOL', 'SURF'))
        kind = next(iter(NEOCLASSICAL_CHANNELS))
        # On raw CDF units the verified helper gives N m cm^-3 * cm^3 / cm^2 * 1e4 = N m m^-2.
        verified = divergence_to_surface_flux(v['AMTR_MOD'], v['DVOL'], v['SURF'], kind)
        independent = independent_flux(ITER)
        np.testing.assert_allclose(verified, independent, rtol=1e-10, atol=1e-12 * np.abs(independent).max())
        self.assertGreater(np.abs(independent).max(), 0.0)


def convert_and_open(cdf, dd_version, outdir):
    result = subprocess.run(
        [sys.executable, str(ROOT / 'convert.py'), '-i', str(cdf), '-b', 'HDF5', '-u', str(outdir),
         '--dd-version', dd_version, '-new', '--ids', 'core_transport'],
        cwd=ROOT, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stdout[-1500:] + result.stderr[-1500:])
    match = re.search(r'ids stored in (\S+)', result.stdout)
    if not match:
        raise RuntimeError('convert.py did not report the stored entry:\n' + result.stdout[-1500:])
    path = Path(match.group(1))
    database, _major, shot, run = path.parts[-4:]
    entry = imas.DBEntry(imas.ids_defs.HDF5_BACKEND, database, int(shot), int(run), str(path.parents[3]),
                         dd_version=dd_version)
    entry.open()
    return entry


def solver_model(path, dd_version):
    with tempfile.TemporaryDirectory() as tmp:
        entry = convert_and_open(path, dd_version, tmp)
        try:
            ids = entry.get('core_transport', occurrence=0)
        finally:
            entry.close()
    models = [m for m in ids.model if str(m.identifier.name) == 'transport_solver']
    return ids, models


@unittest.skipIf(imas is None, 'IMAS-Python is required')
class ConvertedEntryTests(unittest.TestCase):
    def setUp(self):
        if RUN['quick']:
            self.skipTest('end-to-end readback skipped (--quick)')

    def _check_iter(self, dd):
        if not ITER.is_file():
            self.skipTest('ITER example not present')
        ids, models = solver_model(ITER, dd)
        self.assertEqual(len(models), 1)
        name = 'momentum_tor' if dd.startswith('3.') else 'momentum_phi'
        flux = np.array([np.asarray(getattr(p, name).flux) for p in models[0].profiles_1d])
        expected = independent_flux(ITER)
        self.assertEqual(flux.shape, expected.shape)
        sign = np.sign(np.sum(flux * expected))
        self.assertIn(sign, (-1.0, 1.0))
        np.testing.assert_allclose(flux, sign * expected, rtol=1e-6, atol=1e-6 * np.abs(expected).max())
        # zero where the source divergence is zero (ITER's first three times), nonzero after
        active = np.abs(expected).sum(axis=1) > 0
        self.assertTrue(active.any() and not active.all())
        np.testing.assert_array_equal(np.abs(flux).sum(axis=1) > 0, active)
        # a divergence copied into the flux node would not be a cumulative integral over SURF
        divergence = read(ITER, ('AMTR_MOD',))['AMTR_MOD'] * 1.0e6
        self.assertGreater(np.abs(np.abs(flux) - np.abs(divergence)).max(), 1e-3 * np.abs(divergence).max())

    def test_iter_dd4(self):
        self._check_iter('4.1.1')

    def test_iter_dd3(self):
        self._check_iter('3.39.0')

    def test_inactive_amtr_mod_writes_no_momentum_flux(self):
        if not NSTX.is_file():
            self.skipTest('NSTX example not present')
        ids, models = solver_model(NSTX, '4.1.1')
        for model in models:                      # an (all-zero) MOD solver model exists, as before
            for profile in model.profiles_1d:
                self.assertEqual(np.asarray(profile.momentum_phi.flux).size, 0)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='skip the HDF5 end-to-end readback')
    args, rest = parser.parse_known_args(argv)
    RUN['quick'] = args.quick
    result = unittest.main(argv=[sys.argv[0]] + rest, exit=False, verbosity=1).result
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
