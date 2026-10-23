#!/usr/bin/python
"""Verify the B0 used to normalize parallel current densities and rho_tor.

IMAS defines ``j_parallel = average(J.B)/B0`` and ``rho_tor = sqrt(phi/(pi*B0))`` with B0 the
IDS's own ``vacuum_toroidal_field/b0``.  The converter divided by TRANSP's ``BZ`` (the field
at a fixed radius, which differs from ``BZXR/r0`` by up to ~40% in the bundled cases), so the
exported ``j_parallel`` and ``rho_tor`` did not match the exported ``b0``.

1. static policy: ``BZ`` is read only for the three ``beta_tor_norm`` calculations, whose beta
   is TRANSP's own definition (``BZXR**2 * <1/R**2>``) and is a separate, open question,
2. (not with ``--quick``) a real conversion of NSTX 142202U01 (the axis moves 17%) read back
   from HDF5 in DD 4.1.1 (and 3.39.0 for core_sources/core_profiles): the invariants that do
   not depend on the normalization choice hold with the IDS's own ``b0`` -
   ``j_parallel * B0 = <J.B>`` (TRANSP ``PLJB``/``PLJBOH``/``CURB*CONPLJB``) and
   ``pi * B0 * rho_tor**2 = toroidal flux``, and ``r0 * b0`` agrees across IDSs and the TF IDS.

    python examples/verify_parallel_current.py [--quick]
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
CASE = ROOT / 'examples' / 'data' / 'nstx' / '142202U01' / '142202U01.CDF'

try:
    import imas
except ImportError:      # pragma: no cover - environment dependent
    imas = None

RUN = {'quick': False}
# beta_tor_norm keeps TRANSP's BZ: its beta_tor is TRANSP's own (BZXR^2 <1/R^2>) quantity.
BZ_READS_ALLOWED = {'imas_core_profiles.py': 1, 'imas_equilibrium.py': 1, 'imas_summary.py': 1}


class StaticPolicyTests(unittest.TestCase):
    def test_bz_is_only_read_for_beta_tor_norm(self):
        counts = {}
        for path in sorted(ROOT.glob('*.py')):
            text = path.read_text(encoding='utf-8')
            found = len(re.findall(r"""(?:get_ncvar|_get)\(\s*['"]BZ['"]""", text))
            if found:
                counts[path.name] = found
        self.assertEqual(counts, BZ_READS_ALLOWED,
                         'BZ must not normalize j_parallel/rho_tor (use nc2imas.vacuum_b0_magnitude)')


def read(path, names):
    with netcdf_file(path, mmap=False) as nc:
        return {n: np.array(nc.variables[n][:], dtype=float) for n in names}


def convert_and_open(cdf, ids_names, dd_version, outdir):
    result = subprocess.run(
        [sys.executable, str(ROOT / 'convert.py'), '-i', str(cdf), '-b', 'HDF5', '-u', str(outdir),
         '--dd-version', dd_version, '-new', '--ids', ','.join(ids_names)],
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


def b0_of(ids):
    return np.asarray(ids.vacuum_toroidal_field.b0, dtype=float)


def assert_same_up_to_constant_sign(actual, expected, rtol, label):
    """``actual == s * expected`` for one global sign ``s`` (the IMAS toroidal convention)."""
    sign = np.sign(np.sum(actual * expected))
    assert sign in (-1.0, 1.0), label
    scale = np.abs(expected).max()
    np.testing.assert_allclose(actual, sign * expected, rtol=rtol, atol=rtol * scale, err_msg=label)


@unittest.skipIf(imas is None, 'IMAS-Python is required')
class ConvertedEntryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if RUN['quick'] or not CASE.is_file():
            raise unittest.SkipTest('end-to-end readback skipped (--quick or case not present)')
        cls.tmp = tempfile.TemporaryDirectory()
        cls.entry = convert_and_open(CASE, ['tf', 'core_profiles', 'core_sources', 'equilibrium', 'summary'],
                                     '4.1.1', cls.tmp.name)
        cls.cdf = read(CASE, ('PLJB', 'PLJBOH', 'CURB', 'CONPLJB', 'TRFLX', 'X', 'XB', 'BZ', 'BZXR', 'RAXIS'))

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'entry'):
            cls.entry.close()
            cls.tmp.cleanup()

    def test_bz_would_have_been_a_materially_different_normalization(self):
        """On this case TRANSP's BZ differs from the exported |b0| by more than 10%."""
        c = self.cdf
        r0 = 0.5 * (c['RAXIS'][0] + c['RAXIS'][-1]) * 0.01                  # m
        b0 = np.abs(c['BZXR']) * 0.01 / r0                                  # T
        self.assertGreater(float(np.max(np.abs(c['BZ'] / b0 - 1.0))), 0.1)

    def test_core_profiles_parallel_currents_use_the_ids_b0(self):
        ids = self.entry.get('core_profiles')
        b0 = np.abs(b0_of(ids))[:, None]
        for name, source in (('j_total', 'PLJB'), ('j_ohmic', 'PLJBOH')):
            values = np.array([np.asarray(getattr(p, name)) for p in ids.profiles_1d])
            assert_same_up_to_constant_sign(values * b0, self.cdf[source] * 1.0e4, 1e-9, name)   # A T m^-2
        nonind = np.array([np.asarray(p.j_non_inductive) for p in ids.profiles_1d])
        total = np.array([np.asarray(p.j_total) for p in ids.profiles_1d])
        ohmic = np.array([np.asarray(p.j_ohmic) for p in ids.profiles_1d])
        np.testing.assert_allclose(nonind, total - ohmic, rtol=1e-9, atol=1e-9 * np.abs(total).max())

    def test_core_profiles_rho_tor_uses_the_ids_b0(self):
        ids = self.entry.get('core_profiles')
        b0 = np.abs(b0_of(ids))
        rho = np.array([np.asarray(p.grid.rho_tor) for p in ids.profiles_1d])
        flux = np.array([np.interp(self.cdf['X'][i], np.r_[0.0, self.cdf['XB'][i]],
                                   np.r_[0.0, self.cdf['TRFLX'][i]]) for i in range(rho.shape[0])])
        np.testing.assert_allclose(np.pi * b0[:, None] * rho ** 2, np.abs(flux), rtol=2e-3,
                                   atol=2e-3 * np.abs(flux).max())

    def test_core_sources_ohmic_and_nbi_parallel_currents_use_the_ids_b0(self):
        ids = self.entry.get('core_sources')
        b0 = np.abs(b0_of(ids))[:, None]
        by_name = {str(s.identifier.name): s for s in ids.source}
        ohmic = np.array([np.asarray(p.j_parallel) for p in by_name['ohmic'].profiles_1d])
        assert_same_up_to_constant_sign(ohmic * b0, self.cdf['PLJBOH'] * 1.0e4, 1e-9, 'ohmic')
        nbi = np.array([np.asarray(p.j_parallel) for p in by_name['nbi'].profiles_1d])
        assert_same_up_to_constant_sign(nbi * b0, self.cdf['CURB'] * 1.0e4 * self.cdf['CONPLJB'], 1e-9, 'nbi')

    def test_equilibrium_parallel_current_and_rho_tor_use_the_ids_b0(self):
        ids = self.entry.get('equilibrium')
        b0 = np.abs(b0_of(ids))
        rho = np.array([np.asarray(ts.profiles_1d.rho_tor) for ts in ids.time_slice])
        np.testing.assert_allclose(np.pi * b0[:, None] * rho ** 2, np.abs(self.cdf['TRFLX']), rtol=1e-9,
                                   atol=1e-9 * np.abs(self.cdf['TRFLX']).max())
        # j_parallel is built on X and mapped to XB (linear, anchored at the axis value).
        jpar = np.array([np.asarray(ts.profiles_1d.j_parallel) for ts in ids.time_slice])
        jb_si = self.cdf['PLJB'] * 1.0e4
        expected = np.array([np.interp(self.cdf['XB'][i], np.r_[0.0, self.cdf['X'][i]],
                                       np.r_[jb_si[i, 0], jb_si[i]]) for i in range(jb_si.shape[0])])
        assert_same_up_to_constant_sign(jpar * b0[:, None], expected, 1e-6, 'equilibrium j_parallel')

    def test_r0_times_b0_is_one_convention_across_ids_and_the_tf_ids(self):
        tf = self.entry.get('tf')
        tf_t, tf_br = np.asarray(tf.b_field_phi_vacuum_r.time), np.asarray(tf.b_field_phi_vacuum_r.data)
        radii = []
        for name in ('core_profiles', 'core_sources', 'equilibrium'):
            ids = self.entry.get(name)
            r0 = float(ids.vacuum_toroidal_field.r0)
            radii.append(r0)
            np.testing.assert_allclose(r0 * b0_of(ids), np.interp(np.asarray(ids.time), tf_t, tf_br), rtol=1e-9,
                                       err_msg=name)
        summary = self.entry.get('summary')
        gq = summary.global_quantities
        r0 = float(np.asarray(gq.r0.value).ravel()[0])
        radii.append(r0)
        b0 = np.asarray(gq.b0.value, dtype=float)
        np.testing.assert_allclose(r0 * b0, np.interp(np.asarray(summary.time), tf_t, tf_br), rtol=1e-9,
                                   err_msg='summary')
        self.assertEqual(len({round(r, 12) for r in radii}), 1, f'reference radius differs between IDSs: {radii}')


@unittest.skipIf(imas is None, 'IMAS-Python is required')
class Dd3Tests(unittest.TestCase):
    def test_core_sources_and_core_profiles_in_dd3(self):
        if RUN['quick'] or not CASE.is_file():
            self.skipTest('end-to-end readback skipped (--quick or case not present)')
        cdf = read(CASE, ('PLJB', 'PLJBOH'))
        with tempfile.TemporaryDirectory() as tmp:
            entry = convert_and_open(CASE, ['core_profiles', 'core_sources'], '3.39.0', tmp)
            try:
                cp, cs = entry.get('core_profiles'), entry.get('core_sources')
            finally:
                entry.close()
        b0 = np.abs(b0_of(cp))[:, None]
        total = np.array([np.asarray(p.j_total) for p in cp.profiles_1d])
        assert_same_up_to_constant_sign(total * b0, cdf['PLJB'] * 1.0e4, 1e-9, 'DD3 j_total')
        ohmic = [s for s in cs.source if str(s.identifier.name) == 'ohmic'][0]
        values = np.array([np.asarray(p.j_parallel) for p in ohmic.profiles_1d])
        assert_same_up_to_constant_sign(values * np.abs(b0_of(cs))[:, None], cdf['PLJBOH'] * 1.0e4, 1e-9, 'DD3 ohmic')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='skip the HDF5 end-to-end readback')
    args, rest = parser.parse_known_args(argv)
    RUN['quick'] = args.quick
    result = unittest.main(argv=[sys.argv[0]] + rest, exit=False, verbosity=1).result
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
