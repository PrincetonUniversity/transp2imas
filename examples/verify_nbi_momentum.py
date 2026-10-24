#!/usr/bin/python
# registry ids verified here (mapping_registry.yaml): nbi_momentum_source
"""Verify the NBI total toroidal momentum source written to ``core_sources``.

TRANSP assembles its total input torque ``TQIN`` (xplbal.f90, ``TQINP``) from the beam terms
``TQBE + TQBI + TQBTH + TQJXB`` plus ``TQDIFB`` (CDF name ``TQABDF``, fast-ion anomalous
diffusion), ripple, ``TQANOM`` and ``TQNTV``.  The IMAS ``momentum_*`` node is the total source
for the thermal plasma and the ``*_j_cross_b_field`` node is a subset already included in it,
so the NBI total must be the sum of those terms, not ``TQBTH`` alone.

1. source closure on every bundled CDF (SciPy only): ``TQIN`` equals the beam sum, and ITER
   needs ``TQABDF``; the per-beam totals ``TQTOT<nn>`` sum to the four-term total,
2. (not with ``--quick``) real conversions read back from HDF5: the nbi ``momentum`` profile is
   ``ip_sign * 1e6 * TQIN`` and the JxB subset is ``ip_sign * 1e6 * TQJXB``, in DD 4.1.1 and 3.39.0.

    python examples/verify_nbi_momentum.py [--quick]
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
DATA = ROOT / 'examples' / 'data'

try:
    import imas
except ImportError:      # pragma: no cover - environment dependent
    imas = None

RUN = {'quick': False}
BEAM_TERMS = ('TQBE', 'TQBI', 'TQBTH', 'TQJXB')


def beam_cdfs():
    """Bundled TRANSP result CDFs that have beam torque and a total input torque."""
    found = []
    for path in sorted(DATA.rglob('*.CDF')):
        with netcdf_file(path, mmap=False) as nc:
            if {'TQIN', 'TQBE', 'DVOL'} <= set(nc.variables):
                found.append(path)
    return found


def read(path, names):
    with netcdf_file(path, mmap=False) as nc:
        return {n: np.array(nc.variables[n][:], dtype=float) if n in nc.variables else None for n in names}


def rms_rel(residual, reference):
    scale = np.sqrt(np.mean(reference ** 2))
    return float(np.sqrt(np.mean(residual ** 2)) / scale) if scale > 0 else 0.0


class SourceClosureTests(unittest.TestCase):
    def test_tqin_is_the_beam_sum_plus_fast_ion_anomalous_diffusion(self):
        cases = beam_cdfs()
        if not cases:
            self.skipTest('no CDF with beam torques')
        for path in cases:
            with self.subTest(case=path.name):
                v = read(path, BEAM_TERMS + ('TQABDF', 'TQIN', 'DVOL'))
                dv = v['DVOL']
                beam = sum(v[n] for n in BEAM_TERMS)
                abdf = v['TQABDF'] if v['TQABDF'] is not None else 0.0
                self.assertLess(rms_rel((v['TQIN'] - beam - abdf) * dv, v['TQIN'] * dv), 1e-5)

    def test_tqabdf_is_needed_where_present(self):
        iter_cdf = DATA / 'iter' / '38530P94.CDF'
        if not iter_cdf.is_file():
            self.skipTest('ITER example not present')
        v = read(iter_cdf, BEAM_TERMS + ('TQABDF', 'TQIN', 'DVOL'))
        beam = sum(v[n] for n in BEAM_TERMS)
        self.assertGreater(rms_rel((v['TQIN'] - beam) * v['DVOL'], v['TQIN'] * v['DVOL']), 0.02)

    def test_per_beam_totals_close_in_the_volume_integral_except_iter(self):
        """Per-beam ``TQTOT<nn>`` sum to the four-term total when volume integrated (profiles differ
        by ~1e-3), but not for ITER (8.5%), which is why they are not used to build the total."""
        checked = 0
        for path in beam_cdfs():
            with netcdf_file(path, mmap=False) as nc:
                names = sorted(n for n in nc.variables if re.fullmatch(r'TQTOT\d\d', n))
            if not names or path.name.startswith('38530'):
                continue
            v = read(path, BEAM_TERMS + tuple(names) + ('DVOL',))
            beam = (sum(v[n] for n in BEAM_TERMS) * v['DVOL']).sum(axis=1)
            total = (sum(v[n] for n in names) * v['DVOL']).sum(axis=1)
            self.assertLess(rms_rel(total - beam, beam), 1e-6, path.name)
            checked += 1
        if not checked:
            self.skipTest('no CDF with per-beam TQTOT')


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


@unittest.skipIf(imas is None, 'IMAS-Python is required')
class ConvertedEntryTests(unittest.TestCase):
    CASES = (('nstx/141711Y01/141711Y01.CDF', '4.1.1'), ('nstx/141711Y01/141711Y01.CDF', '3.39.0'),
             ('iter/38530P94.CDF', '4.1.1'))

    def _check(self, relpath, dd):
        cdf = DATA / relpath
        if RUN['quick'] or not cdf.is_file():
            self.skipTest('end-to-end readback skipped (--quick or case not present)')
        with tempfile.TemporaryDirectory() as tmp:
            entry = convert_and_open(cdf, ['core_sources'], dd, tmp)
            try:
                ids = entry.get('core_sources')
            finally:
                entry.close()
        nbi = [s for s in ids.source if str(s.identifier.name) == 'nbi']
        self.assertEqual(len(nbi), 1, 'expected exactly one nbi source')
        v = read(cdf, ('TQIN', 'TQJXB', 'TIME3'))
        mom_name = 'momentum_tor' if dd.startswith('3.') else 'momentum_phi'
        jxb_name = mom_name + '_j_cross_b_field'
        profiles = nbi[0].profiles_1d
        self.assertEqual(len(profiles), v['TIME3'].size)
        mom = np.array([np.asarray(getattr(p, mom_name)) for p in profiles])
        jxb = np.array([np.asarray(getattr(p, jxb_name)) for p in profiles])
        self.assertEqual(mom.shape, v['TQIN'].shape)
        self.assertGreater(np.abs(mom).max(), 0.0)
        # TRANSP torque (N m cm^-3, + along Ip) -> N m^-2 with the IMAS toroidal sign: one constant +-1
        sign = np.sign(np.sum(mom * v['TQIN']))
        self.assertIn(sign, (-1.0, 1.0))
        np.testing.assert_allclose(mom, sign * 1.0e6 * v['TQIN'], rtol=1e-5, atol=1e-5 * np.abs(mom).max())
        np.testing.assert_allclose(jxb, sign * 1.0e6 * v['TQJXB'], rtol=1e-5, atol=1e-5 * np.abs(mom).max())

    def test_nstx_dd4(self):
        self._check(*self.CASES[0])

    def test_nstx_dd3(self):
        self._check(*self.CASES[1])

    def test_iter_needs_the_anomalous_diffusion_term(self):
        self._check(*self.CASES[2])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='skip the HDF5 end-to-end readback')
    args, rest = parser.parse_known_args(argv)
    RUN['quick'] = args.quick
    result = unittest.main(argv=[sys.argv[0]] + rest, exit=False, verbosity=1).result
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
