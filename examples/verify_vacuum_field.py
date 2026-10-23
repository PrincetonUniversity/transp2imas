#!/usr/bin/python
"""Verify that ``vacuum_toroidal_field`` is consistent with the TF IDS vacuum field.

IMAS requires ``vacuum_toroidal_field.r0 * b0(t)`` to equal the TF IDS vacuum field times
radius.  The converter used ``b0 = BZXR / RAXIS(t)`` with a constant ``r0``, which breaks
the product when the magnetic axis moves (0.15% for ITER 38530P94 up to 17% for NSTX
142202U01).  This checks the shared ``utilities.vacuum_field_reference`` and, not with
``--quick``, real conversions read back from HDF5.

    python examples/verify_vacuum_field.py [--quick]
"""

import argparse
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from utilities import vacuum_field_reference

try:
    import imas
except ImportError:      # pragma: no cover - environment dependent
    imas = None

CASE = ROOT / 'examples' / 'data' / 'nstx' / '142202U01' / '142202U01.CDF'
RUN = {'quick': False}


class ReferenceFunctionTests(unittest.TestCase):
    def setUp(self):
        self.t = np.linspace(0.0, 1.0, 11)
        self.raxis = 0.8 + 0.25 * np.sin(2.0 * self.t)             # moving axis: 0.8 ... 1.03 m
        self.bzxr = 0.5 * 0.9 * (1.0 + 0.05 * self.t)               # B R in T m

    def test_product_equals_field_times_radius_for_either_sign(self):
        for sign in (1.0, -1.0):
            r0, b0 = vacuum_field_reference(self.raxis, self.bzxr, sign, self.t, self.t)
            np.testing.assert_allclose(r0 * b0, sign * self.bzxr, rtol=1e-14)

    def test_reference_radius_is_one_constant_from_the_source_grid(self):
        r0, _b0 = vacuum_field_reference(self.raxis, self.bzxr, 1.0, self.t[3:6], self.t)
        self.assertAlmostEqual(r0, 0.5 * (self.raxis[0] + self.raxis[-1]), places=14)

    def test_interpolates_the_product_to_another_time_grid(self):
        target = np.linspace(0.05, 0.95, 7)
        r0, b0 = vacuum_field_reference(self.raxis, self.bzxr, -1.0, target, self.t)
        np.testing.assert_allclose(r0 * b0, -np.interp(target, self.t, self.bzxr), rtol=1e-14)

    def test_the_old_formula_violated_the_requirement(self):
        r0 = 0.5 * (self.raxis[0] + self.raxis[-1])
        old_b0 = self.bzxr / self.raxis
        self.assertGreater(np.max(np.abs(r0 * old_b0 / self.bzxr - 1.0)), 0.05)


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
    """Real HDF5 write/readback: every IDS that carries a vacuum field agrees with the TF IDS."""

    WRITERS = ('radiation', 'distribution_sources')

    @classmethod
    def setUpClass(cls):
        if RUN['quick'] or not CASE.is_file():
            raise unittest.SkipTest('end-to-end readback skipped (--quick or case not present)')
        cls.tmp = tempfile.TemporaryDirectory()
        cls.entry = convert_and_open(CASE, ('tf',) + cls.WRITERS, '4.1.1', cls.tmp.name)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'entry'):
            cls.entry.close()
            cls.tmp.cleanup()

    def test_r0_times_b0_equals_tf_vacuum_field_times_radius(self):
        tf = self.entry.get('tf')
        tf_t = np.asarray(tf.b_field_phi_vacuum_r.time)
        tf_br = np.asarray(tf.b_field_phi_vacuum_r.data)
        self.assertGreater(tf_t.size, 1)
        checked = 0
        for name in self.WRITERS:
            ids = self.entry.get(name)
            r0 = float(ids.vacuum_toroidal_field.r0)
            b0 = np.asarray(ids.vacuum_toroidal_field.b0)
            times = np.asarray(ids.time)
            self.assertEqual(b0.size, times.size, name)
            np.testing.assert_allclose(r0 * b0, np.interp(times, tf_t, tf_br), rtol=1e-9, err_msg=name)
            checked += 1
        self.assertEqual(checked, len(self.WRITERS))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='skip the HDF5 end-to-end readback')
    args, rest = parser.parse_known_args(argv)
    RUN['quick'] = args.quick
    result = unittest.main(argv=[sys.argv[0]] + rest, exit=False, verbosity=1).result
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
