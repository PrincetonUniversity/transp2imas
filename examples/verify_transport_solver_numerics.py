#!/usr/bin/python
"""Verify the ``transport_solver_numerics`` IDS built from PT_SOLVER runtime diagnostics.

The slices are on the native CDF ``TIME`` coordinate and carry only quantities that exist in the
CDF: ``NPTRITER`` as integer ``newton_iterations`` and the ``PT_*`` timers (hours) as real
``cumulative_*_time_hours``.  Guards: timers whose units are not hours, non-integer counts and
all-zero channels are not written, and nothing is written where the DD (3.39.0) has unnamed
control parameters.

1. unit tests with a stand-in CDF and real IDS objects (DD 3.39.0, 3.42.2, 4.1.1),
2. (not with ``--quick``) the ITER CDF converted and read back from HDF5 in several DD versions.

    python examples/verify_transport_solver_numerics.py [--quick]
"""

import argparse
import logging
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

from imas_transport_solver_numerics import imas_transport_solver_numerics
from mapping_status import MappingLedger

try:
    import imas
except ImportError:      # pragma: no cover - environment dependent
    imas = None

ITER = ROOT / 'examples' / 'data' / 'iter' / '38530P94.CDF'
RUN = {'quick': False}


class Var:
    def __init__(self, data, dim='TIME', units=''):
        self.data = np.asarray(data, dtype=float)
        self.shape = self.data.shape
        self.dtype = self.data.dtype
        self.dimensions = (dim,)
        self.units = units

    def __getitem__(self, key):
        return self.data[key]


class FakeNC:
    """A CDF whose native TIME (4 samples) differs in length from TIME3 (7)."""

    def __init__(self, newton=(0, 4, 2, 3), cpu=(0.0, 0.1, 0.2, 0.3), wall=(0.0, 1.0, 2.0, 3.0), wall_units=b'HOURS'):
        self.variables = {
            'TIME': Var([10.0, 11.0, 12.5, 13.0], 'TIME', b'SEC'),
            'TIME3': Var(np.linspace(10.0, 13.0, 7), 'TIME3', b'SEC'),
            'NPTRITER': Var(newton),
            'PT_CPTIME': Var(cpu, units=b'HOURS'),
            'PT_WTIME': Var(wall, units=wall_units),
        }


class Shim:
    """Just enough of ``nc2imas`` for this writer, with a real IDS factory."""

    def __init__(self, nc, dd_version):
        self.nc = nc
        self.dd_version = dd_version
        self.logger = logging.getLogger('verify_transport_solver_numerics')
        self.mapping = MappingLedger()
        self.written = None

    def new_ids(self, name):
        return imas.IDSFactory(self.dd_version).new(name)

    def fill_properties(self, ids):
        ids.ids_properties.comment = 'synthetic'
        return ids

    def fill_code_properties(self, ids):
        return ids

    def vacuum_field_reference(self, times):
        return 1.0, np.full(len(times), 2.0)

    def putdb_imas(self, ids):
        self.written = ids


LAST = {}


def write(nc, dd_version='4.1.1'):
    shim = Shim(nc, dd_version)
    imas_transport_solver_numerics(shim)
    LAST['shim'] = shim
    return shim.written


@unittest.skipIf(imas is None, 'IMAS-Python is required')
class UnitTests(unittest.TestCase):
    def test_slices_are_on_native_time_with_typed_named_parameters(self):
        for dd in ('3.42.2', '4.1.1'):
            with self.subTest(dd=dd):
                ids = write(FakeNC(), dd)
                ids.validate()
                np.testing.assert_array_equal(np.asarray(ids.time), [10.0, 11.0, 12.5, 13.0])   # not TIME3 (7)
                self.assertEqual(len(ids.solver_1d), 4)
                self.assertEqual((str(ids.solver.name), int(ids.solver.index)), ('PT_SOLVER', -1))
                it = 2
                parameters = ids.solver_1d[it].control_parameters
                self.assertEqual([(str(s.name), int(s.value)) for s in parameters.integer0d], [('newton_iterations', 2)])
                self.assertEqual([(str(s.name), float(s.value)) for s in parameters.real0d],
                                 [('cumulative_cpu_time_hours', 0.2), ('cumulative_wall_time_hours', 2.0)])
                # no per-equation convergence is invented from a global count
                self.assertEqual(len(ids.solver_1d[it].equation), 0)
                np.testing.assert_allclose(float(ids.vacuum_toroidal_field.r0) * np.asarray(ids.vacuum_toroidal_field.b0), 2.0)

    def test_dd_without_named_control_parameters_writes_nothing_and_records_why(self):
        self.assertIsNone(write(FakeNC(), '3.39.0'))
        (entry,) = LAST['shim'].mapping.entries()
        self.assertEqual((entry['kind'], entry['target']), ('skipped', 'transport_solver_numerics'))
        self.assertIn('unnamed control parameters', entry['reason'])

    def test_inactive_solver_writes_nothing_and_records_why(self):
        self.assertIsNone(write(FakeNC(newton=(0, 0, 0, 0), cpu=(0, 0, 0, 0), wall=(0, 0, 0, 0))))
        (entry,) = LAST['shim'].mapping.entries()
        self.assertIn('no active PT_SOLVER runtime channels', entry['reason'])

    def test_timer_not_in_hours_is_skipped_not_mislabelled(self):
        ids = write(FakeNC(wall_units=b'SEC'))
        names = [str(s.name) for s in ids.solver_1d[1].control_parameters.real0d]
        self.assertEqual(names, ['cumulative_cpu_time_hours'])

    def test_non_integer_iteration_count_is_skipped(self):
        ids = write(FakeNC(newton=(0, 4.5, 2, 3)))
        self.assertEqual(len(ids.solver_1d[1].control_parameters.integer0d), 0)
        self.assertEqual(len(ids.solver_1d[1].control_parameters.real0d), 2)

    def test_channel_that_is_inactive_is_dropped_but_the_others_remain(self):
        ids = write(FakeNC(cpu=(0, 0, 0, 0)))                    # CPU timer all zero: unused
        self.assertEqual([str(s.name) for s in ids.solver_1d[1].control_parameters.real0d],
                         ['cumulative_wall_time_hours'])


def read(path, names):
    with netcdf_file(path, mmap=False) as nc:
        return {n: np.array(nc.variables[n][:], dtype=float) for n in names}


def convert_and_open(cdf, dd_version, outdir):
    result = subprocess.run(
        [sys.executable, str(ROOT / 'convert.py'), '-i', str(cdf), '-b', 'HDF5', '-u', str(outdir),
         '--dd-version', dd_version, '-new', '--ids', 'transport_solver_numerics'],
        cwd=ROOT, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stdout[-1500:] + result.stderr[-1500:])
    path = Path(re.search(r'ids stored in (\S+)', result.stdout).group(1))
    database, _major, shot, run = path.parts[-4:]
    entry = imas.DBEntry(imas.ids_defs.HDF5_BACKEND, database, int(shot), int(run), str(path.parents[3]),
                         dd_version=dd_version)
    entry.open()
    return entry, path


@unittest.skipIf(imas is None, 'IMAS-Python is required')
class ConvertedEntryTests(unittest.TestCase):
    def setUp(self):
        if RUN['quick'] or not ITER.is_file():
            self.skipTest('end-to-end readback skipped (--quick or ITER example not present)')

    def test_iter_values_match_the_cdf_in_each_dd(self):
        cdf = read(ITER, ('TIME', 'NPTRITER', 'PT_CPTIME', 'PT_WTIME'))
        for dd in ('3.42.2', '4.0.0', '4.1.1'):
            with self.subTest(dd=dd), tempfile.TemporaryDirectory() as tmp:
                entry, _path = convert_and_open(ITER, dd, tmp)
                try:
                    ids = entry.get('transport_solver_numerics')
                finally:
                    entry.close()
                ids.validate()
                np.testing.assert_allclose(np.asarray(ids.time), cdf['TIME'], rtol=1e-12)
                self.assertEqual(len(ids.solver_1d), cdf['TIME'].size)
                newton = [int(s.value) for sl in ids.solver_1d for s in sl.control_parameters.integer0d]
                cpu = [float(s.value) for sl in ids.solver_1d for s in sl.control_parameters.real0d
                       if str(s.name) == 'cumulative_cpu_time_hours']
                wall = [float(s.value) for sl in ids.solver_1d for s in sl.control_parameters.real0d
                        if str(s.name) == 'cumulative_wall_time_hours']
                np.testing.assert_array_equal(newton, cdf['NPTRITER'].astype(int))
                np.testing.assert_allclose(cpu, cdf['PT_CPTIME'], rtol=1e-12)
                np.testing.assert_allclose(wall, cdf['PT_WTIME'], rtol=1e-12)

    def test_dd_3_39_writes_no_entry(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = subprocess.run(
                [sys.executable, str(ROOT / 'convert.py'), '-i', str(ITER), '-b', 'HDF5', '-u', tmp,
                 '--dd-version', '3.39.0', '-new', '--ids', 'transport_solver_numerics'],
                cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout[-800:] + result.stderr[-800:])
            self.assertEqual(list(Path(tmp).rglob('transport_solver_numerics.h5')), [])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='skip the HDF5 end-to-end readback')
    args, rest = parser.parse_known_args(argv)
    RUN['quick'] = args.quick
    result = unittest.main(argv=[sys.argv[0]] + rest, exit=False, verbosity=1).result
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
