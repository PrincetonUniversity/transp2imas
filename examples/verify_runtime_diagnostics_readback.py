#!/usr/bin/python
# registry ids verified here (mapping_registry.yaml): pt_runtime_diagnostics_xml
"""Read back the PT_SOLVER runtime diagnostics stored in ``core_transport`` code parameters.

The diagnostic channels (``NPTRITER``, ``PT_CPTIME``, ``PT_WTIME``) live on the native CDF ``TIME``
coordinate, which is not the ``TIME3`` grid of the surrounding IDS, so each channel must carry its
own coordinate values.  Discovery must also work without PT_SOLVER input provenance.  This
converts real CDFs and reads the stored XML back from HDF5:

* ITER with its inputs: three channels, each with the CDF ``TIME`` values and hour units,
* the ITER CDF alone in an empty directory (no TR.DAT or model files): the solver output
  is still found, as the review requires,
* DIII-D 147634T13 (inactive solver, ``TIME`` has 4 samples and ``TIME3`` 7): no channel is invented.

    python examples/verify_runtime_diagnostics_readback.py [--quick]
"""

import argparse
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from pt_solver_diagnostics import KNOWN_CHANNELS, channel_coordinates, decode_channel

DATA = ROOT / 'examples' / 'data'
ITER = DATA / 'iter' / '38530P94.CDF'
DIIID = DATA / 'd3d' / '147634T13' / '147634T13.CDF'

try:
    import imas
except ImportError:      # pragma: no cover - environment dependent
    imas = None

RUN = {'quick': False}


def stored_diagnostics(cdf, dd_version='4.1.1'):
    """Convert ``cdf`` (core_transport only) and return the diagnostics element read from HDF5."""
    with tempfile.TemporaryDirectory() as tmp:
        result = subprocess.run(
            [sys.executable, str(ROOT / 'convert.py'), '-i', str(cdf), '-b', 'HDF5', '-u', tmp,
             '--dd-version', dd_version, '-new', '--ids', 'core_transport'],
            cwd=ROOT, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(result.stdout[-1500:] + result.stderr[-1500:])
        match = re.search(r'ids stored in (\S+)', result.stdout)
        path = Path(match.group(1))
        database, _major, shot, run = path.parts[-4:]
        entry = imas.DBEntry(imas.ids_defs.HDF5_BACKEND, database, int(shot), int(run), str(path.parents[3]),
                             dd_version=dd_version)
        entry.open()
        try:
            ids = entry.get('core_transport', occurrence=0)
            parameters = str(ids.code.parameters)
            n_time3 = len(ids.time)
        finally:
            entry.close()
    root = ET.fromstring(parameters)
    return root.find('pt_solver_runtime_diagnostics'), n_time3


def check_channels(test, diagnostics, cdf):
    with netcdf_file(cdf, mmap=False) as nc:
        time = np.asarray(nc.variables['TIME'][:], dtype='<f8')
        for channel in diagnostics.findall('channel'):
            name = channel.get('name')
            test.assertEqual(channel.get('quantity'), KNOWN_CHANNELS[name][0])
            np.testing.assert_array_equal(decode_channel(channel), np.asarray(nc.variables[name][:], dtype='<f8'))
            np.testing.assert_array_equal(channel_coordinates(diagnostics, channel)['TIME'], time)


@unittest.skipIf(imas is None, 'IMAS-Python is required')
class ReadbackTests(unittest.TestCase):
    def setUp(self):
        if RUN['quick']:
            self.skipTest('end-to-end readback skipped (--quick)')

    def test_iter_channels_keep_native_time_and_hour_units(self):
        if not ITER.is_file():
            self.skipTest('ITER example not present')
        diagnostics, _n = stored_diagnostics(ITER)
        self.assertEqual(sorted(c.get('name') for c in diagnostics.findall('channel')),
                         ['NPTRITER', 'PT_CPTIME', 'PT_WTIME'])
        check_channels(self, diagnostics, ITER)
        self.assertEqual(diagnostics.find("channel[@name='PT_WTIME']").get('units'), 'HOURS')
        self.assertEqual(diagnostics.find("channel[@name='PT_WTIME']").get('canonical_unit'), 'hour')

    def test_solver_output_is_found_without_any_input_provenance(self):
        if not ITER.is_file():
            self.skipTest('ITER example not present')
        with tempfile.TemporaryDirectory() as lone:
            copy = Path(lone) / ITER.name
            shutil.copy(ITER, copy)                       # CDF only: no TR.DAT, no model namelists
            self.assertEqual(sorted(p.name for p in Path(lone).iterdir()), [ITER.name])
            diagnostics, _n = stored_diagnostics(copy)
            self.assertIsNotNone(diagnostics, 'diagnostics dropped because no PT plan was available')
            self.assertEqual(sorted(c.get('name') for c in diagnostics.findall('channel')),
                             ['NPTRITER', 'PT_CPTIME', 'PT_WTIME'])
            check_channels(self, diagnostics, copy)

    def test_inactive_solver_records_nothing_when_time_and_time3_differ(self):
        if not DIIID.is_file():
            self.skipTest('DIII-D example not present')
        with netcdf_file(DIIID, mmap=False) as nc:
            self.assertNotEqual(nc.variables['TIME'].shape, nc.variables['TIME3'].shape)
        diagnostics, n_time3 = stored_diagnostics(DIIID)
        self.assertEqual(n_time3, 7)
        self.assertTrue(diagnostics is None or diagnostics.findall('channel') == [])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='skip the HDF5 end-to-end readback')
    args, rest = parser.parse_known_args(argv)
    RUN['quick'] = args.quick
    result = unittest.main(argv=[sys.argv[0]] + rest, exit=False, verbosity=1).result
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
