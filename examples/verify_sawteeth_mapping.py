#!/usr/bin/python
"""Validate sawteeth crash histories on every bundled CDF and mock IDS writes.

The all-case regression distinguishes observed TLASTSAW changes from an
initial/constant 'last sawtooth' timestamp; it does not fabricate crash types
or inversion/mixing radii from the ad-hoc XSAWBOUND control.
"""

import argparse
import csv
import json
from pathlib import Path
import sys
import types
import unittest
from unittest import mock

import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sawteeth_mapping import (IMAS_MISSING_FLOAT, build_sawteeth_plan,
                              populate_sawteeth_ids)


class FakeVariable:
    def __init__(self, values, dim):
        self.values = np.asarray(values)
        self.dimensions = (dim,)

    def __getitem__(self, key):
        return self.values[key]


class FakeCDF:
    def __init__(self, time, trace, time3=None, xiq1=None, boundary=None):
        self.variables = {'TIME': FakeVariable(time, 'TIME')}
        if trace is not None:
            self.variables['TLASTSAW'] = FakeVariable(trace, 'TIME')
        if time3 is not None:
            self.variables['TIME3'] = FakeVariable(time3, 'TIME3')
        if xiq1 is not None:
            self.variables['XIQ1'] = FakeVariable(xiq1, 'TIME')
        if boundary is not None:
            self.variables['XSAWBOUND'] = FakeVariable(boundary, 'TIME')


def mock_ids():
    return types.SimpleNamespace(
        ids_properties=types.SimpleNamespace(homogeneous_time=None, comment='TRANSP'),
        diagnostics=types.SimpleNamespace(), time=None)


class TestSawteethMapping(unittest.TestCase):
    def test_event_detection_and_censored_initial_crash(self):
        nc = FakeCDF([.10,.12,.14,.16,.18,.20],
                     [.09,.09,.135,.135,.175,.175],
                     xiq1=[-.1,.1,.2,-.2,.3,.4],
                     boundary=[0,0,.25,.26,0,0])
        plan = build_sawteeth_plan(nc)
        self.assertEqual(plan.events, 2)
        np.testing.assert_array_equal(plan.event_indices, [2,4])
        np.testing.assert_allclose(plan.event_times, [.135,.175])
        self.assertEqual(plan.valid_q1_samples, 4)
        self.assertEqual(plan.source_boundary_samples, 2)
        self.assertTrue(np.all(plan.previous_period[:4] == IMAS_MISSING_FLOAT))
        np.testing.assert_allclose(plan.previous_period[4:], [.04,.04])
        ids = mock_ids()
        self.assertEqual(populate_sawteeth_ids(ids, plan), 2)
        self.assertEqual(ids.ids_properties.homogeneous_time, 1)
        np.testing.assert_allclose(ids.time, nc.variables['TIME'].values)
        np.testing.assert_allclose(ids.diagnostics.previous_crash_time, nc.variables['TLASTSAW'].values)
        self.assertEqual(ids.diagnostics.rho_tor_norm_q1[0], IMAS_MISSING_FLOAT)
        self.assertEqual(ids.diagnostics.rho_tor_norm_q1[3], IMAS_MISSING_FLOAT)
        self.assertFalse(hasattr(ids, 'crash_trigger'))
        self.assertFalse(hasattr(ids.diagnostics, 'rho_tor_norm_mixing'))
        self.assertFalse(hasattr(ids.diagnostics, 'rho_tor_norm_inversion'))

    def test_single_new_crash_no_invented_period(self):
        plan = build_sawteeth_plan(FakeCDF([0,1,2], [0,0,1]))
        self.assertEqual(plan.events, 1)
        self.assertIsNone(plan.previous_period)
        self.assertIsNone(plan.rho_tor_norm_q1)
        ids = mock_ids()
        populate_sawteeth_ids(ids, plan)
        self.assertFalse(hasattr(ids.diagnostics, 'previous_period'))

    def test_constant_or_no_trace_is_not_sawtooth_event(self):
        self.assertIsNone(build_sawteeth_plan(FakeCDF([0,1,2], [0,0,0])))
        self.assertIsNone(build_sawteeth_plan(FakeCDF([0,1,2], None,
                                                        xiq1=[.1,.2,.3])))
        self.assertIsNone(build_sawteeth_plan(FakeCDF([550,551], [550,550])))

    def test_time_and_time3_do_not_get_swapped(self):
        plan = build_sawteeth_plan(FakeCDF([1,2,3,4], [.5,.5,2.5,2.5],
                                               time3=[1,1.5,2,2.5,3,3.5,4],
                                               xiq1=[.1,.2,.3,.4]))
        self.assertEqual(plan.time.size, 4)
        self.assertEqual(plan.events, 1)

    def test_invalid_trace_shape_and_native_dimension(self):
        nc = FakeCDF([0,1,2], [0,0,1])
        nc.variables['TLASTSAW'].dimensions = ('TIME3',)
        with self.assertRaisesRegex(ValueError, 'native TIME'):
            build_sawteeth_plan(nc)
        nc.variables['TLASTSAW'].dimensions = ('TIME',)
        nc.variables['TLASTSAW'].values = np.array([0,1])
        with self.assertRaisesRegex(ValueError, 'expected 3'):
            build_sawteeth_plan(nc)

    def test_reset_and_future_crash_times_rejected(self):
        with self.assertRaisesRegex(ValueError, 'decreases'):
            build_sawteeth_plan(FakeCDF([1,2,3], [.5,1.5,.5]))
        with self.assertRaisesRegex(ValueError, 'future'):
            build_sawteeth_plan(FakeCDF([1,2,3], [.5,2.5,2.6]))

    def test_q1_invalid_samples_not_mapped_as_physical_zero(self):
        plan = build_sawteeth_plan(FakeCDF([1,2,3], [.5,1.5,2.5],
                                               xiq1=[-.1,.1,1.1]))
        np.testing.assert_allclose(plan.rho_tor_norm_q1,
                                   [IMAS_MISSING_FLOAT,.1,IMAS_MISSING_FLOAT])

    def test_existing_real_cmod_and_iter(self):
        cmod = ROOT / 'examples/data/cmod/89846.CDF'
        with netcdf_file(str(cmod), 'r', mmap=False) as nc:
            plan = build_sawteeth_plan(nc)
            self.assertEqual(plan.events, 38)
            self.assertEqual(plan.time.size, 397)
            self.assertEqual(plan.valid_q1_samples, 184)
            self.assertEqual(plan.source_boundary_samples, 0)
            np.testing.assert_allclose(plan.event_times,
                                       nc.variables['TLASTSAW'][:][plan.event_indices], rtol=0, atol=0)
        with netcdf_file(str(ROOT / 'examples/data/iter/38530P94.CDF'), 'r', mmap=False) as nc:
            self.assertIsNone(build_sawteeth_plan(nc))

    def test_real_writer_does_not_save_spurious_ids(self):
        # Import under mock external modules, as all other writer tests do.
        with mock.patch.dict(sys.modules, {
                'imas': types.ModuleType('imas'),
                'netCDF4': types.ModuleType('netCDF4')}):
            sys.modules['netCDF4'].Dataset = object
            from imas_sawteeth import imas_sawteeth
        im = mock.Mock()
        im.logger = mock.Mock()
        im.fill_properties.side_effect = lambda ids: ids
        im.fill_code_properties.side_effect = lambda ids: ids
        im.new_ids.side_effect = lambda name: mock_ids()
        im.nc = FakeCDF([.1,.2,.3,.4], [.05,.05,.24,.24], xiq1=[.1,-.1,.2,.3])
        # Real IMAS-Python strings are IDSString0D wrappers, not str.
        wrapped = mock_ids()
        wrapped.ids_properties.comment = types.SimpleNamespace(value='TRANSP')
        im.new_ids.return_value = wrapped
        im.new_ids.side_effect = None
        imas_sawteeth(im)
        im.putdb_imas.assert_called_once()
        ids = im.putdb_imas.call_args[0][0]
        self.assertIn('new crash updates=1', ids.ids_properties.comment)
        self.assertEqual(len(ids.time), 4)
        im.putdb_imas.reset_mock()
        im.nc = FakeCDF([.1,.2,.3], [.05,.05,.05])
        imas_sawteeth(im)
        im.putdb_imas.assert_not_called()


def verify_case(cdf):
    with netcdf_file(str(cdf), 'r', mmap=False) as nc:
        plan = build_sawteeth_plan(nc)
        if plan is None:
            return {'events': 0, 'time_samples': int(nc.variables['TIME'].shape[0]),
                    'q1_valid_samples': 0, 'event_times': [], 'observed_at': []}
        ids = mock_ids()
        assert populate_sawteeth_ids(ids, plan) == plan.events
        assert len(ids.time) == len(ids.diagnostics.previous_crash_time)
        assert plan.events == len(np.unique(plan.event_times))
        return {'events': plan.events, 'time_samples': int(plan.time.size),
                'q1_valid_samples': plan.valid_q1_samples,
                'event_times': [float(v) for v in plan.event_times],
                'observed_at': [float(v) for v in plan.time[plan.event_indices]]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true',
                        help='run synthetic tests (and known C-Mod fixture) only')
    parser.add_argument('--data-root', type=Path, default=ROOT / 'examples/data')
    parser.add_argument('--report-dir', type=Path, default=ROOT / 'build/sawteeth')
    args = parser.parse_args(argv)
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestSawteethMapping)
    if not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful():
        return 1
    if args.quick:
        print('SAWTEETH: PASS — source/IDS contract tests')
        return 0
    files = sorted(args.data_root.rglob('*.CDF'))
    if not files:
        parser.error(f'no CDF fixtures under {args.data_root}')
    args.report_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in files:
        case = verify_case(path)
        case['cdf'] = str(path.relative_to(args.data_root))
        print(f"{case['cdf']}: PASS — observed new crashes={case['events']}, "
              f"q=1 samples={case['q1_valid_samples']}")
        rows.append(case)
    (args.report_dir / 'sawteeth.json').write_text(json.dumps(rows, indent=2) + '\n')
    with (args.report_dir / 'sawteeth.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['cdf', 'observation_time', 'reported_crash_time'])
        for row in rows:
            for observed, crashed in zip(row['observed_at'], row['event_times']):
                writer.writerow([row['cdf'], observed, crashed])
    print(f'SAWTEETH: PASS — all {len(rows)} available CDFs; '
          f"{sum(row['events'] for row in rows)} observed new crash updates")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
