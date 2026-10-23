#!/usr/bin/python
import sys
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pt_solver_diagnostics import (append_runtime_diagnostics_xml, channel_coordinates, decode_channel,
                                   diagnostic_xml)

class Var:
    def __init__(self, data, units='1'):
        self.data = np.asarray(data)
        self.shape = self.data.shape
        self.dtype = self.data.dtype
        self.dimensions = ('TIME3', 'XB')[:self.data.ndim]
        self.units = units
    def __getitem__(self, k): return self.data[k]

class NC:
    def __init__(self):
        self.variables = {
            'NPTRITER': Var([[1,2,3], [4,5,6]]),
            'PT_CPTIME': Var([[0.1, np.nan, 0.3], [0,1,2]]),
            'TE': Var([[1,2,3]]),
            'PT_WTIME': Var([[np.nan]]),
            'PT_PECLET_GUESS': Var([[1,2]]),
        }

class TestDiagnostics(unittest.TestCase):
    def test_lossless_and_selective(self):
        root = ET.fromstring(append_runtime_diagnostics_xml('<base/>', NC()))
        diag = root.find('pt_solver_runtime_diagnostics')
        self.assertEqual(diag.get('count'), '2')
        self.assertEqual([x.get('name') for x in diag.findall('channel')],
                         ['NPTRITER','PT_CPTIME'])
        np.testing.assert_array_equal(decode_channel(diag.findall('channel')[0]), [[1,2,3],[4,5,6]])
        self.assertTrue(np.isnan(decode_channel(diag.findall('channel')[1])[0,1]))
        self.assertEqual(len(diag.findall('skipped')),1)
        self.assertEqual(len(diag.findall('unavailable')),3)
    def test_empty_inventory(self):
        diag = ET.fromstring(append_runtime_diagnostics_xml('<base/>', object())).find('pt_solver_runtime_diagnostics')
        self.assertEqual(diag.get('status'), 'CDF_variable_inventory_unavailable')
    def test_idempotent(self):
        xml = append_runtime_diagnostics_xml('<base/>', NC())
        xml = append_runtime_diagnostics_xml(xml, NC())
        self.assertEqual(len(ET.fromstring(xml).findall('pt_solver_runtime_diagnostics')), 1)


class Var1D:
    """1-D CDF variable on a named dimension (a TIME-type coordinate or a channel on it)."""

    def __init__(self, data, dim='TIME', units='1'):
        self.data = np.asarray(data)
        self.shape = self.data.shape
        self.dtype = self.data.dtype
        self.dimensions = (dim,)
        self.units = units

    def __getitem__(self, k):
        return self.data[k]


class TimeNC:
    """A CDF whose native TIME (4 samples) differs in length from TIME3 (7 samples)."""

    def __init__(self, newton=(3, 4, 5, 6), cpu=(0.1, 0.2, 0.3, 0.4)):
        self.variables = {
            'TIME': Var1D([0.1, 0.3, 0.5, 0.9], 'TIME', b'SEC     '),
            'TIME3': Var1D(np.linspace(0.1, 0.9, 7), 'TIME3', b'SEC     '),
            'NPTRITER': Var1D(newton, 'TIME'),
            'PT_CPTIME': Var1D(cpu, 'TIME', b'HOURS   '),
        }


class TestNativeCoordinates(unittest.TestCase):
    def diag(self, nc, **kw):
        return diagnostic_xml(nc, **kw)

    def test_channels_carry_their_own_time_coordinate_not_time3(self):
        diag = self.diag(TimeNC())
        for channel in diag.findall('channel'):
            self.assertEqual(channel.get('dimensions'), 'TIME')
            self.assertEqual(channel.get('coordinates'), 'TIME')
            coordinates = channel_coordinates(diag, channel)
            np.testing.assert_array_equal(coordinates['TIME'], [0.1, 0.3, 0.5, 0.9])
            self.assertEqual(coordinates['TIME'].size, decode_channel(channel).size)
        self.assertEqual([c.get('name') for c in diag.findall('coordinate')], ['TIME'])   # TIME3 is not attached
        self.assertEqual(diag.find('coordinate').get('units'), 'SEC')

    def test_coordinate_checksum_is_verified(self):
        diag = self.diag(TimeNC())
        coordinate = diag.find('coordinate')
        coordinate.set('sha256', '0' * 64)
        with self.assertRaises(ValueError):
            channel_coordinates(diag, diag.find('channel'))

    def test_coordinate_of_the_wrong_length_is_not_attached(self):
        nc = TimeNC()
        nc.variables['TIME'] = Var1D([0.1, 0.3, 0.5], 'TIME')          # 3 samples for 4-sample channels
        diag = self.diag(nc)
        self.assertEqual({c.get('coordinates') for c in diag.findall('channel')}, {''})
        self.assertEqual(diag.findall('coordinate'), [])

    def test_text_attributes_are_decoded_not_bytes_reprs(self):
        channel = {c.get('name'): c for c in self.diag(TimeNC()).findall('channel')}['PT_CPTIME']
        self.assertEqual(channel.get('units'), 'HOURS')
        self.assertEqual(channel.get('canonical_unit'), 'hour')

    def test_discovery_does_not_need_a_pt_plan(self):
        nc = TimeNC()
        xml = append_runtime_diagnostics_xml('<base/>', nc, record_inactive=False)
        diag = ET.fromstring(xml).find('pt_solver_runtime_diagnostics')
        self.assertEqual(diag.get('count'), '2')

    def test_inactive_channels_are_kept_only_when_a_plan_exists(self):
        zeros = TimeNC(newton=(0, 0, 0, 0), cpu=(0, 0, 0, 0))
        self.assertEqual(append_runtime_diagnostics_xml('<base/>', zeros, record_inactive=False), '<base/>')
        self.assertEqual(append_runtime_diagnostics_xml('', zeros, record_inactive=False), '')
        reason = [x.get('reason') for x in diagnostic_xml(zeros, record_inactive=False).findall('skipped')]
        self.assertTrue(all('zero' in r for r in reason) and len(reason) == 2)
        recorded = diagnostic_xml(zeros, record_inactive=True)
        self.assertEqual(recorded.get('count'), '2')           # plan exists: raw zeros are kept

    def test_mixed_zero_and_nonzero_counts_as_active(self):
        diag = diagnostic_xml(TimeNC(newton=(0, 0, 2, 0)), record_inactive=False)
        self.assertIn('NPTRITER', [c.get('name') for c in diag.findall('channel')])


class TestRealCdfCoordinates(unittest.TestCase):
    DIIID = ROOT / 'examples' / 'data' / 'd3d' / '147634T13' / '147634T13.CDF'
    ITER = ROOT / 'examples' / 'data' / 'iter' / '38530P94.CDF'

    def test_diii_d_time_and_time3_differ_and_the_channel_keeps_time(self):
        if not self.DIIID.is_file():
            self.skipTest('DIII-D example not present')
        from scipy.io import netcdf_file
        with netcdf_file(self.DIIID, mmap=False) as nc:
            self.assertEqual((nc.variables['TIME'].shape, nc.variables['TIME3'].shape), ((4,), (7,)))
            diag = diagnostic_xml(nc, record_inactive=True)
            channel = diag.find("channel[@name='NPTRITER']")
            time = channel_coordinates(diag, channel)['TIME']
            np.testing.assert_array_equal(time, np.asarray(nc.variables['TIME'][:], dtype='<f8'))
            self.assertEqual(time.size, 4)
            self.assertEqual(decode_channel(channel).shape, (4,))
            # not an active solver run: nothing is recorded without a PT plan
            self.assertEqual(diagnostic_xml(nc, record_inactive=False).findall('channel'), [])

    def test_iter_active_channels_are_found_without_a_plan_and_keep_time_and_units(self):
        if not self.ITER.is_file():
            self.skipTest('ITER example not present')
        from scipy.io import netcdf_file
        with netcdf_file(self.ITER, mmap=False) as nc:
            diag = diagnostic_xml(nc, record_inactive=False)
            self.assertEqual(sorted(c.get('name') for c in diag.findall('channel')),
                             ['NPTRITER', 'PT_CPTIME', 'PT_WTIME'])
            for channel in diag.findall('channel'):
                np.testing.assert_array_equal(channel_coordinates(diag, channel)['TIME'],
                                              np.asarray(nc.variables['TIME'][:], dtype='<f8'))
            self.assertEqual(diag.find("channel[@name='PT_WTIME']").get('units'), 'HOURS')


if __name__ == '__main__': unittest.main()
