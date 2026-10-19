import unittest
import xml.etree.ElementTree as ET
import numpy as np
from pt_solver_diagnostics import append_runtime_diagnostics_xml, decode_channel

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

if __name__ == '__main__': unittest.main()
