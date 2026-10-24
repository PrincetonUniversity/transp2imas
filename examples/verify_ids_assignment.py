#!/usr/bin/python
"""Verify the validity checks, array arithmetic and IDS-node assignment primitives (no IMAS needed).

The nodes are small stand-ins that behave like IDS array nodes (``resize`` and slice assignment); the
ledger is the real ``MappingLedger``.  What is pinned:

* ``is_valid``: None, 0-d, too few rows, empty second axis, all-falsy sequences are invalid,
* ``sum_arrays`` / ``operate_arrays`` / ``np_multiply_sum``: results, division by zero giving 0, and bad
  input giving an empty array (logged) instead of an exception,
* ``fill_data``: absent data is a recorded skip, a size mismatch is a broken mapping for every time
  slice (``MappingError`` in strict mode, a recorded skip with the reason in lenient mode) and leaves the
  node untouched, ``ind``/``ind2``/``mult``/``interp``/``resize`` behave as documented,
* ``fill_value``: scalar and array forms.

    python examples/verify_ids_assignment.py
"""

import logging
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import mapping_status
from ids_assignment import (check_1dim, fill_data, fill_value, is_valid, node_path, np_multiply_sum,
                            operate_arrays, sum_arrays)
from mapping_status import MappingError, MappingLedger, SKIPPED


class ArrayNode:
    """Behaves like an IDS float-array node: resizable, assignable by slice."""

    class metadata:
        path = 'profiles_1d/electrons/temperature'

    def __init__(self, n=0):
        self.values = np.zeros(n)

    def resize(self, n, refcheck=False):
        self.values = np.zeros(n)

    def __setitem__(self, key, value):
        self.values[key] = value


class ValueNode:
    def __init__(self):
        self.value = ArrayNode()


class Recorder(logging.Handler):
    def __init__(self):
        super().__init__()
        self.messages = []

    def emit(self, record):
        self.messages.append(record.getMessage())


class IsValidTests(unittest.TestCase):
    def test_invalid_inputs(self):
        for bad in (None, np.array(5.0), np.array([]), np.zeros((0, 3)), np.zeros((3, 0)), [], [0, 0], (), set()):
            self.assertFalse(is_valid(bad), repr(bad))
        self.assertFalse(is_valid(np.ones(3), None))              # one bad argument is enough
        self.assertFalse(is_valid(np.ones(3), ind=3))             # too few rows for the index

    def test_valid_inputs(self):
        for good in (np.ones(3), np.ones((2, 4)), [0, 1], (2.0,), {1}):
            self.assertTrue(is_valid(good), repr(good))
        self.assertTrue(is_valid(np.ones(3), [1], ind=0))
        self.assertTrue(is_valid())

    def test_the_reason_is_logged_only_for_index_one(self):
        recorder, logger = Recorder(), logging.getLogger('assign-test-1')
        logger.addHandler(recorder)
        logger.setLevel(logging.INFO)
        is_valid(None, name='A', ind=0, logger=logger)
        is_valid(None, name='B', ind=1, logger=logger)
        is_valid(np.zeros((3, 0)), name='C', ind=1, logger=logger)
        self.assertEqual(recorder.messages, ['B is None', 'C shape (3, 0)'])

    def test_check_1dim(self):
        self.assertTrue(check_1dim(np.ones(3)) and check_1dim([1, 2]))
        self.assertFalse(check_1dim(np.ones((2, 2))) or check_1dim([[1], [2]]))
        self.assertIsNone(check_1dim(5))


class ArithmeticTests(unittest.TestCase):
    def test_sum_arrays_skips_empty_inputs(self):
        np.testing.assert_array_equal(sum_arrays(np.array([1.0, 2]), np.array([]), np.array([10.0, 20])), [11, 22])
        self.assertEqual(len(sum_arrays(np.array([]), [])), 0)

    def test_sum_arrays_reports_a_shape_error_instead_of_raising(self):
        recorder, logger = Recorder(), logging.getLogger('assign-test-2')
        logger.addHandler(recorder)
        logger.setLevel(logging.INFO)
        sum_arrays(np.ones(3), np.ones(4), logger=logger)
        self.assertTrue(any('sum arrays error' in m for m in recorder.messages))

    def test_operate_arrays(self):
        a, b = np.array([4.0, 9.0, 1.0]), np.array([2.0, 0.0, 4.0])
        np.testing.assert_array_equal(operate_arrays(a, '+', b), [6, 9, 5])
        np.testing.assert_array_equal(operate_arrays(a, '-', b), [2, 9, -3])
        np.testing.assert_array_equal(operate_arrays(a, '*', b), [8, 0, 4])
        np.testing.assert_array_equal(operate_arrays(a, '/', b), [2, 0, 0.25])          # x/0 -> 0
        np.testing.assert_array_equal(operate_arrays(a, 'sqrt'), [2, 3, 1])
        self.assertEqual(operate_arrays([4.0, 9.0], 'sqrt'), [2.0, 3.0])
        self.assertEqual(operate_arrays((4.0,), 'sqrt'), (2.0,))

    def test_operate_arrays_bad_input_gives_an_empty_array(self):
        self.assertEqual(operate_arrays(None, '+', np.ones(2)).size, 0)
        self.assertEqual(operate_arrays(np.ones(2), '+', np.array([])).size, 0)
        self.assertEqual(operate_arrays(np.ones(3), '+', np.ones(4)).size, 0)             # shape error
        self.assertEqual(operate_arrays(np.ones(2), '%', np.ones(2)).size, 0)             # unsupported
        self.assertEqual(operate_arrays(np.array([-1.0]), 'sqrt').size, 0)                # invalid value is an error here

    def test_np_multiply_sum(self):
        a, b = np.arange(6.0).reshape(2, 3), np.ones((2, 3))
        np.testing.assert_array_equal(np_multiply_sum(a, b, axis=1), [3, 12])
        np.testing.assert_array_equal(np_multiply_sum(a, b, axis=1, scale=2.0), [6, 24])
        self.assertEqual(np_multiply_sum(a, np.array([]), axis=1).size, 0)
        self.assertEqual(np_multiply_sum(a, np.ones((3, 2)), axis=1).size, 0)


class FillDataTests(unittest.TestCase):
    def setUp(self):
        self.ledger = MappingLedger()
        self.ledger.writer = 'core_profiles'
        self.logger = logging.getLogger('assign-test-3')
        self.times = lambda: (np.array([0.0, 1.0, 2.0, 3.0]), np.array([0.5, 1.5]))

    def tearDown(self):
        mapping_status.set_strict(True)

    def fill(self, node, dat, **kw):
        return fill_data(node, dat, self.ledger, self.logger, times=self.times, **kw)

    def test_a_valid_array_resizes_and_fills(self):
        node = ArrayNode(2)
        self.assertTrue(self.fill(node, np.arange(5.0)))
        np.testing.assert_array_equal(node.values, np.arange(5.0))
        self.assertEqual(self.ledger.entries(), [])

    def test_time_slice_range_scale_and_no_resize(self):
        data = np.arange(12.0).reshape(3, 4)
        node = ArrayNode()
        self.fill(node, data, ind=1)                                  # one time slice of 2-D data
        np.testing.assert_array_equal(node.values, [4, 5, 6, 7])
        self.fill(node, np.arange(10.0), ind=2, ind2=5, mult=2.0)     # a range of 1-D data, scaled
        np.testing.assert_array_equal(node.values, [4, 6, 8])
        fixed = ArrayNode(4)
        self.fill(fixed, np.arange(4.0), resize=False)
        np.testing.assert_array_equal(fixed.values, [0, 1, 2, 3])
        with self.assertRaises(MappingError):                         # no resize: a different length cannot be assigned
            self.fill(ArrayNode(3), np.arange(4.0), resize=False)

    def test_interpolation_to_the_profile_time_base(self):
        node = ArrayNode()
        self.fill(node, np.array([0.0, 10.0, 20.0, 30.0]), interp=True)
        np.testing.assert_allclose(node.values, [5.0, 15.0])
        same = lambda: (np.array([0.0, 1.0]), np.array([0.0, 1.0]))
        node = ArrayNode()
        fill_data(node, np.array([3.0, 4.0]), self.ledger, self.logger, interp=True, times=same)
        np.testing.assert_array_equal(node.values, [3.0, 4.0])

    def test_absent_data_is_a_recorded_skip(self):
        node = ArrayNode(3)
        for absent in (None, np.array([]), np.array(2.0)):
            self.assertFalse(self.fill(node, absent))
        rows = self.ledger.entries('core_profiles')
        self.assertEqual({(r['kind'], r['reason'], r['target']) for r in rows},
                         {(SKIPPED, 'no source data', 'profiles_1d/electrons/temperature')})
        self.assertEqual(node.values.size, 3)

    def test_a_size_mismatch_is_a_broken_mapping_for_every_time_slice(self):
        for ind in (0, 1, 2):
            node = ArrayNode(2)
            with self.subTest(ind=ind), self.assertRaises(MappingError) as caught:
                self.fill(node, np.ones((3, 5)), ind=ind, chk=4, name='TE')
            self.assertIn('size mismatch 5 != 4', str(caught.exception))
            self.assertIn('TE', str(caught.exception))
            self.assertEqual(node.values.size, 2)                                      # untouched

    def test_lenient_policy_records_the_reason_and_leaves_the_node_alone(self):
        mapping_status.set_strict(False)
        node = ArrayNode(2)
        self.assertFalse(self.fill(node, np.arange(5.0), chk=4, name='TE'))
        (row,) = self.ledger.entries('core_profiles')
        self.assertEqual((row['kind'], row['source']), (SKIPPED, 'TE'))
        self.assertIn('size mismatch 5 != 4', row['reason'])
        self.assertEqual(node.values.size, 2)

    def test_a_scalar_for_an_array_node_and_a_schema_error_are_broken_mappings(self):
        with self.assertRaises(MappingError) as caught:
            self.fill(ArrayNode(), [5.0], ind=0)                      # element 0 of the data is a scalar
        self.assertIn('cannot fill a scalar value', str(caught.exception))
        with self.assertRaises(MappingError):
            self.fill(object(), np.arange(3.0))                       # not a node: AttributeError on resize
        with self.assertRaises(MappingError):
            self.fill(ArrayNode(), np.ones((3, 2)))                   # 2-D data into a 1-D node

    def test_node_path_never_raises(self):
        self.assertEqual(node_path(ArrayNode()), 'profiles_1d/electrons/temperature')
        self.assertEqual(node_path(object()), 'object')


class FillValueTests(unittest.TestCase):
    def test_scalar_and_array_forms(self):
        node = ValueNode()
        fill_value(node, 3.5)                                         # a scalar goes straight to .value
        self.assertEqual(node.value, 3.5)
        node = ValueNode()
        fill_value(node, np.arange(4.0), mult=2.0)                    # an array resizes .value to its length
        np.testing.assert_array_equal(node.value.values, [0, 2, 4, 6])
        node = ValueNode()
        fill_value(node, 7.0, rsz=0, mult=2.0)                        # rsz=0: never resize (scalar summary values)
        self.assertEqual(node.value, 14.0)

    def test_invalid_data_is_ignored(self):
        node = ValueNode()
        fill_value(node, np.array([]))
        fill_value(node, None)
        self.assertEqual(node.value.values.size, 0)


if __name__ == '__main__':
    unittest.main(verbosity=1)
