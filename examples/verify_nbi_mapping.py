#!/usr/bin/python
# registry ids verified here (mapping_registry.yaml): nbi_native_time
"""Verify the nbi writer's time axis and source identities with real ``nbi`` IDS objects.

Cases the review of the t2i2 snapshot asked for, each built from a small fake CDF:

* a beam with only a full-energy component keeps the native ``TIME`` (it used to be zeroed, which wrote an
  all-zero time axis for both ITER beams and all 14 DIII-D 163303U09 beams);
* pulsed and off intervals are expressed by the power values, not by the clock;
* non-consecutive source numbers keep their own numbers and their own channels;
* no species-total channel, several species, and a species-only file never index past the list or guess;
* a trace whose length differs from ``TIME`` is a broken mapping, not silently written.

    python examples/verify_nbi_mapping.py
"""

import re
from pathlib import Path
import sys
import types
import unittest
import logging

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

try:
    import imas
except ImportError:      # pragma: no cover - environment dependent
    imas = None

if imas is not None:
    import ids_assignment
    import ids_identity
    from imas_nbi import imas_nbi
    from mapping_status import MappingError, MappingLedger

DD = '4.1.1'
TIME = np.array([0.1, 0.2, 0.3, 0.4])


class FakeConverter:
    """The part of the converter the nbi writer uses, over a dict of CDF arrays and a real IDS factory."""

    def __init__(self, variables):
        self.variables = {k: np.asarray(v, dtype=float) for k, v in variables.items()}
        self.logger = logging.getLogger('verify-nbi-mapping')
        self.mapping = MappingLedger()
        self.mapping.writer = 'nbi'
        self.dd_version = DD
        self.factory = imas.IDSFactory(DD)

    # CDF
    def get_ncvar(self, name, scale=1.0, **kwargs):
        return self.variables[name] * scale if name in self.variables else np.array([])

    def get_var_names_re(self, expr, nd1=0, nd2=0, length=False, **kwargs):
        regex = re.escape(expr[:-1]) + rf'\d{{{nd1},{nd2}}}$'
        return sorted(n for n in self.variables if re.fullmatch(regex, n))

    def get_var_names(self, expr, cut='', **kwargs):
        regex = re.escape(expr[:-1]) + r'[A-Z0-9]$'
        return sorted(n[len(cut):] for n in self.variables if re.fullmatch(regex, n))

    # DD and filling
    def dd_at_least(self, major, minor=0, patch=0):
        return (4, 1, 1) >= (major, minor, patch)

    def new_ids(self, name):
        return self.factory.new(name)

    def is_valid(self, *values, ind=0, name=''):
        return ids_assignment.is_valid(*values, ind=ind, name=name, logger=self.logger)

    def checkErrExit(self, *values, name=''):
        pass

    def set_nbi_unit_identity(self, unit, identifier, description=None):
        ids_identity.set_nbi_unit_identity(self.dd_at_least, unit, identifier, description)

    def set_nuclear_charge(self, node, value):
        ids_identity.set_nuclear_charge(self.dd_at_least, node, value)

    def set_species_name(self, node, value):
        ids_identity.set_species_name(self.dd_at_least, node, value)

    def fill_properties(self, ids):
        return ids

    def fill_code_properties(self, ids):
        return ids

    def fill_external_model_inputs(self, ids, family):
        return []


def build(variables):
    converter = FakeConverter({'TIME': TIME, **variables})
    ids = imas_nbi(converter).build()
    return converter, (ids[0] if ids else None)


def beam(number, power, e1=1.0e5, e2=0.0, e3=0.0):
    suffix = f'{number:02d}'
    return {f'PINJ{suffix}': power, f'EINJ{suffix}_E1': [e1] * len(TIME),
            f'EINJ{suffix}_E2': [e2] * len(TIME), f'EINJ{suffix}_E3': [e3] * len(TIME)}


@unittest.skipIf(imas is None, 'IMAS-Python is required (real IDS objects)')
class NbiMappingTests(unittest.TestCase):
    def test_a_full_energy_only_beam_keeps_the_native_time(self):
        _converter, nbi = build({**beam(1, [0.0, 1.0e6, 1.0e6, 1.0e6]), 'PINJ_D': [0.0, 1.0e6, 1.0e6, 1.0e6]})
        unit = nbi.unit[0]
        np.testing.assert_array_equal(np.asarray(unit.energy.time), TIME)
        np.testing.assert_array_equal(np.asarray(unit.power_launched.time), TIME)
        np.testing.assert_array_equal(np.asarray(nbi.time), TIME)

    def test_off_and_pulsed_intervals_are_in_the_power_not_the_clock(self):
        power = [0.0, 2.0e6, 0.0, 2.0e6]
        _converter, nbi = build({**beam(1, power, e2=3.0e4, e3=2.0e4), 'PINJ_D': power})
        unit = nbi.unit[0]
        np.testing.assert_array_equal(np.asarray(unit.power_launched.data), power)
        np.testing.assert_array_equal(np.asarray(unit.power_launched.time), TIME)
        self.assertTrue(np.all(np.diff(np.asarray(unit.power_launched.time)) > 0))

    def test_a_single_time_sample(self):
        converter = FakeConverter({'TIME': [0.5], 'PINJ01': [1.0e6], 'EINJ01_E1': [8.0e4]})
        nbi = imas_nbi(converter).build()[0]
        np.testing.assert_array_equal(np.asarray(nbi.unit[0].power_launched.time), [0.5])

    def test_non_consecutive_sources_keep_their_numbers_and_channels(self):
        _converter, nbi = build({**beam(3, [1.0, 1.0, 1.0, 1.0], e1=7.0e4), **beam(7, [5.0, 5.0, 5.0, 5.0], e1=9.0e4),
                                 'PINJ_D': [6.0] * 4})
        self.assertEqual([str(u.name) for u in nbi.unit], ['unit_3', 'unit_7'])
        self.assertEqual(np.asarray(nbi.unit[0].power_launched.data).tolist(), [1.0] * 4)
        self.assertEqual(np.asarray(nbi.unit[1].power_launched.data).tolist(), [5.0] * 4)
        self.assertEqual(np.asarray(nbi.unit[1].energy.data).tolist(), [9.0e4] * 4)

    def test_sources_are_ordered_by_number_not_by_name(self):
        _converter, nbi = build({**beam(10, [1.0] * 4), **beam(2, [2.0] * 4)})
        self.assertEqual([str(u.name) for u in nbi.unit], ['unit_2', 'unit_10'])

    def test_no_species_total_channel_leaves_the_species_unset_without_an_index_error(self):
        converter, nbi = build({**beam(1, [1.0] * 4)})
        self.assertEqual(len(nbi.unit), 1)
        self.assertFalse(nbi.unit[0].species.a.has_value)

    def test_one_species_total_applies_to_every_unit(self):
        _converter, nbi = build({**beam(1, [1.0] * 4), **beam(2, [1.0] * 4), 'PINJ_D': [2.0] * 4})
        for unit in nbi.unit:
            self.assertAlmostEqual(float(unit.species.a), 2.0, delta=0.1)

    def test_several_species_are_not_matched_by_list_position(self):
        converter, nbi = build({**beam(1, [1.0] * 4), **beam(2, [1.0] * 4), 'PINJ_D': [1.0] * 4, 'PINJ_T': [1.0] * 4})
        for unit in nbi.unit:
            self.assertFalse(unit.species.a.has_value)
        reasons = [row['reason'] for row in converter.mapping.entries('nbi')]
        self.assertTrue(any('cannot be matched' in reason for reason in reasons), reasons)

    def test_a_species_only_file_builds_one_unit_per_species_with_its_own_channel(self):
        _converter, nbi = build({'PINJ_D': [1.0, 2.0, 3.0, 4.0], 'EINJ': [8.0e4] * 4})
        self.assertEqual(len(nbi.unit), 1)
        self.assertEqual(np.asarray(nbi.unit[0].power_launched.data).tolist(), [1.0, 2.0, 3.0, 4.0])
        self.assertAlmostEqual(float(nbi.unit[0].species.a), 2.0, delta=0.1)

    def test_a_trace_of_another_length_is_a_broken_mapping(self):
        variables = {**beam(1, [1.0, 1.0, 1.0])}                    # 3 samples on a 4-sample TIME
        with self.assertRaises(MappingError):
            build(variables)

    def test_no_nbi_channels_writes_nothing(self):
        self.assertEqual(imas_nbi(FakeConverter({'TIME': TIME})).build(), [])


if __name__ == '__main__':
    unittest.main(verbosity=1)
