#!/usr/bin/python
"""Verify the required-versus-optional mapping policy.

A CDF may legitimately lack a variable, and the converter then skips the matching IMAS fields.  A
variable that *is* in the CDF but cannot be mapped (wrong shape, unreadable, schema mismatch) is a
different case and must not look like absent optional data.  Policy:

* absent optional inputs and skipped fields are recorded, with reasons, in the ``MappingLedger`` and
  reported per writer in the conversion journal,
* a broken mapping raises ``MappingError`` (strict, the default); with ``--lenient-mapping`` the
  failure is recorded as a skipped field with its reason instead, so nothing is hidden either way.

1. unit tests of the ledger, the journal and the policy switch (no IMAS),
2. real conversions of bundled CDFs through the actual converter, checking the journal.

    python examples/verify_mapping_status.py [--quick]
"""

import argparse
import json
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

from conversion_plan import ConversionJournal
import mapping_status
from mapping_status import (ADJUSTED, KINDS, MISSING_OPTIONAL, SKIPPED, UNAVAILABLE_IN_DD, MappingError, MappingLedger,
                            RecordingDataset)

try:
    import imas
    import netCDF4
    from nc2imas import nc2imas
    from netcdf_class import netcdf_connection
except ImportError:      # pragma: no cover - environment dependent
    imas = netCDF4 = nc2imas = netcdf_connection = None

RUN = {'quick': False}
CASE = ROOT / 'examples' / 'data' / 'nstx' / '141711Y05' / '141711Y05.CDF'


class LedgerTests(unittest.TestCase):
    def test_entries_are_per_writer_and_deduplicated_with_counts(self):
        ledger = MappingLedger()
        ledger.writer = 'core_profiles'
        for _ in range(3):
            ledger.missing_optional('NE_ERR')
        ledger.skipped('no source data', target='core_profiles/profiles_1d/zeff')
        ledger.writer = 'summary'
        ledger.missing_optional('NE_ERR')
        mine = ledger.entries('core_profiles')
        self.assertEqual([(e['kind'], e['source'], e['count']) for e in mine if e['kind'] == MISSING_OPTIONAL],
                         [(MISSING_OPTIONAL, 'NE_ERR', 3)])
        self.assertEqual(len(ledger.entries('summary')), 1)
        self.assertEqual(len(ledger.entries()), 3)

    def test_unknown_kind_is_rejected(self):
        with self.assertRaises(ValueError):
            MappingLedger().record('maybe', source='X')
        self.assertEqual(KINDS, (MISSING_OPTIONAL, SKIPPED, UNAVAILABLE_IN_DD, ADJUSTED))

    def test_summary_is_compact_json_with_complete_counts(self):
        ledger = MappingLedger()
        ledger.writer = 'w'
        for i in range(250):
            ledger.missing_optional(f'VAR{i:03d}')
        ledger.skipped('shape mismatch', source='A', target='x/y')
        ledger.unavailable_in_dd('workflow', 'IDS not in DD 3.39.0')
        summary = ledger.summary('w', max_listed=10)
        self.assertEqual(summary[MISSING_OPTIONAL]['count'], 250)
        self.assertEqual(len(summary[MISSING_OPTIONAL]['sources']), 10)
        self.assertTrue(summary[MISSING_OPTIONAL]['truncated'])
        self.assertEqual(summary[SKIPPED]['fields'], [{'source': 'A', 'target': 'x/y', 'reason': 'shape mismatch',
                                                       'count': 1}])
        self.assertEqual(summary[UNAVAILABLE_IN_DD]['fields'][0]['target'], 'workflow')
        json.dumps(summary)                                   # JSON-serializable

    def test_mapping_error_names_the_source_target_and_writer(self):
        error = MappingError('size mismatch 19 != 20', source='PBE', target='core_sources/source/profiles_1d',
                             writer='core_sources')
        for part in ('size mismatch 19 != 20', 'PBE', 'core_sources/source/profiles_1d', 'core_sources'):
            self.assertIn(part, str(error))
        self.assertIsInstance(error, RuntimeError)


class JournalTests(unittest.TestCase):
    def test_journal_carries_the_mapping_summary_of_each_writer(self):
        ledger = MappingLedger()
        ledger.writer = 'summary'
        ledger.missing_optional('Q95_ERR')
        journal = ConversionJournal('case.CDF', ('summary',))
        journal.record_writer('summary', [], 0.1, mapping=ledger.summary('summary'))
        record = journal.data['writers'][0]
        self.assertEqual(record['mapping'][MISSING_OPTIONAL]['sources'], ['Q95_ERR'])
        json.dumps(journal.data)

    def test_journal_without_mapping_information_is_unchanged(self):
        journal = ConversionJournal('case.CDF', ('workflow',))
        journal.record_writer('workflow', [], 0.1)
        self.assertNotIn('mapping', journal.data['writers'][0])


class ConsumptionTests(unittest.TestCase):
    class Dataset:
        variables = {'A': 1, 'B': 2, 'C': 3}
        dimensions = {'X': 3}

        def ncattrs(self):
            return ['shot']

    def test_only_item_access_and_get_count_as_reads_and_are_attributed_to_the_writer(self):
        ledger = MappingLedger()
        nc = RecordingDataset(self.Dataset(), ledger)
        _ = nc.variables['A']                       # before any writer: setup
        ledger.writer = 'summary'
        self.assertIn('B', nc.variables)            # probing is not a read
        self.assertEqual(sorted(nc.variables), ['A', 'B', 'C'])
        self.assertEqual(list(nc.variables.keys()), ['A', 'B', 'C'])
        self.assertEqual(nc.variables.get('B'), 2)
        self.assertIsNone(nc.variables.get('MISSING'))
        with self.assertRaises(KeyError):
            nc.variables['MISSING']
        self.assertEqual(ledger.consumed(None), ['A'])
        self.assertEqual(ledger.consumed('summary'), ['B'])
        self.assertEqual(nc.ncattrs(), ['shot'])    # everything else is delegated
        self.assertEqual(nc.dimensions, {'X': 3})

    def test_journal_lists_what_a_writer_read(self):
        journal = ConversionJournal('case.CDF', ('summary',))
        journal.record_writer('summary', [], 0.1, consumed=['TE', 'NE'])
        journal.record_writer('wall', [], 0.1)
        self.assertEqual(journal.data['writers'][0]['consumed'], ['NE', 'TE'])
        self.assertNotIn('consumed', journal.data['writers'][1])


class PolicyTests(unittest.TestCase):
    def tearDown(self):
        mapping_status.set_strict(True)

    def test_broken_raises_when_strict_and_records_a_reason_when_lenient(self):
        ledger = MappingLedger()
        ledger.writer = 'core_profiles'
        with self.assertRaises(MappingError) as caught:
            ledger.broken(ValueError('size mismatch 19 != 20'), source='PBE', target='source/profiles_1d')
        self.assertIn('size mismatch 19 != 20', str(caught.exception))
        self.assertIn('core_profiles', str(caught.exception))
        self.assertEqual(ledger.entries(), [])                    # strict records nothing, it fails
        mapping_status.set_strict(False)
        ledger.broken(ValueError('size mismatch 19 != 20'), source='PBE', target='source/profiles_1d')
        (entry,) = ledger.entries('core_profiles')
        self.assertEqual((entry['kind'], entry['source'], entry['target']), (SKIPPED, 'PBE', 'source/profiles_1d'))
        self.assertIn('ValueError: size mismatch 19 != 20', entry['reason'])


class FakeCdf:
    """The slice of a netCDF dataset that ``get_ncvar`` reads."""

    def __init__(self, **variables):
        self.variables = {k: np.asarray(v) for k, v in variables.items()}


@unittest.skipIf(netcdf_connection is None, 'netCDF4 and IMAS-Python are required')
class GetNcvarTests(unittest.TestCase):
    def setUp(self):
        self.conn = netcdf_connection.__new__(netcdf_connection)
        self.conn.mapping = MappingLedger()
        self.conn.mapping.writer = 'summary'
        self.conn.nc = FakeCdf(PRESENT=np.arange(4, dtype=np.float32), SHORT=np.arange(3.0))
        self.conn.timesc = np.array([0.0, 1.0, 2.0, 3.0])
        self.conn.timep = np.array([0.5, 1.5])

    def tearDown(self):
        mapping_status.set_strict(True)

    def test_an_absent_variable_is_optional_input_recorded_by_name(self):
        self.assertEqual(self.conn.get_ncvar('NOT_IN_CDF').size, 0)
        self.assertEqual([(e['kind'], e['source']) for e in self.conn.mapping.entries('summary')],
                         [(MISSING_OPTIONAL, 'NOT_IN_CDF')])

    def test_a_present_variable_is_returned_as_float64(self):
        values = self.conn.get_ncvar('PRESENT', scale=2.0)
        self.assertEqual(values.dtype, np.float64)
        np.testing.assert_array_equal(values, [0, 2, 4, 6])
        self.assertEqual(self.conn.mapping.entries(), [])

    def test_a_present_variable_that_cannot_be_transformed_is_a_broken_mapping(self):
        with self.assertRaises(MappingError) as caught:                     # 3 samples cannot be on the 4-sample TIME
            self.conn.get_ncvar('SHORT', interp=True)
        self.assertIn('SHORT', str(caught.exception))
        with self.assertRaises(MappingError):
            self.conn.get_ncvar('PRESENT', ind=10)                          # time index out of range

    def test_lenient_policy_records_the_reason_and_returns_empty(self):
        mapping_status.set_strict(False)
        self.assertEqual(self.conn.get_ncvar('SHORT', interp=True).size, 0)
        (entry,) = self.conn.mapping.entries('summary')
        self.assertEqual((entry['kind'], entry['source']), (SKIPPED, 'SHORT'))
        self.assertIn('ValueError', entry['reason'])


@unittest.skipIf(nc2imas is None, 'netCDF4 and IMAS-Python are required')
class FillDataTests(unittest.TestCase):
    def setUp(self):
        import logging
        self.im = nc2imas.__new__(nc2imas)
        self.im.mapping = MappingLedger()
        self.im.mapping.writer = 'core_sources'
        self.im.logger = logging.getLogger('verify_mapping_status')
        self.im.timesc = self.im.timep = np.array([0.0, 1.0])
        self.ids = imas.IDSFactory('4.1.1').new('core_sources')

    def tearDown(self):
        mapping_status.set_strict(True)

    def test_filling_valid_data_returns_true(self):
        self.assertTrue(self.im.fill_data(self.ids.time, np.arange(5.0)))
        np.testing.assert_array_equal(np.asarray(self.ids.time), np.arange(5.0))

    def test_absent_data_is_a_recorded_skip_not_an_error(self):
        self.assertFalse(self.im.fill_data(self.ids.time, np.array([])))
        (entry,) = self.im.mapping.entries('core_sources')
        self.assertEqual((entry['kind'], entry['reason'], entry['target']), (SKIPPED, 'no source data', 'time'))

    def test_a_size_mismatch_is_a_broken_mapping_for_every_time_slice(self):
        for ind in (0, 1, 2):                                   # it used to be reported only for slice 0
            with self.subTest(ind=ind), self.assertRaises(MappingError) as caught:
                self.im.fill_data(self.ids.time, np.ones((3, 5)), ind=ind, chk=4)
            self.assertIn('size mismatch 5 != 4', str(caught.exception))
            self.assertIn('core_sources', str(caught.exception))

    def test_lenient_size_mismatch_is_skipped_with_the_reason(self):
        mapping_status.set_strict(False)
        self.assertFalse(self.im.fill_data(self.ids.time, np.arange(5.0), chk=4))
        (entry,) = self.im.mapping.entries('core_sources')
        self.assertIn('size mismatch 5 != 4', entry['reason'])
        self.assertEqual(len(np.asarray(self.ids.time)), 0)    # nothing half-written


def corrupt_radial_length(src, dst, variable):
    """Copy a CDF, giving ``variable`` (TIME3 x X) a radial length one shorter than X."""
    with netCDF4.Dataset(src) as source, netCDF4.Dataset(dst, 'w', format=source.data_model) as target:
        target.setncatts({k: source.getncattr(k) for k in source.ncattrs()})
        for name, dim in source.dimensions.items():
            target.createDimension(name, None if dim.isunlimited() else len(dim))
        short = len(source.dimensions[source.variables[variable].dimensions[1]]) - 1
        target.createDimension('X_SHORT', short)
        for name, var in source.variables.items():
            dims = (var.dimensions[0], 'X_SHORT') if name == variable else var.dimensions
            fill = var.getncattr('_FillValue') if '_FillValue' in var.ncattrs() else None
            out = target.createVariable(name, var.dtype, dims, fill_value=fill)
            out.setncatts({k: var.getncattr(k) for k in var.ncattrs() if k != '_FillValue'})
            data = var[:]
            out[:] = data[:, :short] if name == variable else data


def convert(cdf, outdir, *extra, ids='core_sources', report=None):
    command = [sys.executable, str(ROOT / 'convert.py'), '-i', str(cdf), '-b', 'HDF5', '-u', str(outdir),
               '--dd-version', '4.1.1', '-new', '--ids', ids, *extra]
    if report:
        command += ['--report-json', str(report)]
    return subprocess.run(command, cwd=ROOT, capture_output=True, text=True)


@unittest.skipIf(imas is None or netCDF4 is None, 'netCDF4 and IMAS-Python are required')
class ConvertedJournalTests(unittest.TestCase):
    def setUp(self):
        if RUN['quick'] or not CASE.is_file():
            self.skipTest('real conversions skipped (--quick or case not present)')

    def test_a_normal_conversion_reports_absent_inputs_and_reasoned_skips(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / 'journal.json'
            result = convert(CASE, tmp, ids='core_profiles,lh_antennas,summary', report=report)
            self.assertEqual(result.returncode, 0, result.stdout[-600:] + result.stderr[-600:])
            journal = json.loads(report.read_text())
        self.assertEqual((journal['status'], journal['mapping_policy']), ('PASS', 'strict'))
        writers = {w['ids']: w for w in journal['writers']}
        absent = writers['core_profiles']['mapping'][MISSING_OPTIONAL]
        self.assertGreater(absent['count'], 0)
        self.assertIn('BDENS_3', absent['sources'])             # a species fast-ion density this run does not have
        self.assertEqual(writers['lh_antennas']['status'], 'SKIP')
        reasons = [f['reason'] for f in writers['lh_antennas']['mapping'][SKIPPED]['fields']]
        self.assertTrue(reasons and all(reasons))               # a skipped writer always carries a reason

    def test_a_present_variable_with_the_wrong_radial_length_fails_strict_and_is_recorded_lenient(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / CASE.name
            corrupt_radial_length(CASE, bad, 'PBE')
            for name in ('TR.DAT',):
                source = CASE.with_name(CASE.stem + name)
                if source.is_file():
                    (Path(tmp) / source.name).write_bytes(source.read_bytes())
            strict_report = Path(tmp) / 'strict.json'
            strict = convert(bad, Path(tmp) / 'db_strict', report=strict_report)
            self.assertNotEqual(strict.returncode, 0, 'a broken mapping must not be silently dropped')
            failed = json.loads(strict_report.read_text())
            self.assertEqual(failed['status'], 'FAIL')
            self.assertIn('size mismatch 19 != 20', failed['error'])
            self.assertIn('PBE', failed['error'])                    # the offending CDF variable is named
            self.assertIn('core_sources', failed['error'])
            lenient_report = Path(tmp) / 'lenient.json'
            lenient = convert(bad, Path(tmp) / 'db_lenient', '--lenient-mapping', report=lenient_report)
            self.assertEqual(lenient.returncode, 0, lenient.stdout[-600:] + lenient.stderr[-600:])
            journal = json.loads(lenient_report.read_text())
        self.assertEqual((journal['status'], journal['mapping_policy']), ('PASS', 'lenient'))
        skipped = journal['writers'][0]['mapping'][SKIPPED]['fields']
        mismatches = [f for f in skipped if 'size mismatch' in f['reason']]
        self.assertTrue(mismatches, skipped)
        self.assertIn('PBE', {f.get('source') for f in mismatches})


D3D_CASE = ROOT / 'examples' / 'data' / 'd3d' / '147634B13' / '147634B13.CDF'
# (writer, a profile variable it reads on TIME3 x X, a small bundled case that has it): one per reader that used
# to drop a radial-size mismatch silently
READER_VARIABLES = (('waves', 'PEECH1', D3D_CASE, 'PEECH1'), ('edge_sources', 'SEGF_D', D3D_CASE, 'SEGF_D'),
                    ('radiation', 'PRAD_BR', D3D_CASE, 'PRAD_BR'),
                    ('distribution_sources', 'BDEP01_TOT', CASE, 'BDEP01_TOT'),
                    # core_profiles and equilibrium read radial profiles through IdsWriter.profile, which names the CDF
                    # variable; core_transport and summary through the same check on their own grid
                    ('core_profiles', 'ZEFFP', D3D_CASE, 'ZEFFP'), ('core_profiles', 'UTHRM', D3D_CASE, 'UTHRM'),
                    ('core_profiles', 'CUROH', D3D_CASE, 'CUROH'), ('core_profiles', 'PTOWB', D3D_CASE, 'PTOWB'),
                    ('equilibrium', 'ELONG', D3D_CASE, 'ELONG'), ('equilibrium', 'GXI', D3D_CASE, 'GXI'),
                    ('equilibrium', 'TRIANG', D3D_CASE, 'TRIANG'), ('equilibrium', 'DRAV', D3D_CASE, 'DRAV'),
                    ('core_transport', 'CONDE', D3D_CASE, 'CONDE'), ('core_transport', 'AMTR_OBS', D3D_CASE, 'AMTR_OBS'),
                    ('summary', 'UTHRM', D3D_CASE, 'UTHRM'), ('summary', 'PTOWB', D3D_CASE, 'PTOWB'),
                    ('summary', 'CUROH', D3D_CASE, 'CUROH'))


@unittest.skipIf(imas is None or netCDF4 is None, 'netCDF4 and IMAS-Python are required')
class EveryReaderShapeTests(unittest.TestCase):
    """A present profile with the wrong radial size is a broken mapping in every writer's reader."""

    def setUp(self):
        if RUN['quick']:
            self.skipTest('real conversions skipped (--quick)')

    def conversions(self, ids, variable, case):
        if not case.is_file():
            self.skipTest(f'{case.name} not present')
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / case.name
            corrupt_radial_length(case, bad, variable)
            for source in case.parent.iterdir():
                if source.suffix != '.CDF' and source.is_file():
                    (Path(tmp) / source.name).write_bytes(source.read_bytes())
            strict_report, lenient_report = Path(tmp) / 'strict.json', Path(tmp) / 'lenient.json'
            strict = convert(bad, Path(tmp) / 'db_strict', ids=ids, report=strict_report)
            lenient = convert(bad, Path(tmp) / 'db_lenient', '--lenient-mapping', ids=ids, report=lenient_report)
            yield strict, json.loads(strict_report.read_text()), lenient, json.loads(lenient_report.read_text())

    def test_strict_fails_naming_the_variable_and_lenient_records_it(self):
        for ids, variable, case, named in READER_VARIABLES:
            with self.subTest(writer=ids, variable=variable):
                for strict, failed, lenient, journal in self.conversions(ids, variable, case):
                    self.assertNotEqual(strict.returncode, 0, f'{ids}: {variable} was silently dropped')
                    self.assertEqual(failed['status'], 'FAIL')
                    self.assertIn(named, failed['error'])
                    self.assertRegex(failed['error'], 'radial size mismatch|size mismatch')
                    self.assertEqual(lenient.returncode, 0, lenient.stdout[-500:] + lenient.stderr[-500:])
                    self.assertEqual(journal['mapping_policy'], 'lenient')
                    skipped = journal['writers'][0]['mapping'][SKIPPED]['fields']
                    hits = [f for f in skipped if f.get('source') == named and 'size mismatch' in f['reason']]
                    self.assertTrue(hits, skipped)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='skip the real conversions')
    args, rest = parser.parse_known_args(argv)
    RUN['quick'] = args.quick
    result = unittest.main(argv=[sys.argv[0]] + rest, exit=False, verbosity=1).result
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
