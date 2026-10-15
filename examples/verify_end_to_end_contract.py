#!/usr/bin/python
"""Dependency-free tests for --ids, write journals, and IMAS readback contracts.

The real IMAS conversion is deliberately NOT mocked here.  It is executed by
verify_end_to_end.py on a workstation with an IMAS-Python installation.
"""

import ast
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from conversion_plan import ConversionJournal, WRITER_ORDER, select_ids
import verify_end_to_end as e2e


class TestEndToEndContract(unittest.TestCase):
    def test_every_writer_registered_and_selection_in_canonical_order(self):
        module = ast.parse((ROOT / 'convert.py').read_text())
        func = next(n for n in module.body if isinstance(n, ast.FunctionDef)
                    and n.name == 'make_writers')
        dictionary = next(n.value for n in ast.walk(func) if isinstance(n, ast.Return))
        self.assertIsInstance(dictionary, ast.Dict)
        self.assertEqual(tuple(node.value for node in dictionary.keys), WRITER_ORDER)
        self.assertEqual(select_ids(None), WRITER_ORDER)
        self.assertEqual(select_ids('ALL'), WRITER_ORDER)
        self.assertEqual(select_ids('ic_antenna,lh_antenna'),
                         ('lh_antennas', 'ic_antennas'))
        self.assertEqual(select_ids('equilibrium,workflow,equilibrium,core_profiles'),
                         ('workflow', 'core_profiles', 'equilibrium'))
        # Journal names must be exact IMAS DBEntry.get() names, not progress labels.
        self.assertIn("new_ids('lh_antennas')", (ROOT / 'imas_lh_antennas.py').read_text())
        self.assertIn("new_ids('ic_antennas')", (ROOT / 'imas_ic_antennas.py').read_text())
        with self.assertRaisesRegex(ValueError, 'Unknown --ids'):
            select_ids('pf_alternative')
        with self.assertRaisesRegex(ValueError, 'requires'):
            select_ids('equilibrium,')

    def test_convert_parser_exposes_report_and_ids(self):
        source = (ROOT / 'utilities.py').read_text()
        self.assertIn("parser.add_argument('--ids'", source)
        self.assertIn("parser.add_argument('--report-json'", source)
        source = (ROOT / 'convert.py').read_text()
        self.assertIn('ConversionJournal(', source)
        self.assertIn('nci._written_ids[before:]', source)

    def test_journal_does_not_claim_attempts_as_writes(self):
        with tempfile.TemporaryDirectory() as td:
            journal = ConversionJournal('case.CDF', select_ids('workflow,waves,equilibrium'))
            journal.record_writer('workflow', [{'ids': 'workflow', 'occurrence': 0}], 0.1)
            journal.record_writer('waves', [], 0.1)
            journal.record_writer('equilibrium', [], 0.1, ValueError('bad TIME3'))
            journal.fail(ValueError('bad TIME3'))
            path = Path(td) / 'report.json'
            journal.save(path)
            import json
            obj = json.loads(path.read_text())
            self.assertEqual(obj['status'], 'FAIL')
            self.assertEqual([r['status'] for r in obj['writers']],
                             ['WROTE', 'SKIP', 'FAIL'])
            self.assertEqual(obj['written_ids'], [{'ids':'workflow', 'occurrence':0}])

    def test_numpy_cdf_metadata_and_events_survive_journal_json(self):
        # Real CDF examples expose NumPy int32 shot/run metadata.  Without
        # normalization every converted IDS can succeed, yet the JSON journal
        # fails at finalization and the entire end-to-end suite reports FAIL.
        converter = types.SimpleNamespace(
            imas_dbname=np.str_('D3D'), shot=np.int32(147634),
            run=np.int64(2013), occur=np.int32(2), dd_version=np.str_('4.1.1'),
            outdir=Path('/tmp/imasdb/D3D/4/147634/2013'),
            backend=np.int32(9))
        journal = ConversionJournal('147634T13.CDF', ('equilibrium',))
        journal.attach_converter(converter)
        journal.record_writer('equilibrium', [
            {'ids': np.str_('equilibrium'), 'occurrence': np.int32(2)}],
            np.float32(0.4))
        journal.finish()
        with tempfile.TemporaryDirectory() as td:
            output = Path(td) / 'e2e' / 'case.json'
            journal.save(output)
            import json
            saved = json.loads(output.read_text())
        self.assertEqual(saved['status'], 'PASS')
        self.assertEqual(saved['machine'], 'D3D')
        self.assertEqual(saved['shot'], 147634)
        self.assertEqual(saved['run'], 2013)
        self.assertEqual(saved['occurrence'], 2)
        self.assertEqual(saved['dd_version'], '4.1.1')
        self.assertEqual(saved['backend'], '9')
        self.assertEqual(saved['output_directory'], str(converter.outdir))
        self.assertEqual(saved['written_ids'],
                         [{'ids': 'equilibrium', 'occurrence': 2}])
        self.assertEqual(saved['writers'][0]['occurrences'], [2])
        self.assertIs(type(saved['shot']), int)
        self.assertIs(type(saved['written_ids'][0]['occurrence']), int)

    def test_numpy_nested_payload_and_failure_journal(self):
        # The defensive JSON encoder supports unexpected nested NumPy values,
        # without stringifying numbers or inventing successful DB writes.
        journal = ConversionJournal('run.CDF', ('workflow', 'equilibrium'))
        journal.record_writer('workflow', [
            {'ids': 'workflow', 'occurrence': np.int64(0)}], 0.1)
        journal.record_writer('equilibrium', [], 0.2,
                              ValueError('invalid coordinates'))
        journal.fail(ValueError('invalid coordinates'))
        journal.data['diagnostic'] = {
            'n': np.int32(7), 'timing': np.float32(1.25),
            'valid': np.bool_(True), 'counts': np.array([3, 5], dtype=np.int32)}
        with tempfile.TemporaryDirectory() as td:
            output = Path(td) / 'partial.json'
            journal.save(output)
            import json
            saved = json.loads(output.read_text())
        self.assertEqual(saved['status'], 'FAIL')
        self.assertEqual(saved['diagnostic'],
                         {'n': 7, 'timing': 1.25, 'valid': True,
                          'counts': [3, 5]})
        self.assertEqual(saved['written_ids'],
                         [{'ids': 'workflow', 'occurrence': 0}])
        self.assertEqual([row['status'] for row in saved['writers']],
                         ['WROTE', 'FAIL'])

    def test_put_journal_only_records_success(self):
        import importlib.util
        module_name = '_imas_class_for_test'
        spec = importlib.util.spec_from_file_location(module_name, ROOT / 'imas_class.py')
        module = importlib.util.module_from_spec(spec)
        fake_imas = types.ModuleType('imas')
        with mock.patch.dict(sys.modules, {'imas': fake_imas}):
            spec.loader.exec_module(module)
        conn = object.__new__(module.imas_connection)
        conn.occur = 0
        conn._active_writer = 'equilibrium'
        conn._written_ids = []
        conn.ids = mock.Mock()
        conn.putdb_imas(object())
        self.assertEqual(conn._written_ids, [{'ids': 'equilibrium', 'occurrence': 0}])
        conn.ids.put.side_effect = ValueError('IMAS invalid coordinates')
        with self.assertRaisesRegex(ValueError, 'invalid coordinates'):
            conn.putdb_imas(object())
        self.assertEqual(len(conn._written_ids), 1)

    def test_time_dimension_and_values_are_both_required(self):
        original = np.array([0.1, 0.3, 0.6])
        self.assertEqual(e2e.check_time_vector('equilibrium', original, original, 3), 3)
        with self.assertRaisesRegex(AssertionError, 'slice count'):
            e2e.check_time_vector('equilibrium', original, original, 2)
        with self.assertRaisesRegex(AssertionError, 'differs'):
            e2e.check_time_vector('equilibrium', [0.1, 0.2, 0.6], original, 3)

    def test_readback_schema_time_and_radial_checks(self):
        times = np.array([0.1, 0.3, 0.6])
        outline = types.SimpleNamespace(r=[1, 2], z=[3, 4])
        eq = mock.Mock()
        eq.time = times
        eq.time_slice = [types.SimpleNamespace(boundary=types.SimpleNamespace(outline=outline))] * 3
        self.assertEqual(e2e.validate_readback('equilibrium', eq, times)['native_TIME3_samples'], 3)
        eq.time = np.array([0.1, 0.2, 0.6])
        with self.assertRaises(AssertionError):
            e2e.validate_readback('equilibrium', eq, times)
        profile = types.SimpleNamespace(grid=types.SimpleNamespace(rho_tor_norm=[0., 1.]),
                                        electrons=types.SimpleNamespace(temperature=[1., 2.]))
        cp = mock.Mock()
        cp.time = times
        cp.profiles_1d = [profile] * 3
        self.assertEqual(e2e.validate_readback('core_profiles', cp, times)['native_TIME3_samples'], 3)
        profile.electrons.temperature = [1., 2., 3.]
        with self.assertRaisesRegex(AssertionError, 'radial grid'):
            e2e.validate_readback('core_profiles', cp, times)

    def test_sawteeth_native_time_and_crash_readback(self):
        from sawteeth_mapping import build_sawteeth_plan, IMAS_MISSING_FLOAT
        from verify_sawteeth_mapping import FakeCDF
        plan = build_sawteeth_plan(FakeCDF([.1,.2,.3,.4], [.05,.15,.15,.35],
                                                time3=[.1,.15,.2,.25,.3,.35,.4],
                                                xiq1=[-.1,.2,.3,-.2]))
        self.assertEqual(plan.events, 2)
        ids = mock.Mock()
        ids.time = plan.time.copy()
        ids.diagnostics.previous_crash_time = plan.previous_crash_time.copy()
        ids.diagnostics.previous_period = plan.previous_period.copy()
        ids.diagnostics.rho_tor_norm_q1 = plan.rho_tor_norm_q1.copy()
        ids.crash_trigger = np.array([], dtype=int)
        checks = e2e.validate_readback('sawteeth', ids, np.arange(7), sawteeth_source=plan)
        self.assertEqual(checks['observed_new_crashes'], 2)
        self.assertEqual(checks['native_TIME_samples'], 4)
        self.assertEqual(ids.diagnostics.rho_tor_norm_q1[0], IMAS_MISSING_FLOAT)
        ids.diagnostics.previous_crash_time[2] = .14
        with self.assertRaisesRegex(AssertionError, 'previous_crash_time'):
            e2e.validate_readback('sawteeth', ids, np.arange(7), sawteeth_source=plan)
        ids.diagnostics.previous_crash_time = plan.previous_crash_time.copy()
        ids.crash_trigger = np.array([1, 1])
        with self.assertRaisesRegex(AssertionError, 'invented crash mechanism'):
            e2e.validate_readback('sawteeth', ids, np.arange(7), sawteeth_source=plan)

    def test_core_transport_mmm_native_grid_readback(self):
        times = np.array([0.1, 0.3])
        x = np.array([[0.2, 0.5], [0.25, 0.55]])
        xb = np.array([[0.3, 0.6], [0.35, 0.65]])
        source = {'X': x, 'XB': xb, 'mode_models': {'mmm07_mtm'}}

        def build_model(name, use_x):
            profiles = []
            for i, time in enumerate(times):
                profiles.append(types.SimpleNamespace(
                    time=float(time),
                    grid_d=types.SimpleNamespace(rho_tor_norm=xb[i].copy()),
                    grid_flux=types.SimpleNamespace(
                        rho_tor_norm=(x[i] if use_x else xb[i]).copy())))
            return types.SimpleNamespace(
                identifier=types.SimpleNamespace(name=name),
                profiles_1d=profiles)

        ct = mock.Mock()
        ct.time = times
        ct.model = [build_model('combined', False),
                    build_model('mmm07_total', False),
                    build_model('mmm07_mtm', True)]
        check = e2e.validate_readback('core_transport', ct, times, source)
        self.assertEqual(check['mmm_native_coordinate_models'], 2)
        self.assertEqual(check['native_TIME3_samples'], 2)
        ct.model[-1].profiles_1d[0].grid_flux.rho_tor_norm = xb[0]
        with self.assertRaisesRegex(AssertionError, 'grid_flux/X differs'):
            e2e.validate_readback('core_transport', ct, times, source)

    def test_discovery_matches_every_present_case_including_nested(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / 'd3d').mkdir()
            (root / 'd3d' / '1.CDF').write_bytes(b'CDF')
            (root / 'nstx.CDF').write_bytes(b'CDF')
            self.assertEqual(len(e2e.discover_cases(root, [])), 2)
            self.assertEqual([x.stem for x in e2e.discover_cases(root, ['1'])], ['1'])
            with self.assertRaises(FileNotFoundError):
                e2e.discover_cases(root, ['missing'])

    def test_no_runtime_success_claimed_for_plan_only(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / 'test.CDF').write_bytes(b'not a real CDF')
            result = root / 'reports'
            self.assertEqual(e2e.main(['--plan-only', '--data-root', td, '--report-dir', str(result)]), 0)
            self.assertIn('"status": "NOT_RUN"', (result / 'summary.json').read_text())

    def test_mock_process_readback_and_partial_failure(self):
        import json
        cdf = ROOT / 'examples/data/d3d/147634T13/147634T13.CDF'
        _t, t3 = e2e.source_time3(cdf)
        self.assertNotEqual(len(_t), len(t3), 'T13 must exercise unequal scalar/profile grids')
        outline = types.SimpleNamespace(r=[1, 2], z=[3, 4])
        eq = mock.Mock()
        eq.time = t3
        eq.time_slice = [types.SimpleNamespace(boundary=types.SimpleNamespace(outline=outline))] * len(t3)
        db = mock.Mock()
        db.get.return_value = eq
        fake_imas = types.SimpleNamespace(ids_defs=types.SimpleNamespace(HDF5_BACKEND=9),
                                          DBEntry=mock.Mock(return_value=db))
        with tempfile.TemporaryDirectory() as td:
            args = e2e.parse_arguments(['--case', cdf.stem, '--ids', 'equilibrium',
                                        '--report-dir', td])
            args.data_root = (ROOT / 'examples/data').resolve()
            args.report_dir = Path(td)
            payload = {'format_version': 1, 'status': 'PASS', 'error': None,
                       'machine': 'D3D', 'shot': 147634, 'run': 2013,
                       'dd_version': '4.1.1', 'occurrence': 0,
                       'writers': [{'ids': 'equilibrium', 'status': 'WROTE', 'write_count': 1}],
                       'written_ids': [{'ids': 'equilibrium', 'occurrence': 0}]}

            def simulate(args_command, **kwargs):
                self.assertIn('--ids', args_command)
                self.assertIn('--report-json', args_command)
                target = Path(args_command[args_command.index('--report-json') + 1])
                target.write_text(json.dumps(payload))
                return types.SimpleNamespace(returncode=0)

            with mock.patch.object(e2e.subprocess, 'run', side_effect=simulate):
                output = e2e.run_one(cdf, args, Path(td) / 'db', fake_imas)
            self.assertEqual(output['status'], 'PASS')
            self.assertEqual(len(output['readback']), 1)
            self.assertEqual(output['readback'][0]['checks']['native_TIME3_samples'], 7)
            self.assertEqual(db.get.call_args.args[0], 'equilibrium')
            self.assertTrue(db.close.called)
            payload['status'] = 'FAIL'
            payload['error'] = 'CoordinateError: mismatch'
            with mock.patch.object(e2e.subprocess, 'run', side_effect=simulate):
                result = e2e.run_one(cdf, args, Path(td) / 'db', fake_imas)
            self.assertEqual(result['status'], 'FAIL')
            self.assertEqual(result['readback'][0]['status'], 'PASS',
                             'readback of earlier successful puts must survive later failures')

    def test_validation_must_remain_enabled(self):
        with mock.patch.dict('os.environ', {'IMAS_AL_DISABLE_VALIDATE': '1'}):
            with self.assertRaisesRegex(RuntimeError, 'must NOT be enabled'):
                e2e.require_imas_runtime()


if __name__ == '__main__':
    unittest.main(verbosity=2)
