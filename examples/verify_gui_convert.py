#!/usr/bin/python
"""Tests of the conversion interface: option schema, planning and the overwrite guard, the job runner, and the POST API.

The planner and the job runner are tested without IMAS (a stub converter stands in for ``convert.py`` where a process is
needed).  Real conversions and real plotting are in ``verify_gui_plots.py``.

    python examples/verify_gui_convert.py
"""

import http.client
import json
import os
from pathlib import Path
import shlex
import sys
import tempfile
import textwrap
import threading
import time
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import conversion_plan
import utilities
from gui import __main__ as launcher
from gui import browse, convert_plan, jobs, server

DATA = ROOT / 'examples' / 'data'
NSTX = DATA / 'nstx' / '141711Y05' / '141711Y05.CDF'


def write_cdf(path, device='SYNTHETIC', shot=1):
    """A small classic CDF with the global attributes the converter (and the planner) read: device and shot."""
    from scipy.io import netcdf_file
    from verify_curation import variables
    with netcdf_file(path, 'w') as nc:
        nc.createDimension('TIME', 3)
        nc.createDimension('TIME3', 2)
        nc.createDimension('X', 3)
        nc.device = device
        nc.shot = shot
        for name, var in variables().items():
            output = nc.createVariable(name, 'd', var.dimensions)
            output[:] = var.values
            output.units = var.units


def make_world(tmp, **overrides):
    """Input root with a generated CDF (device ``SYNTHETIC``, shot 1) and an output root."""
    tmp = Path(tmp)
    (tmp / 'inputs').mkdir(exist_ok=True)
    (tmp / 'out').mkdir(exist_ok=True)
    cdf = tmp / 'inputs' / 'run1Z10.CDF'
    if not cdf.exists():
        write_cdf(cdf)
    (tmp / 'inputs' / 'geom.json').write_text('{}')
    config = convert_plan.Config([tmp / 'inputs'], [tmp / 'out'], sys.executable, True, False, home=tmp / 'home')
    for key, value in overrides.items():
        setattr(config, key, value)
    return config, cdf


class FakeCdfIdentity:
    """``identify`` reads the CDF; tests of the planner use the generated one, whose device and shot are fixed below."""


class SchemaTests(unittest.TestCase):
    def test_every_converter_option_is_decided_by_the_gui(self):
        """A new option in ``convert.py`` fails this test until someone exposes it or marks it managed."""
        dests = {a.dest for a in utilities.build_parser()._actions if a.dest != 'help'}
        decided = set(convert_plan.EXPOSED) | set(convert_plan.MANAGED)
        self.assertEqual(dests - decided, set(), 'the converter gained an option the GUI has not decided about')
        self.assertEqual(decided - dests, set(), 'the GUI names an option the converter no longer has')
        self.assertFalse(set(convert_plan.EXPOSED) & set(convert_plan.MANAGED))

    def test_the_schema_comes_from_the_converter_parser_with_its_help_text(self):
        schema = convert_plan.option_schema()
        by_dest = {o['dest']: o for o in schema['options']}
        self.assertEqual(by_dest['backend']['choices'], ['MDS', 'MDSPLUS', 'HDF', 'HDF5'])
        self.assertIn('Data Dictionary', by_dest['dd_version']['help'])
        self.assertEqual(by_dest['new_db']['kind'], 'flag')
        self.assertEqual(by_dest['verbose']['kind'], 'count')
        self.assertEqual(by_dest['report_json']['managed_by_gui'], 'always written into the job directory')
        self.assertEqual(schema['ids'], list(conversion_plan.WRITER_ORDER))
        self.assertEqual(schema['default_backend'], 'HDF5')

    def test_the_converter_default_backend_is_not_silently_inherited(self):
        """The parser's own default is MDS; the plan always sends the backend explicitly."""
        with tempfile.TemporaryDirectory() as tmp:
            config, cdf = make_world(tmp)
            plan = convert_plan.validate({'cdf': str(cdf), 'options': {}}, config)
            self.assertEqual(plan['argv'][plan['argv'].index('-b') + 1], 'HDF5')


class PlanTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.config, self.cdf = make_world(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def plan(self, options=None, **request):
        return convert_plan.validate({'cdf': str(self.cdf), 'options': options or {}, **request}, self.config)

    def errors(self, options=None, **request):
        with self.assertRaises(convert_plan.PlanError) as caught:
            self.plan(options, **request)
        return {e['field']: e['message'] for e in caught.exception.errors}

    def test_the_target_is_computed_the_way_the_converter_computes_it(self):
        plan = self.plan({'dd_version': '4.1.1'})
        target = plan['target']
        self.assertEqual((target['machine'], target['shot'], target['run'], target['backend']), ('SYNTHETIC', 1, 2610, 'HDF5'))
        expected = self.config.scratch_database / 'SYNTHETIC' / '4' / '1' / '2610'
        self.assertEqual(target['primary'], str(expected))
        argv = plan['argv']
        self.assertEqual(argv[argv.index('-i') + 1], str(self.cdf.resolve()))
        self.assertEqual(argv[argv.index('--input-format') + 1], 'netcdf')
        self.assertEqual(argv[argv.index('-b') + 1], 'HDF5')
        self.assertEqual(argv[argv.index('-u') + 1], str(self.config.scratch_database))

    def test_a_machine_override_and_run_number_change_the_target(self):
        target = self.plan({'dd_version': '3.39.0', 'machine': 'NSTX', 'run': 77})['target']
        self.assertEqual(target['primary'], str(self.config.scratch_database / 'NSTX' / '3' / '1' / '77'))

    def test_auto_dd_version_checks_both_majors_and_says_so(self):
        plan = self.plan({'dd_version': 'auto'})
        self.assertEqual(len(plan['target']['entry_paths']), 2)
        self.assertTrue(any('DD version is auto' in w for w in plan['warnings']))

    def test_an_existing_entry_needs_the_exact_path_typed_back(self):
        entry = self.config.scratch_database / 'SYNTHETIC' / '4' / '1' / '2610'
        entry.mkdir(parents=True)
        (entry / 'master.h5').write_text('x')
        plan = self.plan({'dd_version': '4.1.1'})
        self.assertTrue(plan['target']['exists'] and plan['needs_confirmation'])
        self.assertEqual(plan['confirm_phrase'], str(entry))
        self.assertFalse(convert_plan.confirm(plan, None))
        self.assertFalse(convert_plan.confirm(plan, str(entry) + '/'))
        self.assertFalse(convert_plan.confirm(plan, 'yes'))
        self.assertTrue(convert_plan.confirm(plan, str(entry)))
        self.assertTrue(any('replaces it' in w for w in plan['warnings']))

    def test_extending_an_existing_entry_without_new_also_needs_confirmation(self):
        entry = self.config.scratch_database / 'SYNTHETIC' / '4' / '1' / '2610'
        entry.mkdir(parents=True)
        (entry / 'wall.h5').write_text('x')
        plan = self.plan({'dd_version': '4.1.1', 'new_db': False})
        self.assertTrue(plan['needs_confirmation'])
        self.assertNotIn('-new', plan['argv'])
        self.assertTrue(any('existing entry' in w for w in plan['warnings']))

    def test_a_missing_entry_without_create_new_is_an_error(self):
        self.assertIn('new_db', self.errors({'dd_version': '4.1.1', 'new_db': False}))

    def test_the_next_free_run_number_skips_every_existing_entry(self):
        base = self.config.scratch_database / 'SYNTHETIC' / '4' / '1'
        for run in (2610, 2611, 2613):
            (base / str(run)).mkdir(parents=True)
            (base / str(run) / 'master.h5').write_text('x')
        target = self.plan({'dd_version': '4.1.1'})['target']
        self.assertEqual((target['run'], target['exists'], target['next_free_run']), (2610, True, 2612))
        self.assertEqual(self.plan({'dd_version': '4.1.1', 'run': 2612})['target']['next_free_run'], None)      # nothing exists there
        taken = self.plan({'dd_version': '4.1.1', 'run': 2613})['target']
        self.assertEqual(taken['next_free_run'], 2614)

    def test_no_run_number_is_suggested_when_the_entry_cannot_be_checked(self):
        self.assertIsNone(self.plan({'backend': 'MDS', 'dd_version': '4.1.1'})['target']['next_free_run'])
        self.assertIsNone(convert_plan.next_free_run(self.config.scratch_database, 'M', 1, 5, 'MDS', '4.1.1'))

    def test_an_empty_directory_is_not_an_existing_entry(self):
        (self.config.scratch_database / 'SYNTHETIC' / '4' / '1' / '2610').mkdir(parents=True)
        self.assertFalse(self.plan({'dd_version': '4.1.1'})['needs_confirmation'])

    def test_mds_entries_cannot_be_checked_so_they_always_need_confirmation(self):
        plan = self.plan({'backend': 'MDS', 'dd_version': '4.1.1'})
        self.assertTrue(plan['needs_confirmation'] and plan['target']['exists_unknown'])
        self.assertTrue(any('cannot tell' in w for w in plan['warnings']))

    def test_database_modes_scratch_custom_and_default(self):
        scratch = self.plan({'dd_version': '4.1.1'}, database={'mode': 'scratch'})
        self.assertEqual(scratch['target']['database_root'], str(self.config.scratch_database))
        custom = self.plan({'dd_version': '4.1.1'}, database={'mode': 'custom', 'path': str(self.tmp.name) + '/out/mydb'})
        self.assertTrue(custom['target']['primary'].startswith(str(Path(self.tmp.name).resolve() / 'out' / 'mydb')))
        self.assertIn('-u', custom['argv'])
        relative = self.plan({'dd_version': '4.1.1'}, database={'mode': 'custom', 'path': 'relative-db'})
        self.assertTrue(relative['target']['database_root'].endswith('/out/relative-db'))

    def test_a_custom_database_outside_the_output_roots_is_refused(self):
        for path in ('/tmp/elsewhere', str(self.tmp.name) + '/inputs', str(self.tmp.name) + '/out/../inputs', '', None):
            errors = self.errors({'dd_version': '4.1.1'}, database={'mode': 'custom', 'path': path})
            self.assertIn('database.path', errors, path)

    def test_the_default_database_needs_server_permission_and_then_passes_no_u(self):
        self.assertIn('database.mode', self.errors({'dd_version': '4.1.1'}, database={'mode': 'default'}))
        self.config.allow_default_database = True
        plan = self.plan({'dd_version': '4.1.1'}, database={'mode': 'default'})
        self.assertNotIn('-u', plan['argv'])                      # the converter's own default location
        self.assertEqual(plan['target']['database_root'], str(self.config.home / 'public' / 'imasdb'))

    def test_unknown_options_and_values_are_rejected_with_field_names(self):
        errors = self.errors({'dd_version': '4.1.1', 'sneaky': 1, 'user': '/etc', 'report_json': '/tmp/x'})
        self.assertIn('options', errors)
        self.assertIn('sneaky', errors['options'])
        self.assertIn('user', errors['options'])
        for options, field in (({'backend': 'ORACLE'}, 'backend'), ({'dd_version': '4.1'}, 'dd_version'), ({'dd_version': 'latest; rm -rf /'}, 'dd_version'),
                               ({'run': -1}, 'run'), ({'run': True}, 'run'), ({'run': '5'}, 'run'), ({'occur': -2}, 'occur'), ({'nzones': 10 ** 6}, 'nzones'),
                               ({'verbose': 9}, 'verbose'), ({'ip_direction': 'sideways'}, 'ip_direction'), ({'log': 'loud'}, 'log'),
                               ({'machine': 'NSTX; ls'}, 'machine'), ({'machine': 'a' * 40}, 'machine'), ({'lenient_mapping': 'yes'}, 'lenient_mapping'),
                               ({'ids': ['nope']}, 'ids'), ({'ids': []}, 'ids'), ({'ids': 'core_profiles'}, 'ids'), ({'moments': 'maybe'}, 'moments')):
            self.assertIn(field, self.errors({'dd_version': '4.1.1', **options}), options)

    def test_the_cdf_must_be_an_existing_cdf_below_an_input_root(self):
        (self.config.input_roots[0] / 'notes.txt').write_text('x')
        for requested in ('/etc/passwd', '../out/x.CDF', str(Path(self.tmp.name) / 'out' / 'x.CDF'), 'inputs', '', 'missing.CDF', 'notes.txt'):
            with self.assertRaises(convert_plan.PlanError, msg=requested) as caught:
                convert_plan.validate({'cdf': requested, 'options': {}}, self.config)
            self.assertIn('cdf', {e['field'] for e in caught.exception.errors}, requested)

    def test_a_symlink_escaping_the_input_root_is_refused(self):
        outside = Path(self.tmp.name) / 'out' / 'secret.CDF'
        outside.write_bytes(self.cdf.read_bytes())
        os.symlink(outside, self.config.input_roots[0] / 'link.CDF')
        with self.assertRaises(convert_plan.PlanError):
            convert_plan.validate({'cdf': str(self.config.input_roots[0] / 'link.CDF'), 'options': {}}, self.config)

    def test_geometry_files_must_be_json_below_the_input_roots(self):
        good = self.plan({'dd_version': '4.1.1', 'pf_passive_geometry': 'geom.json'})
        self.assertEqual(good['argv'][good['argv'].index('--pf-passive-geometry') + 1], str((self.config.input_roots[0] / 'geom.json').resolve()))
        for path in ('/etc/passwd', '../out/geom.json', 'run1Z10.CDF'):
            self.assertIn('pf_active_geometry', self.errors({'dd_version': '4.1.1', 'pf_active_geometry': path}), path)

    def test_ids_are_ordered_by_the_converter_and_all_selected_means_default(self):
        order = list(conversion_plan.WRITER_ORDER)
        plan = self.plan({'dd_version': '4.1.1', 'ids': ['summary', 'wall', 'core_profiles', 'wall']})
        self.assertEqual(plan['selected_ids'], [i for i in order if i in {'summary', 'wall', 'core_profiles'}])
        self.assertEqual(plan['argv'][plan['argv'].index('--ids') + 1], ','.join(plan['selected_ids']))
        self.assertTrue(any('partial' in w for w in plan['warnings']))
        full = self.plan({'dd_version': '4.1.1', 'ids': order})
        self.assertNotIn('--ids', full['argv'])
        self.assertEqual(full['selected_ids'], order)

    def test_every_option_reaches_the_command_in_the_converters_own_spelling(self):
        plan = self.plan({'backend': 'HDF5', 'dd_version': '4.1.1', 'run': 5, 'occur': 2, 'machine': 'NSTX', 'ip_direction': 'cw', 'bt_direction': 'ccw',
                          'pf_passive_geometry': 'geom.json', 'ids': ['wall'], 'lenient_mapping': True, 'moments': 'yes', 'nzones': 101,
                          'ec_nml_launch': 'yes', 'debug': True, 'debug_max': True, 'verbose': 2, 'log': 'info', 'new_db': True})
        argv = plan['argv']
        for flag in ('-r', '-o', '-m', '--ip-direction', '--bt-direction', '--pf-passive-geometry', '--ids', '--lenient-mapping', '-mom', '-zones',
                     '-ecnml', '--debug', '--debug-max', '-l', '-new', '--report-json'):
            self.assertIn(flag, argv, flag)
        self.assertEqual(argv.count('-v'), 2)
        parsed = utilities.build_parser().parse_args(argv)           # the converter's own parser accepts exactly this command
        self.assertEqual((parsed.run, parsed.occur, parsed.machine, parsed.ip_direction, parsed.bt_direction, parsed.moments, parsed.nzones),
                         (5, 2, 'NSTX', 'cw', 'ccw', 'yes', 101))
        self.assertEqual((parsed.ids, parsed.lenient_mapping, parsed.new_db, parsed.verbose, parsed.log), ('wall', True, True, 2, 'info'))

    def test_the_default_command_is_minimal_and_parses(self):
        argv = self.plan({'dd_version': '4.1.1'})['argv']
        self.assertEqual([a for a in argv if a.startswith('-') and not a.startswith('--')], ['-i', '-b', '-u', '-new'])
        utilities.build_parser().parse_args(argv)

    def test_the_command_is_an_argument_list_and_awkward_paths_stay_single_arguments(self):
        awkward_dir = self.config.input_roots[0] / "my dir; echo 'x' $(id)"
        awkward_dir.mkdir()
        awkward = awkward_dir / 'run1Z10.CDF'
        awkward.write_bytes(self.cdf.read_bytes())
        plan = convert_plan.validate({'cdf': str(awkward), 'options': {'dd_version': '4.1.1'}}, self.config)
        self.assertEqual(plan['argv'][1], str(awkward.resolve()))
        self.assertEqual(shlex.split(plan['display'])[2:], plan['argv'])           # displayed quoting round-trips exactly
        self.assertTrue(plan['display'].startswith('python convert.py '))

    def test_a_file_name_without_a_shot_and_run_is_a_field_error_not_a_crash(self):
        (self.config.input_roots[0] / 'geom.CDF').write_bytes(self.cdf.read_bytes())
        with self.assertRaises(convert_plan.PlanError) as caught:
            convert_plan.validate({'cdf': 'geom.CDF', 'options': {}}, self.config)
        self.assertEqual([e['field'] for e in caught.exception.errors], ['cdf'])
        self.assertIn('shot and run id', caught.exception.errors[0]['message'])

    def test_the_report_goes_to_the_given_job_path_and_a_failed_validation_returns_every_error(self):
        plan = convert_plan.validate({'cdf': str(self.cdf), 'options': {'dd_version': '4.1.1'}}, self.config, report_json='/jobs/ab/r.json')
        self.assertEqual(plan['argv'][plan['argv'].index('--report-json') + 1], '/jobs/ab/r.json')
        errors = self.errors({'dd_version': 'x', 'backend': 'y', 'run': -3})
        self.assertGreaterEqual(len(errors), 3)

    def test_a_large_cdf_is_never_refused_for_its_size_but_a_large_report_still_is(self):
        """TRANSP output is routinely hundreds of MB: the cap is for reports the viewer reads into memory."""
        with mock.patch.object(browse, 'MAX_FILE_BYTES', 1):                      # every file here is now "large"
            plan = self.plan({'dd_version': '4.1.1'})
            self.assertEqual(plan['cdf'], str(self.cdf.resolve()))
            self.assertEqual(convert_plan.scan_inputs(self.config.input_roots)[0]['path'], str(self.cdf))
            report = self.config.input_roots[0] / 'r.json'
            report.write_text('{}')
            with self.assertRaises(browse.AccessError) as caught:
                browse.resolve(self.config.input_roots, str(report))
            self.assertIn('size limit for report files', str(caught.exception))
            self.assertRegex(str(caught.exception), r'file is \d+ MB, over the \d+ MB')
            with self.assertRaises(convert_plan.PlanError):                       # geometry JSON is read, so it keeps the cap
                self.plan({'dd_version': '4.1.1', 'pf_passive_geometry': 'geom.json'})

    def test_a_real_size_over_the_report_cap_is_still_a_cdf_the_planner_accepts(self):
        big = self.config.input_roots[0] / 'big1Z10.CDF'
        write_cdf(big)
        with big.open('ab') as stream:
            stream.truncate(browse.MAX_FILE_BYTES + 5 * 1024 * 1024)               # a sparse file larger than the report cap
        self.assertGreater(big.stat().st_size, browse.MAX_FILE_BYTES)
        try:
            plan = convert_plan.validate({'cdf': str(big), 'options': {'dd_version': '4.1.1'}}, self.config)
        except convert_plan.PlanError as exc:
            self.assertNotIn('size limit', str(exc))                              # reading the padded file may fail, never for size
        else:
            self.assertEqual(plan['target']['shot'], 1)

    def test_the_converter_option_is_worded_for_what_it_does(self):
        plan = self.plan({'dd_version': '4.1.1', 'lenient_mapping': True})
        text = ' '.join(plan['warnings'])
        self.assertIn('cannot be mapped', text)
        self.assertIn('--lenient-mapping', text)
        self.assertNotIn('Lenient mapping', text)
        self.assertIn('--lenient-mapping', plan['argv'])                          # the converter's own flag is unchanged

    def test_the_cdf_identity_is_read_from_the_file_not_trusted_from_the_request(self):
        machine, shot, run, rid = convert_plan.identify(self.cdf)
        self.assertEqual((machine, shot, run), ('SYNTHETIC', 1, 2610))
        with self.assertRaises(convert_plan.PlanError):
            convert_plan.identify(self.config.input_roots[0] / 'geom.json')              # no shot/run in the name
        unreadable = self.config.input_roots[0] / '5Y01.CDF'
        unreadable.write_text('not a netCDF file')
        with self.assertRaises(convert_plan.PlanError) as caught:
            convert_plan.identify(unreadable)
        self.assertIn('cannot read the CDF', caught.exception.errors[0]['message'])

    @unittest.skipUnless(NSTX.is_file(), 'bundled NSTX fixture required')
    def test_the_bundled_nstx_case_plans_the_entry_seen_in_real_conversions(self):
        config = convert_plan.Config([DATA], [Path(self.tmp.name) / 'out'], sys.executable, True, False)
        plan = convert_plan.validate({'cdf': str(NSTX), 'options': {'dd_version': '4.1.1'}}, config)
        self.assertEqual((plan['target']['machine'], plan['target']['shot'], plan['target']['run']), ('NSTX', 141711, 2505))
        self.assertTrue(plan['target']['primary'].endswith('/imasdb/NSTX/4/141711/2505'))

    def test_only_the_run_cdf_is_offered_not_its_companions(self):
        for name in ('141711Y05.CDF', '141711Y05.cdf', '38530P94.CDF', '89846.CDF', '120967E03.CDF', 'run1Z10.CDF'):
            self.assertTrue(convert_plan.is_run_cdf(name), name)
        for name in ('141711Y05PH.CDF', '141711Y05ph.cdf', '120967E03_ps.cdf', '120967E03_fi_1.cdf', '120967E03_state.CDF', 'plasma_state.cdf',
                     '141711Y05.CDF.bak', '141711Y05.DATA1', 'notes.txt', 'CDF', '.CDF', '141711Y05_TAR.GZ1'):
            self.assertFalse(convert_plan.is_run_cdf(name), name)
        root = self.config.input_roots[0]
        for name in ('77Y01.CDF', '77Y01PH.CDF', '77Y01_ps.cdf', '77Y01_fi_2.cdf'):
            (root / name).write_bytes(b'x')
        names = [Path(f['path']).name for f in convert_plan.scan_inputs(self.config.input_roots)]
        self.assertIn('77Y01.CDF', names)
        self.assertEqual([n for n in names if n.startswith('77Y01')], ['77Y01.CDF'])
        for name in ('77Y01.CDF', '77Y01PH.CDF', '77Y01_ps.cdf', '77Y01_fi_2.cdf'):
            (root / name).unlink()
        plan = self.plan({'dd_version': '4.1.1'})
        self.assertTrue(plan['cdf'])                                                    # typing a path of another name is not blocked here

    def test_scans_list_cdfs_and_entries(self):
        found = convert_plan.scan_inputs(self.config.input_roots)
        self.assertEqual([Path(f['path']).name for f in found], ['run1Z10.CDF'])
        entry = self.config.scratch_database / 'M' / '4' / '7' / '3'
        entry.mkdir(parents=True)
        for name in ('master.h5', 'wall.h5', 'core_transport_1.h5', 'core_transport.h5'):
            (entry / name).write_bytes(b'12345')
        entries = convert_plan.scan_entries([self.config.scratch_database])
        self.assertEqual((entries[0]['machine'], entries[0]['dd_major'], entries[0]['shot'], entries[0]['run']), ('M', '4', '7', '3'))
        self.assertEqual(entries[0]['ids'], ['core_transport', 'wall'])


STUB = textwrap.dedent('''\
    #!/usr/bin/python
    import json, os, sys, time
    args = sys.argv[1:]
    mode = os.environ.get('STUB_MODE', 'ok')
    print('stub converter', flush=True)
    print('PYTHONPATH=' + os.environ.get('PYTHONPATH', '<unset>'), flush=True)
    print('cwd=' + os.getcwd(), flush=True)
    print('args=' + json.dumps(args), flush=True)
    report = args[args.index('--report-json') + 1]
    if mode == 'sleep':
        print('sleeping', flush=True)
        time.sleep(60)
    if mode == 'fail':
        print('boom', file=sys.stderr, flush=True)
        sys.exit(3)
    json.dump({'format_version': 1, 'status': 'PASS', 'writers': [], 'selected_ids': [], 'publication': {'status': 'complete'}}, open(report, 'w'))
    print('done', flush=True)
''')


class JobTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.stub = Path(self.tmp.name) / 'convert_stub.py'
        self.stub.write_text(STUB)
        self.manager = jobs.JobManager(Path(self.tmp.name) / 'work', sys.executable, ROOT, self.stub)

    def tearDown(self):
        for job_id in list(self.manager._processes):
            self.manager.cancel(job_id)
        end = time.time() + 20
        while self.manager._processes and time.time() < end:      # the watcher thread writes the final record
            time.sleep(0.05)
        self.tmp.cleanup()

    def start(self, mode='ok'):
        os.environ['STUB_MODE'] = mode
        job_id, directory = self.manager.reserve()
        report = directory / 'conversion-report.json'
        plan = {'argv': ['-i', 'x.CDF', '--report-json', str(report)], 'display': 'python convert.py -i x.CDF'}
        return self.manager.start(job_id, plan)

    def test_a_job_runs_logs_writes_its_report_and_records_the_outcome(self):
        job = self.start()
        self.assertEqual(job['state'], 'running')
        final = self.manager.wait(job['id'], 20)
        self.assertEqual((final['state'], final['returncode'], final['report_written']), ('succeeded', 0, True))
        log = self.manager.log(job['id'])
        self.assertIn('stub converter', log['text'])
        self.assertIn('done', log['text'])
        self.assertTrue(log['done'])
        self.assertEqual(json.loads(Path(final['report']).read_text())['status'], 'PASS')
        self.assertEqual(json.loads((self.manager.jobs_dir / job['id'] / 'job.json').read_text())['state'], 'succeeded')

    def test_the_child_runs_in_its_job_directory_without_pythonpath_and_with_the_given_arguments(self):
        os.environ['PYTHONPATH'] = '/old/access/layer'
        try:
            job = self.start()
            self.manager.wait(job['id'], 20)
        finally:
            os.environ.pop('PYTHONPATH')
        text = self.manager.log(job['id'])['text']
        self.assertIn('PYTHONPATH=<unset>', text)
        self.assertIn(f"cwd={self.manager.jobs_dir / job['id']}", text)
        self.assertIn('"-i", "x.CDF"', text)

    def test_a_failing_converter_is_a_failed_job_with_its_exit_code_and_stderr_in_the_log(self):
        job = self.start('fail')
        final = self.manager.wait(job['id'], 20)
        self.assertEqual((final['state'], final['returncode'], final['report_written']), ('failed', 3, False))
        self.assertIn('boom', self.manager.log(job['id'])['text'])

    def test_the_log_can_be_followed_by_offset(self):
        job = self.start()
        self.manager.wait(job['id'], 20)
        first = self.manager.log(job['id'], 0)
        part = self.manager.log(job['id'], 10)
        self.assertEqual(part['text'], first['text'][10:])
        self.assertEqual(self.manager.log(job['id'], first['next'])['text'], '')

    def test_only_one_conversion_runs_at_a_time(self):
        job = self.start('sleep')
        with self.assertRaises(jobs.JobBusy):
            self.start('ok')
        self.manager.cancel(job['id'])
        self.manager.wait(job['id'], 20)
        self.assertEqual(self.start('ok')['state'], 'running')

    def test_cancel_stops_the_whole_process_group_and_marks_the_job_cancelled(self):
        job = self.start('sleep')
        time.sleep(0.5)
        self.assertTrue(self.manager.cancel(job['id']))
        final = self.manager.wait(job['id'], 20)
        self.assertEqual(final['state'], 'cancelled')
        self.assertFalse(self.manager.cancel(job['id']))                # nothing left to cancel

    def test_a_job_left_running_by_a_dead_server_is_lost_not_succeeded(self):
        job = self.start('sleep')
        record = json.loads((self.manager.jobs_dir / job['id'] / 'job.json').read_text())
        self.assertEqual(record['state'], 'running')
        fresh = jobs.JobManager(self.manager.work_dir, sys.executable, ROOT, self.stub)       # a restarted server
        self.assertEqual(fresh.get(job['id'])['state'], 'lost')
        self.assertIn('outcome is unknown', fresh.get(job['id'])['note'])

    def test_job_ids_cannot_escape_the_jobs_directory(self):
        for bad in ('../x', '..', '', 'ZZ', 'a/b', '0' * 40, None):
            with self.assertRaises(KeyError, msg=str(bad)):
                self.manager.get(bad)

    def test_the_job_list_is_newest_first_and_survives_a_new_manager(self):
        first = self.start()
        self.manager.wait(first['id'], 20)
        second = self.start()
        self.manager.wait(second['id'], 20)
        listed = jobs.JobManager(self.manager.work_dir, sys.executable, ROOT, self.stub).list()
        self.assertEqual([j['id'] for j in listed][:2], [second['id'], first['id']])

    def test_a_missing_interpreter_is_a_failed_job_not_an_exception(self):
        manager = jobs.JobManager(Path(self.tmp.name) / 'w2', '/nonexistent/python', ROOT, self.stub)
        job_id, directory = manager.reserve()
        job = manager.start(job_id, {'argv': ['--report-json', str(directory / 'r.json')], 'display': 'x'})
        self.assertEqual(job['state'], 'failed')
        self.assertIn('could not start', job['note'])


class InterpreterOptionTests(unittest.TestCase):
    """``--imas-python`` takes the interpreter or the virtual-environment directory that has it."""

    def test_a_virtual_environment_directory_resolves_to_its_interpreter(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / 'bin').mkdir()
            os.symlink(sys.executable, Path(tmp) / 'bin' / 'python')
            self.assertEqual(launcher.resolve_python(tmp), str(Path(tmp) / 'bin' / 'python'))
            self.assertEqual(launcher.resolve_python(tmp + '/'), str(Path(tmp) / 'bin' / 'python'))
            self.assertEqual(launcher.resolve_python(str(Path(tmp) / 'bin' / 'python')), str(Path(tmp) / 'bin' / 'python'))

    def test_the_interpreter_itself_and_a_home_relative_path_are_accepted(self):
        self.assertEqual(launcher.resolve_python(sys.executable), sys.executable)
        with mock.patch.dict(os.environ, {'HOME': str(Path(sys.executable).parent.parent)}):
            self.assertEqual(launcher.resolve_python('~/bin/' + Path(sys.executable).name), str(Path(sys.executable).parent.parent / 'bin' / Path(sys.executable).name))

    def test_bad_paths_give_a_reason_instead_of_a_silent_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError) as caught:
                launcher.resolve_python(tmp)
            self.assertIn('directory without bin/python', str(caught.exception))
            (Path(tmp) / 'plain.txt').write_text('x')
            with self.assertRaises(ValueError) as caught:
                launcher.resolve_python(str(Path(tmp) / 'plain.txt'))
            self.assertIn('not executable', str(caught.exception))
            with self.assertRaises(ValueError) as caught:
                launcher.resolve_python(str(Path(tmp) / 'missing'))
            self.assertIn('does not exist', str(caught.exception))

    def test_the_launcher_exits_2_with_the_reason_for_a_bad_interpreter(self):
        import subprocess
        done = subprocess.run([sys.executable, '-m', 'gui', '--imas-python', '/nonexistent/python'], cwd=ROOT, capture_output=True, text=True, timeout=60)
        self.assertEqual(done.returncode, 2)
        self.assertIn('--imas-python: /nonexistent/python does not exist', done.stderr)

    def launch(self, *arguments):
        """Start ``python -m gui`` for real and return ``(process, port, token, startup text)``."""
        import re
        import subprocess
        proc = subprocess.Popen([sys.executable, '-m', 'gui', *arguments], cwd=ROOT, stdout=subprocess.PIPE, text=True,
                                env={**os.environ, 'PYTHONPATH': ''})
        lines = []
        for line in proc.stdout:
            lines.append(line)
            match = re.match(r'Open: http://127\.0\.0\.1:(\d+)/\?token=(\S+)', line)
            if match:
                return proc, int(match.group(1)), match.group(2), ''.join(lines)
        raise AssertionError('the server did not start: ' + ''.join(lines))

    @staticmethod
    def stop_server(proc):
        """Reap a launched test server and close its captured stdout pipe."""
        import subprocess
        if proc.poll() is None:
            proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
        if proc.stdout is not None:
            proc.stdout.close()

    def config_of(self, port, token):
        connection = http.client.HTTPConnection('127.0.0.1', port, timeout=30)
        connection.request('GET', '/api/config', headers={'Host': f'127.0.0.1:{port}', 'X-GUI-Token': token})
        body = json.loads(connection.getresponse().read())
        connection.close()
        return body

    def test_a_virtual_environment_directory_with_enable_run_turns_running_on_when_it_has_imas(self):
        try:
            import imas  # noqa: F401
        except ImportError:
            self.skipTest('IMAS-Python is required to show running turning on (skipped, not passed)')
        prefix = Path(sys.executable).parent.parent
        if not (prefix / 'bin' / 'python').is_file():
            self.skipTest('this interpreter does not live in a bin/ directory of an environment')
        with tempfile.TemporaryDirectory() as tmp:
            proc, port, token, text = self.launch('--root', tmp, '--imas-python', str(prefix) + '/', '--enable-run', '--output-root', str(Path(tmp) / 'o'))
            try:
                config = self.config_of(port, token)
                self.assertTrue(config['run_enabled'], text)
                self.assertTrue(config['imas_python'].endswith('/bin/python'))
                self.assertIn('conversion enabled', text)
            finally:
                self.stop_server(proc)

    def test_without_enable_run_the_launcher_and_the_config_say_running_is_off_even_with_imas(self):
        try:
            import imas  # noqa: F401
        except ImportError:
            self.skipTest('IMAS-Python is required (skipped, not passed)')
        with tempfile.TemporaryDirectory() as tmp:
            proc, port, token, text = self.launch('--root', tmp, '--imas-python', sys.executable, '--output-root', str(Path(tmp) / 'o'))
            try:
                config = self.config_of(port, token)
                self.assertEqual((config['run_enabled'], config['run_requested']), (False, False))
                self.assertIn('add --enable-run to allow conversions', text)
                self.assertTrue(config['imas_python'])                                  # plotting still works
            finally:
                self.stop_server(proc)

    def test_a_directory_without_imas_starts_the_server_and_says_why_running_is_off(self):
        import re
        import subprocess
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / 'bin').mkdir()
            os.symlink(sys.executable, Path(tmp) / 'bin' / 'python')
            env = {**os.environ, 'PYTHONPATH': ''}
            proc = subprocess.Popen([sys.executable, '-m', 'gui', '--root', tmp, '--imas-python', tmp, '--enable-run', '--output-root', str(Path(tmp) / 'o')],
                                    cwd=ROOT, stdout=subprocess.PIPE, text=True, env=env)
            try:
                lines = []
                for line in proc.stdout:
                    lines.append(line)
                    if line.startswith('Open:'):
                        break
                output = ''.join(lines)
                has_imas = 'imas' in subprocess.run([sys.executable, '-c', 'import importlib.util as u; print(u.find_spec("imas") is not None and "imas")'],
                                                    capture_output=True, text=True).stdout
                if not has_imas:
                    self.assertIn('conversion off', output)
                    self.assertIn('cannot import imas', output)
                self.assertIn('Open: http://127.0.0.1:', output)
            finally:
                self.stop_server(proc)


class Client:
    def __init__(self, srv):
        self.srv, self.port = srv, srv.server_address[1]

    def call(self, method, path, body=None, headers=None, token=True, raw=None, content_type='application/json'):
        connection = http.client.HTTPConnection('127.0.0.1', self.port, timeout=60)
        merged = {'Host': f'127.0.0.1:{self.port}', **(headers or {})}
        if token:
            merged['X-GUI-Token'] = self.srv.token
        data = raw
        if body is not None:
            data = json.dumps(body).encode()
        if data is not None and content_type:
            merged['Content-Type'] = content_type
        connection.request(method, path, body=data, headers=merged)
        response = connection.getresponse()
        payload = response.read()
        connection.close()
        try:
            return response.status, json.loads(payload)
        except ValueError:
            return response.status, payload


class PostApiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.config, cls.cdf = make_world(cls.tmp.name)
        cls.stub = Path(cls.tmp.name) / 'convert_stub.py'
        cls.stub.write_text(STUB)
        cls.readonly_config = convert_plan.Config([Path(cls.tmp.name) / 'inputs'], [Path(cls.tmp.name) / 'out'], sys.executable, False, False)
        cls.srv = server.make_server([Path(cls.tmp.name) / 'inputs'], config=cls.config, repo_root=ROOT, convert_script=cls.stub)
        cls.ro = server.make_server([Path(cls.tmp.name) / 'inputs'], config=cls.readonly_config, repo_root=ROOT)
        for s in (cls.srv, cls.ro):
            threading.Thread(target=s.serve_forever, daemon=True).start()
        cls.client, cls.ro_client = Client(cls.srv), Client(cls.ro)
        cls.request = {'cdf': str(cls.cdf), 'options': {'dd_version': '4.1.1'}}

    @classmethod
    def tearDownClass(cls):
        for job in cls.srv.jobs.list():
            cls.srv.jobs.cancel(job['id'])
        end = time.time() + 20
        while cls.srv.jobs._processes and time.time() < end:
            time.sleep(0.05)
        for s in (cls.srv, cls.ro):
            s.shutdown()
            s.server_close()
        cls.tmp.cleanup()

    def wait_idle(self):
        end = time.time() + 30
        while self.srv.jobs._processes and time.time() < end:
            time.sleep(0.1)

    def test_planning_works_on_a_read_only_server_but_running_is_refused(self):
        status, body = self.ro_client.call('POST', '/api/convert/plan', self.request)
        self.assertEqual(status, 200)
        self.assertFalse(body['run_enabled'])
        status, body = self.ro_client.call('POST', '/api/convert/start', self.request)
        self.assertEqual(status, 403)
        self.assertIn('--enable-run', body['error'])
        self.assertEqual(self.ro_client.call('POST', '/api/jobs/abc/cancel', {})[0], 403)
        self.assertFalse(self.ro.config.enable_run)

    def test_config_reports_what_is_enabled_and_where(self):
        status, body = self.client.call('GET', '/api/config')
        self.assertEqual(status, 200)
        self.assertTrue(body['run_enabled'])
        self.assertEqual(body['scratch_database'], str(self.config.scratch_database))
        self.assertFalse(body['allow_default_database'])
        self.assertEqual(self.ro_client.call('GET', '/api/config')[1]['run_enabled'], False)

    def test_a_post_needs_the_header_token_a_cookie_or_query_token_is_not_enough(self):
        self.assertEqual(self.client.call('POST', '/api/convert/plan', self.request, token=False)[0], 401)
        self.assertEqual(self.client.call('POST', '/api/convert/plan', self.request, token=False,
                                          headers={'Cookie': f'gui_token={self.srv.token}'})[0], 401)
        self.assertEqual(self.client.call('POST', f'/api/convert/plan?token={self.srv.token}', self.request, token=False)[0], 401)
        self.assertEqual(self.client.call('POST', '/api/convert/plan', self.request, headers={'X-GUI-Token': 'wrong'}, token=False)[0], 401)

    def test_a_cross_origin_post_is_refused_even_with_the_token(self):
        for origin in ('http://evil.example', 'https://127.0.0.1', f'http://127.0.0.1:{self.client.port + 1}', 'null'):
            self.assertEqual(self.client.call('POST', '/api/convert/plan', self.request, headers={'Origin': origin})[0], 403, origin)
        own = f'http://127.0.0.1:{self.client.port}'
        self.assertEqual(self.client.call('POST', '/api/convert/plan', self.request, headers={'Origin': own})[0], 200)

    def test_a_foreign_host_header_is_refused_on_posts(self):
        self.assertEqual(self.client.call('POST', '/api/convert/plan', self.request, headers={'Host': 'evil.example'})[0], 421)

    def test_content_type_size_and_json_are_enforced(self):
        self.assertEqual(self.client.call('POST', '/api/convert/plan', raw=b'{}', content_type='text/plain')[0], 415)
        self.assertEqual(self.client.call('POST', '/api/convert/plan', raw=b'{}', content_type=None)[0], 415)
        self.assertEqual(self.client.call('POST', '/api/convert/plan', raw=b'{' + b'"a":1,' * 20000 + b'"b":2}')[0], 413)
        self.assertEqual(self.client.call('POST', '/api/convert/plan', raw=b'not json')[0], 400)
        self.assertEqual(self.client.call('POST', '/api/convert/plan', raw=b'[1, 2]')[0], 400)

    def test_only_the_listed_post_routes_exist(self):
        for path in ('/api/config', '/api/curation', '/api/jobs', '/api/convert/schema', '/api/imas/data', '/', '/api/convert/nope', '/api/jobs/zz/cancel'):
            self.assertEqual(self.client.call('POST', path, {})[0], 405, path)
        for method in ('PUT', 'DELETE', 'PATCH'):
            self.assertEqual(self.client.call(method, '/api/convert/start', {})[0], 405, method)

    def test_plan_errors_are_a_400_with_every_field_error(self):
        status, body = self.client.call('POST', '/api/convert/plan', {'cdf': str(self.cdf), 'options': {'dd_version': 'x', 'run': -1}})
        self.assertEqual(status, 400)
        self.assertEqual({e['field'] for e in body['errors']}, {'dd_version', 'run'})
        self.assertEqual(self.client.call('POST', '/api/convert/plan', {'cdf': '/etc/passwd'})[0], 400)

    def test_start_needs_confirmation_for_an_existing_entry_and_then_runs(self):
        self.wait_idle()
        entry = self.config.scratch_database / 'SYNTHETIC' / '4' / '1' / '2610'
        entry.mkdir(parents=True, exist_ok=True)
        (entry / 'master.h5').write_text('x')
        status, body = self.client.call('POST', '/api/convert/start', self.request)
        self.assertEqual(status, 409)
        self.assertEqual(body['confirm_phrase'], str(entry))
        self.assertEqual(self.client.call('POST', '/api/convert/start', {**self.request, 'confirm_overwrite': 'nope'})[0], 409)
        self.assertFalse(self.srv.jobs._processes)                        # nothing started without the exact phrase
        os.environ['STUB_MODE'] = 'ok'
        status, body = self.client.call('POST', '/api/convert/start', {**self.request, 'confirm_overwrite': str(entry)})
        self.assertEqual(status, 202, body)
        final = self.srv.jobs.wait(body['job']['id'], 30)
        self.assertEqual(final['state'], 'succeeded')
        entry.joinpath('master.h5').unlink()
        entry.rmdir()

    def test_a_job_runs_end_to_end_through_the_api_with_log_and_report(self):
        self.wait_idle()
        os.environ['STUB_MODE'] = 'ok'
        status, body = self.client.call('POST', '/api/convert/start', {'cdf': str(self.cdf), 'options': {'dd_version': '3.39.0'}})
        self.assertEqual(status, 202, body)
        job = body['job']
        self.assertEqual(job['plan']['argv'][-1], '-new')
        self.srv.jobs.wait(job['id'], 30)
        status, record = self.client.call('GET', f"/api/jobs/{job['id']}")
        self.assertEqual((status, record['job']['state'], record['job']['report_written']), (200, 'succeeded', True))
        status, chunk = self.client.call('GET', f"/api/jobs/{job['id']}/log?offset=0")
        self.assertIn('stub converter', chunk['text'])
        status, listed = self.client.call('GET', '/api/jobs')
        self.assertIn(job['id'], [j['id'] for j in listed['jobs']])
        status, report = self.client.call('GET', f"/api/conversion?path={record['job']['report']}")      # the viewer opens the job's report
        self.assertEqual((status, report['view']['header']['status']), (200, 'PASS'))

    def test_a_second_start_while_one_runs_is_a_409_and_cancel_works(self):
        self.wait_idle()
        os.environ['STUB_MODE'] = 'sleep'
        status, body = self.client.call('POST', '/api/convert/start', {'cdf': str(self.cdf), 'options': {'dd_version': '3.39.0', 'run': 41}})
        self.assertEqual(status, 202, body)
        job_id = body['job']['id']
        status, busy = self.client.call('POST', '/api/convert/start', {'cdf': str(self.cdf), 'options': {'dd_version': '3.39.0', 'run': 42}})
        self.assertEqual(status, 409)
        self.assertIn('already running', busy['error'])
        status, cancelled = self.client.call('POST', f'/api/jobs/{job_id}/cancel', {})
        self.assertEqual((status, cancelled['cancelled']), (200, True))
        self.assertEqual(self.srv.jobs.wait(job_id, 30)['state'], 'cancelled')
        self.assertEqual(self.client.call('POST', '/api/jobs/ffffffff/cancel', {})[0], 404)

    def test_job_endpoints_reject_bad_ids(self):
        for path in ('/api/jobs/../../etc', '/api/jobs/zzzz', '/api/jobs/ffff/log?offset=0', '/api/jobs/%2e%2e'):
            status, _ = self.client.call('GET', path)
            self.assertIn(status, (404, 403), path)

    def test_imas_endpoints_refuse_entries_outside_the_database_roots(self):
        outside = Path(self.tmp.name) / 'inputs'
        for entry in (str(outside), '/etc', '', '../../..'):
            status, body = self.client.call('GET', f'/api/imas/describe?entry={entry}')
            self.assertEqual(status, 403, entry)

    def test_an_entry_directory_must_hold_master_h5(self):
        entry = self.config.scratch_database / 'X' / '4' / '1' / 'empty'
        entry.mkdir(parents=True)
        status, body = self.client.call('GET', f'/api/imas/describe?entry={entry}')
        self.assertEqual(status, 403)
        self.assertIn('master.h5', body['error'])

    def test_the_config_explains_why_running_is_off(self):
        config = convert_plan.Config([Path(self.tmp.name) / 'inputs'], [Path(self.tmp.name) / 'out'], None, True, False, imas_note='/x/python cannot import imas')
        srv = server.make_server([Path(self.tmp.name) / 'inputs'], config=config, repo_root=ROOT)
        threading.Thread(target=srv.serve_forever, daemon=True).start()
        try:
            body = Client(srv).call('GET', '/api/config')[1]
            self.assertEqual((body['run_enabled'], body['run_requested'], body['imas_note']), (False, True, '/x/python cannot import imas'))
        finally:
            srv.shutdown()
            srv.server_close()

    def test_a_server_without_an_imas_interpreter_says_why_plotting_is_off(self):
        config = convert_plan.Config([Path(self.tmp.name) / 'inputs'], [Path(self.tmp.name) / 'out'], None, True, False)
        srv = server.make_server([Path(self.tmp.name) / 'inputs'], config=config, repo_root=ROOT)
        threading.Thread(target=srv.serve_forever, daemon=True).start()
        try:
            client = Client(srv)
            self.assertFalse(client.call('GET', '/api/config')[1]['run_enabled'])
            status, body = client.call('POST', '/api/convert/start', self.request)
            self.assertEqual(status, 409)
            self.assertIn('IMAS interpreter', body['error'])
            entry = config.scratch_database / 'X' / '4' / '1' / '1'
            entry.mkdir(parents=True)
            (entry / 'master.h5').write_text('x')
            status, body = client.call('GET', f'/api/imas/describe?entry={entry}')
            self.assertEqual(status, 502)
            self.assertIn('IMAS interpreter', body['error'])
        finally:
            srv.shutdown()
            srv.server_close()

    def test_schema_inputs_and_entries_endpoints(self):
        status, schema = self.client.call('GET', '/api/convert/schema')
        self.assertEqual((status, schema['default_backend']), (200, 'HDF5'))
        status, inputs = self.client.call('GET', '/api/inputs')
        self.assertEqual([Path(c['path']).name for c in inputs['cdfs']], ['run1Z10.CDF'])
        self.assertEqual(self.client.call('GET', '/api/entries')[0], 200)


if __name__ == '__main__':
    unittest.main(verbosity=1)
