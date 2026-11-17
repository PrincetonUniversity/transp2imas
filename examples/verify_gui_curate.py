#!/usr/bin/python
"""Tests of the curation interface of the workbench: the option schema against ``curate.py``'s own parser, validation and
the command line, running the REAL ``curate.py`` through the job manager (exit codes mean different things), and the guarded
HTTP endpoints.

    python examples/verify_gui_curate.py

Needs the curation requirements (jsonschema, referencing, f90nml) but no IMAS.  The real runs use the bundled ITER fixture.
"""

import json
from pathlib import Path
import sys
import tempfile
import textwrap
import threading
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import curate
from gui import convert_plan, curate_plan, jobs, server
from verify_gui_convert import Client

EXAMPLES = ROOT / 'examples'
ITER_CDF = 'data/iter/38530P94.CDF'
EVIDENCE = 'curation/evidence-template.json'
PILOT = 'curation/pilot-template.json'

SLEEP_STUB = textwrap.dedent('''\
    #!/usr/bin/python
    import sys, time
    print('stub curation', sys.argv[1:], flush=True)
    time.sleep(60)
''')


def config(tmp, enable_run=False):
    return convert_plan.Config([EXAMPLES], [Path(tmp) / 'out'], None, enable_run, False)


class SchemaDriftTests(unittest.TestCase):
    def test_every_curate_option_is_exposed_or_managed_by_the_gui(self):
        dests = {a.dest for a in curate.build_parser()._actions if a.dest != 'help'}
        decided = set(curate_plan.EXPOSED) | set(curate_plan.MANAGED)
        self.assertEqual(dests, decided, 'curate.py gained or lost an option: decide whether the GUI exposes it')
        self.assertFalse(set(curate_plan.EXPOSED) & set(curate_plan.MANAGED))

    def test_the_schema_is_read_from_the_parser(self):
        schema = curate_plan.option_schema()
        by = {o['dest']: o for o in schema['options']}
        self.assertEqual(by['fail_on']['choices'], ['error', 'fail', 'warning'])
        self.assertEqual(set(schema['fail_on']), set(by['fail_on']['choices']))
        self.assertTrue(by['plugin']['managed_by_gui'] and not by['plugin']['exposed'])
        self.assertTrue(by['validate_only']['exposed'])
        self.assertTrue(schema['available'])

    def test_the_form_fields_are_the_exposed_options_plus_the_two_choices(self):
        self.assertEqual(set(curate_plan.FORM_FIELDS), {'trdat_mode', 'trdat', 'evidence_kind', 'evidence', 'policy', 'fail_on',
                                                         'require_complete', 'validate_only'})

    def test_curate_main_still_uses_the_shared_parser(self):
        self.assertEqual(curate.build_parser().parse_args(['-i', 'x.CDF']).fail_on, 'fail')


class ClassifyTests(unittest.TestCase):
    def test_files_are_recognised_by_their_schema_version(self):
        names = curate_plan.schema_names()
        self.assertEqual(curate_plan.classify({'schema_version': names['policy'][0]}), 'policy')
        self.assertEqual(curate_plan.classify({'schema_version': names['pilot'][0]}), 'pilot')
        for version in names['evidence']:
            self.assertEqual(curate_plan.classify({'schema_version': version}), 'evidence')
        for other in ({}, [], 'text', {'schema_version': 'nope'}, {'schema_version': None}, 3):
            self.assertIsNone(curate_plan.classify(other))

    def test_the_templates_are_found_and_sorted_into_kinds(self):
        found = curate_plan.scan_json([EXAMPLES])
        self.assertIn('evidence-template.json', [f['name'] for f in found['evidence']])
        self.assertIn('pilot-template.json', [f['name'] for f in found['pilot']])
        self.assertEqual([f for f in found['policy'] if 'template' in f['name'] and 'policy' not in f['name']], [])

    def test_scan_skips_unreadable_and_oversized_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / 'broken.json').write_text('{nope')
            (Path(tmp) / 'binary.json').write_bytes(b'\xff\xfe\x00')
            (Path(tmp) / 'big.json').write_text(json.dumps({'schema_version': curate_plan.schema_names()['policy'][0], 'pad': 'x' * (5 * 1024 * 1024)}))
            self.assertEqual(curate_plan.scan_json([tmp]), {'evidence': [], 'pilot': [], 'policy': []})

    def test_the_tr_dat_beside_a_cdf_is_found_in_any_letter_case(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / '77Y01.CDF').write_bytes(b'x')
            self.assertIsNone(curate_plan.trdat_beside(Path(tmp) / '77Y01.CDF'))
            (Path(tmp) / '77Y01tr.dat').write_text('x')
            self.assertEqual(curate_plan.trdat_beside(Path(tmp) / '77Y01.CDF').name, '77Y01tr.dat')


class PlanTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.config = config(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def plan(self, options=None, cdf=ITER_CDF, **extra):
        return curate_plan.validate({'cdf': cdf, 'options': options or {}, **extra}, self.config, output_dir='/out/curation')

    def errors(self, options=None, cdf=ITER_CDF):
        with self.assertRaises(curate_plan.PlanError) as caught:
            self.plan(options, cdf)
        return {e['field']: e['message'] for e in caught.exception.errors}

    def test_the_default_plan_screens_the_cdf_with_the_tr_dat_beside_it(self):
        plan = self.plan()
        self.assertEqual(plan['argv'][:4], ['-i', str((EXAMPLES / ITER_CDF).resolve()), '--trdat', str((EXAMPLES / 'data/iter/38530P94TR.DAT').resolve())])
        self.assertEqual(plan['argv'][-4:], ['--fail-on', 'fail', '--output-dir', '/out/curation'])
        self.assertEqual((plan['report'], plan['validate_only'], plan['fail_on']), ('/out/curation/curation.json', False, 'fail'))
        self.assertTrue(plan['display'].startswith('python curate.py -i '))
        self.assertTrue(any('not_assessed' in w for w in plan['warnings']))
        self.assertTrue(any('not a scientific certification' in w for w in plan['warnings']))

    def test_the_tr_dat_choices(self):
        none = self.plan({'trdat_mode': 'none'})
        self.assertNotIn('--trdat', none['argv'])
        self.assertTrue(any('timeline is reported as unknown' in w for w in none['warnings']))
        chosen = self.plan({'trdat_mode': 'path', 'trdat': 'data/iter/38530P94TR.DAT'})
        self.assertIn('--trdat', chosen['argv'])
        self.assertIn('trdat', self.errors({'trdat_mode': 'path', 'trdat': 'data/iter/38530P94.CDF'}))          # wrong suffix
        self.assertIn('trdat', self.errors({'trdat_mode': 'path'}))
        self.assertIn('trdat_mode', self.errors({'trdat_mode': 'guess'}))

    def test_evidence_pilot_and_policy_must_be_what_they_say(self):
        plan = self.plan({'evidence_kind': 'evidence', 'evidence': EVIDENCE})
        self.assertEqual(plan['argv'][plan['argv'].index('--evidence') + 1], str((EXAMPLES / EVIDENCE).resolve()))
        self.assertEqual(plan['evidence']['kind'], 'evidence')
        self.assertFalse(any('Without' in w or 'No evidence' in w for w in plan['warnings']))
        pilot = self.plan({'evidence_kind': 'pilot', 'evidence': PILOT})
        self.assertIn('--pilot', pilot['argv'])
        self.assertNotIn('--evidence', pilot['argv'])
        self.assertIn('pilot facts', self.errors({'evidence_kind': 'evidence', 'evidence': PILOT})['evidence'])
        self.assertIn('evidence sidecar', self.errors({'evidence_kind': 'pilot', 'evidence': EVIDENCE})['evidence'])
        self.assertIn('not a curation file', self.errors({'evidence_kind': 'evidence', 'evidence': 'pf_active_geometry.synthetic.json'})['evidence'])
        self.assertIn('not a curation policy', self.errors({'policy': EVIDENCE})['policy'])
        self.assertIn('evidence', self.errors({'evidence_kind': 'evidence'}))
        self.assertIn('evidence_kind', self.errors({'evidence_kind': 'both'}))

    def test_paths_outside_the_input_roots_and_bad_files_are_refused(self):
        self.assertIn('outside', self.errors(cdf='/etc/passwd')['cdf'])
        self.assertIn('outside', self.errors(cdf='../../etc/hosts')['cdf'])
        self.assertIn('cdf', self.errors(cdf='data/iter/missing.CDF'))
        self.assertIn('cdf', self.errors(cdf='data/iter/38530P94TR.DAT'))
        self.assertIn('cdf', self.errors(cdf=''))
        self.assertIn('outside', self.errors({'evidence_kind': 'evidence', 'evidence': '/etc/hostname'})['evidence'])
        with tempfile.TemporaryDirectory() as tmp:
            broken = Path(tmp) / 'x.json'
            broken.write_text('{nope')
            self.config.input_roots.append(Path(tmp))
            self.assertIn('not readable JSON', self.errors({'evidence_kind': 'evidence', 'evidence': str(broken)})['evidence'])

    def test_a_file_name_is_one_argument_never_shell_syntax(self):
        with tempfile.TemporaryDirectory() as tmp:
            odd = Path(tmp) / 'a b; touch pwned $(id).json'
            odd.write_text(json.dumps({'schema_version': curate_plan.schema_names()['policy'][0]}))
            self.config.input_roots.append(Path(tmp))
            plan = self.plan({'policy': str(odd)})
            self.assertEqual(plan['argv'][plan['argv'].index('--policy') + 1], str(odd.resolve()))
            self.assertIn("'", plan['display'])                                               # quoted for display only

    def test_validate_only_needs_no_cdf_and_writes_no_report(self):
        plan = curate_plan.validate({'options': {'validate_only': True, 'evidence_kind': 'evidence', 'evidence': EVIDENCE,
                                                'require_complete': True}}, self.config)
        self.assertNotIn('-i', plan['argv'])
        self.assertNotIn('--output-dir', plan['argv'])
        self.assertNotIn('--require-complete', plan['argv'])                                  # meaningless when only validating
        self.assertIn('--validate-only', plan['argv'])
        self.assertIsNone(plan['report'])
        self.assertTrue(any('not run' in w for w in plan['warnings']))
        bare = curate_plan.validate({'options': {'validate_only': True}}, self.config)
        self.assertTrue(any('built-in defaults' in w for w in bare['warnings']))
        self.assertIn('-i', self.plan({'validate_only': True})['argv'])                       # a CDF given while validating is passed on

    def test_options_the_form_does_not_have_are_refused(self):
        for key in ('plugin', 'output_dir', 'pilot', 'input', 'version', 'anything'):
            self.assertIn('options', self.errors({key: 'x'}), key)
        self.assertIn('fail_on', self.errors({'fail_on': 'never'}))
        self.assertIn('require_complete', self.errors({'require_complete': 'yes'}))
        self.assertIn('validate_only', self.errors({'validate_only': 1.5}))

    def test_strictness_options_reach_the_command(self):
        plan = self.plan({'fail_on': 'warning', 'require_complete': True})
        self.assertIn('--require-complete', plan['argv'])
        self.assertEqual(plan['argv'][plan['argv'].index('--fail-on') + 1], 'warning')

    def test_the_request_is_not_modified(self):
        request = {'cdf': ITER_CDF, 'options': {'fail_on': 'error'}}
        before = json.dumps(request)
        curate_plan.validate(request, self.config)
        self.assertEqual(json.dumps(request), before)


class JobTests(unittest.TestCase):
    """The real curate.py, started the way the server starts it."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.config = config(cls.tmp.name)
        cls.manager = jobs.JobManager(cls.config.output_roots[0], None, ROOT)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_job(self, options, cdf=ITER_CDF, timeout=120):
        job_id, directory = self.manager.reserve()
        plan = curate_plan.validate({'cdf': cdf, 'options': options}, self.config, output_dir=str(directory / 'curation'))
        self.manager.start(job_id, plan, kind='curate')
        return self.manager.wait(job_id, timeout), plan

    def test_a_clean_screening_succeeds_and_its_counts_come_from_the_report(self):
        job, plan = self.run_job({})
        self.assertEqual((job['state'], job['returncode'], job['kind']), ('succeeded', 0, 'curate'), self.manager.log(job['id'])['text'][-400:])
        self.assertTrue(job['report_written'])
        report = json.loads(Path(job['report']).read_text())
        self.assertEqual(job['outcome']['counts'], report['summary']['counts'])
        self.assertFalse(job['outcome']['scientifically_certified'])
        self.assertTrue((Path(job['report']).parent / 'curation.npz').is_file())
        self.assertIn('Scientific certification: NOT performed', self.manager.log(job['id'])['text'])

    def test_findings_at_the_threshold_are_a_completed_screening_not_a_failure(self):
        job, _plan = self.run_job({'fail_on': 'warning'})                  # the ITER fixture has a warning
        self.assertEqual((job['state'], job['returncode']), ('findings', 1), self.manager.log(job['id'])['text'][-400:])
        self.assertTrue(job['report_written'])
        self.assertGreater(job['outcome']['counts']['warning'], 0)

    def test_a_screening_that_cannot_complete_is_a_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / 'bad.json'
            bad.write_text(json.dumps({'schema_version': curate_plan.schema_names()['evidence'][-1], 'nonsense': True}))
            self.config.input_roots.append(Path(tmp))
            try:
                job, _plan = self.run_job({'evidence_kind': 'evidence', 'evidence': str(bad)})
            finally:
                self.config.input_roots.pop()
        self.assertEqual((job['state'], job['returncode']), ('failed', 2))
        self.assertFalse(job['report_written'])
        self.assertIn('could not complete', self.manager.log(job['id'])['text'])

    def test_validate_only_runs_without_a_cdf_and_records_no_report(self):
        job_id, directory = self.manager.reserve()
        plan = curate_plan.validate({'options': {'validate_only': True, 'evidence_kind': 'evidence', 'evidence': EVIDENCE}}, self.config)
        self.manager.start(job_id, plan, kind='curate')
        job = self.manager.wait(job_id, 60)
        log = self.manager.log(job_id)['text']
        self.assertEqual(job['state'], 'succeeded', log)
        self.assertFalse(job['report_written'])
        self.assertIsNone(job.get('outcome'))
        self.assertIn('Evidence valid', log)
        self.assertIn('Not checked without a source', log)

    def test_the_input_files_are_not_modified(self):
        import hashlib
        files = [EXAMPLES / ITER_CDF, EXAMPLES / 'data/iter/38530P94TR.DAT', EXAMPLES / EVIDENCE]
        before = [hashlib.sha256(f.read_bytes()).hexdigest() for f in files]
        self.run_job({'evidence_kind': 'evidence', 'evidence': EVIDENCE})
        self.assertEqual(before, [hashlib.sha256(f.read_bytes()).hexdigest() for f in files])

    def test_the_final_state_mapping(self):
        f = jobs.JobManager.final_state
        self.assertEqual([f('curate', c, False) for c in (0, 1, 2, 137)], ['succeeded', 'findings', 'failed', 'failed'])
        self.assertEqual([f('convert', c, False) for c in (0, 1, 2)], ['succeeded', 'failed', 'failed'])
        self.assertEqual(f('curate', 1, True), 'cancelled')
        self.assertIsNone(jobs.JobManager.curation_outcome('/nonexistent/curation.json'))

    def test_one_job_per_kind_and_cancel(self):
        stub = Path(self.tmp.name) / 'curate_stub.py'
        stub.write_text(SLEEP_STUB)
        manager = jobs.JobManager(Path(self.tmp.name) / 'busy', None, ROOT, curate_script=stub)
        plan = {'argv': ['--fail-on', 'fail'], 'display': 'stub', 'report': None}
        first, _ = manager.reserve()
        manager.start(first, plan, kind='curate')
        second, _ = manager.reserve()
        with self.assertRaises(jobs.JobBusy):
            manager.start(second, plan, kind='curate')
        convert_stub = Path(self.tmp.name) / 'convert_stub.py'
        convert_stub.write_text(SLEEP_STUB)
        manager.convert_script = convert_stub
        manager.python = sys.executable
        third, _ = manager.reserve()
        manager.start(third, {'argv': ['--report-json', '/nonexistent/r.json'], 'display': 'stub'}, kind='convert')   # a conversion may run beside it
        self.assertTrue(manager.cancel(first))
        self.assertTrue(manager.cancel(third))
        self.assertEqual(manager.wait(first, 30)['state'], 'cancelled')
        self.assertEqual(manager.wait(third, 30)['state'], 'cancelled')
        with self.assertRaises(ValueError):
            manager.start(manager.reserve()[0], plan, kind='shell')


class ApiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.config = config(cls.tmp.name)                     # no --enable-run and no IMAS interpreter: curation must still work
        cls.stub = Path(cls.tmp.name) / 'curate_stub.py'
        cls.stub.write_text(SLEEP_STUB)
        cls.srv = server.make_server([EXAMPLES], config=cls.config, repo_root=ROOT)
        threading.Thread(target=cls.srv.serve_forever, daemon=True).start()
        cls.client = Client(cls.srv)

    @classmethod
    def tearDownClass(cls):
        cls.srv.shutdown()
        cls.srv.server_close()
        cls.tmp.cleanup()

    def wait_for(self, job_id, timeout=120):
        end = time.time() + timeout
        while time.time() < end:
            status, body = self.client.call('GET', f'/api/jobs/{job_id}')
            if body['job']['state'] != 'running':
                return body['job']
            time.sleep(0.2)
        self.fail('job did not finish')

    def test_curation_is_offered_without_enable_run(self):
        status, config_body = self.client.call('GET', '/api/config')
        self.assertEqual((status, config_body['run_enabled'], config_body['curate_enabled']), (200, False, True))
        status, schema = self.client.call('GET', '/api/curate/schema')
        self.assertEqual((status, schema['available']), (200, True))
        status, inputs = self.client.call('GET', '/api/curate/inputs')
        self.assertEqual(status, 200)
        self.assertIn('38530P94.CDF', [Path(c['path']).name for c in inputs['cdfs']])
        self.assertTrue(inputs['files']['evidence'] and inputs['files']['pilot'])
        self.assertNotIn('38530P94PH.CDF', [Path(c['path']).name for c in inputs['cdfs']])

    def test_plan_start_poll_and_open_the_report(self):
        body = {'cdf': ITER_CDF, 'options': {'evidence_kind': 'evidence', 'evidence': EVIDENCE}}
        status, planned = self.client.call('POST', '/api/curate/plan', body)
        self.assertEqual(status, 200)
        self.assertIn('<job directory>/curation', planned['plan']['display'])
        status, started = self.client.call('POST', '/api/curate/start', body)
        self.assertEqual(status, 202, started)
        job = self.wait_for(started['job']['id'])
        self.assertIn(job['state'], ('succeeded', 'findings'))
        self.assertTrue(job['report_written'])
        self.assertEqual(job['kind'], 'curate')
        status, log = self.client.call('GET', f'/api/jobs/{job["id"]}/log?offset=0')
        self.assertEqual(status, 200)
        self.assertIn('Curation screening:', log['text'])
        status, view = self.client.call('GET', '/api/curation?path=' + job['report'])         # the viewer can open what the job wrote
        self.assertEqual(status, 200, view)
        status, listed = self.client.call('GET', '/api/jobs')
        self.assertIn(job['id'], [j['id'] for j in listed['jobs']])
        self.assertEqual([j['outcome']['counts'] for j in listed['jobs'] if j['id'] == job['id']], [job['outcome']['counts']])
        status, scanned = self.client.call('GET', '/api/scan')
        self.assertIn(job['report'], [r['path'] for r in scanned['reports']])

    def test_invalid_requests_say_which_field(self):
        status, body = self.client.call('POST', '/api/curate/plan', {'cdf': ITER_CDF, 'options': {'evidence_kind': 'pilot', 'evidence': EVIDENCE}})
        self.assertEqual(status, 400)
        self.assertEqual(body['errors'][0]['field'], 'evidence')
        status, body = self.client.call('POST', '/api/curate/start', {'cdf': '/etc/passwd', 'options': {}})
        self.assertEqual(status, 400)
        status, body = self.client.call('POST', '/api/curate/start', {'cdf': ITER_CDF, 'options': {'plugin': 'os'}})
        self.assertEqual(status, 400)                                                           # a plugin would run code: never accepted
        self.assertIn('plugin', body['errors'][0]['message'])

    def test_a_start_that_is_refused_starts_nothing(self):
        before = len(self.client.call('GET', '/api/jobs')[1]['jobs'])
        self.client.call('POST', '/api/curate/start', {'cdf': '/etc/passwd'})
        self.client.call('POST', '/api/curate/start', {'cdf': ITER_CDF}, token=False)
        self.assertEqual(len(self.client.call('GET', '/api/jobs')[1]['jobs']), before)

    def test_the_post_guards_apply_to_curation_too(self):
        body = {'cdf': ITER_CDF}
        self.assertEqual(self.client.call('POST', '/api/curate/start', body, token=False)[0], 401)
        self.assertEqual(self.client.call('POST', '/api/curate/plan', body, headers={'Origin': 'http://evil.example'})[0], 403)
        self.assertEqual(self.client.call('POST', '/api/curate/plan', body, content_type='text/plain')[0], 415)
        self.assertEqual(self.client.call('POST', '/api/curate/plan', raw=b'[1]')[0], 400)
        self.assertEqual(self.client.call('POST', '/api/curate/plan', headers={'Host': 'evil.example'})[0], 421)
        self.assertEqual(self.client.call('POST', '/api/curate/plan', raw=b'x' * 70000)[0], 413)
        self.assertEqual(self.client.call('GET', '/api/curate/schema', token=False)[0], 401)

    def test_conversion_is_still_off_without_enable_run_and_cannot_be_cancelled(self):
        status, body = self.client.call('POST', '/api/convert/start', {'cdf': ITER_CDF})
        self.assertEqual(status, 403)
        status, body = self.client.call('POST', '/api/jobs/' + 'a' * 12 + '/cancel', {})
        self.assertEqual(status, 403)                          # unknown ids are refused like any cancel: no hint which jobs exist

    def test_the_page_serves_the_curate_view(self):
        connection = __import__('http.client').client.HTTPConnection('127.0.0.1', self.srv.server_address[1])
        connection.request('GET', '/static/curate.js', headers={'Host': f'127.0.0.1:{self.srv.server_address[1]}', 'X-GUI-Token': self.srv.token})
        response = connection.getresponse()
        text = response.read().decode()
        self.assertEqual(response.status, 200)
        self.assertIn('views.curate', text)
        self.assertNotIn('innerHTML', text)
        self.assertNotIn('eval(', text)


class CancelAndBusyApiTests(unittest.TestCase):
    """A second curation while one runs is refused; a curation can be cancelled without --enable-run."""

    def test_busy_and_cancel(self):
        with tempfile.TemporaryDirectory() as tmp:
            stub = Path(tmp) / 'curate_stub.py'
            stub.write_text(SLEEP_STUB)
            srv = server.make_server([EXAMPLES], config=config(tmp), repo_root=ROOT, curate_script=stub)
            threading.Thread(target=srv.serve_forever, daemon=True).start()
            try:
                client = Client(srv)
                status, first = client.call('POST', '/api/curate/start', {'cdf': ITER_CDF})
                self.assertEqual(status, 202, first)
                status, second = client.call('POST', '/api/curate/start', {'cdf': ITER_CDF})
                self.assertEqual(status, 409)
                self.assertIn('already running', second['error'])
                status, body = client.call('POST', f'/api/jobs/{first["job"]["id"]}/cancel', {})
                self.assertEqual((status, body['cancelled']), (200, True))
                end = time.time() + 30
                state = 'running'
                while time.time() < end and state == 'running':
                    state = client.call('GET', f'/api/jobs/{first["job"]["id"]}')[1]['job']['state']
                    time.sleep(0.2)
                self.assertEqual(state, 'cancelled')
            finally:
                srv.shutdown()
                srv.server_close()


if __name__ == '__main__':
    unittest.main(verbosity=1)
