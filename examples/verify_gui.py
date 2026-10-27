#!/usr/bin/python
"""Tests of the read-only report viewer (``gui``): view models, masks, path rules, the HTTP server and the static page.

No IMAS and no browser: the server is exercised over real HTTP on an ephemeral loopback port.  A real-browser render
check is separate (``verify_gui_browser.py``) because it needs Playwright.

    python examples/verify_gui.py
"""

import hashlib
import http.client
import json
import os
from pathlib import Path
import re
import sys
import tempfile
import threading
import unittest
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from curation.reporting import write_report
from gui import browse, masks, reports, server
import verify_curation as base
from verify_curation import run_screen, variables

NSTX_CDF = ROOT / 'examples' / 'data' / 'nstx' / '141711Y05' / '141711Y05.CDF'


def write_curation_report(directory, manifest=None, data=None, timeline=None):
    report, _found = run_screen(data, manifest, timeline)
    write_report(report, directory)
    return Path(directory) / 'curation.json'


CONVERSION = {
    'format_version': 1, 'status': 'PASS', 'error': None, 'machine': 'NSTX', 'shot': 142202, 'run': 2601, 'occurrence': 0,
    'backend': 'HDF5', 'dd_version': '4.1.1', 'mapping_policy': 'strict', 'source_cdf': '/data/142202U01.CDF',
    'output_directory': '/tmp/x/NSTX/4/142202/2601', 'selected_ids': ['wall', 'lh_antennas'],
    'written_ids': [{'ids': 'wall', 'occurrence': 0}], 'setup_reads': ['TIME'],
    'publication': {'accepted': [{'ids': 'wall', 'occurrence': 0}], 'failed': None, 'kind': None, 'status': 'complete', 'unattempted': []},
    'writers': [
        {'ids': 'wall', 'occurrences': [0], 'seconds': 0.07, 'status': 'WROTE', 'write_count': 1, 'write_seconds': 0.003,
         'consumed': ['RLIM', 'TIME', 'YLIM'], 'reads': {'total': 3, 'variables': 3, 'most_repeated': {}}},
        {'ids': 'lh_antennas', 'occurrences': [0], 'seconds': 0.01, 'status': 'SKIPPED', 'write_count': 0, 'write_seconds': 0.0,
         'mapping': {'skipped': {'count': 1, 'truncated': False, 'fields': [
             {'count': 1, 'reason': 'no IDS written; the writer logged its reason', 'target': 'lh_antennas'}]}}}]}


class GraphLayoutTests(unittest.TestCase):
    def test_parents_are_left_of_children_and_edges_follow_parentage(self):
        nodes = {'a': {'parents': []}, 'b': {'parents': ['a']}, 'c': {'parents': ['b', 'a']}, 'd': {'parents': []}}
        layout = reports.graph_layout(nodes)
        at = {n['id']: n for n in layout['nodes']}
        self.assertEqual((at['a']['x'], at['b']['x'], at['c']['x'], at['d']['x']), (0, 240, 480, 0))
        self.assertEqual({(e['from'], e['to']) for e in layout['edges']}, {('a', 'b'), ('b', 'c'), ('a', 'c')})
        self.assertEqual(sorted(n['y'] for n in layout['nodes'] if n['x'] == 0), [0, 84])

    def test_a_cycle_cannot_loop_forever(self):
        layout = reports.graph_layout({'a': {'parents': ['b']}, 'b': {'parents': ['a']}})
        self.assertEqual(len(layout['nodes']), 2)

    def test_unknown_parents_are_ignored_and_scope_is_kept(self):
        layout = reports.graph_layout({'a': {'parents': ['ghost'], 'scope': False, 'usage': ['constraint']}})
        self.assertEqual((layout['edges'], layout['nodes'][0]['scope'], layout['nodes'][0]['x']), ([], False, 0))

    def test_an_empty_graph_has_a_size(self):
        self.assertEqual(reports.graph_layout({}), {'nodes': [], 'edges': [], 'width': 240, 'height': 1 * 84})


class CurationViewTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = write_curation_report(self.tmp.name)
        self.document = json.loads(self.path.read_text())
        self.view = reports.curation_view(self.document, masks.mask_prefixes(Path(self.tmp.name) / 'curation.npz'))

    def tearDown(self):
        self.tmp.cleanup()

    def test_every_finding_and_status_count_is_carried_unchanged(self):
        self.assertEqual(len(self.view['findings']), len(self.document['findings']))
        self.assertEqual(self.view['summary']['counts'], self.document['summary']['counts'])
        for shown, raw in zip(self.view['findings'], self.document['findings']):
            self.assertEqual((shown['id'], shown['status'], shown['message']), (raw['id'], raw['status'], raw['message']))

    def test_the_certification_notice_is_always_present_and_nothing_is_certified(self):
        self.assertIn('not a scientific certification', self.view['notices'][0])
        self.assertFalse(self.view['summary']['scientifically_certified'])

    def test_the_ledger_policy_and_fingerprints_come_through(self):
        self.assertEqual(self.view['ledger'], self.document['execution_ledger'])
        self.assertEqual(self.view['header']['fingerprints'], self.document['fingerprints'])
        self.assertEqual(self.view['summary']['rule_coverage'], self.document['summary']['rule_coverage'])

    def test_findings_with_stored_masks_are_marked_and_others_are_not(self):
        flagged = {f['id']: f['has_masks'] for f in self.view['findings']}
        self.assertTrue(flagged['GRID.X'])
        self.assertFalse(flagged['EVIDENCE.rotation'])
        self.assertEqual(set(self.view['mask_prefixes']), set(masks.mask_prefixes(Path(self.tmp.name) / 'curation.npz')))

    def test_an_older_report_gets_a_notice_and_is_not_rewritten(self):
        old = {**self.document, 'schema_version': 'transp2imas.curation.report.v2'}
        for key in ('execution_ledger', 'fingerprints', 'policy', 'registry'):
            old.pop(key)
        view = reports.curation_view(old)
        self.assertTrue(any('report.v2' in n for n in view['notices']))
        self.assertEqual((view['ledger'], view['header']['fingerprints']), ([], None))

    def test_the_assignment_list_is_bounded(self):
        document = {**self.document, 'declared_input_timeline': {
            'available': True, 'assignments': [{'order': i} for i in range(1000)], 'states': [], 'gaps': [], 'unresolved': []}}
        timeline = reports.curation_view(document)['timeline']
        self.assertEqual((timeline['assignments_total'], len(timeline['assignments'])), (1000, reports.MAX_ASSIGNMENTS_SHOWN))

    def test_the_view_is_json_serializable_without_nonfinite_numbers(self):
        json.dumps(self.view, allow_nan=False)

    def test_report_kinds_are_recognized_by_content(self):
        self.assertTrue(reports.is_curation_report(self.document))
        self.assertFalse(reports.is_curation_report(CONVERSION))
        self.assertTrue(reports.is_conversion_report(CONVERSION))
        self.assertFalse(reports.is_conversion_report(self.document))
        self.assertFalse(reports.is_curation_report([]))


class ConversionViewTests(unittest.TestCase):
    def test_header_publication_writers_and_journal(self):
        view = reports.conversion_view(CONVERSION)
        self.assertEqual((view['header']['status'], view['header']['machine'], view['header']['shot']), ('PASS', 'NSTX', 142202))
        self.assertEqual((view['publication']['status'], len(view['publication']['accepted'])), ('complete', 1))
        rows = {r['ids']: r for r in view['writers']}
        self.assertEqual((rows['wall']['consumed'], rows['wall']['variables_read'], rows['lh_antennas']['mapping_entries']), (3, 3, 1))
        self.assertEqual(view['mapping_journal'], [{'ids': 'lh_antennas', 'kind': 'skipped', 'target': 'lh_antennas',
                                                    'reason': 'no IDS written; the writer logged its reason', 'count': 1, 'truncated': False}])
        self.assertIn('not a scientific validation', view['notices'][0])

    def test_a_partial_publication_is_shown_as_it_is(self):
        partial = json.loads(json.dumps(CONVERSION))
        partial['status'] = 'FAIL'
        partial['error'] = 'backend rejected a write'
        partial['publication'] = {'accepted': [{'ids': 'wall', 'occurrence': 0}], 'status': 'partial', 'kind': 'backend_failure',
                                  'failed': {'ids': 'core_profiles', 'error': 'x'}, 'unattempted': [{'ids': 'equilibrium', 'occurrence': 0}]}
        view = reports.conversion_view(partial)
        self.assertEqual((view['publication']['status'], view['publication']['kind'], len(view['publication']['unattempted'])),
                         ('partial', 'backend_failure', 1))
        self.assertEqual(view['header']['error'], 'backend rejected a write')

    def test_a_real_conversion_report_from_a_verification_run_is_readable_when_present(self):
        found = sorted((ROOT / 'build' / 'end-to-end').glob('*-conversion.json'))
        if not found:
            self.skipTest('no end-to-end conversion report in build/end-to-end (run make verify-end-to-end)')
        for path in found:
            view = reports.conversion_view(json.loads(path.read_text()))
            self.assertTrue(view['writers'], path.name)


class MaskTests(unittest.TestCase):
    def test_categories_are_exact_and_flagged_outranks_everything(self):
        evaluated = np.array([[0, 1, 1, 1, 1]], dtype=bool)
        flagged = np.array([[0, 0, 1, 0, 1]], dtype=bool)
        excluded = np.array([[0, 0, 0, 1, 1]], dtype=bool)
        self.assertEqual(masks.categories(evaluated, flagged, excluded).tolist(), [[0, 1, 2, 3, 4]])

    def test_a_single_flagged_cell_survives_block_reduction(self):
        code = np.ones((400, 400), dtype=np.int8)
        code[123, 321] = 2
        grid, steps = masks.reduce_grid(code, limit=100)
        self.assertEqual((grid.shape, steps), ((100, 100), (4, 4)))
        self.assertEqual(int((grid == 2).sum()), 1)
        self.assertEqual(int(grid[123 // 4, 321 // 4]), 2)

    def test_small_grids_are_not_reduced_and_one_dimensional_masks_are_drawn_as_one_row(self):
        grid, steps = masks.reduce_grid(np.array([1, 2, 1]))
        self.assertEqual((grid.tolist(), steps), ([[1, 2, 1]], (1, 1)))
        with self.assertRaises(ValueError):
            masks.reduce_grid(np.zeros((2, 2, 2)))

    def test_excluded_outranks_unflagged_and_unevaluated_is_lowest(self):
        grid, _ = masks.reduce_grid(np.array([[0, 3], [1, 3]], dtype=np.int8), limit=1)
        self.assertEqual(grid.tolist(), [[3]])

    def test_a_stored_mask_reads_back_with_its_axes_and_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = {'schema_version': base.MANIFEST_SCHEMA, 'events': [
                {'id': 'p', 'type': 'pellet', 'start': 0.9, 'end': 1.1, 'rules': ['PROFILE.TE.FINITE']}]}
            write_curation_report(tmp, manifest)
            view = masks.mask_view(Path(tmp) / 'curation.npz', 'PROFILE_TE_FINITE', ['TIME3', 'X'])
        self.assertEqual((view['shape'], view['dimensions']), ([2, 3], ['TIME3', 'X']))
        self.assertEqual(view['counts'], {'not evaluated': 0, 'evaluated': 3, 'flagged': 0, 'excluded': 3, 'flagged and excluded': 0})
        self.assertEqual([a['name'] for a in view['axes']], ['TIME3', 'X'])
        self.assertEqual(view['axes'][0]['values'], [0.0, 1.0])
        self.assertEqual(view['axes'][1]['values'], [0.1, 0.5, 0.9])

    def test_a_missing_prefix_is_a_key_error_and_a_pickled_file_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_curation_report(tmp)
            with self.assertRaises(KeyError):
                masks.mask_view(Path(tmp) / 'curation.npz', 'NO_SUCH', ['TIME3', 'X'])
            bad = Path(tmp) / 'pickled.npz'
            np.savez(bad, BAD__evaluated=np.array([{'a': 1}], dtype=object))
            self.assertEqual(masks.mask_prefixes(bad), ['BAD'])                  # listing names reads no array
            with self.assertRaises(ValueError):
                masks.mask_view(bad, 'BAD', ['TIME3', 'X'])                      # reading one would unpickle: refused


class PathRuleTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name, 'root')
        self.root.mkdir()
        self.outside = Path(self.tmp.name, 'outside')
        self.outside.mkdir()
        (self.root / 'a.json').write_text('{}')
        (self.outside / 'secret.json').write_text('{"x": 1}')

    def tearDown(self):
        self.tmp.cleanup()

    def test_files_inside_a_root_resolve_relative_or_absolute(self):
        self.assertEqual(browse.resolve([self.root], 'a.json'), (self.root / 'a.json').resolve())
        self.assertEqual(browse.resolve([self.root], str(self.root / 'a.json')), (self.root / 'a.json').resolve())

    def test_traversal_absolute_outside_and_escaping_symlinks_are_refused(self):
        os.symlink(self.outside / 'secret.json', self.root / 'link.json')
        os.symlink(self.outside, self.root / 'dirlink')
        for requested in ('../outside/secret.json', str(self.outside / 'secret.json'), 'link.json', 'dirlink/secret.json',
                          '/etc/passwd', 'a.json/../../outside/secret.json'):
            with self.assertRaises(browse.AccessError, msg=requested):
                browse.resolve([self.root], requested)

    def test_only_json_and_npz_files_below_the_size_cap_are_readable(self):
        (self.root / 'notes.txt').write_text('x')
        (self.root / 'big.json').write_text('{}')
        for requested, fragment in (('notes.txt', 'only .json and .npz'), ('', 'no path'), ('missing.json', 'not a file'), ('.', 'only .json')):
            with self.assertRaises(browse.AccessError, msg=requested) as caught:
                browse.resolve([self.root], requested)
            self.assertIn(fragment, str(caught.exception))
        with mock.patch.object(browse, 'MAX_FILE_BYTES', 1):
            with self.assertRaises(browse.AccessError) as caught:
                browse.resolve([self.root], 'big.json')
            self.assertIn('size limit', str(caught.exception))

    def test_a_null_byte_is_refused(self):
        with self.assertRaises(browse.AccessError):
            browse.resolve([self.root], 'a.json\0.txt')

    def test_a_non_json_content_is_an_access_error_not_a_crash(self):
        (self.root / 'bad.json').write_text('not json')
        with self.assertRaises(browse.AccessError):
            browse.load_json([self.root], 'bad.json')

    def test_scan_recognizes_reports_by_content_ignores_junk_and_respects_depth_and_hidden_dirs(self):
        write_curation_report(self.root / 'one')
        (self.root / 'conv.json').write_text(json.dumps(CONVERSION))
        (self.root / 'junk.json').write_text('{"a": 1}')
        (self.root / 'broken.json').write_text('{')
        write_curation_report(self.root / '.hidden')
        deep = self.root
        for part in 'abcdef':
            deep = deep / part
        write_curation_report(deep)
        found = browse.scan([self.root], max_depth=3)
        kinds = {(Path(f['path']).relative_to(self.root.resolve()).as_posix(), f['kind']) for f in found}
        self.assertEqual(kinds, {('one/curation.json', 'curation'), ('conv.json', 'conversion')})


class Client:
    def __init__(self, server_obj, host=None):
        self.port = server_obj.server_address[1]
        self.token = server_obj.token
        self.host = host or f'127.0.0.1:{self.port}'

    def request(self, method, path, headers=None, token=True, query_token=False):
        connection = http.client.HTTPConnection('127.0.0.1', self.port, timeout=10)
        merged = {'Host': self.host, **(headers or {})}
        if token and not query_token:
            merged['X-GUI-Token'] = self.token
        if query_token:
            path += ('&' if '?' in path else '?') + 'token=' + self.token
        connection.request(method, path, headers=merged)
        response = connection.getresponse()
        body = response.read()
        result = (response.status, dict(response.getheaders()), body)
        connection.close()
        return result

    def api(self, name, **params):
        from urllib.parse import urlencode
        status, headers, body = self.request('GET', f'/api/{name}?{urlencode(params)}')
        return status, json.loads(body)


class ServerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.tmp.name, 'root')
        cls.root.mkdir()
        Path(cls.tmp.name, 'outside').mkdir()
        (Path(cls.tmp.name) / 'outside' / 'secret.json').write_text('{"secret": true}')
        manifest = {'schema_version': base.MANIFEST_SCHEMA, 'events': [
            {'id': 'p', 'type': 'pellet', 'start': 0.9, 'end': 1.1, 'rules': ['PROFILE.TE.FINITE']}]}
        cls.curation = write_curation_report(cls.root / 'run', manifest)
        (cls.root / 'conversion.json').write_text(json.dumps(CONVERSION))
        cls.server = server.make_server([cls.root])
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()
        cls.client = Client(cls.server)

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.tmp.cleanup()

    def test_it_listens_on_loopback_only_and_refuses_another_host_without_opt_in(self):
        self.assertEqual(self.server.server_address[0], '127.0.0.1')
        with self.assertRaises(ValueError) as caught:
            server.make_server([self.root], host='0.0.0.0')
        self.assertIn('--allow-remote', str(caught.exception))
        remote = server.make_server([self.root], host='0.0.0.0', allow_remote=True)
        remote.server_close()
        with self.assertRaises(ValueError):
            server.make_server([self.root / 'nope'])

    def test_a_request_without_the_token_is_unauthorized_and_a_wrong_token_too(self):
        status, _h, body = self.client.request('GET', '/', token=False)
        self.assertEqual(status, 401)
        self.assertNotIn(b'<html', body)
        self.assertEqual(self.client.request('GET', '/api/config', token=False)[0], 401)
        self.assertEqual(self.client.request('GET', '/api/config', headers={'X-GUI-Token': 'wrong'}, token=False)[0], 401)
        self.assertEqual(self.client.request('GET', '/?token=wrong', token=False)[0], 401)

    def test_the_query_token_serves_the_page_and_sets_an_httponly_samesite_cookie_that_then_works(self):
        status, headers, body = self.client.request('GET', '/', query_token=True)
        self.assertEqual(status, 200)
        self.assertIn(b'<title>TRANSP2IMAS workbench</title>', body)
        cookie = headers['Set-Cookie']
        self.assertIn('HttpOnly', cookie)
        self.assertIn('SameSite=Strict', cookie)
        value = re.match(r'gui_token=([^;]+)', cookie).group(1)
        status, _h, _b = self.client.request('GET', '/api/config', headers={'Cookie': f'gui_token={value}'}, token=False)
        self.assertEqual(status, 200)

    def test_a_foreign_host_header_is_refused_even_with_the_token(self):
        for host in ('evil.example', f'localhost:{self.client.port + 1}', 'attacker.test:80'):
            self.assertEqual(self.client.request('GET', '/api/config', headers={'Host': host})[0], 421, host)
        self.assertEqual(self.client.request('GET', '/api/config', headers={'Host': f'localhost:{self.client.port}'})[0], 200)

    def test_only_get_and_head_are_served_and_head_has_no_body(self):
        for method in ('POST', 'PUT', 'DELETE', 'PATCH', 'OPTIONS'):
            status, headers, _b = self.client.request(method, '/api/config')
            self.assertEqual((status, headers.get('Allow')), (405, 'GET, HEAD'), method)
        status, headers, body = self.client.request('HEAD', '/')
        self.assertEqual((status, body), (200, b''))
        self.assertGreater(int(headers['Content-Length']), 0)

    def test_security_headers_are_sent_on_pages_and_api_responses(self):
        for path in ('/api/config', '/static/app.js', '/nonexistent'):
            _s, headers, _b = self.client.request('GET', path)
            self.assertIn("default-src 'self'", headers['Content-Security-Policy'], path)
            self.assertIn("script-src 'self'", headers['Content-Security-Policy'], path)
            self.assertIn("frame-ancestors 'none'", headers['Content-Security-Policy'], path)
            self.assertEqual(headers['X-Content-Type-Options'], 'nosniff')
            self.assertEqual(headers['Cache-Control'], 'no-store')

    def test_static_files_are_an_allowlist(self):
        self.assertEqual(self.client.request('GET', '/static/app.js')[0], 200)
        self.assertEqual(self.client.request('GET', '/static/style.css')[0], 200)
        for path in ('/static/../server.py', '/static/%2e%2e/server.py', '/static/__init__.py', '/static/', '/static/../../curate.py'):
            self.assertEqual(self.client.request('GET', path)[0], 404, path)

    def test_scan_lists_both_report_kinds(self):
        status, body = self.client.api('scan')
        self.assertEqual(status, 200)
        self.assertEqual({r['kind'] for r in body['reports']}, {'curation', 'conversion'})

    def test_the_curation_endpoint_returns_the_view_of_the_report_on_disk(self):
        status, body = self.client.api('curation', path=str(self.curation))
        self.assertEqual(status, 200)
        document = json.loads(self.curation.read_text())
        self.assertEqual(len(body['view']['findings']), len(document['findings']))
        self.assertEqual(body['view']['summary']['counts'], document['summary']['counts'])
        self.assertIn('PROFILE_TE_FINITE', body['view']['mask_prefixes'])
        self.assertEqual(self.client.api('curation', path='run/curation.json')[0], 200)        # relative to the first root

    def test_the_conversion_endpoint_and_the_kind_check(self):
        status, body = self.client.api('conversion', path='conversion.json')
        self.assertEqual((status, body['view']['header']['machine']), (200, 'NSTX'))
        self.assertEqual(self.client.api('conversion', path=str(self.curation))[0], 400)
        self.assertEqual(self.client.api('curation', path='conversion.json')[0], 400)

    def test_the_mask_endpoint_returns_a_bounded_grid_and_404_for_unknown_masks(self):
        status, body = self.client.api('masks', path=str(self.curation), prefix='PROFILE_TE_FINITE')
        self.assertEqual(status, 200)
        self.assertEqual(body['counts']['excluded'], 3)
        self.assertLessEqual(len(body['grid']), masks.MAX_CELLS)
        self.assertEqual(self.client.api('masks', path=str(self.curation), prefix='NOPE')[0], 404)
        self.assertEqual(self.client.api('masks', path=str(self.curation), prefix='')[0], 404)

    def test_paths_outside_the_roots_and_traversal_are_refused_with_403(self):
        outside = Path(self.tmp.name, 'outside', 'secret.json')
        for requested in (str(outside), '../outside/secret.json', '/etc/passwd', 'run/../../outside/secret.json'):
            for endpoint in ('curation', 'conversion'):
                status, body = self.client.api(endpoint, path=requested)
                self.assertEqual(status, 403, (endpoint, requested))
                self.assertNotIn('secret', json.dumps(body))
        self.assertEqual(self.client.api('masks', path='../outside/secret.json', prefix='x')[0], 403)

    def test_a_report_whose_mask_store_name_tries_to_leave_its_directory_cannot_read_elsewhere(self):
        crafted = json.loads(self.curation.read_text())
        crafted['mask_storage']['path'] = '../../outside/evil.npz'
        (self.root / 'run' / 'crafted.json').write_text(json.dumps(crafted))
        status, body = self.client.api('masks', path='run/crafted.json', prefix='PROFILE_TE_FINITE')
        self.assertEqual(status, 404)                                    # only the file name is used, beside the report

    def test_unknown_endpoints_and_errors_are_json_without_a_traceback(self):
        status, body = self.client.api('nope')
        self.assertEqual((status, list(body)), (404, ['error']))
        status, body = self.client.api('curation')
        self.assertEqual(status, 403)
        self.assertNotIn('Traceback', json.dumps(body))

    def test_exercising_every_endpoint_writes_nothing(self):
        def snapshot():
            return {str(p): (p.stat().st_mtime_ns, hashlib.sha256(p.read_bytes()).hexdigest())
                    for p in sorted(Path(self.tmp.name).rglob('*')) if p.is_file()}
        before = snapshot()
        for name, params in (('config', {}), ('scan', {}), ('curation', {'path': str(self.curation)}),
                             ('conversion', {'path': 'conversion.json'}),
                             ('masks', {'path': str(self.curation), 'prefix': 'PROFILE_TE_FINITE'})):
            self.assertEqual(self.client.api(name, **params)[0], 200, name)
        for method in ('POST', 'PUT', 'DELETE'):
            self.client.request(method, '/api/curation')
        self.assertEqual(snapshot(), before)


class StaticPageSafetyTests(unittest.TestCase):
    """The page builds its DOM from report text, which analysts write: it must never reach an HTML parser or eval."""

    SCRIPTS = ('core.js', 'charts.js', 'app.js', 'convert.js', 'explore.js')

    def setUp(self):
        static = ROOT / 'gui' / 'static'
        self.js = {name: (static / name).read_text() for name in self.SCRIPTS}
        self.html = (static / 'index.html').read_text()
        self.css = (static / 'style.css').read_text()

    def test_every_served_script_is_linted_here(self):
        served = {name for name in server.STATIC_FILES if name.endswith('.js')}
        self.assertEqual(served, set(self.SCRIPTS))
        self.assertEqual({p.name for p in (ROOT / 'gui' / 'static').glob('*.js')}, set(self.SCRIPTS))

    def test_no_script_uses_an_html_or_code_sink(self):
        for name, content in self.js.items():
            for forbidden in ('innerHTML', 'outerHTML', 'insertAdjacentHTML', 'document.write', 'eval(', 'new Function',
                              "setTimeout('", 'setTimeout("', 'setInterval("', "setInterval('", 'srcdoc', 'javascript:',
                              'location.href =', 'window.open', 'localStorage', 'sessionStorage', 'document.cookie'):
                self.assertNotIn(forbidden, content, f'{name}: {forbidden}')

    def test_no_remote_resource_is_referenced_anywhere(self):
        pages = {**self.js, 'index.html': self.html, 'style.css': self.css}
        for name, content in pages.items():
            for url in re.findall(r'https?://[^\s\'")]+', content):
                self.assertEqual(url, 'http://www.w3.org/2000/svg', (name, url))
        self.assertNotIn('@import', self.css)
        self.assertNotIn('url(', self.css)

    def test_the_html_has_no_inline_script_style_or_event_handler(self):
        self.assertEqual(re.findall(r'<script(?![^>]*\bsrc=)[^>]*>', self.html), [])
        self.assertEqual(re.findall(r'\son[a-z]+\s*=', self.html), [])
        self.assertNotIn('style=', self.html)
        sources = re.findall(r'<script src="/static/([^"]+)"', self.html)
        self.assertEqual(sources, ['core.js', 'charts.js', 'app.js', 'convert.js', 'explore.js'])      # core first: the others use it

    def test_the_views_come_in_the_order_convert_explore_reports_and_the_page_lands_on_convert(self):
        app = self.js['app.js']
        self.assertIn("const VIEW_NAMES = [['convert', 'Convert'], ['explore', 'Explore'], ['reports', 'Reports']];", app)
        self.assertIn(": wanted ? 'reports' : 'convert'", app)

    def test_the_skip_option_is_worded_for_what_it_does_not_with_the_jargon_word(self):
        convert = self.js['convert.js']
        self.assertIn('Skip variables that cannot be mapped', convert)
        self.assertNotIn("'Lenient mapping'", convert)
        shown = re.findall(r"check\('[a-z_]+', '([^']+)'", convert)
        self.assertFalse([label for label in shown if 'enient' in label], shown)

    def test_the_page_says_why_running_is_off_next_to_the_button(self):
        convert = self.js['convert.js']
        self.assertIn('--enable-run', convert)
        self.assertIn('--imas-python', convert)
        self.assertIn("class: 'notice'", convert[convert.index('const reasons'):convert.index('const why')+400])

    def test_the_page_says_it_is_read_only_until_running_is_enabled(self):
        self.assertIn('read-only', self.html)
        self.assertIn("config.run_enabled ? 'conversion enabled' : 'read-only'", self.js['app.js'])

    def test_report_and_form_text_goes_through_textcontent_or_text_nodes(self):
        self.assertIn('textContent', self.js['core.js'])
        self.assertIn('createTextNode', self.js['core.js'])

    def test_a_post_carries_the_token_in_a_header_and_never_in_the_url(self):
        core = self.js['core.js']
        self.assertIn("'X-GUI-Token': state.config.token", core)
        self.assertIn("method: 'POST'", core)
        for name, content in self.js.items():
            self.assertNotIn('?token=', content.replace('/?token=', ''), name)
            self.assertNotIn('token=${', content, name)

    def test_only_the_core_script_issues_requests(self):
        for name, content in self.js.items():
            if name != 'core.js':
                self.assertNotIn('fetch(', content, name)
                self.assertNotIn('XMLHttpRequest', content, name)
                self.assertNotIn('WebSocket', content, name)


if __name__ == '__main__':
    unittest.main(verbosity=1)
