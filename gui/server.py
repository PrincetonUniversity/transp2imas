#!/usr/bin/python
"""The standard-library HTTP server behind the viewer and the conversion interface.

Security posture (the viewer reads files, and with ``--enable-run`` starts conversions, on a possibly shared machine):

* it binds to the loopback interface unless ``allow_remote`` is set, and a request must carry the one-time token
  (query string on the first visit, then an HttpOnly SameSite cookie, or an ``X-GUI-Token`` header);
* the ``Host`` header must name the address the server listens on, which defeats DNS-rebinding from a web page;
* GET and HEAD read; the only POST endpoints are ``/api/convert/plan`` and ``/api/curate/plan`` (dry runs),
  ``/api/convert/start``, ``/api/curate/start`` and ``/api/jobs/<id>/cancel``.  A POST needs the token in the ``X-GUI-Token`` *header* (not the cookie, so another site
  cannot make the browser send it), ``Content-Type: application/json``, a body of at most 64 KiB and, when an ``Origin``
  header is present, the viewer's own origin.  Starting or cancelling a *conversion* also needs ``--enable-run``; curation
  does not (it needs no IMAS, reads its inputs and writes only a report into its job directory).  The only programs ever
  started are this repository's ``convert.py`` and ``curate.py`` from a validated argument list, never a shell string;
* every file goes through ``browse.resolve`` (roots, ``.json``/``.npz``, no escaping symlink, size cap);
* a strict Content-Security-Policy and ``nosniff`` are sent, and the page builds its DOM with ``textContent`` only, so
  analyst-written text in a report cannot run as script.
"""

import hmac
import re
import http.cookies
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import secrets
from urllib.parse import parse_qs, urlsplit

from gui import __version__, browse, convert_plan, curate_plan, imas_client, jobs as jobs_module, masks, reports

STATIC_DIR = Path(__file__).resolve().parent / 'static'
STATIC_FILES = {'index.html': 'text/html; charset=utf-8', 'style.css': 'text/css; charset=utf-8',
                **{name: 'text/javascript; charset=utf-8' for name in ('core.js', 'charts.js', 'app.js', 'convert.js', 'curate.js', 'explore.js')}}
CSP = ("default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; connect-src 'self'; "
       "object-src 'none'; base-uri 'none'; frame-ancestors 'none'; form-action 'none'")
COOKIE = 'gui_token'
MAX_BODY = 64 * 1024
POST_ROUTES = ('convert/plan', 'convert/start', 'curate/plan', 'curate/start')
LOOPBACK = ('127.0.0.1', 'localhost', '::1')


def sibling_npz(document, json_path):
    """The mask store named by a curation report, required to sit in the same directory as the report."""
    name = Path((document.get('mask_storage') or {}).get('path') or 'curation.npz').name
    return Path(json_path).parent / name


class Handler(BaseHTTPRequestHandler):
    server_version = 'transp2imas-gui/' + __version__
    protocol_version = 'HTTP/1.0'

    def log_message(self, format, *args):                       # quiet by default; the console shows the URL only
        if getattr(self.server, 'verbose', False):
            super().log_message(format, *args)

    # -- helpers
    def _send(self, status, body, content_type, extra=()):
        data = body if isinstance(body, bytes) else body.encode('utf-8')
        self.send_response(status)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(len(data)))
        self.send_header('Content-Security-Policy', CSP)
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.send_header('Referrer-Policy', 'no-referrer')
        self.send_header('Cache-Control', 'no-store')
        for key, value in extra:
            self.send_header(key, value)
        self.end_headers()
        if self.command != 'HEAD':
            self.wfile.write(data)

    def _json(self, status, document):
        self._send(status, json.dumps(document, allow_nan=False, default=str), 'application/json; charset=utf-8')

    def _error(self, status, message):
        self._json(status, {'error': message})

    def _host_ok(self):
        host = (self.headers.get('Host') or '').strip().lower()
        return host in self.server.allowed_hosts

    def _token_ok(self, query):
        supplied = [self.headers.get('X-GUI-Token')] + query.get('token', [])
        cookie = http.cookies.SimpleCookie(self.headers.get('Cookie') or '')
        if COOKIE in cookie:
            supplied.append(cookie[COOKIE].value)
        return any(value and hmac.compare_digest(value, self.server.token) for value in supplied)

    # -- dispatch
    def do_GET(self):
        url = urlsplit(self.path)
        query = parse_qs(url.query, keep_blank_values=True)
        if not self._host_ok():
            return self._error(421, 'unexpected Host header')
        if not self._token_ok(query):
            return self._send(401, 'unauthorized: open the URL printed by the server (it carries the token)\n',
                              'text/plain; charset=utf-8')
        route = url.path
        if route == '/' or route == '/index.html':
            extra = [('Set-Cookie', f'{COOKIE}={self.server.token}; HttpOnly; SameSite=Strict; Path=/')]
            return self._static('index.html', extra)
        if route.startswith('/static/'):
            return self._static(route[len('/static/'):])
        if route.startswith('/api/'):
            return self._api(route[len('/api/'):], query)
        return self._error(404, 'not found')

    do_HEAD = do_GET

    def _method_not_allowed(self):
        self._send(405, 'only GET and HEAD are served here (the guarded POST endpoints are listed in gui/server.py)\n',
                   'text/plain; charset=utf-8', [('Allow', 'GET, HEAD')])

    do_PUT = do_DELETE = do_PATCH = do_OPTIONS = _method_not_allowed

    def do_POST(self):
        url = urlsplit(self.path)
        route = url.path[len('/api/'):] if url.path.startswith('/api/') else ''
        cancel = re.fullmatch(r'jobs/([0-9a-f]{1,32})/cancel', route)
        if not (route in POST_ROUTES or cancel):
            return self._method_not_allowed()
        if not self._host_ok():
            return self._error(421, 'unexpected Host header')
        supplied = self.headers.get('X-GUI-Token') or ''
        if not hmac.compare_digest(supplied, self.server.token):
            return self._error(401, 'a POST needs the token in the X-GUI-Token header')
        origin = self.headers.get('Origin')
        if origin is not None and origin != f'http://{self.headers.get("Host")}':
            return self._error(403, 'cross-origin request refused')
        if not (self.headers.get('Content-Type') or '').lower().startswith('application/json'):
            return self._error(415, 'send application/json')
        try:
            length = int(self.headers.get('Content-Length') or 0)
        except ValueError:
            return self._error(400, 'bad Content-Length')
        if length > MAX_BODY:
            return self._error(413, f'request body over {MAX_BODY} bytes')
        try:
            body = json.loads(self.rfile.read(length) or b'{}')
        except (ValueError, UnicodeDecodeError):
            return self._error(400, 'the body is not valid JSON')
        if not isinstance(body, dict):
            return self._error(400, 'the body must be a JSON object')
        if cancel:
            return self._cancel(cancel.group(1))
        if route.startswith('curate/'):
            return self._curate(body, start=route == 'curate/start')
        if route == 'convert/plan':
            return self._plan(body, start=False)
        return self._plan(body, start=True)

    # -- conversion
    def _plan(self, body, start):
        server = self.server
        if start and not server.config.enable_run:
            return self._error(403, 'running conversions is disabled: start the viewer with --enable-run')
        if start and not server.config.imas_python:
            return self._error(409, 'no IMAS interpreter is configured (--imas-python); conversion needs IMAS')
        job_id = server.jobs.reserve()[0] if server.jobs else 'planned'
        report = server.jobs._dir(job_id) / 'conversion-report.json' if server.jobs and start else '<job directory>/conversion-report.json'
        try:
            plan = convert_plan.validate(body, server.config, report_json=str(report))
        except convert_plan.PlanError as exc:
            return self._json(400, {'error': str(exc), 'errors': exc.errors})
        if not start:
            return self._json(200, {'plan': plan, 'run_enabled': server.config.enable_run})
        if not convert_plan.confirm(plan, body.get('confirm_overwrite')):
            return self._json(409, {'error': 'this would replace or extend an existing entry: type its path to confirm',
                                    'confirm_phrase': plan['confirm_phrase'], 'plan': plan})
        try:
            job = server.jobs.start(job_id, plan)
        except jobs_module.JobBusy as exc:
            return self._error(409, str(exc))
        return self._json(202, {'job': job})

    # -- curation: it reads the inputs and writes only into its job directory, so it is not behind --enable-run
    def _curate(self, body, start):
        server = self.server
        job_id = server.jobs.reserve()[0]
        output = server.jobs._dir(job_id) / 'curation' if start else '<job directory>/curation'
        try:
            plan = curate_plan.validate(body, server.config, output_dir=str(output))
        except curate_plan.PlanError as exc:
            return self._json(400, {'error': str(exc), 'errors': exc.errors})
        if not start:
            return self._json(200, {'plan': plan})
        try:
            job = server.jobs.start(job_id, plan, kind='curate')
        except jobs_module.JobBusy as exc:
            return self._error(409, str(exc))
        return self._json(202, {'job': job})

    def _cancel(self, job_id):
        if not self.server.jobs:
            return self._error(403, 'no jobs are available')
        try:
            job = self.server.jobs.get(job_id)
        except KeyError:
            # with conversion off, an unknown id is refused like any cancel, so the answer does not reveal which jobs exist
            return self._error(404 if self.server.config.enable_run else 403, 'no such job' if self.server.config.enable_run else 'running conversions is disabled')
        if job.get('kind') != 'curate' and not self.server.config.enable_run:
            return self._error(403, 'running conversions is disabled')
        return self._json(200, {'cancelled': self.server.jobs.cancel(job_id)})

    def _static(self, name, extra=()):
        if name not in STATIC_FILES:
            return self._error(404, 'not found')
        self._send(200, (STATIC_DIR / name).read_bytes(), STATIC_FILES[name], extra)

    # -- API
    def _api(self, name, query):
        roots = self.server.roots
        arg = lambda key: (query.get(key) or [''])[0]
        try:
            server = self.server
            if name == 'config':
                config = server.config
                return self._json(200, {'version': __version__, 'roots': [str(r) for r in roots], 'token': server.token,
                                        'run_enabled': bool(config.enable_run and config.imas_python),
                                        'run_requested': config.enable_run, 'imas_python': config.imas_python,
                                        'imas_note': config.imas_note, 'curate_enabled': curate_plan.availability()[0],
                                        'curate_reason': curate_plan.availability()[1],
                                        'input_roots': [str(r) for r in config.input_roots],
                                        'output_roots': [str(r) for r in config.output_roots],
                                        'scratch_database': str(config.scratch_database),
                                        'allow_default_database': config.allow_default_database,
                                        'default_database': str(config.default_database),
                                        'read_database_roots': [str(r) for r in server.entry_roots()]})
            if name == 'convert/schema':
                return self._json(200, convert_plan.option_schema())
            if name == 'curate/schema':
                return self._json(200, curate_plan.option_schema())
            if name == 'curate/inputs':
                available, _reason = curate_plan.availability()
                return self._json(200, {'cdfs': convert_plan.scan_inputs(server.config.input_roots),
                                        'files': curate_plan.scan_json(server.config.input_roots) if available else {}})
            if name == 'inputs':
                return self._json(200, {'cdfs': convert_plan.scan_inputs(server.config.input_roots)})
            if name == 'entries':
                return self._json(200, {'entries': convert_plan.scan_entries(server.entry_roots())})
            if name == 'jobs':
                return self._json(200, {'jobs': [self._job_view(j) for j in (server.jobs.list() if server.jobs else [])]})
            job_match = re.fullmatch(r'jobs/([0-9a-f]{1,32})(/log)?', name)
            if job_match:
                if not server.jobs:
                    return self._error(404, 'no such job')
                if job_match.group(2):
                    return self._json(200, server.jobs.log(job_match.group(1), int((query.get('offset') or ['0'])[0] or 0)))
                return self._json(200, {'job': self._job_view(server.jobs.get(job_match.group(1)))})
            if name.startswith('imas/'):
                return self._imas(name[len('imas/'):], query)
            if name == 'scan':
                return self._json(200, {'reports': browse.scan(roots)})
            if name == 'curation':
                path, document = browse.load_json(roots, arg('path'))
                if not reports.is_curation_report(document):
                    return self._error(400, 'not a curation report')
                npz = sibling_npz(document, path)
                prefixes = masks.mask_prefixes(npz) if npz.is_file() else []
                return self._json(200, {'path': str(path), 'view': reports.curation_view(document, prefixes)})
            if name == 'conversion':
                path, document = browse.load_json(roots, arg('path'))
                if not reports.is_conversion_report(document):
                    return self._error(400, 'not a conversion report')
                return self._json(200, {'path': str(path), 'view': reports.conversion_view(document)})
            if name == 'masks':
                path, document = browse.load_json(roots, arg('path'))
                prefix = arg('prefix')
                finding = next((f for f in document.get('findings', [])
                                if (f.get('scope') or {}).get('mask_prefix') == prefix), None)
                npz = sibling_npz(document, path)
                if finding is None or not npz.is_file():
                    return self._error(404, 'no such mask')
                browse.resolve(roots, str(npz))
                return self._json(200, masks.mask_view(npz, prefix, finding['scope'].get('source_dimensions', [])))
        except browse.AccessError as exc:
            return self._error(403, str(exc))
        except KeyError:
            return self._error(404, 'not found')
        except imas_client.WorkerFailure as exc:
            return self._error(502, str(exc))
        except (ValueError, OSError) as exc:
            return self._error(400, str(exc))
        return self._error(404, 'unknown endpoint')

    def _job_view(self, job):
        """A job record for the page: the command, plan, outcome and where its report is, without the full argv paths twice."""
        keep = ('id', 'kind', 'state', 'started', 'finished', 'returncode', 'display', 'plan', 'report', 'report_written', 'note', 'outcome')
        return {k: job.get(k) for k in keep}

    def _imas(self, name, query):
        server = self.server
        arg = lambda key: (query.get(key) or [''])[0]
        entry = imas_client.resolve_entry(server.entry_roots(), arg('entry'))
        request = {'entry': str(entry), 'dd_version': arg('dd_version') or None}
        if name == 'describe':
            request['action'] = 'describe'
        elif name == 'nodes':
            request.update(action='nodes', ids=arg('ids'), occurrence=int(arg('occurrence') or 0))
        elif name == 'data':
            paths = json.loads(arg('paths') or '[]')
            if not isinstance(paths, list) or not all(isinstance(p, str) for p in paths):
                raise ValueError('paths must be a JSON list of node paths')
            request.update(action='data', recipe=arg('recipe') or None, ids=arg('ids') or None, paths=paths,
                           occurrence=int(arg('occurrence') or 0), time_index=int(arg('time_index') or 0),
                           all_times=arg('all_times') == '1', frames=arg('frames') == '1')
        elif name == 'df':
            if arg('what') not in ('info', 'ep', 'energy', 'pitch', 'zones'):
                raise ValueError('what must be info, ep, energy, pitch or zones')
            request.update(action='df', what=arg('what'), dist=int(arg('dist') or 0), time_index=int(arg('time_index') or 0),
                           occurrence=int(arg('occurrence') or 0), zone=int(arg('zone')) if arg('zone') != '' else None)
            for key in ('e_min', 'e_max', 'p_min', 'p_max'):
                if arg(key) != '':
                    value = float(arg(key))
                    if value != value or value in (float('inf'), float('-inf')):
                        raise ValueError(f'{key} must be a finite number')
                    request[key] = value
        else:
            return self._error(404, 'unknown endpoint')
        return self._json(200, imas_client.call(server.config.imas_python, server.repo_root, request))


class ViewerServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, address, roots, token, allowed_hosts, verbose=False, config=None, jobs=None, repo_root=None,
                 database_roots=()):
        super().__init__(address, Handler)
        self.roots, self.token, self.allowed_hosts, self.verbose = roots, token, allowed_hosts, verbose
        self.config, self.jobs, self.repo_root = config, jobs, repo_root
        self.database_roots = [Path(r) for r in database_roots]

    def entry_roots(self):
        """Database directories whose entries may be listed and plotted (read only)."""
        roots = [self.config.scratch_database, *self.database_roots]
        if self.config.allow_default_database:
            roots.append(self.config.default_database)
        return [r for r in browse.stat_entry_roots(roots)]


def make_server(roots, host='127.0.0.1', port=0, token=None, allow_remote=False, verbose=False, config=None,
                database_roots=(), repo_root=None, convert_script=None, curate_script=None):
    """A configured, not yet serving, server.  Non-loopback hosts need ``allow_remote`` and are logged as such.

    ``config`` (``convert_plan.Config``) says where inputs may come from, where entries and jobs may go, whether
    running is enabled and which interpreter has IMAS; without it the server is a read-only report viewer.
    """
    if host not in LOOPBACK and not allow_remote:
        raise ValueError(f'refusing to listen on {host!r}: the viewer binds to the loopback interface unless '
                         '--allow-remote is given (use an SSH tunnel instead where you can)')
    roots = [Path(os.path.realpath(r)) for r in roots] or [Path.cwd()]
    for root in roots:
        if not root.is_dir():
            raise ValueError(f'root {str(root)!r} is not a directory')
    repo_root = Path(repo_root or Path(__file__).resolve().parent.parent)
    if config is None:
        config = convert_plan.Config(roots, [Path(os.path.realpath(roots[0]))], None, False, False)
    if config.enable_run and config.imas_python:
        config.output_roots[0].mkdir(parents=True, exist_ok=True)
    # a manager always exists: curation is not behind --enable-run, and past jobs are listed either way
    manager = jobs_module.JobManager(config.output_roots[0], config.imas_python, repo_root, convert_script, curate_script)
    report_roots = list(roots) + [r for r in config.output_roots if r not in roots]    # curation reports land in the output root
    token = token or secrets.token_urlsafe(18)
    server = ViewerServer((host, port), report_roots, token, set(), verbose, config, manager, repo_root, database_roots)
    bound = server.server_address[1]
    names = {host, 'localhost', '127.0.0.1'} if host in LOOPBACK else {host}
    server.allowed_hosts = {f'{name}:{bound}'.lower() for name in names}
    return server
