#!/usr/bin/python
"""Running ``convert.py`` or ``curate.py`` as a subprocess, with a live log, cancellation and a persistent record.

The only programs this module can start are ``convert.py`` (in the configured IMAS interpreter) and ``curate.py`` (in the
viewer's own interpreter, which needs no IMAS) of this repository, from an argument *list* produced by ``convert_plan`` or
``curate_plan`` (never a shell string).  One job of each kind runs at a time.  Each job gets a directory
``<work>/jobs/<id>/`` holding ``job.json`` (the command, plan and outcome), ``log.txt`` and the report the program writes;
nothing there is ever deleted by the GUI.  A job found ``running`` after a restart whose process is gone is marked
``lost``, not guessed to have succeeded.

Exit codes mean different things.  A conversion either succeeds (0) or fails.  Curation exits 0 when no finding reaches the
``--fail-on`` threshold, 1 when some does (the screening *completed*; this is the state ``findings``, not a failure) and 2
when it could not complete; exit 0 never means scientific acceptance.
"""

import json
import os
from pathlib import Path
import secrets
import signal
import subprocess
import sys
import threading
import time

STATES = ('queued', 'running', 'succeeded', 'findings', 'failed', 'cancelled', 'lost')
MAX_LOG_CHUNK = 256 * 1024


class JobBusy(RuntimeError):
    """Another job of the same kind is still running."""


def clean_environment(base=None):
    """The child environment: an old Access Layer on ``PYTHONPATH`` breaks IMAS-Python, and output must not be buffered."""
    env = dict(os.environ if base is None else base)
    env.pop('PYTHONPATH', None)
    env['PYTHONUNBUFFERED'] = '1'
    return env


class JobManager:
    def __init__(self, work_dir, python, repo_root, convert_script=None, curate_script=None, curate_python=None):
        self.work_dir = Path(work_dir)
        self.jobs_dir = self.work_dir / 'jobs'
        self.python = python
        self.repo_root = Path(repo_root)
        self.convert_script = Path(convert_script) if convert_script else self.repo_root / 'convert.py'   # tests substitute a stub
        self.curate_script = Path(curate_script) if curate_script else self.repo_root / 'curate.py'
        self.curate_python = curate_python or sys.executable                  # curation needs no IMAS
        self._lock = threading.Lock()
        self._processes = {}
        self._kinds = {}
        self._cancelled = set()

    # -- records
    def _dir(self, job_id):
        if not job_id or not all(c in '0123456789abcdef' for c in job_id) or len(job_id) > 32:
            raise KeyError(job_id)
        return self.jobs_dir / job_id

    def _write(self, job):
        path = self._dir(job['id']) / 'job.json'
        temporary = path.with_suffix('.tmp')
        temporary.write_text(json.dumps(job, indent=2, allow_nan=False, default=str))
        os.replace(temporary, path)

    def get(self, job_id):
        path = self._dir(job_id) / 'job.json'
        if not path.is_file():
            raise KeyError(job_id)
        with self._lock:
            job = json.loads(path.read_text())
            if job['state'] == 'running' and job_id not in self._processes:
                job.update(state='lost', finished=time.time(),
                           note='the server restarted while this job was running; its outcome is unknown')
                self._write(job)
        return job

    def list(self, limit=50):
        jobs = []
        if self.jobs_dir.is_dir():
            for path in sorted(self.jobs_dir.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True)[:limit]:
                try:
                    jobs.append(self.get(path.name))
                except (KeyError, OSError, ValueError):
                    continue
        return jobs

    def reserve(self):
        """A fresh job id and its (not yet created) directory, so the report path can go into the command."""
        job_id = secrets.token_hex(6)
        return job_id, self._dir(job_id)

    def running(self):
        return [j for j in self._processes]

    # -- running
    def start(self, job_id, plan, kind='convert'):
        """Start ``convert.py`` or ``curate.py`` for a validated ``plan``; ``JobBusy`` if a job of that kind is running."""
        if kind not in ('convert', 'curate'):
            raise ValueError(f'unknown job kind {kind!r}')
        with self._lock:
            if any(self._kinds.get(j) == kind for j in self._processes):
                raise JobBusy(f'a {"conversion" if kind == "convert" else "curation"} is already running; wait for it or cancel it')
            directory = self._dir(job_id)
            directory.mkdir(parents=True, exist_ok=False)
            python, script = (self.python, self.convert_script) if kind == 'convert' else (self.curate_python, self.curate_script)
            argv = [python, str(script), *plan['argv']]
            report = plan['argv'][plan['argv'].index('--report-json') + 1] if kind == 'convert' else plan.get('report')
            job = {'id': job_id, 'kind': kind, 'state': 'running', 'started': time.time(), 'finished': None,
                   'returncode': None, 'argv': argv, 'display': plan['display'], 'plan': plan,
                   'log': str(directory / 'log.txt'), 'report': report}
            self._write(job)
            log = open(directory / 'log.txt', 'wb')
            try:
                process = subprocess.Popen(argv, cwd=directory, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                           env=clean_environment(), start_new_session=True)
            except OSError as exc:
                log.close()
                job.update(state='failed', finished=time.time(), note=f'could not start {script.name}: {exc}')
                self._write(job)
                return job
            self._processes[job_id] = process
            self._kinds[job_id] = kind
            threading.Thread(target=self._watch, args=(job_id, process, log), daemon=True).start()
            return job

    @staticmethod
    def final_state(kind, code, cancelled):
        """The state a finished process maps to (curation exit 1 is "findings at the threshold", not a failure)."""
        if cancelled:
            return 'cancelled'
        if code == 0:
            return 'succeeded'
        return 'findings' if kind == 'curate' and code == 1 else 'failed'

    @staticmethod
    def curation_outcome(report_path):
        """Counts and flags read back from the report the screening wrote (``None`` if there is none or it is unreadable)."""
        try:
            document = json.loads(Path(report_path).read_text(encoding='utf-8'))
            summary = document['summary']
            return {'counts': summary['counts'], 'scientifically_certified': bool(summary.get('scientifically_certified')),
                    'policy': (document.get('policy') or {}).get('name')}
        except (OSError, ValueError, KeyError, TypeError):
            return None

    def _watch(self, job_id, process, log):
        code = process.wait()
        log.close()
        with self._lock:
            job = json.loads((self._dir(job_id) / 'job.json').read_text())
            cancelled = job_id in self._cancelled
            self._cancelled.discard(job_id)
            job.update(state=self.final_state(job['kind'], code, cancelled), returncode=code, finished=time.time())
            job['report_written'] = bool(job.get('report')) and Path(job['report']).is_file()
            if job['kind'] == 'curate' and job['report_written']:
                job['outcome'] = self.curation_outcome(job['report'])
            self._write(job)                        # the final record first: ``get`` must never see a finished process as running
            self._processes.pop(job_id, None)
            self._kinds.pop(job_id, None)

    def cancel(self, job_id):
        """Ask the converter's process group to stop (SIGTERM, then SIGKILL after a grace period)."""
        process = self._processes.get(job_id)
        if process is None:
            return False
        self._cancelled.add(job_id)
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            return False

        def kill_later():
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
        threading.Thread(target=kill_later, daemon=True).start()
        return True

    def log(self, job_id, offset=0):
        """A chunk of the log from byte ``offset``: ``{'text', 'next', 'done'}`` (bounded, undecodable bytes replaced)."""
        job = self.get(job_id)
        path = self._dir(job_id) / 'log.txt'
        data = b''
        if path.is_file():
            with open(path, 'rb') as stream:
                stream.seek(max(0, int(offset)))
                data = stream.read(MAX_LOG_CHUNK)
        return {'text': data.decode('utf-8', errors='replace'), 'next': max(0, int(offset)) + len(data),
                'done': job['state'] not in ('running', 'queued') and len(data) < MAX_LOG_CHUNK}

    def wait(self, job_id, timeout=60.0):
        """Block until the job leaves ``running`` (tests and scripts); returns the final record."""
        end = time.time() + timeout
        while time.time() < end:
            job = self.get(job_id)
            if job['state'] != 'running':
                return job
            time.sleep(0.05)
        raise TimeoutError(job_id)
