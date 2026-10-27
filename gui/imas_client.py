#!/usr/bin/python
"""Calling the IMAS plot worker as a bounded subprocess, and checking which entries it may be pointed at.

The viewer process may live in an environment without IMAS; reading an entry happens in ``gui.imas_worker`` run by the
configured IMAS interpreter.  The worker only reads (``'r'``).  An entry must be an HDF5 entry directory (it holds
``master.h5``) below one of the readable database roots; the request is JSON on the command line, the answer one
``GUI_RESULT:`` line, bounded in size and time.
"""

import json
import os
from pathlib import Path
import subprocess
import threading

from gui import browse
from gui.imas_worker import PREFIX
from gui.jobs import clean_environment

TIMEOUT_SECONDS = 90
MAX_RESULT_BYTES = 16 * 1024 * 1024
_slots = threading.BoundedSemaphore(2)               # at most two readers at once; a third request waits
_entry_locks = {}
_registry_lock = threading.Lock()


def _entry_lock(entry):
    """One reader at a time per entry: HDF5 refuses to open a file another process is reading, so two requests for the
    same entry must not overlap (different entries, as when comparing, still run side by side)."""
    with _registry_lock:
        return _entry_locks.setdefault(str(entry), threading.Lock())


class WorkerFailure(RuntimeError):
    """The worker could not produce a result (its own message is kept)."""


def resolve_entry(roots, requested):
    """The real directory of an HDF5 entry below ``roots``; ``AccessError`` for anything else."""
    if not requested or '\0' in str(requested):
        raise browse.AccessError('no entry given')
    real = Path(os.path.realpath(requested))
    if not any(real == Path(os.path.realpath(r)) or Path(os.path.realpath(r)) in real.parents for r in roots):
        raise browse.AccessError('entry is outside the readable database roots')
    if not (real / 'master.h5').is_file():
        raise browse.AccessError('not an HDF5 IMAS entry (no master.h5); only HDF5 entries can be plotted')
    return real


def call(python, repo_root, request):
    """Run one worker request; returns the decoded result or raises ``WorkerFailure``."""
    if not python:
        raise WorkerFailure('no IMAS interpreter is configured (start the viewer with --imas-python)')
    with _entry_lock(request.get('entry')), _slots:
        try:
            done = subprocess.run([python, '-m', 'gui.imas_worker', json.dumps(request)], cwd=repo_root, capture_output=True,
                                  env=clean_environment(), timeout=TIMEOUT_SECONDS, stdin=subprocess.DEVNULL)
        except subprocess.TimeoutExpired as exc:
            raise WorkerFailure(f'reading the entry took longer than {TIMEOUT_SECONDS} s') from exc
        except OSError as exc:
            raise WorkerFailure(f'cannot start the IMAS interpreter: {exc}') from exc
    lines = [line for line in done.stdout.decode('utf-8', errors='replace').splitlines() if line.startswith(PREFIX)]
    if not lines:
        tail = done.stderr.decode('utf-8', errors='replace').strip().splitlines()[-3:]
        raise WorkerFailure('the IMAS worker produced no result' + (': ' + ' | '.join(tail) if tail else ''))
    if len(lines[-1]) > MAX_RESULT_BYTES:
        raise WorkerFailure('the result is too large to show')
    result = json.loads(lines[-1][len(PREFIX):])
    if 'error' in result and done.returncode != 0:
        raise WorkerFailure(result['error'])
    return result
