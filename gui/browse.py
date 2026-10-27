#!/usr/bin/python
"""Finding reports under the allowed roots, and the one place that decides which paths may be read.

Only files below a root, only ``.json`` and ``.npz`` names, no symlink that leaves the roots, and a size cap.  Nothing
here writes.  ``scan`` recognises a curation report or a conversion report by content, so any file name works.
"""

import json
import os
from pathlib import Path

from gui import reports

MAX_FILE_BYTES = 64 * 1024 * 1024
MAX_SCAN_BYTES = 32 * 1024 * 1024
ALLOWED_SUFFIXES = ('.json', '.npz')
MAX_RESULTS = 500


class AccessError(ValueError):
    """The requested path is outside the roots or is not a readable report file."""


REPORT_CAP = object()           # default for ``max_bytes``: the report cap in force when the call is made


def resolve(roots, requested, suffixes=ALLOWED_SUFFIXES, max_bytes=REPORT_CAP):
    """The real path of ``requested`` if it is an allowed file below one of ``roots``, else ``AccessError``.

    ``max_bytes`` bounds files the viewer *reads into memory* (reports); pass ``None`` for inputs that are only handed to the
    converter, such as a TRANSP CDF, which is routinely hundreds of megabytes.
    """
    if not requested or '\0' in requested:
        raise AccessError('no path given')
    path = Path(requested)
    if not path.is_absolute():
        path = Path(roots[0]) / path
    real = Path(os.path.realpath(path))
    if not any(real == Path(os.path.realpath(r)) or Path(os.path.realpath(r)) in real.parents for r in roots):
        raise AccessError('path is outside the allowed roots')
    if real.suffix.lower() not in suffixes:
        raise AccessError('only ' + ' and '.join(suffixes) + ' files can be opened' if suffixes != ALLOWED_SUFFIXES
                          else 'only .json and .npz files can be opened')
    if not real.is_file():
        raise AccessError('not a file')
    size = real.stat().st_size
    if max_bytes is REPORT_CAP:
        max_bytes = MAX_FILE_BYTES
    if max_bytes is not None and size > max_bytes:
        raise AccessError(f'file is {size / 1e6:.0f} MB, over the {max_bytes / 1e6:.0f} MB size limit for report files')
    return real


def load_json(roots, requested):
    path = resolve(roots, requested)
    if path.suffix.lower() != '.json':
        raise AccessError('not a JSON file')
    try:
        return path, json.loads(path.read_text(encoding='utf-8'))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise AccessError(f'not a readable JSON report: {exc}') from exc


def _headline(kind, document):
    if kind == 'curation':
        source = document.get('source') or {}
        counts = (document.get('summary') or {}).get('counts', {})
        return {'source': Path(source.get('path') or 'in-memory').name, 'counts': counts,
                'policy': (document.get('policy') or {}).get('name'), 'generated': document.get('generated_at_utc')}
    return {'source': Path(document.get('source_cdf') or '').name, 'status': document.get('status'),
            'machine': document.get('machine'), 'shot': document.get('shot'), 'run': document.get('run'),
            'backend': document.get('backend')}


def scan(roots, max_depth=4):
    """Reports found under the roots: ``[{'kind', 'path', 'mtime', 'size', 'headline'}]``, newest first."""
    found = []
    for root in roots:
        root = Path(os.path.realpath(root))
        base_depth = len(root.parts)
        for directory, subdirs, files in os.walk(root):
            if len(Path(directory).parts) - base_depth >= max_depth:
                subdirs[:] = []
            subdirs[:] = [d for d in subdirs if not d.startswith('.') and d != '__pycache__']
            for name in sorted(files):
                if not name.lower().endswith('.json'):
                    continue
                path = Path(directory) / name
                try:
                    size = path.stat().st_size
                    if size > MAX_SCAN_BYTES or not path.is_file():
                        continue
                    document = json.loads(path.read_text(encoding='utf-8'))
                except (OSError, UnicodeDecodeError, json.JSONDecodeError):
                    continue
                kind = ('curation' if reports.is_curation_report(document)
                        else 'conversion' if reports.is_conversion_report(document) else None)
                if kind:
                    found.append({'kind': kind, 'path': str(path), 'mtime': path.stat().st_mtime, 'size': size,
                                  'headline': _headline(kind, document)})
                if len(found) >= MAX_RESULTS:
                    break
    return sorted(found, key=lambda item: -item['mtime'])


def stat_entry_roots(roots):
    """Existing, real database root directories from ``roots`` (non-directories are dropped)."""
    return [Path(os.path.realpath(r)) for r in roots if Path(r).is_dir()]
