#!/usr/bin/python
"""Resolved-input evidence: what the run actually used, compared with what TR.DAT declares.

``resolved_inputs`` is supplied by the analyst (from the TRANSP log, the final namelist dump or the run record).
It is evidence about the run, not something this package derives.  A control marked ``resolved: false`` is a
declaration, not a runtime value, and is reported as such.  Nothing is repaired: a differing declaration is
reported next to the resolved value.
"""

import hashlib
from pathlib import Path


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def segment_value(control, time):
    """The value of a resolved control at ``time`` (the start of a declared interval), or ``(None, None)``.

    Segments are half-open ``[start, end)``; ``None`` bounds are unbounded.  A ``None`` time (an interval of
    unknown start) is covered only by a segment that is itself unbounded below.
    """
    for number, segment in enumerate(control['segments']):
        lo, hi = segment.get('start'), segment.get('end')
        if lo is not None and (time is None or time < lo):
            continue
        if hi is not None and time is not None and time >= hi:
            continue
        return segment['value'], number
    return None, None


def same_value(declared, resolved):
    """Compare a declared TR.DAT scalar with a resolved one without coercing text to numbers or booleans."""
    if isinstance(declared, bool) or isinstance(resolved, bool):
        return isinstance(declared, bool) and isinstance(resolved, bool) and declared == resolved
    if isinstance(declared, (int, float)) and isinstance(resolved, (int, float)):
        return abs(declared - resolved) <= 1e-12 * max(1.0, abs(resolved))
    return declared == resolved


def check_files(inputs):
    """``[{'role', 'path', 'status', ...}]`` where status is ``verified``, ``mismatch``, ``unreadable`` or ``unhashed``."""
    rows = []
    for item in inputs:
        row = {'role': item['role'], 'path': item['path'], 'declared_sha256': item.get('sha256')}
        path = Path(item['path'])
        if 'sha256' not in item:
            row['status'] = 'unhashed'
        else:
            try:
                row['actual_sha256'] = sha256_of(path)
                row['status'] = 'verified' if row['actual_sha256'].lower() == item['sha256'].lower() else 'mismatch'
            except OSError as exc:
                row['status'] = 'unreadable'
                row['error'] = str(exc)
        rows.append(row)
    return rows
