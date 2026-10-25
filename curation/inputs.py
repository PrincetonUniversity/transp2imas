#!/usr/bin/python
"""Conservative declared-input timeline, not a replacement for TRANSP preprocessing.

Only explicit scalar literals are interpreted. Expressions, includes, inline
multiple assignments and invalid update times create parser gaps, never defaults.
An omitted switch is unknown; later explicit declarations persist in snapshots.
"""

from pathlib import Path
import re

from input_provenance import strip_fortran_comment
from curation.source import digest_file

NUMBER = re.compile(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[dDeE][+-]?\d+)?\Z')
KEY = re.compile(r'^([A-Za-z][A-Za-z0-9_%]*(?:\(\s*\d+\s*\))?)\s*=\s*(.*)$')


def normalize_name(name):
    """TR.DAT keys and reference names compare without whitespace, upper-cased (``EINJA ( 1 )`` is ``EINJA(1)``)."""
    return re.sub(r'\s+', '', name).upper()


def literal(raw):
    value = raw.strip().rstrip(',').strip()
    logical = value.upper()
    if logical in ('.TRUE.', '.T.', 'T', 'TRUE'):
        return True
    if logical in ('.FALSE.', '.F.', 'F', 'FALSE'):
        return False
    if NUMBER.fullmatch(value):
        result = float(value.replace('d', 'e').replace('D', 'E'))
        import math
        return result if math.isfinite(result) else None
    if len(value) >= 2 and value[0] in "\"'" and value[-1] == value[0]:
        return value[1:-1].replace(value[0] * 2, value[0])
    return None


def declared_timeline(path=None):
    if path is None:
        return {'available': False, 'states': [], 'gaps': [], 'source': None,
                'semantics': 'explicit declarations only; no runtime/default resolution'}
    path = Path(path)
    if path.stat().st_size > 10 * 1024 * 1024:
        raise ValueError('TR.DAT exceeds the curation 10 MiB input limit')
    body = path.read_bytes().decode('utf-8', errors='replace')
    states, gaps, values, locations = [], [], {}, {}
    start = None
    temporal_valid = True
    for line_number, line in enumerate(body.splitlines(), 1):
        clean = strip_fortran_comment(line).strip()
        if clean.startswith('&'):
            clean = re.sub(r'^&\w+\s*', '', clean)     # a group header may share its line with an assignment
        if not clean or clean == '/':
            continue
        update = re.match(r'^~\s*update_time\s*=\s*(.*)$', clean, re.I)
        if update:
            time = literal(update[1])
            if isinstance(time, bool) or not isinstance(time, (int, float)) or (start is not None and time <= start):
                gaps.append({'line': line_number, 'reason': 'unparsed or nonincreasing update time'})
                temporal_valid = False
                continue
            states.append({'start_time': start, 'end_time': time, 'values': dict(values),
                           'source_lines': dict(locations)})
            start = time
            continue
        match = KEY.match(clean)
        if not match:
            gaps.append({'line': line_number, 'reason': 'unsupported statement; not interpreted'})
            continue
        name, raw = match.groups()
        name = normalize_name(name)
        value = literal(raw)
        if value is None:
            values.pop(name, None)  # do not reuse a superseded, now-unknown value
            locations.pop(name, None)
            gaps.append({'line': line_number, 'name': name, 'reason': 'not a single scalar literal'})
        else:
            values[name] = value
            locations[name] = line_number
    states.append({'start_time': start, 'end_time': None, 'values': dict(values),
                   'source_lines': dict(locations)})
    return {'available': True, 'states': states, 'gaps': gaps, 'temporal_valid': temporal_valid,
            'source': {'path': str(path.resolve()), 'sha256': digest_file(path)},
            'semantics': 'explicit declarations only; no runtime/default resolution'}
