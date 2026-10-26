#!/usr/bin/python
"""Declared-input timeline: stage 3 (state resolution) over the ``curation.namelist`` adapter.

Stages stay separate.  ``curation.namelist`` collects sources and parses statement syntax (f90nml); this module only
applies the ordered updates.  Only explicit declarations are resolved: an omitted switch is unknown, there are no TRANSP
defaults, no declared array bounds and no model execution semantics.  A later explicit declaration persists into later
update blocks; a partial update keeps the earlier known elements it does not touch (element-wise namelist semantics);
a null entry changes nothing; an unreadable or unsupported update invalidates everything it may have changed
(``invalidate``), never reusing a superseded value.

If f90nml cannot be imported the old scalar reader runs instead and the timeline says ``parser.kind ==
'scalar_fallback'`` with its limitations; it is never presented as the full parser.
"""

from pathlib import Path
import re

from input_provenance import strip_fortran_comment
from curation import namelist

NUMBER = re.compile(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[dDeE][+-]?\d+)?\Z')
KEY = re.compile(r'^([A-Za-z][A-Za-z0-9_%]*)\s*(?:\(([^()]*)\))?\s*=\s*(.*)$')
INTEGER_INDEX = re.compile(r'\s*([+-]?\d+)\s*\Z')
SLICE_INDEX = re.compile(r'\s*([+-]?\d+)\s*:\s*([+-]?\d+)\s*\Z')
SEMANTICS = 'explicit declarations only; no runtime/default resolution'
LIMIT_BYTES = namelist.MAX_BYTES
FALLBACK_LIMITATIONS = [
    'f90nml is not importable: only single scalar assignments, one per line, are read',
    'repeat notation (n*v), arrays, array sections and multidimensional elements are not resolved',
    'derived-type paths are read as flat names; null entries and several values per line are not interpreted',
    'dotted logicals such as .T are not read',
    'referenced files (PT templates) are not followed or hashed, and includes are not recorded',
]


def normalize_name(name):
    """TR.DAT keys and reference names compare without whitespace, upper-cased (``EINJA ( 1 )`` is ``EINJA(1)``)."""
    return re.sub(r'\s+', '', name).upper()


def literal(raw):
    """One supported scalar (``bool``, number as ``float``, ``str``) or ``None``; the contract is ``fortran_literals``."""
    return namelist.scalar_value(raw, dotted=False)


def invalidate(values, locations, base, first=None, last=None):
    """Drop every stored element of array ``base`` that an unresolved update may have changed; return their keys.

    ``first``/``last`` give an established inclusive index range; without one the whole array is unresolved.  An
    element outside an established range is a known, unaffected value and is kept; a multidimensional element is
    always dropped, since the range is only a one-subscript statement.  The whole-name key ``base`` itself is removed
    too, since it no longer describes the array.
    """
    removed = []
    for key in list(values):
        name, _, rest = key.partition('(')
        if name != base:
            continue
        if rest and first is not None:
            indices = rest.rstrip(')').split(',')
            if len(indices) == 1 and not first <= int(indices[0]) <= last:
                continue
        removed.append(key)
        del values[key]
        locations.pop(key, None)
    return sorted(removed)


_MISSING = object()


def _drop(values, locations, key):
    locations.pop(key, None)
    return [key] if values.pop(key, _MISSING) is not _MISSING else []


def _invalidate_region(values, locations, root, region):
    """Drop what an unresolved statement may have changed (see ``namelist.affected_region``); the dropped keys."""
    removed = []
    if region[0] == 'keys':
        for key in region[1] + [root]:
            removed += _drop(values, locations, key)
    elif region[0] == 'from':
        removed = invalidate(values, locations, root, region[1], 10 ** 9)
    else:
        removed = invalidate(values, locations, root)
    return sorted(set(removed))


def _apply(item, values, locations):
    """Apply one parsed assignment to the running state; return ``(removed keys, gap reason or None)``."""
    parsed, line = item['parsed'], item['line']
    root = parsed['root']
    removed = []
    if parsed['status'] == 'unresolved':
        removed = _invalidate_region(values, locations, root, parsed['affected'])
        return removed, parsed['problems'][0]['reason']
    reason = None
    if parsed['scalar_name']:
        element = parsed['elements'][0]
        if element['state'] == 'null':
            return [], None
        removed = invalidate(values, locations, root)    # a bare scalar for a name that holds elements is ambiguous
        if element['state'] == 'value':
            values[root] = element['value']
            locations[root] = line
            return [key for key in removed if key != root], None
        return removed, 'not a single scalar literal'
    if parsed['subs'] is None:
        removed += _drop(values, locations, root)          # NAME=a,b is NAME(1),NAME(2); the scalar NAME is superseded
    for element in parsed['elements']:
        key = element['key']
        if element['state'] == 'value':
            values[key] = element['value']
            locations[key] = line
        elif element['state'] == 'unknown':
            removed += _drop(values, locations, key)          # do not reuse a superseded, now-unknown value
            reason = 'not a single scalar literal'
    if reason and parsed['section']:                # an unreadable value: the token count is not trusted, so
        for key in parsed['section'][len(parsed['elements']):]:     # the rest of the section is unknown too
            removed += _drop(values, locations, key)
    return sorted(set(removed)), reason


def _finish(states, gaps, temporal_valid, path, sha, extra):
    return {'available': True, 'states': states, 'gaps': gaps, 'temporal_valid': temporal_valid,
            'source': {'path': str(Path(path).resolve()), 'sha256': sha}, 'semantics': SEMANTICS, **extra}


def _resolve(items, path):
    """Stage 3: ordered updates over the update blocks; returns ``(states, gaps, temporal_valid, assignments)``."""
    states, gaps, assignments = [], [], []
    values, locations = {}, {}
    start, temporal_valid, order = None, True, 0
    for item in items:
        kind = item['type']
        if kind == 'update_time':
            time = literal(item['raw'])
            if isinstance(time, bool) or not isinstance(time, (int, float)) or (start is not None and time <= start):
                gaps.append({'line': item['line'], 'kind': 'update_time',
                             'reason': 'unparsed or nonincreasing update time'})
                temporal_valid = False
                continue
            states.append({'start_time': start, 'end_time': time, 'values': dict(values),
                           'source_lines': dict(locations)})
            start = time
        elif kind == 'include':
            gaps.append({'line': item['line'], 'kind': 'include_not_executed',
                         'reason': 'include directive is recorded, not executed'})
        elif kind == 'unsupported':
            removed = []
            for name in item['names']:                  # each target on the line, as precisely as its subscripts allow
                target = namelist.parse_target(name)
                removed += _invalidate_region(values, locations, target['root'], namelist.affected_region(target, None))
            gaps.append({'line': item['line'], 'kind': item['kind'], 'invalidated': sorted(set(removed)),
                         'names': item['names'], 'reason': item['reason']})
        elif kind == 'assignment':
            order += 1
            removed, reason = _apply(item, values, locations)
            record = namelist.assignment_record(order, len(states), path, item)
            record['invalidated'] = removed
            assignments.append(record)
            if reason:
                parsed = item['parsed']
                problem = parsed['problems'][0]['kind'] if parsed['problems'] else 'unsupported_literal'
                gaps.append({'line': item['line'], 'name': parsed['name'], 'invalidated': removed,
                             'kind': problem, 'reason': reason})
    states.append({'start_time': start, 'end_time': None, 'values': dict(values), 'source_lines': dict(locations)})
    return states, gaps, temporal_valid, assignments


def declared_timeline(path=None, allowed_roots=None, max_depth=namelist.MAX_DEPTH,
                      max_assignments=namelist.MAX_ASSIGNMENTS):
    """Declared state per update block from a TR.DAT, with parser, assignments, files and unresolved constructs.

    ``allowed_roots`` (default: the TR.DAT directory) bound where referenced PT templates may be read.
    """
    if path is None:
        return {'available': False, 'states': [], 'gaps': [], 'source': None, 'semantics': SEMANTICS,
                'parser': None, 'assignments': [], 'files': None, 'unresolved': []}
    path = Path(path)
    if path.stat().st_size > LIMIT_BYTES:
        raise ValueError('TR.DAT exceeds the curation 10 MiB input limit')
    data = path.read_bytes()
    body = data.decode('utf-8', errors='replace')
    sha = namelist.sha256_bytes(data)
    f90nml = namelist.load_f90nml()
    if f90nml is None:
        return _fallback(body, path, sha)
    context = namelist.Context(f90nml)
    document = namelist.parse_document(context, body, 'transp_dat')
    states, gaps, temporal_valid, assignments = _resolve(document['items'], path.resolve())
    files, documents, unresolved, reference_gaps = namelist.load_references(
        context, path, body, sha, allowed_roots, max_depth, max_assignments)
    files['entries'][0]['syntax_check'] = document['syntax_check']
    files['entries'][0]['groups'] = document['groups']
    gaps = sorted(gaps + reference_gaps, key=lambda gap: gap['line'] if gap['line'] is not None else 10 ** 9)
    unresolved = namelist.unresolved_records(path.resolve(), document['items']) + unresolved
    return _finish(states, gaps, temporal_valid, path, sha, {
        'parser': context.record, 'assignments': assignments[:max_assignments],
        'assignments_truncated': len(assignments) > max_assignments, 'assignment_limit': max_assignments,
        'files': files, 'documents': documents, 'unresolved': unresolved})


def _fallback(body, path, sha):
    states, gaps, temporal_valid = scalar_reader(body)
    gaps.insert(0, {'line': None, 'kind': 'parser_fallback',
                    'reason': 'f90nml unavailable: the limited scalar reader was used (see parser.limitations)'})
    return _finish(states, gaps, temporal_valid, path, sha, {
        'parser': {'kind': 'scalar_fallback', 'version': None, 'config': {}, 'limitations': list(FALLBACK_LIMITATIONS)},
        'assignments': [], 'assignments_truncated': False, 'files': None, 'documents': [],
        'unresolved': [{'kind': 'parser_fallback', 'file': str(Path(path).resolve()), 'line': None,
                        'reason': 'f90nml is not importable; arrays, repeats and references are unresolved'}]})


def scalar_reader(body):
    """The limited fallback: one scalar assignment per line (the pre-f90nml reader); ``(states, gaps, temporal_valid)``."""
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
        base, index, raw = match.groups()
        base = normalize_name(base)
        value = literal(raw)
        single = INTEGER_INDEX.match(index) if index is not None else None
        name = f'{base}({int(single[1])})' if single else base
        if index is not None and not single:
            # a slice, a bound-less range or a multidimensional index: the assignment is not resolved here, so
            # whatever it may have changed is unknown (an established slice range keeps its unaffected elements)
            span = SLICE_INDEX.match(index)
            lo, hi = (int(span[1]), int(span[2])) if span else (None, None)
            removed = invalidate(values, locations, base, lo, hi)
            gaps.append({'line': line_number, 'name': f'{base}({index.strip()})', 'invalidated': removed,
                         'reason': 'array slice or multidimensional assignment is not resolved by the scalar reader'})
        elif value is None:
            removed = [name] if name in values else []
            values.pop(name, None)  # do not reuse a superseded, now-unknown value
            locations.pop(name, None)
            if not single:                                  # a whole-array update that cannot be read
                removed = invalidate(values, locations, base) + removed
            gaps.append({'line': line_number, 'name': name, 'invalidated': sorted(set(removed)),
                         'reason': 'not a single scalar literal'})
        else:
            if not single:                                  # a bare scalar for a name that holds elements is ambiguous
                invalidate(values, locations, base)
            values[name] = value
            locations[name] = line_number
    states.append({'start_time': start, 'end_time': None, 'values': dict(values),
                   'source_lines': dict(locations)})
    return states, gaps, temporal_valid
