#!/usr/bin/python
"""TRANSP namelist adapter: source collection (stage 1) and syntax parsing with f90nml (stage 2).

This module never resolves a declared state: that is stage 3, ``curation.inputs``.  It reads text, splits it into
statements with original line/column spans, collects file references, and parses each statement's syntax.

Stage 1 handles what is not a standard namelist: a TR.DAT holds bare ``NAME=value`` assignments without a ``&group``
wrapper, ``~update_time`` lines and ``PT_TEMPLATE='file.nml'`` references.  Only the reference forms that
``input_provenance`` recognises are collected; ``include`` directives are recorded and never executed.  There is no
``eval``, no shell and no expression evaluation.  Referenced files resolve relative to the including document, inside
the allowed roots, with a depth limit and cycle detection.

Stage 2 hands every statement, as its own one-assignment namelist, to ``f90nml``; an element is trusted only when the
raw token (``fortran_literals``) and f90nml agree.  f90nml alone is not trusted: it reads ``1+2`` as 100.0 (a Fortran
exponent form), reads a ``$word`` as a null, silently drops surplus values and strides, and is lenient about logicals.
An element that is not a supported scalar is UNKNOWN, never 0, ``False`` or an earlier value.  A null entry
(``A=1,,3``) is not an unknown: under standard namelist semantics it leaves that element unchanged.

Limits: syntax only.  TRANSP defaults, declared array bounds, preprocessor behaviour and model execution are not known
here; a whole-array assignment without a subscript assumes lower bound 1 and says so.
"""

import contextlib
import hashlib
import io
import itertools
import math
from pathlib import Path
import re
import warnings

import fortran_literals
from input_provenance import (assignments as provenance_assignments, references_from_text, scalar_file_value,
                              strip_fortran_comment)

MAX_BYTES = 10 * 1024 * 1024
MAX_DEPTH = 4
MAX_FILES = 32
MAX_ELEMENTS = 100000           # one statement; guards ``1000000000*0``
MAX_ASSIGNMENTS = 5000
MAX_RECORD_VALUES = 256
KINDS = ('transp_dat', 'namelist')
INCLUDE = re.compile(r'^#?\s*include\s+(?:[\'\"]([^\'\"]+)[\'\"]|(\S+))\s*$', re.I)
UPDATE = re.compile(r'^\s*~\s*update_time\s*=\s*(.*)$', re.I)
START = re.compile(r'(?<![\w%.])([A-Za-z_]\w*(?:\s*\([^()]*\))?(?:\s*%\s*[A-Za-z_]\w*(?:\s*\([^()]*\))?)*)\s*=(?!=)')
COMPONENT = re.compile(r'\s*([A-Za-z_]\w*)\s*(?:\(([^()]*)\))?\s*\Z')
INTEGER = re.compile(r'\s*([+-]?\d+)\s*\Z')
SECTION = re.compile(r'\s*([+-]?\d+)\s*:\s*([+-]?\d+)\s*\Z')
REPEAT = re.compile(r'(\d+)\*')
# Fortran reads ``.T`` as true, but ``fortran_literals`` (the shared contract) lists only ``.T.``; these four forms
# are the whole adapter-level extension, needed by the bundled TR.DAT files (``NLMDIF=.T``).
DOTTED_LOGICALS = {'.T': True, '.F': False, '.TRUE': True, '.FALSE': False}
OPTIONS = {'comment_tokens': '!', 'default_start_index': 1, 'global_start_index': None, 'row_major': False,
           'sparse_arrays': False, 'strict_logical': False}


def load_f90nml():
    """The ``f90nml`` module, or ``None`` when it cannot be imported (the caller then reports a fallback)."""
    try:
        import f90nml
    except ImportError:
        return None
    return f90nml


def scalar_value(raw, dotted=True):
    """One supported scalar (``bool``, number as ``float``, ``str``) or ``None``: ``fortran_literals``, plus ``.T``/``.F``."""
    value = fortran_literals.parse_scalar(raw)
    if value is None and dotted:
        value = DOTTED_LOGICALS.get(str(raw).strip().rstrip(',').strip().upper())
    return float(value) if isinstance(value, int) and not isinstance(value, bool) else value


class Context:
    """One f90nml parser with the explicit configuration recorded in every document."""

    def __init__(self, f90nml):
        self.f90nml = f90nml
        self.parser = f90nml.Parser()
        applied = {}
        for name, value in OPTIONS.items():
            if hasattr(self.parser, name):
                setattr(self.parser, name, value)
                applied[name] = value
        self.record = {
            'kind': 'f90nml', 'version': getattr(f90nml, '__version__', None),
            'config': {
                'parser_options': applied,
                'unsupported_parser_options': sorted(set(OPTIONS) - set(applied)),
                'comments': "'!' outside quotes, removed before parsing; quoted '!' and doubled quotes are kept",
                'start_index': 'subscripts are kept as written; a whole-array assignment without a subscript '
                               'assumes lower bound 1 and is reported as unresolved (assumed_lower_bound)',
                'array_order': 'Fortran column-major; flat keys NAME(i) or NAME(i,j,...); no positional lists',
                'dtype': 'derived-type components a%b%c are flat keys A%B%C; a subscript inside the path is unsupported',
                'repeat_and_null': 'n*value is expanded by the adapter token layer and cross-checked with f90nml; a '
                                   'null entry leaves its element unchanged and is not an unknown',
                'scalars': 'fortran_literals.parse_scalar plus the dotted logicals .T/.F/.TRUE/.FALSE; anything else '
                           'is unknown (so 1+2 is unknown although f90nml reads 100.0)',
                'statements': 'one f90nml call per statement; f90nml warnings or errors make the statement unresolved',
                'multiple_assignments_per_line': {'transp_dat': False, 'namelist': True},
                'continuation': 'namelist documents: a non-assignment line continues the previous statement; TR.DAT: only '
                                'after a trailing comma',
                'template_header': 'a bare integer before the first group of a namelist document is the PT template '
                                   'group count; it is recorded and checked against the group count, not interpreted',
            }}


def parser_record():
    """``{'kind', 'version', 'config'}`` for the f90nml adapter, or ``None`` when f90nml is unavailable."""
    f90nml = load_f90nml()
    return Context(f90nml).record if f90nml else None


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def mask_quotes(text):
    """``text`` with every quoted region (quotes included) replaced by ``x``; an unclosed quote runs to the end."""
    out, quote = [], None
    for ch in text:
        if quote:
            out.append('x')
            if ch == quote:
                quote = None
        elif ch in '\'"':
            quote = ch
            out.append('x')
        else:
            out.append(ch)
    return ''.join(out)


# stage 1: statements ----------------------------------------------------------------------------------------------

def collect(text, kind='transp_dat'):
    """Ordered items of a document: assignments (with line/column spans), groups, update times, includes, unsupported.

    ``kind='transp_dat'`` allows ``~update_time`` and takes one assignment per line; ``'namelist'`` allows several.
    """
    multi = kind == 'namelist'
    items, last = [], None
    group, pending = None, None
    for number, raw in enumerate(text.splitlines(), 1):
        if pending:                                                     # the terminator that ended the previous line
            items.append(pending)
            pending, last = None, None
        clean = strip_fortran_comment(raw).rstrip()
        pos = len(clean) - len(clean.lstrip())
        if not clean.strip():
            continue
        header = re.match(r'[&$](\w*)', clean[pos:])
        if header:                                                      # &name, $name, &end, $end
            pos += header.end()
            closing = header.group(1).upper() == 'END'
            group = None if closing else header.group(1)
            items.append({'type': 'group', 'line': number, 'name': header.group(1), 'closing': closing})
            last = None
        rest = clean[pos:]
        if multi and mask_quotes(rest).rstrip().endswith('/') and not rest.strip().startswith('/'):
            rest = rest.rstrip()[:-1]                                   # 'a=1 /' ends the group on the same line
            pending = {'type': 'group', 'line': number, 'name': group, 'closing': True}
            group = None
        if rest.strip().startswith('/'):                                # namelist terminator
            items.append({'type': 'group', 'line': number, 'name': group, 'closing': True})
            group, last = None, None
            continue
        if not rest.strip():
            continue
        if kind == 'transp_dat':
            update = UPDATE.match(rest)
            if update:
                items.append({'type': 'update_time', 'line': number, 'raw': update.group(1).strip()})
                last = None
                continue
        include = INCLUDE.match(rest.strip())
        if include:
            filename = scalar_file_value(include.group(1) or include.group(2))
            items.append({'type': 'include', 'line': number, 'filename': filename, 'text': rest.strip()})
            last = None
            continue
        starts = list(START.finditer(mask_quotes(rest)))
        prefix = rest[:starts[0].start()].strip() if starts else rest.strip()
        if prefix:
            if kind == 'namelist' and not starts and not items and re.fullmatch(r'\d+', prefix):
                items.append({'type': 'header', 'line': number, 'count': int(prefix)})     # PT template group count
                continue
            if last is not None and (multi or last['value'].rstrip().endswith(',')):
                last['value'] = (last['value'] + ' ' + prefix).strip()
                last['end_line'] = number
                last['text'] += ' ' + prefix
            else:
                items.append({'type': 'unsupported', 'line': number, 'text': rest.strip(), 'kind': 'unsupported_statement',
                              'reason': 'unsupported statement; not interpreted', 'names': []})
                last = None
                continue
            if not starts:
                continue
        if len(starts) > 1 and not multi:
            names = [normalize_target(m.group(1)) for m in starts]
            items.append({'type': 'unsupported', 'line': number, 'text': rest.strip(), 'names': names,
                          'kind': 'multiple_assignments_per_line',
                          'reason': 'several assignments on one line are not interpreted in a TR.DAT'})
            last = None
            continue
        for index, match in enumerate(starts):
            stop = starts[index + 1].start() if index + 1 < len(starts) else len(rest)
            item = {'type': 'assignment', 'line': number, 'end_line': number, 'column': pos + match.start() + 1,
                    'lhs': match.group(1), 'value': rest[match.end():stop].strip(), 'group': group,
                    'text': rest[match.start():stop].strip()}
            items.append(item)
            last = item
    if pending:
        items.append(pending)
    return items


def normalize_target(lhs):
    return re.sub(r'\s+', '', lhs).upper()


def synthetic_namelist(items):
    """A standard namelist text for the assignments (one group per update block) and a line map to the source.

    ``line_map[k]`` is ``{'line', 'end_line', 'column'}`` of synthetic line ``k`` (1-based) in the original file; the
    group header and terminator lines map to ``None``.  Statements whose value list cannot be tokenised (unbalanced
    quote, absurd repeat count) are left out, so that f90nml never receives them.
    """
    lines, line_map, block, open_group = [], {}, 0, False
    for item in items + [{'type': 'update_time', 'line': None}]:
        if item['type'] == 'update_time':
            if open_group:
                lines.append('/')
                line_map[len(lines)] = None
                open_group = False
                block += 1
            continue
        if item['type'] != 'assignment' or item.get('parsed', {}).get('count') is None:
            continue                                    # a statement whose values cannot be tokenised is not passed on
        if not open_group:
            lines.append(f'&tr_block_{block}')
            line_map[len(lines)] = None
            open_group = True
        lines.append(f"{item['lhs']}={item['value']}")
        line_map[len(lines)] = {'line': item['line'], 'end_line': item['end_line'], 'column': item['column']}
    return {'text': '\n'.join(lines) + '\n', 'line_map': line_map}


# stage 2: statement syntax ----------------------------------------------------------------------------------------

def parse_subscripts(text):
    """``([('i', k) | ('r', lo, hi), ...], None)`` or ``(None, (kind, reason))`` for a subscript list."""
    subs = []
    for piece in text.split(','):
        single, span = INTEGER.match(piece), SECTION.match(piece)
        if single:
            subs.append(('i', int(single.group(1))))
        elif span:
            lo, hi = int(span.group(1)), int(span.group(2))
            if lo > hi:
                return None, ('empty_section', 'empty array section')
            subs.append(('r', lo, hi))
        elif ':' in piece and piece.count(':') == 1:
            return None, ('bounds_less_range', 'array section without both bounds; the array bounds are not declared here')
        elif ':' in piece:
            return None, ('strided_section', 'strided array section is not supported')
        else:
            return None, ('expression_subscript', 'subscript is not an integer literal; expressions are not evaluated')
    return subs, None


def parse_target(lhs):
    """``{'root', 'subs', 'name'}`` or ``{'error': (kind, reason), 'name'}`` for the left-hand side of an assignment."""
    name = normalize_target(lhs)
    parts = lhs.split('%')
    names, subs = [], None
    for position, part in enumerate(parts):
        match = COMPONENT.match(part)
        if not match:
            return {'error': ('syntax_error', 'unreadable assignment target'), 'name': name, 'root': name}
        names.append(match.group(1).upper())
        if match.group(2) is not None:
            if position < len(parts) - 1:
                return {'error': ('subscripted_component', 'a subscript inside a derived-type path is not supported'),
                        'name': name, 'root': '%'.join(names)}
            subs, error = parse_subscripts(match.group(2))
            if error:
                return {'error': error, 'name': name, 'root': '%'.join(names)}
    return {'root': '%'.join(names), 'subs': subs, 'name': name}


def split_values(text):
    """``([(raw | None, repeat), ...], None)`` or ``(None, reason)``; ``None`` raw is a null entry.

    Values are separated by commas and/or blanks outside quotes and parentheses.  An empty slot between separators is
    null; one trailing comma is only a terminator.
    """
    tokens, current, expect = [], '', True
    index, depth = 0, 0

    def flush():
        nonlocal current, expect
        if current:
            tokens.append(current)
            current, expect = '', False

    while index < len(text):
        ch = text[index]
        if ch in '\'"':
            end = index + 1
            while True:
                end = text.find(ch, end)
                if end < 0:
                    return None, 'unterminated quoted string'
                if text[end + 1:end + 2] == ch:
                    end += 2
                    continue
                break
            current += text[index:end + 1]
            index = end + 1
            continue
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
            if depth < 0:
                return None, 'unbalanced parenthesis'
        if depth == 0 and ch == ',':
            if current:
                flush()
            elif expect:
                tokens.append(None)
            expect = True
        elif depth == 0 and ch.isspace():
            flush()
        else:
            current += ch
        index += 1
    if depth:
        return None, 'unbalanced parenthesis'
    flush()
    if not tokens:
        return None, 'no value'
    out, total = [], 0
    for token in tokens:
        repeat, raw = 1, token
        if token is not None:
            match = REPEAT.match(token)
            if match:
                repeat = int(match.group(1))
                raw = token[match.end():] or None
                if repeat < 1:
                    return None, 'repeat count must be at least 1'
        total += repeat
        if total > MAX_ELEMENTS:
            return None, f'more than {MAX_ELEMENTS} values in one statement'
        out.append((raw, repeat))
    return out, None


def section_keys(root, subs):
    """Flat keys of an array section in Fortran column-major order (the first subscript varies fastest)."""
    ranges = [range(sub[1], sub[2] + 1) if sub[0] == 'r' else (sub[1],) for sub in subs]
    size = 1
    for item in ranges:
        size *= len(item)
    if size > MAX_ELEMENTS:
        return None
    return [element_key(root, combo[::-1]) for combo in itertools.product(*reversed(ranges))]


def element_key(root, index):
    return f"{root}({','.join(str(i) for i in index)})"


def flatten(value):
    if isinstance(value, list):
        return [leaf for item in value for leaf in flatten(item)]
    return [value]


def same_value(mine, theirs):
    if isinstance(mine, bool) or isinstance(theirs, bool):
        return mine is theirs
    if isinstance(mine, (int, float)) and isinstance(theirs, (int, float)):
        return math.isclose(mine, theirs, rel_tol=1e-12, abs_tol=0.0)
    return mine == theirs


def f90nml_parse(context, target, value_text, count):
    """``(flat values, start_index, None)`` from f90nml for one statement, or ``(None, None, reason)``.

    f90nml drops every value after the first of ``A(k)=v1,v2`` (with a warning), although a Fortran namelist assigns
    ``A(k)``, ``A(k+1)``, ...; for one subscript the equivalent section ``A(k:k+n-1)`` is given to it instead.
    """
    name, subs = target['name'], target['subs']
    if subs is not None and len(subs) == 1 and subs[0][0] == 'i' and count > 1:
        name = f"{target['root']}({subs[0][1]}:{subs[0][1] + count - 1})"
    text = f"&tr\n{name}={value_text}\n/\n"
    with warnings.catch_warnings(record=True) as caught, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('always')
        try:
            group = context.parser.reads(text)['tr']
            node, names = group, target['root'].lower().split('%')
            for name in names[:-1]:
                node = node[name]
            leaf = node[names[-1]]
        except Exception as error:               # f90nml raises assorted types (AssertionError, ValueError, KeyError)
            return None, None, f'f90nml could not parse the statement ({type(error).__name__})'
    if caught:
        return None, None, 'f90nml warned: ' + str(caught[0].message)
    return flatten(leaf), (getattr(node, 'start_index', None) or {}).get(names[-1]), None


def affected_region(target, count):
    """What an unresolved statement may have changed: ``('keys', [...])``, ``('from', first)`` or ``('base',)``."""
    subs = target.get('subs')
    if 'error' in target and target['error'][0] == 'subscripted_component':
        return ('keys', [])
    if 'error' in target or subs is None:
        return ('base',)
    if any(sub[0] == 'r' for sub in subs):
        keys = section_keys(target['root'], subs)
        return ('keys', keys) if keys is not None else ('base',)
    if len(subs) == 1:
        first = subs[0][1]
        if count is None:
            return ('from', first)
        return ('keys', [element_key(target['root'], (first + k,)) for k in range(count)])
    if count == 1:
        return ('keys', [element_key(target['root'], [sub[1] for sub in subs])])
    return ('base',)


def placement(target, count):
    """``(keys, assumption, None)`` for the first ``count`` values, or ``(None, None, (kind, reason))``."""
    root, subs = target['root'], target['subs']
    if subs is None:
        if count == 1:
            return [root], None, None
        return [element_key(root, (k,)) for k in range(1, count + 1)], 'assumed_lower_bound', None
    if any(sub[0] == 'r' for sub in subs):
        keys = section_keys(root, subs)
        if keys is None or count > len(keys):
            return None, None, ('too_many_values', f'{count} values for a section of {len(keys or [])} elements')
        return keys[:count], None, None
    if len(subs) == 1:
        return [element_key(root, (subs[0][1] + k,)) for k in range(count)], None, None
    if count == 1:
        return [element_key(root, [sub[1] for sub in subs])], None, None
    return None, None, ('unbounded_continuation',
                        'several values from a multidimensional element: the extent of the first dimension is not declared')


def parse_statement(context, item):
    """Annotate an assignment item with ``item['parsed']``: placement, per-element state, problems, affected region."""
    target = parse_target(item['lhs'])
    tokens, token_problem = split_values(item['value'])
    count = sum(repeat for _raw, repeat in tokens) if tokens else None
    parsed = {'name': target['name'], 'root': target.get('root'), 'subs': target.get('subs'), 'count': count,
              'status': 'unresolved', 'problems': [], 'elements': [], 'assumption': None,
              'scalar_name': False, 'affected': None, 'section': None}
    item['parsed'] = parsed

    def fail(kind, reason):
        parsed['problems'].append({'kind': kind, 'reason': reason})
        parsed['affected'] = affected_region(target, count)

    if 'error' in target:
        return fail(*target['error'])
    if token_problem:
        return fail('syntax_error', token_problem)
    keys, assumption, problem = placement(target, count)
    if problem:
        return fail(*problem)
    flat, start, problem = f90nml_parse(context, target, item['value'], count)
    if problem:
        return fail('f90nml_rejected', problem)
    if len(flat) != count:
        return fail('parser_disagreement', f'f90nml found {len(flat)} values, the token layer {count}')
    subs = target['subs']
    if subs is not None and start is not None and list(start)[:len(subs)] != [s[1] for s in subs]:
        return fail('parser_disagreement', f'f90nml start index {start} differs from the subscripts')
    parsed['assumption'] = assumption
    parsed['scalar_name'] = subs is None and count == 1
    if subs is not None and any(sub[0] == 'r' for sub in subs):
        parsed['section'] = section_keys(target['root'], subs)
    expanded = [(raw, repeat) for raw, repeat in tokens for _ in range(repeat)]
    for key, (raw, _repeat), theirs in zip(keys, expanded, flat):
        if raw is None:
            element = {'key': key, 'state': 'null'}
        else:
            value = scalar_value(raw)
            if value is None:
                element = {'key': key, 'state': 'unknown', 'raw': raw, 'reason': 'not a supported scalar literal'}
            elif not same_value(value, theirs):
                element = {'key': key, 'state': 'unknown', 'raw': raw,
                           'reason': f'parser disagreement (f90nml read {theirs!r})'}
            else:
                element = {'key': key, 'state': 'value', 'value': value, 'raw': raw}
        parsed['elements'].append(element)
    parsed['status'] = 'partial' if any(e['state'] == 'unknown' for e in parsed['elements']) else 'resolved'
    return None


def json_value(element):
    if element['state'] == 'value':
        return element['value']
    if element['state'] == 'null':
        return {'null': True}
    return {'unknown': element['raw']}


def assignment_record(order, block, path, item):
    """JSON-safe record of one parsed assignment, with its source span."""
    parsed = item['parsed']
    elements = parsed['elements']
    values = [json_value(e) for e in elements[:MAX_RECORD_VALUES]]
    subs = parsed['subs']
    return {'order': order, 'block': block, 'file': str(path), 'line': item['line'], 'end_line': item['end_line'],
            'column': item['column'], 'name': parsed['root'] or normalize_target(item['lhs']),
            'indices': None if subs is None else [str(s[1]) if s[0] == 'i' else f'{s[1]}:{s[2]}' for s in subs],
            'keys': [e['key'] for e in elements[:MAX_RECORD_VALUES]],
            'value': values[0] if len(values) == 1 and parsed['scalar_name'] else values,
            'values_truncated': len(elements) > MAX_RECORD_VALUES,
            'raw': item['value'], 'status': parsed['status'], 'problems': list(parsed['problems']),
            'assumption': parsed['assumption']}


def unresolved_records(path, items):
    """Constructs the syntax stage could not resolve (statement problems, directives), one record each."""
    records, assumed = [], set()
    for item in items:
        if item['type'] == 'include':
            continue                                    # reported once, with the reference, by ``load_references``
        if item['type'] == 'unsupported':
            records.append({'kind': item['kind'], 'file': str(path), 'line': item['line'], 'reason': item['reason'],
                            'names': item['names']})
        if item['type'] != 'assignment':
            continue
        parsed = item['parsed']
        for problem in parsed['problems']:
            records.append({'kind': problem['kind'], 'file': str(path), 'line': item['line'], 'name': parsed['name'],
                            'reason': problem['reason']})
        for element in parsed['elements']:
            if element['state'] == 'unknown':
                records.append({'kind': 'unsupported_literal', 'file': str(path), 'line': item['line'],
                                'name': element['key'], 'raw': element['raw'], 'reason': element['reason']})
        if parsed['assumption'] and parsed['root'] not in assumed:
            assumed.add(parsed['root'])
            records.append({'kind': parsed['assumption'], 'file': str(path), 'line': item['line'],
                            'name': parsed['root'],
                            'reason': 'array lower bound is not declared in the input; 1 is assumed'})
    return records


def whole_text_check(context, text):
    """``(group names, status)`` of f90nml reading a whole text; a failure does not stop the per-statement results."""
    with warnings.catch_warnings(record=True) as caught, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('always')
        try:
            names = list(context.parser.reads(text))
        except Exception as error:               # f90nml raises assorted types
            return [], f'f90nml rejected the text ({type(error).__name__})'
    return names, 'ok with f90nml warnings' if caught else 'ok'


def parse_document(context, text, kind='transp_dat'):
    """Stage 1 + 2 for one document: items (assignments annotated), synthetic namelist and its whole-text syntax check."""
    items = collect(text, kind)
    for item in items:
        if item['type'] == 'assignment':
            parse_statement(context, item)
    synthetic = synthetic_namelist(items)
    check = {'synthetic_namelist': whole_text_check(context, synthetic['text'])[1]}
    if kind == 'namelist':                       # the document itself, as f90nml reads it (repeated groups merge)
        names, status = whole_text_check(context, text)
        check.update({'document': status, 'document_groups': names})
    groups = [{'name': i['name'], 'line': i['line'], 'closing': i['closing']} for i in items if i['type'] == 'group']
    header = [i for i in items if i['type'] == 'header']
    if header:
        opened = sum(1 for g in groups if not g['closing'])
        check['template_count'] = {'declared': header[0]['count'], 'groups_found': opened,
                                   'matches': header[0]['count'] == opened}
    return {'items': items, 'synthetic': synthetic, 'syntax_check': check, 'groups': groups}


# stage 1: files ----------------------------------------------------------------------------------------------------

def references(text):
    """``[{'key', 'filename', 'line', 'follow'}]`` for the reference forms ``input_provenance`` recognises."""
    candidates = [(line, key, scalar_file_value(value)) for line, key, value in provenance_assignments(text)]
    for number, raw in enumerate(text.splitlines(), 1):
        match = INCLUDE.match(strip_fortran_comment(raw).strip())
        if match:
            candidates.append((number, 'INCLUDE', scalar_file_value(match.group(1) or match.group(2))))
    found = []
    for key, filename in references_from_text(text):
        line = None
        for position, candidate in enumerate(candidates):
            if candidate[1] == key and candidate[2] == filename:
                line = candidates.pop(position)[0]
                break
        follow = key != 'INCLUDE' and ('TEMPLATE' in key.upper() or filename.lower().endswith(('.nml', '.namelist')))
        found.append({'key': key, 'filename': filename, 'line': line, 'follow': follow})
    return sorted(found, key=lambda ref: ref['line'] if ref['line'] is not None else 10 ** 9)


def read_bounded(path):
    """``(text, sha256, size)``; raises ``ValueError`` over the size limit, ``OSError`` if unreadable."""
    data = Path(path).read_bytes() if Path(path).stat().st_size <= MAX_BYTES else None
    if data is None:
        raise ValueError('input file exceeds the curation 10 MiB limit')
    return data.decode('utf-8', errors='replace'), sha256_bytes(data), len(data)


def load_references(context, top_path, top_text, top_sha, allowed_roots=None, max_depth=MAX_DEPTH,
                    max_assignments=MAX_ASSIGNMENTS):
    """Follow PT-template style references; return ``(files, documents, unresolved, gaps)``.

    A referenced file is parsed as its own namelist document and its ``pt_update_time`` is preserved, not interpreted.
    """
    top = Path(top_path).resolve()
    roots = [Path(r).resolve() for r in (allowed_roots if allowed_roots else [top.parent])]
    entries = [{'path': str(top), 'sha256': top_sha, 'size': len(top_text.encode('utf-8')), 'role': 'top_level',
                'depth': 0, 'status': 'parsed', 'referenced_by': None}]
    documents, unresolved, gaps = [], [], []

    def problem(entry, kind, reason, ref, parent):
        entry.update({'status': kind, 'sha256': None})
        record = {'kind': kind, 'file': str(parent), 'line': ref['line'], 'name': ref['key'],
                  'reference': ref['filename'], 'reason': reason}
        unresolved.append(record)
        gaps.append({'line': ref['line'], 'name': ref['key'], 'kind': kind, 'reason': reason + ': ' + ref['filename']})

    def visit(parent, text, depth, chain):
        for ref in references(text):
            entry = {'path': ref['filename'], 'sha256': None, 'size': None, 'role': 'reference', 'depth': depth + 1,
                     'status': 'recognised', 'referenced_by': str(parent), 'key': ref['key'], 'line': ref['line']}
            entries.append(entry)
            if ref['key'] == 'INCLUDE':
                entry['status'] = 'include_not_executed'
                unresolved.append({'kind': 'include_not_executed', 'file': str(parent), 'line': ref['line'],
                                   'name': 'INCLUDE', 'reference': ref['filename'],
                                   'reason': 'include directives are recorded, never executed'})
                continue
            if not ref['follow']:
                entry['status'] = 'recognised_not_followed'
                continue
            candidate = Path(ref['filename'])
            resolved = (candidate if candidate.is_absolute() else parent.parent / candidate).resolve()
            entry['path'] = str(resolved)
            if not any(resolved == root or root in resolved.parents for root in roots):
                problem(entry, 'outside_allowed_roots', 'reference resolves outside the allowed roots', ref, parent)
            elif resolved in chain:
                problem(entry, 'cycle', 'reference cycle', ref, parent)
            elif depth + 1 > max_depth:
                problem(entry, 'depth_limit', f'reference depth exceeds {max_depth}', ref, parent)
            elif len(entries) > MAX_FILES:
                problem(entry, 'file_limit', f'more than {MAX_FILES} referenced files', ref, parent)
            elif not resolved.is_file():
                problem(entry, 'missing', 'referenced file not found', ref, parent)
            else:
                try:
                    child_text, sha, size = read_bounded(resolved)
                except (OSError, ValueError) as error:
                    problem(entry, 'unreadable', f'referenced file cannot be read ({type(error).__name__})', ref, parent)
                    continue
                entry.update({'sha256': sha, 'size': size, 'status': 'parsed'})
                document = parse_document(context, child_text, 'namelist')
                records = [assignment_record(order, 0, resolved, i)
                           for order, i in enumerate((i for i in document['items'] if i['type'] == 'assignment'), 1)]
                documents.append({
                    'path': str(resolved), 'sha256': sha, 'parser': context.record,
                    'syntax_check': document['syntax_check'], 'groups': document['groups'],
                    'assignments': records[:max_assignments], 'assignments_truncated': len(records) > max_assignments,
                    'unresolved': unresolved_records(resolved, document['items']),
                    'semantics': 'syntax only; pt_update_time and model semantics are preserved, not interpreted'})
                visit(resolved, child_text, depth + 1, chain | {resolved})

    visit(top, top_text, 0, {top})
    files = {'allowed_roots': [str(r) for r in roots], 'max_depth': max_depth, 'max_files': MAX_FILES,
             'entries': entries}
    return files, documents, unresolved, gaps
