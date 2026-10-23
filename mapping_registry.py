#!/usr/bin/python
"""Load, validate and document the declarative TRANSP -> IMAS mapping registry.

``mapping_registry.yaml`` is the single place that says what is mapped, deferred or only
recognised.  This module has no IMAS or NetCDF dependency; checks that need the Data Dictionary or
a CDF are in ``examples/verify_mapping_registry.py`` and receive those as arguments.
"""

import fnmatch
from pathlib import Path

import yaml

REGISTRY_PATH = Path(__file__).resolve().parent / 'mapping_registry.yaml'

KINDS = ('direct_typed', 'derived_typed', 'xml_only', 'input_consumed', 'deferred')
TYPED_KINDS = ('direct_typed', 'derived_typed')
VERIFICATIONS = ('readback_verified', 'writer_verified', 'unverified')
# Registry destination keys are DD major versions; these are the versions checked for each.
DD_VERSIONS = {'3': '3.39.0', '4': '4.1.1'}
REQUIRED = {
    'direct_typed': ('writer', 'ids', 'source', 'destination', 'units', 'verification'),
    'derived_typed': ('writer', 'ids', 'source', 'destination', 'units', 'transform', 'verification'),
    'xml_only': ('writer', 'ids', 'destination', 'verification'),
    'input_consumed': ('source', 'reason'),
    'deferred': ('source', 'reason'),
}
KNOWN_FIELDS = {'id', 'kind', 'writer', 'ids', 'source', 'input_files', 'destination', 'units',
                'transform', 'sign', 'reason', 'verification', 'tests', 'read_by'}


def load_registry(path=REGISTRY_PATH):
    with open(path, encoding='utf-8') as handle:
        entries = yaml.safe_load(handle)
    if not isinstance(entries, list):
        raise ValueError(f'{path}: the registry must be a list of entries')
    return entries


def matches(pattern, name):
    """Shell-style match where '?' is one character, as in the TRANSP names used by the writers."""
    return fnmatch.fnmatchcase(name, pattern)


def entry_matches_variable(entry, name):
    return any(matches(pattern, name) for pattern in entry.get('source') or ())


def validate_schema(entries, *, writers=None, root=None):
    """Return a list of problems with the registry structure (empty when valid)."""
    problems = []
    seen = set()
    for index, entry in enumerate(entries):
        label = entry.get('id', f'#{index}') if isinstance(entry, dict) else f'#{index}'
        if not isinstance(entry, dict):
            problems.append(f'{label}: entry is not a mapping')
            continue
        if label in seen:
            problems.append(f'{label}: duplicate id')
        seen.add(label)
        kind = entry.get('kind')
        if kind not in KINDS:
            problems.append(f'{label}: kind {kind!r} not in {KINDS}')
            continue
        unknown = set(entry) - KNOWN_FIELDS
        if unknown:
            problems.append(f'{label}: unknown field(s) {sorted(unknown)}')
        for field in REQUIRED[kind]:
            if entry.get(field) in (None, '', {}):
                if field == 'source' and kind == 'xml_only':
                    continue
                problems.append(f'{label}: {kind} requires {field}')
        if 'source' in entry and not isinstance(entry['source'], list):
            problems.append(f'{label}: source must be a list')
        if kind in TYPED_KINDS and not entry.get('source'):
            problems.append(f'{label}: a typed mapping needs at least one source variable')
        if 'verification' in entry and entry['verification'] not in VERIFICATIONS:
            problems.append(f'{label}: verification {entry["verification"]!r} not in {VERIFICATIONS}')
        destination = entry.get('destination')
        if destination is not None:
            bad = set(destination) - set(DD_VERSIONS) if isinstance(destination, dict) else {'<not a mapping>'}
            if bad:
                problems.append(f'{label}: destination keys must be among {sorted(DD_VERSIONS)}, got {sorted(bad)}')
        if entry.get('verification') in ('readback_verified', 'writer_verified') and not entry.get('tests'):
            problems.append(f'{label}: {entry["verification"]} requires tests')
        for test in entry.get('tests') or ():
            if root is not None and not (Path(root) / test).is_file():
                problems.append(f'{label}: test {test} does not exist')
        if writers is not None and entry.get('writer') is not None and entry['writer'] not in writers:
            problems.append(f'{label}: writer {entry["writer"]!r} is not a conversion writer')
    return problems


def render_markdown(entries):
    """Deterministic Markdown table of the registry, checked into docs/mapping-registry.md."""
    lines = ['# TRANSP to IMAS mapping registry', '',
             'Generated from `mapping_registry.yaml` by `python examples/verify_mapping_registry.py --write-docs`;',
             'do not edit by hand.  Recognising a variable name is not an export: only typed entries are exports.', '']
    for kind in KINDS:
        group = [e for e in entries if e['kind'] == kind]
        if not group:
            continue
        lines += [f'## {kind} ({len(group)})', '']
        if kind in ('deferred', 'input_consumed'):
            lines += ['| id | source | reason |', '|---|---|---|']
            for e in group:
                lines.append(f'| {e["id"]} | {", ".join(f"`{s}`" for s in e["source"])} | {e["reason"]} |')
        else:
            lines += ['| id | writer | source | destination (DD4 / DD3) | units | verification | tests |',
                      '|---|---|---|---|---|---|---|']
            for e in group:
                dest = e.get('destination', {})
                d4, d3 = dest.get('4', '-'), dest.get('3', '-')
                dest_text = f'`{d4}`' if d4 == d3 else f'`{d4}` / `{d3}`'
                sources = ', '.join(f'`{s}`' for s in e.get('source') or ()) or ', '.join(e.get('input_files') or ()) or '-'
                tests = ', '.join(f'`{Path(t).name}`' for t in e.get('tests') or ()) or '-'
                lines.append(f'| {e["id"]} | {e["writer"]} | {sources} | {dest_text} | '
                             f'{e.get("units", "-")} | {e["verification"]} | {tests} |')
        lines.append('')
    return '\n'.join(lines)
