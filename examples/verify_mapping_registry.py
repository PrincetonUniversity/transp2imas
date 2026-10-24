#!/usr/bin/python
"""Validate mapping_registry.yaml against the schema, tests, Data Dictionary, CDFs and measured reads.

Checks (the first two need nothing but PyYAML, the rest are skipped with a notice when their
input is unavailable):

* schema, unique ids, existing test scripts, known writers,
* ``docs/mapping-registry.md`` is the current rendering (``--write-docs`` regenerates it),
* every destination path exists in the Data Dictionary version it is registered for (IMAS-Python),
* every source pattern of a typed entry matches a variable in at least one bundled CDF,
* measured reads (``--journals DIR``): a writer that wrote an IDS read at least one of the present
  sources of each of its typed entries, and no variable registered as ``deferred`` was read by a writer not listed in its ``read_by``.
"""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import mapping_registry as reg
from conversion_plan import WRITER_ORDER

DOCS = ROOT / 'docs' / 'mapping-registry.md'


def check_docs(entries, write):
    text = reg.render_markdown(entries) + '\n'
    if write:
        DOCS.write_text(text, encoding='utf-8')
        print(f'wrote {DOCS.relative_to(ROOT)}')
        return []
    if not DOCS.is_file() or DOCS.read_text(encoding='utf-8') != text:
        return [f'{DOCS.relative_to(ROOT)} is stale; run: python examples/verify_mapping_registry.py --write-docs']
    return []


def dd_path_exists(ids_metadata, path):
    node = ids_metadata
    for part in path.split('/'):
        try:
            node = node[part]
        except (KeyError, IndexError):
            return False
    return True


def check_dd(entries):
    try:
        import imas
    except ImportError as exc:
        print(f'SKIP Data Dictionary check: {type(exc).__name__}: {exc}')
        return []
    problems = []
    factories = {}
    for entry in entries:
        for major, path in (entry.get('destination') or {}).items():
            version = reg.DD_VERSIONS[major]
            factory = factories.setdefault(version, imas.IDSFactory(version))
            for ids_name in entry['ids']:
                try:
                    metadata = factory.new(ids_name).metadata
                except (KeyError, ValueError, AttributeError):
                    problems.append(f'{entry["id"]}: IDS {ids_name} not in DD {version}')
                    continue
                if not dd_path_exists(metadata, path):
                    problems.append(f'{entry["id"]}: {ids_name}/{path} not in DD {version}')
    return problems


def cdf_variables(path):
    from scipy.io import netcdf_file
    with netcdf_file(str(path), 'r', mmap=False) as handle:
        return set(handle.variables)


def check_cdfs(entries, data_root):
    cdfs = sorted(Path(data_root).rglob('*.CDF'))
    if not cdfs:
        print(f'SKIP CDF source check: no CDF under {data_root}')
        return []
    present = set()
    for cdf in cdfs:
        present |= cdf_variables(cdf)
    problems = []
    for entry in entries:
        patterns = entry.get('source') or ()
        unmatched = [p for p in patterns if not any(reg.matches(p, name) for name in present)]
        if entry['kind'] in reg.TYPED_KINDS and patterns and len(unmatched) == len(patterns):
            problems.append(f'{entry["id"]}: no source matches a variable in any bundled CDF ({", ".join(patterns)})')
        elif unmatched and entry['kind'] in reg.TYPED_KINDS:
            print(f'note {entry["id"]}: not in the bundled CDFs (species-dependent): {", ".join(unmatched)}')
    print(f'checked sources against {len(cdfs)} CDF file(s)')
    return problems


def check_journals(entries, journals_dir):
    journals = sorted(Path(journals_dir).glob('*-conversion.json'))
    if not journals:
        print(f'SKIP measured-read check: no journals in {journals_dir}')
        return []
    problems = []
    for journal_path in journals:
        journal = json.loads(journal_path.read_text(encoding='utf-8'))
        cdf = Path(journal['source_cdf'])
        if not cdf.is_file():
            continue
        present = cdf_variables(cdf)
        read_by = {w['ids']: set(w.get('consumed', ())) for w in journal['writers']}
        wrote = {w['ids'] for w in journal['writers'] if w['status'] == 'WROTE'}
        for entry in entries:
            here = {name for name in present if reg.entry_matches_variable(entry, name)}
            if entry['kind'] == 'deferred':
                allowed = set(entry.get('read_by') or ())
                for name in sorted(here):
                    readers = {w for w, names in read_by.items() if name in names}
                    if readers - allowed:
                        problems.append(f'{journal_path.name}: {entry["id"]} is deferred but {name} was read by '
                                        f'{sorted(readers - allowed)}')
            elif entry['kind'] in reg.TYPED_KINDS and entry['writer'] in wrote and here:
                if not here & read_by.get(entry['writer'], set()):
                    problems.append(f'{journal_path.name}: {entry["id"]}: writer {entry["writer"]} wrote an IDS '
                                    f'but read none of {sorted(here)}')
    print(f'checked measured reads in {len(journals)} journal(s)')
    return problems


def self_test():
    """The validator must reject a malformed registry (no data needed)."""
    good = {'id': 'g', 'kind': 'direct_typed', 'writer': 'tf', 'ids': ['tf'], 'source': ['BZXR'],
            'destination': {'4': 'x'}, 'units': 'T', 'verification': 'unverified'}
    assert reg.validate_schema([good], writers=WRITER_ORDER) == []
    bad_cases = {
        'duplicate id': [good, dict(good)],
        'unknown kind': [dict(good, kind='guess')],
        'typed without source': [dict(good, source=[])],
        'verified without tests': [dict(good, verification='readback_verified')],
        'missing test file': [dict(good, verification='writer_verified', tests=['examples/nope.py'])],
        'unknown writer': [dict(good, writer='nope')],
        'deferred without reason': [{'id': 'd', 'kind': 'deferred', 'source': ['X']}],
        'unknown field': [dict(good, colour='red')],
        'bad dd key': [dict(good, destination={'5': 'x'})],
    }
    for label, entries in bad_cases.items():
        assert reg.validate_schema(entries, writers=WRITER_ORDER, root=ROOT), f'accepted: {label}'
    assert reg.matches('PICF??N??', 'PICF01N02') and not reg.matches('VPOL?_NC', 'VPOLD2_NC')
    print(f'self-test: {len(bad_cases)} malformed registries rejected')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--write-docs', action='store_true', help='regenerate docs/mapping-registry.md')
    parser.add_argument('--data-root', type=Path, default=ROOT / 'examples' / 'data')
    parser.add_argument('--journals', type=Path, default=None,
                        help='directory of *-conversion.json journals (measured-read check)')
    parser.add_argument('--no-dd', action='store_true', help='skip the Data Dictionary check')
    parser.add_argument('--quick', action='store_true', help='schema and docs only')
    args = parser.parse_args(argv)

    self_test()
    entries = reg.load_registry()
    problems = reg.validate_schema(entries, writers=WRITER_ORDER, root=ROOT)
    problems += check_docs(entries, args.write_docs)
    if not args.quick:
        if not args.no_dd:
            problems += check_dd(entries)
        problems += check_cdfs(entries, args.data_root)
        if args.journals:
            problems += check_journals(entries, args.journals)
    kinds = {k: sum(e['kind'] == k for e in entries) for k in reg.KINDS}
    print(f'{len(entries)} registry entries: {kinds}')
    if problems:
        print('\n'.join(f'FAIL {p}' for p in problems))
        return 1
    print('PASS mapping registry')
    return 0


if __name__ == '__main__':
    sys.exit(main())
