#!/usr/bin/python
"""Snapshot conversions and compare two snapshots IDS by IDS, to prove a refactor changed no output.

``snapshot`` converts bundled CDFs through ``convert.py`` into one IMAS HDF5 tree per case;
``compare`` reads two such trees back with IMAS-Python and reports every IDS node whose value differs
(``imas.util.idsdiffgen``), ignoring only the fields that legitimately change between two runs of the
same code (creation date, run directory and similar, see ``VOLATILE``).  A refactor is accepted when
``compare`` of a snapshot taken before and one taken after reports no difference.

    python examples/compare_conversions.py snapshot build/snap-before [--cases ITER-38530P94 ...] [--dd-version 4.1.1]
    python examples/compare_conversions.py compare build/snap-before build/snap-after

``compare`` exits 1 on any difference or any IDS present on one side only.
"""

import argparse
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'examples' / 'data'

# Small and varied: ITER (multi-impurity, all IDSs), DIII-D with recombination neutrals, NSTX with PT_SOLVER regions,
# C-Mod with sawtooth crashes.
DEFAULT_CASES = {
    'ITER-38530P94': DATA / 'iter' / '38530P94.CDF',
    'D3D-147634T16': DATA / 'd3d' / '147634T16' / '147634T16.CDF',
    'NSTX-141711Y05': DATA / 'nstx' / '141711Y05' / '141711Y05.CDF',
    'CMOD-89846': DATA / 'cmod' / '89846.CDF',
}

# Paths (regex over the IDSDiff path text) that differ between two runs of unchanged code.  Measured: with
# nothing ignored, two runs of the same code differ only in the creation date; provenance, the code
# parameters XML and the properties comment are reproducible and are compared.
VOLATILE = [
    r'ids_properties/creation_date',
]


def snapshot(out_dir, cases, dd_version):
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    failed = []
    for label, cdf in cases.items():
        if not Path(cdf).is_file():
            print(f'SKIP {label}: {cdf} is not installed')
            continue
        target = out_dir / label
        target.mkdir(parents=True, exist_ok=True)
        result = subprocess.run([sys.executable, str(ROOT / 'convert.py'), '-i', str(cdf), '-b', 'HDF5',
                                 '-u', str(target), '--dd-version', dd_version, '-new'],
                                cwd=ROOT, capture_output=True, text=True)
        status = 'ok' if result.returncode == 0 else f'FAILED ({result.returncode})'
        print(f'{label}: {status}')
        if result.returncode:
            failed.append(label)
            (out_dir / f'{label}.log').write_text(result.stdout[-4000:] + result.stderr[-4000:])
    return 1 if failed else 0


def flat(value, width=70):
    """One-line, shortened text of a diff value (arrays print on many lines)."""
    text = ' '.join(str(value).split())
    return text if len(text) <= width else text[:width] + '...'


def entries(root):
    """{(ids name, occurrence): file} of the single HDF5 entry below ``root``."""
    files = {}
    for path in sorted(Path(root).glob('*/*/*/*/*.h5')):
        if path.name == 'master.h5':
            continue
        match = re.fullmatch(r'(.+?)(?:_(\d+))?\.h5', path.name)
        files[(match.group(1), int(match.group(2) or 0))] = path
    return files


def open_entry(file, dd_version):
    import imas
    return imas.DBEntry(f'imas:hdf5?path={file.parent}', 'r', dd_version=dd_version)


def differences(a, b, ignore):
    import imas
    out = []
    for item in imas.util.idsdiffgen(a, b):
        path = str(item[0])
        if any(pattern.search(path) for pattern in ignore):
            continue
        out.append(item)
    return out


def compare(left, right, dd_version='4.1.1', limit=8):
    import imas
    ignore = [re.compile(p) for p in VOLATILE]
    bad = total = 0
    for case in sorted(p.name for p in Path(left).iterdir() if p.is_dir()):
        if not (Path(right) / case).is_dir():
            print(f'{case}: missing from {right}')
            bad += 1
            continue
        a_files, b_files = entries(Path(left) / case), entries(Path(right) / case)
        for key in sorted(set(a_files) ^ set(b_files)):
            print(f'{case} {key[0]}[{key[1]}]: only in {"left" if key in a_files else "right"}')
            bad += 1
        compared = changed = 0
        for key in sorted(set(a_files) & set(b_files)):
            a_entry, b_entry = open_entry(a_files[key], dd_version), open_entry(b_files[key], dd_version)
            try:
                name, occurrence = key
                diffs = differences(a_entry.get(name, occurrence=occurrence), b_entry.get(name, occurrence=occurrence),
                                    ignore)
            finally:
                a_entry.close()
                b_entry.close()
            compared += 1
            if diffs:
                changed += 1
                bad += 1
                print(f'{case} {name}[{occurrence}]: {len(diffs)} differing node(s)')
                for item in diffs[:limit]:
                    print('    ', str(item[0]), '|', flat(item[1]), '|', flat(item[2]))
        total += compared
        print(f'{case}: {compared} IDS compared, {changed} changed')
    if not total:
        print('nothing was compared: no IDS found on both sides')
        bad += 1
    print('IDENTICAL' if not bad else f'DIFFERENT ({bad} problem(s))')
    return 1 if bad else 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)
    snap = sub.add_parser('snapshot')
    snap.add_argument('out_dir')
    snap.add_argument('--cases', nargs='*', choices=sorted(DEFAULT_CASES), default=sorted(DEFAULT_CASES))
    snap.add_argument('--dd-version', default='4.1.1')
    cmp_ = sub.add_parser('compare')
    cmp_.add_argument('left')
    cmp_.add_argument('right')
    cmp_.add_argument('--dd-version', default='4.1.1', help='DD the snapshots were written with (not the environment default)')
    args = parser.parse_args(argv)
    if args.command == 'snapshot':
        return snapshot(args.out_dir, {c: DEFAULT_CASES[c] for c in args.cases}, args.dd_version)
    return compare(args.left, args.right, args.dd_version)


if __name__ == '__main__':
    sys.exit(main())
