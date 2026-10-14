#!/usr/bin/python
"""Convert every supplied TRANSP CDF, then validate and read back written IDSs.

Requires a functioning IMAS-Python installation; the normal source-only CI
intentionally does not run this script.  All data are written into a fresh
isolated HDF5 area, never the user's default ~/public/imasdb directory.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import traceback

import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = ROOT / 'examples' / 'data'


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, default=DEFAULT_DATA,
                        help='recursively discover every *.CDF under this directory')
    parser.add_argument('--report-dir', type=Path, default=ROOT / 'build' / 'end-to-end',
                        help='persist JSON summary, per-case JSON and conversion logs')
    parser.add_argument('--database-dir', type=Path, default=None,
                        help='persist HDF5 test databases here (default: temporary and deleted)')
    parser.add_argument('--dd-version', default='4.1.1', help='target IMAS DD version')
    parser.add_argument('--backend', default='HDF5', choices=['HDF5'],
                        help='isolated end-to-end regression currently requires HDF5')
    parser.add_argument('--ids', default=None, help='optional comma-separated IDS subset')
    parser.add_argument('--case', action='append', default=[], metavar='RUNID',
                        help='restrict to this run ID; may be repeated; default all CDFs')
    parser.add_argument('--timeout', type=int, default=1800,
                        help='per-CDF conversion timeout, seconds (default: 1800)')
    parser.add_argument('--python', default=sys.executable,
                        help='Python interpreter used to run the converter')
    parser.add_argument('--plan-only', action='store_true',
                        help='discover cases without attempting conversion or claiming validation')
    parser.add_argument('--fail-fast', action='store_true',
                        help='stop after first failing case; default processes all')
    return parser.parse_args(argv)


def discover_cases(data_root, cases):
    root = Path(data_root)
    if not root.is_dir():
        raise FileNotFoundError(f'Data root is not a directory: {root}')
    files = sorted(root.rglob('*.CDF'))
    if cases:
        requested = set(cases)
        found = {p.stem for p in files}
        missing = requested - found
        if missing:
            raise FileNotFoundError('Requested example(s) not found: ' + ', '.join(sorted(missing)))
        files = [p for p in files if p.stem in requested]
    if not files:
        raise FileNotFoundError(f'No CDF examples found under {root}')
    ids = [str(p.relative_to(root)) for p in files]
    if len(ids) != len(set(ids)):
        raise AssertionError('Duplicate CDF example paths')
    return files


def require_imas_runtime():
    if os.environ.get('IMAS_AL_DISABLE_VALIDATE', '').lower() in ('1', 'true', 'yes'):
        raise RuntimeError('IMAS_AL_DISABLE_VALIDATE must NOT be enabled during end-to-end verification')
    try:
        import imas
    except ImportError as exc:
        raise RuntimeError('IMAS-Python is required for real conversion/readback. '
                           'Activate your IMAS environment; do not use --plan-only '
                           'as a replacement for runtime verification.') from exc
    if not hasattr(imas, 'DBEntry'):
        raise RuntimeError('Installed imas does not expose DBEntry')
    return imas


def source_time3(cdf):
    # TIME3 is the native equilibrium/core_profiles time axis.  TIME is a
    # distinct scalar-history axis; do not substitute the two or compare sizes.
    with netcdf_file(str(cdf), 'r', mmap=False) as nc:
        t3 = np.asarray(nc.variables['TIME3'][:], dtype=np.float64).copy()
        t = np.asarray(nc.variables['TIME'][:], dtype=np.float64).copy()
    if t3.ndim != 1 or not t3.size or t.ndim != 1 or not t.size:
        raise AssertionError(f'{cdf}: invalid TIME or TIME3')
    return t, t3


def check_time_vector(name, actual, expected, count):
    arr = np.asarray(actual, dtype=np.float64)
    if arr.ndim != 1 or len(arr) != len(expected) or count != len(expected):
        raise AssertionError(
            f'{name}: stored time={arr.shape}, slice count={count}; '
            f'CDF TIME3={len(expected)}')
    if not np.all(np.isfinite(arr)):
        raise AssertionError(f'{name}: nonfinite stored time')
    # Exact CDF values are directly assigned as float64 by the converter;
    # allow only floating-point read/write roundoff, NOT grid interpolation.
    if not np.allclose(arr, expected, atol=1.0e-9, rtol=1.0e-9):
        raise AssertionError(f'{name}: stored time differs from native CDF TIME3 '
                             f'(max residual {np.max(np.abs(arr-expected)):.3g} s)')
    return len(arr)


def validate_readback(ids_name, ids, time3):
    """Validate schema and known time/radial coordinate invariants."""
    ids.validate()  # IMAS DD coordinate and mandatory-field validator.
    checks = {'schema': 'PASS'}
    if ids_name == 'equilibrium':
        n = check_time_vector('equilibrium', ids.time, time3, len(ids.time_slice))
        checks['native_TIME3_samples'] = n
        for i in (0, n - 1):
            eq = ids.time_slice[i]
            if len(eq.boundary.outline.r) and len(eq.boundary.outline.r) != len(eq.boundary.outline.z):
                raise AssertionError(f'equilibrium/time_slice[{i}]: LCFS outline R/Z mismatch')
    elif ids_name == 'core_profiles':
        n = check_time_vector('core_profiles', ids.time, time3, len(ids.profiles_1d))
        checks['native_TIME3_samples'] = n
        for i in (0, n - 1):
            profile = ids.profiles_1d[i]
            rho = np.asarray(profile.grid.rho_tor_norm, dtype=np.float64)
            te = np.asarray(profile.electrons.temperature, dtype=np.float64)
            if rho.size and te.size and (rho.ndim != 1 or te.ndim != 1 or rho.size != te.size):
                raise AssertionError(f'core_profiles[{i}]: Te does not match radial grid')
    return checks


def readback_accepted_ids(imas, manifest, database_root, cdf):
    """Read only the IDSs confirmed by a successful DBEntry.put call."""
    accepted = manifest.get('written_ids', [])
    results = []
    if not accepted:
        return results
    machine, shot, run = (manifest.get(key) for key in ('machine', 'shot', 'run'))
    if not machine or shot is None or run is None:
        raise AssertionError('Conversion write journal lacks IMAS database identity')
    backend = imas.ids_defs.HDF5_BACKEND
    db = imas.DBEntry(backend, str(machine), int(shot), int(run), str(database_root),
                      dd_version=str(manifest['dd_version']))
    _, t3 = source_time3(cdf)
    opened = False
    try:
        db.open()
        opened = True
        seen = set()
        for event in accepted:
            name, occurrence = event['ids'], int(event['occurrence'])
            key = (name, occurrence)
            if key in seen:
                continue
            seen.add(key)
            record = {'ids': name, 'occurrence': occurrence, 'status': 'PASS'}
            try:
                ids = db.get(name, occurrence=occurrence)
                record['checks'] = validate_readback(name, ids, t3)
            except Exception as exc:
                record['status'] = 'FAIL'
                record['error'] = f'{type(exc).__name__}: {exc}'
            results.append(record)
    finally:
        if opened:
            db.close()
    return results


def run_one(cdf, args, database_root, imas):
    relative = str(cdf.relative_to(args.data_root))
    runid = cdf.stem
    # Keep distinct reports even for identically named runs on two machines.
    key = '__'.join(cdf.relative_to(args.data_root).with_suffix('').parts)
    manifest_file = args.report_dir / (key + '-conversion.json')
    log_file = args.report_dir / (key + '-conversion.log')
    command = [args.python, str(ROOT / 'convert.py'), '-i', str(cdf),
               '-b', 'HDF5', '-new', '-u', str(database_root),
               '--dd-version', args.dd_version, '--report-json', str(manifest_file)]
    if args.ids:
        command += ['--ids', args.ids]
    result = {'case': runid, 'cdf': relative, 'status': 'FAIL',
              'conversion_exit_code': None, 'writers': [], 'readback': [],
              'warnings': [], 'errors': [], 'log': str(log_file),
              'journal': str(manifest_file)}
    manifest_file.unlink(missing_ok=True)
    started = time.monotonic()
    try:
        with log_file.open('w') as output:
            completed = subprocess.run(command, cwd=str(ROOT), stdout=output,
                                       stderr=subprocess.STDOUT, timeout=args.timeout,
                                       env={**os.environ, 'MPLBACKEND': 'Agg'},
                                       check=False)
        result['conversion_exit_code'] = completed.returncode
    except subprocess.TimeoutExpired:
        result['errors'].append(f'conversion timeout after {args.timeout} seconds')
    except Exception as exc:
        result['errors'].append(f'could not launch converter: {type(exc).__name__}: {exc}')
    result['seconds'] = round(time.monotonic() - started, 2)
    if manifest_file.is_file():
        try:
            manifest = json.loads(manifest_file.read_text())
            result['writers'] = manifest.get('writers', [])
            if manifest.get('status') != 'PASS':
                result['errors'].append('converter: ' + str(manifest.get('error')))
            result['readback'] = readback_accepted_ids(imas, manifest, database_root, cdf)
            for item in result['readback']:
                if item['status'] == 'FAIL':
                    result['errors'].append(item['ids'] + ': ' + item['error'])
            if result['conversion_exit_code'] == 0 and manifest.get('status') == 'PASS' and not result['errors']:
                result['status'] = 'PASS' if manifest.get('written_ids') else 'SKIP'
        except Exception as exc:
            result['errors'].append(f'readback error: {type(exc).__name__}: {exc}')
    else:
        result['errors'].append('converter did not write a conversion journal; inspect log')
    if result['conversion_exit_code'] not in (0, None):
        result['errors'].append(f'converter exit code {result["conversion_exit_code"]}')
    return result


def main(argv=None):
    args = parse_arguments(argv)
    args.data_root = args.data_root.resolve()
    args.report_dir = args.report_dir.resolve()
    files = discover_cases(args.data_root, args.case)
    if args.ids is not None:
        # Fail before any conversion if the selected IDS list is misspelled.
        from conversion_plan import select_ids
        select_ids(args.ids)
    args.report_dir.mkdir(parents=True, exist_ok=True)
    if args.plan_only:
        plan = {'status': 'NOT_RUN', 'reason': 'plan-only does not test IMAS',
                'cases': [str(p.relative_to(args.data_root)) for p in files]}
        (args.report_dir / 'summary.json').write_text(json.dumps(plan, indent=2) + '\n')
        print(f'PLAN ONLY: {len(files)} CDFs discovered; no conversion or readback performed')
        return 0
    imas = require_imas_runtime()
    if args.database_dir:
        database_root = args.database_dir.resolve()
        database_root.mkdir(parents=True, exist_ok=True)
        temporary = None
    else:
        temporary = tempfile.TemporaryDirectory(prefix='transp2imas-imas-e2e-')
        database_root = Path(temporary.name)
    results = []
    try:
        for idx, cdf in enumerate(files, 1):
            print(f'[{idx}/{len(files)}] {cdf.relative_to(args.data_root)} ... ',
                  end='', flush=True)
            one = run_one(cdf, args, database_root, imas)
            results.append(one)
            print(one['status'] + (f' ({len(one["errors"])} errors)' if one['errors'] else ''))
            if args.fail_fast and one['status'] == 'FAIL':
                break
    finally:
        if temporary:
            temporary.cleanup()
    failures = sum(row['status'] == 'FAIL' for row in results)
    summary = {'status': 'FAIL' if failures or len(results) != len(files) else 'PASS',
               'dd_version': args.dd_version, 'backend': args.backend,
               'examples_discovered': len(files), 'examples_executed': len(results),
               'examples_passed': sum(row['status'] == 'PASS' for row in results),
               'examples_skipped': sum(row['status'] == 'SKIP' for row in results),
               'examples_failed': failures, 'databases_retained': bool(args.database_dir),
               'cases': results}
    (args.report_dir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    with (args.report_dir / 'summary.md').open('w') as report:
        report.write('# TRANSP2IMAS end-to-end IMAS regression\n\n')
        report.write(f'DD {args.dd_version}: {summary["examples_passed"]}/{len(files)} '
                     f'passed; skipped={summary["examples_skipped"]}; failures={failures}\n\n')
        report.write('| Run | Conversion | IDS readbacks | Errors |\n|---|---|---|---|\n')
        for one in results:
            count = sum(row['status'] == 'PASS' for row in one['readback'])
            report.write(f'| {one["case"]} | {one["status"]} | '
                         f'{count}/{len(one["readback"])} | '
                         f'{"; ".join(one["errors"]).replace("|", "/")} |\n')
    print(f'IMAS end-to-end: {summary["status"]} '
          f'({summary["examples_passed"]}/{len(files)} cases); '
          f'logs and reports in {args.report_dir}')
    return 1 if summary['status'] == 'FAIL' else 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f'END-TO-END CONFIGURATION ERROR: {exc}', file=sys.stderr)
        raise SystemExit(2) from exc
