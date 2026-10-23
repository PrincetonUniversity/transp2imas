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
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from theory_transport_mapping import MODEL_SPECS, AUXILIARY_SPECS, build_mmm_plan
from sawteeth_mapping import build_sawteeth_plan
from input_provenance import InputInventory
from pt_solver_regions import build_pt_region_plan, region_plan_xml
import xml.etree.ElementTree as ET
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


def is_companion_cdf(path):
    """True for a readable NetCDF file without TIME3/X/XB, e.g. ``<run>PH.CDF``.

    Such companion files are not time-dependent TRANSP results (the converter refuses them:
    "machine is not found in CDF header"), so they are not end-to-end cases.  A file that
    cannot be opened as NetCDF is *not* classified here and is left to fail in the converter.
    """
    from scipy.io import netcdf_file
    with open(path, 'rb') as handle:
        magic = handle.read(4)
    if len(magic) < 4 or magic[:3] != b'CDF':         # not a NetCDF classic header
        return False
    try:
        with netcdf_file(path, mmap=False) as nc:
            return not {'TIME3', 'X', 'XB'} <= set(nc.variables)
    except (OSError, ValueError, TypeError, IndexError):
        return False


def discover_cases(data_root, cases):
    root = Path(data_root)
    if not root.is_dir():
        raise FileNotFoundError(f'Data root is not a directory: {root}')
    files = [p for p in sorted(root.rglob('*.CDF')) if not is_companion_cdf(p)]
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


# transport_solver_numerics control parameter -> CDF channel (all on the native TIME coordinate)
RUNTIME_SOURCES = {'newton_iterations': 'NPTRITER', 'cumulative_cpu_time_hours': 'PT_CPTIME',
                   'cumulative_wall_time_hours': 'PT_WTIME'}


def source_runtime_diagnostics(cdf):
    """Native TIME and the PT_SOLVER runtime channels that exist in the CDF."""
    with netcdf_file(str(cdf), 'r', mmap=False) as nc:
        source = {'TIME': np.asarray(nc.variables['TIME'][:], dtype=np.float64).copy()}
        for parameter, variable in RUNTIME_SOURCES.items():
            if variable in nc.variables:
                source[parameter] = np.asarray(nc.variables[variable][:], dtype=np.float64).copy()
    return source


def check_solver_numerics_readback(ids, source):
    """The solver slices are on the native TIME coordinate (never compared with TIME3)."""
    time = source['TIME']
    if not np.allclose(np.asarray(ids.time, dtype=np.float64), time, rtol=1.0e-9, atol=1.0e-9):
        raise AssertionError('transport_solver_numerics.time is not the native CDF TIME')
    if len(ids.solver_1d) != time.size:
        raise AssertionError(f'transport_solver_numerics: {len(ids.solver_1d)} slices for '
                             f'{time.size} native TIME samples')
    names = set()
    for it, solver_slice in enumerate(ids.solver_1d):
        parameters = solver_slice.control_parameters
        for slot in list(parameters.integer0d) + list(parameters.real0d):
            name = str(slot.name)
            if name not in RUNTIME_SOURCES or name not in source:
                raise AssertionError(f'transport_solver_numerics: unexpected control parameter {name!r}')
            if not np.isclose(float(slot.value), source[name][it], rtol=1.0e-9, atol=1.0e-12):
                raise AssertionError(f'transport_solver_numerics.{name}[{it}] does not match the CDF')
            names.add(name)
    if not names:
        raise AssertionError('transport_solver_numerics has no control parameters')
    return {'native_TIME_samples': int(time.size), 'control_parameters': sorted(names)}


def source_mmm_coordinates(cdf):
    """Retrieve native source grids and MMM component identities for readback."""
    with netcdf_file(str(cdf), 'r', mmap=False) as nc:
        plan = build_mmm_plan(nc)
    if not plan['active']:
        return None
    mode_models = {
        spec['name'] for spec in MODEL_SPECS
        if any(source in plan['fields'] and AUXILIARY_SPECS[source][1] == 'X'
               for source in spec['auxiliary'])
    }
    return {'X': plan['x'], 'XB': plan['xb'], 'mode_models': mode_models}


def check_mmm_native_grids(ids, time3, source):
    """Enforce native TIME3/X and TIME3/XB after an actual IMAS readback.

    Coefficients always use grid_d (XB).  Only reduced MMM components with
    mode quantities on X have grid_flux (X); their actual flux nodes remain
    unset.  Comparing full native grids catches accidental regridding.
    """
    count = 0
    for model in ids.model:
        name = str(model.identifier.name)
        if not name.startswith('mmm07_'):
            continue
        if len(model.profiles_1d) != len(time3):
            raise AssertionError(f'{name}: expected {len(time3)} TIME3 profiles, '
                                 f'found {len(model.profiles_1d)}')
        radial = source['X'] if name in source['mode_models'] else source['XB']
        for i, profile in enumerate(model.profiles_1d):
            if not np.isclose(float(profile.time), time3[i], rtol=1e-9, atol=1e-9):
                raise AssertionError(f'{name}[{i}]: profile time != CDF TIME3')
            for node, expected, label in ((profile.grid_d.rho_tor_norm,
                                           source['XB'][i], 'grid_d/XB'),
                                          (profile.grid_flux.rho_tor_norm,
                                           radial[i], 'grid_flux/' +
                                           ('X' if name in source['mode_models'] else 'XB'))):
                actual = np.asarray(node, dtype=np.float64)
                if actual.shape != expected.shape or not np.allclose(
                        actual, expected, rtol=1e-9, atol=1e-9):
                    raise AssertionError(f'{name}[{i}]: {label} differs from CDF native grid')
        count += 1
    return count


def check_time_vector(name, actual, expected, count, coordinate='TIME3'):
    arr = np.asarray(actual, dtype=np.float64)
    if arr.ndim != 1 or len(arr) != len(expected) or count != len(expected):
        raise AssertionError(
            f'{name}: stored time={arr.shape}, slice count={count}; '
            f'CDF {coordinate}={len(expected)}')
    if not np.all(np.isfinite(arr)):
        raise AssertionError(f'{name}: nonfinite stored time')
    # Exact CDF values are directly assigned as float64 by the converter;
    # allow only floating-point read/write roundoff, NOT grid interpolation.
    if not np.allclose(arr, expected, atol=1.0e-9, rtol=1.0e-9):
        raise AssertionError(f'{name}: stored time differs from native CDF {coordinate} '
                             f'(max residual {np.max(np.abs(arr-expected)):.3g} s)')
    return len(arr)


def source_sawteeth(cdf):
    """Build native TIME crash history from CDF, if there are actual updates."""
    with netcdf_file(str(cdf), 'r', mmap=False) as nc:
        return build_sawteeth_plan(nc)


def check_sawteeth_readback(ids, plan):
    """Fail when a saved native crash trace, q1 radius, or period is changed."""
    if plan is None:
        raise AssertionError('sawteeth was written despite no observed source events')
    n = check_time_vector('sawteeth', ids.time, plan.time, len(plan.time), coordinate='TIME')
    fields = [('previous_crash_time', plan.previous_crash_time)]
    if plan.previous_period is not None:
        fields.append(('previous_period', plan.previous_period))
    if plan.rho_tor_norm_q1 is not None:
        fields.append(('rho_tor_norm_q1', plan.rho_tor_norm_q1))
    for name, expected in fields:
        stored = np.asarray(getattr(ids.diagnostics, name), dtype=np.float64)
        if stored.shape != expected.shape or not np.allclose(stored, expected,
                                                                rtol=1.0e-9, atol=1.0e-9):
            raise AssertionError(f'sawteeth.diagnostics.{name}: does not match native TIME source')
    if np.asarray(ids.crash_trigger).size:
        raise AssertionError('sawteeth.crash_trigger has an invented crash mechanism')
    return {'native_TIME_samples': n, 'observed_new_crashes': plan.events,
            'q1_valid_samples': plan.valid_q1_samples}


def source_pt_solver_regions(cdf):
    trdat = cdf.with_name(cdf.stem + 'TR.DAT')
    return build_pt_region_plan(InputInventory(trdat, cdf.stem))


def check_pt_solver_region_readback(ids, expected_plan):
    """Check regional *inputs* after actual DD4 IDS write and readback."""
    wrapped = ids.code.parameters
    stored = getattr(wrapped, 'value', wrapped)
    if not isinstance(stored, str) or not stored.strip():
        if expected_plan is None:
            return {'pt_solver_region_configuration': 'NOT_APPLICABLE'}
        raise AssertionError('core_transport.code.parameters lacks PT_SOLVER input XML')
    root = ET.fromstring(stored)
    actual = root.find('pt_solver_regions')
    if expected_plan is None:
        if actual is not None:
            raise AssertionError('unexpected PT_SOLVER region metadata')
        return {'pt_solver_region_configuration': 'NOT_APPLICABLE'}
    if actual is None:
        raise AssertionError('PT_SOLVER region configuration missing after IMAS readback')
    expected = region_plan_xml(expected_plan)
    if ET.tostring(actual, encoding='unicode') != ET.tostring(expected, encoding='unicode'):
        raise AssertionError('PT_SOLVER regional inputs changed after IMAS readback')
    return {'pt_solver_region_configuration': 'PASS',
            'pt_solver_boundary_states': len(expected_plan['boundary_states']),
            'pt_solver_model_configs': sum(len(x['blocks']) for x in expected_plan['templates'])}


def validate_readback(ids_name, ids, time3, mmm_source=None,
                      sawteeth_source=None, pt_solver_source=None, solver_source=None):
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
    elif ids_name == 'sawteeth':
        checks.update(check_sawteeth_readback(ids, sawteeth_source))
    elif ids_name == 'transport_solver_numerics':
        checks.update(check_solver_numerics_readback(ids, solver_source))
    elif ids_name == 'core_transport':
        checks.update(check_pt_solver_region_readback(ids, pt_solver_source))
        # Regional settings-only models intentionally have no profile array, so the TIME3 check
        # applies to the models that carry profiles; the combined model must always have them.
        for model in ids.model:
            if len(model.profiles_1d):
                n = check_time_vector('core_transport', ids.time, time3, len(model.profiles_1d))
                checks['native_TIME3_samples'] = n
            elif str(model.identifier.name) == 'combined':
                raise AssertionError('core_transport: combined model has no profiles_1d')
        if mmm_source is not None:
            checks['mmm_native_coordinate_models'] = check_mmm_native_grids(
                ids, time3, mmm_source)
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
    mmm_source = (source_mmm_coordinates(cdf)
                  if any(row['ids'] == 'core_transport' for row in accepted)
                  else None)
    sawteeth_source = (source_sawteeth(cdf)
                       if any(row['ids'] == 'sawteeth' for row in accepted)
                       else None)
    pt_solver_source = (source_pt_solver_regions(cdf)
                        if any(row['ids'] == 'core_transport' for row in accepted)
                        else None)
    solver_source = (source_runtime_diagnostics(cdf)
                     if any(row['ids'] == 'transport_solver_numerics' for row in accepted)
                     else None)
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
                record['checks'] = validate_readback(
                    name, ids, t3, mmm_source, sawteeth_source,
                    pt_solver_source if name == 'core_transport' else None,
                    solver_source if name == 'transport_solver_numerics' else None)
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
