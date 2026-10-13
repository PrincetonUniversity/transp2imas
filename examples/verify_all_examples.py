#!/usr/bin/python
"""Discover and verify EVERY bundled TRANSP CDF, with a per-case report.

Source-side only: confirms CDF/time grids, every available EC/IC deposition,
active/passive PF mapping, raw input/XML provenance, and existing source
sanity checks. Optional physics are recorded as N/A, not imaginary data.

No IMAS installation is necessary; a real IDS write/readback remains separate.
Missing TR.DAT/referenced inputs are WARN by default and FAIL with
--strict-inputs. A known ITER ECCUR current closure at 0.543% is WARN unless
--strict-current is requested. Every CDF is processed before returning an
error status so a single broken example cannot hide coverage of other runs.
"""

import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import traceback
import xml.etree.ElementTree as ET

import numpy as np
from scipy.io import netcdf_file

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from input_provenance import InputInventory
from pf_active_mapping import build_active_plan, populate_active_ids
from pf_passive_mapping import build_passive_plan, populate_passive_ids
from verify_pf_active_mapping import ids_node as active_ids_node
from verify_pf_passive_mapping import ids_node as passive_ids_node
import verify_waves_mapping as waves
from verify_core_profiles_timegrids import verify_case as verify_core_profiles_case, verify_writer_contract
from verify_equilibrium_timegrids import verify_case as verify_equilibrium_case, writer_contract as equilibrium_writer_contract, synthetic as equilibrium_synthetic


DEEP_SOURCE_CHECKS = (
    ('cdf-physics', 'verify_cdf_extension.py'),
    ('nbi-sources', 'verify_nbi_distribution_sources.py'),
    ('sources-transport', 'verify_sources_transport.py'),
)
DEEP_INPUT_CHECKS = (
    ('namelist-beams', 'verify_distribution_sources_namelist.py'),
    ('cocos-signs', 'verify_cocos17.py'),
    ('ids-inventory', 'verify_next_ids_inventory.py'),
    ('metadata-provenance', 'verify_metadata_provenance_inputs.py'),
)


def sha256_file(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def source_grids(nc, cdf):
    names = nc.variables
    require('TIME' in names and 'TIME3' in names, f'{cdf}: TIME/TIME3 missing')
    result = {}
    for key in ('TIME', 'TIME3'):
        t = np.asarray(names[key][:], dtype=np.float64)
        require(t.ndim == 1 and t.size > 1 and np.all(np.isfinite(t)),
                f'{cdf}: invalid {key} time base')
        require(np.all(np.diff(t) >= 0), f'{cdf}: nonmonotonic {key}')
        result[key] = int(t.size)
    require('DVOL' in names and 'X' in names, f'{cdf}: DVOL/X missing')
    dvol = names['DVOL']
    x = names['X']
    require(len(dvol.shape) == 2 and len(x.shape) == 2,
            f'{cdf}: expected 2D finite-volume radial grids')
    require(dvol.shape[0] == result['TIME3'] and x.shape == dvol.shape,
            f'{cdf}: DVOL/X do not match TIME3')
    result['radial_points'] = int(dvol.shape[1])
    result['variables'] = len(names)
    return result


def process_case(cdf, report_dir, python, deep, strict_inputs, strict_current):
    runid = cdf.stem
    rel = cdf.relative_to(HERE / 'data')
    trdat = cdf.with_name(f'{runid}TR.DAT')
    result = {'case':runid, 'cdf':str(rel), 'trdat':trdat.exists(),
              'status':'PASS', 'checks':{}, 'warnings':[], 'failures':[],
              'wave_status':'N/A', 'models':[], 'input_files':[],
              'pf_active_coils':0, 'pf_passive_loops':0}
    run_log = report_dir / (runid + '.log')
    lines = [f'Case {runid}', f'CDF: {cdf}', f'TR.DAT: {trdat if trdat.exists() else "MISSING"}']

    def check(label, action):
        try:
            value = action()
            result['checks'][label] = 'PASS'
            return value
        except Exception as exc:
            msg = f'{label}: {type(exc).__name__}: {exc}'
            result['checks'][label] = 'FAIL'
            result['failures'].append(msg)
            lines.append('FAIL ' + msg)
            lines.append(traceback.format_exc())
            return None

    def source():
        with netcdf_file(str(cdf), 'r', mmap=False) as nc:
            grid = source_grids(nc, cdf)
            active = build_active_plan(nc)
            passive = build_passive_plan(nc)
            active_result = populate_active_ids(active_ids_node(), active, 4)
            passive_result = populate_passive_ids(passive_ids_node(), passive, 4)
            require(active_result['coil_current'] == sum(int(x.current is not None)
                     for x in active.coils), 'pf_active coil current count mismatch')
            require(passive_result == len(passive.loops),
                    'pf_passive current count mismatch')
            result['pf_active_coils'] = len(active.coils)
            result['pf_passive_loops'] = len(passive.loops)
            result.update(grid)
            lines.append(f"Grid: {grid}; active coils={len(active.coils)}; "
                         f"passive loops={len(passive.loops)}")
        return True

    check('cdf-and-pf-mapping', source)

    def core_profiles_time_check():
        verify_writer_contract()
        with netcdf_file(str(cdf), 'r', mmap=False) as nc:
            result['core_profiles_peaking'] = verify_core_profiles_case(nc, cdf.name)
        return True

    check('core-profiles-time-grid', core_profiles_time_check)

    def equilibrium_time_check():
        equilibrium_writer_contract()
        equilibrium_synthetic()
        with netcdf_file(str(cdf), 'r', mmap=False) as nc:
            result['equilibrium_timegrid'] = verify_equilibrium_case(nc, cdf.name)
        return True

    check('equilibrium-time-grid', equilibrium_time_check)

    def inputs():
        inv = InputInventory(trdat, runid)
        require(inv.records and inv.records[0].role == 'TRANSP_INPUT',
                'main TRANSP input marker missing')
        result['models'] = inv.models()
        for item in inv.records:
            result['input_files'].append({'name':item.path.name, 'status':item.status,
                'model':item.model, 'role':item.role})
            if item.status == 'embedded':
                require(item.sha256 == sha256_file(item.path),
                        f'input SHA256 mismatch: {item.path}')
                require(item.text, f'empty input unexpectedly embedded: {item.path}')
        # The same XML used for model-specific code.library nodes must be parseable.
        for model in inv.models():
            root = ET.fromstring(inv.scoped_xml(model))
            require(root.attrib['model'] == model, f'wrong XML scope: {model}')
            require(root.find('input_files') is not None, 'missing XML files group')
        for item in inv.missing:
            result['warnings'].append(f'input not archived: {item.path.name} [{item.status}]')
        for item in inv.records:
            if item.role == 'ADJACENT_UNREFERENCED':
                result['warnings'].append(f'adjacent file not proven used: {item.path.name}')
        if inv.truncated:
            result['warnings'].append('input reference traversal truncated')
        lines.append(f"Input inventory: {len(inv.records)} file(s), "
                     f"models={inv.models()}, missing={len(inv.missing)}")
        if strict_inputs:
            require(not inv.missing and not inv.truncated,
                    f'{runid}: input inventory incomplete')
        return True

    check('input-xml', inputs)

    def wave_check():
        if strict_current:
            waves.CURRENT_WARN_RTOL = waves.CURRENT_RTOL
        else:
            waves.CURRENT_WARN_RTOL = 6.0e-3
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            errors = waves.verify_file(cdf)
        output = buf.getvalue()
        lines.append(output)
        if 'WARN ' in output:
            for line in output.splitlines():
                if 'WARN ' in line:
                    result['warnings'].append(line.strip())
        result['wave_status'] = 'N/A' if 'SKIP no supported' in output else (
                                 'FAIL' if errors else ('WARN' if 'WARN ' in output else 'PASS'))
        require(not errors, '; '.join(errors))
        return True

    check('wave-source-closures', wave_check)

    def external(label, script, args):
        cmd = [python, str(HERE / script), *[str(a) for a in args]]
        run = subprocess.run(cmd, cwd=str(ROOT), stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, errors='replace')
        lines.append(f'$ {" ".join(cmd)}\n{run.stdout}')
        if run.returncode:
            raise RuntimeError(f'exit={run.returncode}; see {run_log}')

    if deep:
        for label, script in DEEP_SOURCE_CHECKS:
            check(label, lambda label=label, script=script: external(label, script, [cdf] +
                  ([trdat] if label == 'cdf-physics' else [])))
        for label, script in DEEP_INPUT_CHECKS:
            if not trdat.is_file():
                result['checks'][label] = 'N/A (missing TR.DAT)'
                continue
            args = [trdat] if label == 'namelist-beams' else (
                   [cdf, '--dat', trdat] if label == 'cocos-signs' else [cdf, trdat])
            check(label, lambda label=label, script=script, args=args:
                  external(label, script, args))

    result['status'] = 'FAIL' if result['failures'] else (
                       'WARN' if result['warnings'] else 'PASS')
    for issue in result['warnings']:
        lines.append('WARN ' + issue)
    for issue in result['failures']:
        lines.append('FAIL ' + issue)
    lines.append('RESULT: ' + result['status'])
    run_log.write_text('\n\n'.join(lines) + '\n', encoding='utf-8')
    return result


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root', type=Path, default=HERE / 'data')
    p.add_argument('--report-dir', type=Path, default=ROOT / 'build' / 'verification')
    p.add_argument('--python', default=sys.executable)
    p.add_argument('--quick', action='store_true', help='skip standalone deep source checks')
    p.add_argument('--strict-inputs', action='store_true',
                   help='fail if TR.DAT or a referenced model file is missing')
    p.add_argument('--strict-current', action='store_true',
                   help='fail on EC current residual greater than 0.5%%')
    ns = p.parse_args(argv)
    ns.report_dir.mkdir(parents=True, exist_ok=True)
    cdfs = sorted(ns.data_root.rglob('*.CDF'))
    if not cdfs:
        p.error('no CDF fixtures in ' + str(ns.data_root))
    results = []
    for cdf in cdfs:
        print(f'[{len(results)+1}/{len(cdfs)}] {cdf.relative_to(ns.data_root)} ... ',
              end='', flush=True)
        r = process_case(cdf, ns.report_dir, ns.python, not ns.quick,
                         ns.strict_inputs, ns.strict_current)
        results.append(r)
        print(r['status'], f"(waves={r['wave_status']}, "
                           f"PF coils={r['pf_active_coils']}, "
                           f"passive loops={r['pf_passive_loops']}, "
                           f"warnings={len(r['warnings'])})", flush=True)
    report = {'cases_total':len(cdfs),'cases_checked':len(results),
              'failures':sum(r['status']=='FAIL' for r in results),
              'warnings':sum(r['status']=='WARN' for r in results),
              'cases':results}
    (ns.report_dir / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    markdown = ["# TRANSP source-verification inventory", '',
                f"Checked {len(results)} of {len(cdfs)} discovered CDFs.", '',
                '| CDF | Inputs | Waves | PF active | PF passive | Result |',
                '|---|---|---|---:|---:|---|']
    for r in results:
        markdown.append(f"| `{r['cdf']}` | {'present' if r['trdat'] else '**missing**'} "
                        f"| {r['wave_status']} | {r['pf_active_coils']} "
                        f"| {r['pf_passive_loops']} | **{r['status']}** |")
    markdown.extend(['', '## Outstanding issues', ''])
    for r in results:
        for issue in r['warnings'] + r['failures']:
            markdown.append(f"- `{r['case']}`: {issue}")
    (ns.report_dir / 'summary.md').write_text('\n'.join(markdown) + '\n')
    print(f"Verified {len(results)}/{len(cdfs)}; FAIL={report['failures']}; "
          f"WARN={report['warnings']}; reports: {ns.report_dir}")
    return 1 if report['failures'] else 0


if __name__ == '__main__':
    sys.exit(main())
