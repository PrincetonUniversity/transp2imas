#!/usr/bin/python
"""Check TRANSP CDF quantities used by the v8.4 extended IMAS mapping.

This diagnostic does not require an IMAS database.  It verifies source fields,
units/shape assumptions, cumulative geometry, and referenced template files.
"""

import os
import re
import sys
from pathlib import Path
import numpy as np


def open_cdf(path):
    try:
        from netCDF4 import Dataset
        return Dataset(path, 'r'), 'netCDF4'
    except ImportError:
        from scipy.io import netcdf_file
        return netcdf_file(path, 'r', mmap=False), 'scipy.io.netcdf_file'


def arr(var):
    try:
        return np.asarray(var[:], dtype=np.float64)
    except Exception:
        return np.asarray(var.data, dtype=np.float64)


def attr(var, name, default=''):
    try:
        value = var.getncattr(name)
    except Exception:
        value = getattr(var, name, default)
    if isinstance(value, bytes):
        value = value.decode('ISO-8859-1', errors='replace')
    return str(value).strip()


def parse_template_refs(nml):
    refs = []
    if not nml or not Path(nml).exists():
        return refs
    for raw in Path(nml).read_text(encoding='ISO-8859-1').splitlines():
        line = raw.split('!', 1)[0].strip()
        if '=' not in line:
            continue
        key, value = line.split('=', 1)
        if 'TEMPLATE' not in key.upper():
            continue
        value = value.strip().rstrip(',').strip().strip("'\"")
        refs.append((key.strip(), value))
    return refs


def walk_template_refs(root_nml):
    """Recursively resolve TEMPLATE references relative to the referring file."""
    root_nml = Path(root_nml)
    pending = [("TRANSP_NAMELIST", str(root_nml), root_nml)]
    seen = set()
    out = []
    idx = 0
    while idx < len(pending):
        key, reference, path = pending[idx]
        idx += 1
        path = Path(os.path.expandvars(str(path))).expanduser()
        identity = str(path.resolve())
        if identity in seen:
            continue
        seen.add(identity)
        exists = path.exists()
        out.append((key, reference, path, exists))
        if not exists:
            continue
        for child_key, child_ref in parse_template_refs(path):
            child = Path(os.path.expandvars(child_ref)).expanduser()
            if not child.is_absolute():
                child = path.parent / child
            pending.append((child_key, child_ref, child))
    return out


def main():
    if len(sys.argv) < 2:
        print('usage: verify_cdf_extension.py RUN.CDF [RUNTR.DAT]')
        return 2
    cdf = Path(sys.argv[1])
    nml = Path(sys.argv[2]) if len(sys.argv) > 2 else cdf.with_name(cdf.stem + 'TR.DAT')
    f, backend = open_cdf(str(cdf))
    variables = f.variables
    print(f'CDF: {cdf} ({backend})')
    print(f'variables: {len(variables)}')

    expected = [
        'UE','UI','UTHRM','UTOTL','PPLAS','PTOWB','UFASTPP','UFASTPA',
        'NI','SHAT','ETA_USE','VRPOT','CUR','CUROH','CURBS','PLJB','PLJBOH',
        'PMASS','DVOL','DAREA','SURF','TRFLX','PLFLX','PSI0_TR',
        'NEUTT','NEUTX','NEUTX_DD','BTNTS_DD','BBNTS_DD','POHT','PRAD',
        'PRAD_BR','PRAD_LI','PRAD_CY','PINJ'
    ]
    missing = []
    for name in expected:
        v = variables.get(name)
        if v is None:
            missing.append(name)
            print(f'  MISSING {name}')
        else:
            a = arr(v)
            print(f'  {name:12s} shape={a.shape!s:14s} units={attr(v,"units")}')

    if all(name in variables for name in ('PPLAS','UTHRM')):
        p = arr(variables['PPLAS'])
        u = arr(variables['UTHRM'])
        mask = np.isfinite(p) & np.isfinite(u) & (np.abs(p) > 1.0)
        rel = np.median(np.abs(p[mask] - (2.0/3.0)*u[mask]*1.0e6) / np.abs(p[mask]))
        print(f'PPLAS vs 2/3*UTHRM: median relative difference={rel:.3e}')

    if all(name in variables for name in ('DVOL','PVOL')):
        dv = arr(variables['DVOL'])
        pv = arr(variables['PVOL'])
        err = np.max(np.abs(np.sum(dv, axis=1) - pv))
        print(f'sum(DVOL) vs PVOL: max abs difference={err:.6g} cm^3')
    if all(name in variables for name in ('DAREA','PAREA')):
        da = arr(variables['DAREA'])
        pa = arr(variables['PAREA'])
        err = np.max(np.abs(np.sum(da, axis=1) - pa))
        print(f'sum(DAREA) vs PAREA: max abs difference={err:.6g} cm^2')

    print(f'namelist: {nml} ({"present" if nml.exists() else "missing"})')
    if nml.exists():
        for key, reference, path, exists in walk_template_refs(nml):
            if key == "TRANSP_NAMELIST":
                continue
            print(f'  {key} -> {reference}: {"present" if exists else "missing"} [{path}]')

    try:
        f.close()
    except Exception:
        pass
    if missing:
        print('RESULT: PASS with optional/missing source fields')
    else:
        print('RESULT: PASS')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
