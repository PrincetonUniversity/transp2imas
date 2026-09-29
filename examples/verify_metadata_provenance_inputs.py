#!/usr/bin/python
"""Check source-file metadata used by v8.6.2 provenance without IMAS."""

from pathlib import Path
import hashlib
import re
import sys
from datetime import datetime, timezone
from scipy.io import netcdf_file

CASE_DIR = Path(__file__).resolve().parent / 'data' / 'nstx' / '141711Y01'
cdf = Path(sys.argv[1]) if len(sys.argv) > 1 else CASE_DIR / '141711Y01.CDF'
trdat = Path(sys.argv[2]) if len(sys.argv) > 2 else CASE_DIR / '141711Y01TR.DAT'


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def timestamp(path):
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')

print('CDF:', cdf.name, 'sha256=', sha256(cdf), 'timestamp=', timestamp(cdf))
print('TR.DAT:', trdat.name, 'sha256=', sha256(trdat), 'timestamp=', timestamp(trdat))

text = trdat.read_text(errors='replace')
templates = []
for match in re.finditer(r'(?im)^\s*([A-Z0-9_]*TEMPLATE)\s*=\s*[\'\"]?([^\'\"!,\n]+)', text):
    key, raw = match.group(1), match.group(2).strip()
    p = Path(raw).expanduser()
    if not p.is_absolute():
        p = trdat.parent / p
    templates.append((key, raw, p.exists()))
print('template references:')
for key, raw, exists in templates:
    print(f'  {key}={raw!r}:', 'present' if exists else 'missing')

f = netcdf_file(str(cdf), 'r', mmap=False)
for key in ('TRANSP_version', 'TRANSP_DOI', 'TRANSP_BUILD_DATE', 'TGLF_hash', 'GLF23_hash'):
    value = getattr(f, key, '')
    if isinstance(value, bytes):
        value = value.decode(errors='replace').strip()
    print(f'CDF attribute {key}: {value!r}')
for name in ('PPHI', 'AMOI', 'OMEGA', 'AMOM'):
    if name in f.variables:
        v = f.variables[name]
        print(f'{name}: dimensions={v.dimensions}, units={getattr(v, "units", b"")!r}')
f.close()
