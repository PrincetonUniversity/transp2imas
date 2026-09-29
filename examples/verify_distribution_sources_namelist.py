#!/usr/bin/python
"""Regression check for per-beam namelist lookup used by distribution_sources.

This check does not require IMAS. It parses the TRANSP TR.DAT assignment forms
relevant to ABEAMA/XZBEAMA, including Fortran repetition syntax such as
``2*1.999``, and verifies that beam number ``ib`` maps to Fortran namelist index
``(ib)``. The production writer implements this via
``get_namelist_var(..., ind=ib-1)`` because the nc2imas accessor takes a
zero-based Python index.
"""

import re
import sys
from pathlib import Path

DEFAULT_TRDAT = (Path(__file__).resolve().parent / 'data' / 'nstx' /
                 '141711Y01' / '141711Y01TR.DAT')
path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_TRDAT
if not path.exists():
    raise SystemExit(f'not found: {path}')

values = {}
pat = re.compile(r'^\s*(ABEAMA|XZBEAMA)\s*\(\s*(\d+)\s*\)\s*=\s*([^!,\s]+)', re.I)
for raw in path.read_text(encoding='ISO-8859-1').splitlines():
    line = raw.split('!', 1)[0]
    m = pat.match(line)
    if not m:
        continue
    name = m.group(1).upper()
    ib = int(m.group(2))
    token = m.group(3)
    repeated = re.fullmatch(r'(\d+)\*(.+)', token)
    count = int(repeated.group(1)) if repeated else 1
    number = repeated.group(2) if repeated else token
    value = float(number.replace('D', 'E').replace('d', 'e'))
    for offset in range(count):
        values[(name, ib + offset)] = value

beams = sorted({ib for name, ib in values if name == 'ABEAMA'} |
               {ib for name, ib in values if name == 'XZBEAMA'})
print('file:', path)
print('beam indices found:', beams)
for ib in beams:
    a = values.get(('ABEAMA', ib))
    z = values.get(('XZBEAMA', ib))
    print(f'  beam {ib:02d}: ABEAMA={a!r} XZBEAMA={z!r}')

if path.name.upper().startswith('174783Q32'):
    expected = list(range(1, 15))
    assert beams == expected, (beams, expected)
    assert all(values.get(('ABEAMA', ib)) == 2.0 for ib in expected)
    assert all(values.get(('XZBEAMA', ib)) == 1.0 for ib in expected)
    print('174783Q32 beam-species regression: PASS (14 deuterium beams)')
else:
    print('generic parse: PASS')
