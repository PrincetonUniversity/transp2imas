#!/usr/bin/python
"""Static regression check for writer calls into nc2imas.

This check does not import IMAS.  It parses the converter writers and verifies
that every ``self.im.some_method(...)`` call names a method implemented by
``nc2imas`` or one of its local base classes.
"""

import ast
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
WRITERS = sorted(ROOT.glob('imas_*.py'))
IMPLEMENTATIONS = [ROOT / 'nc2imas.py', ROOT / 'imas_class.py', ROOT / 'netcdf_class.py']

methods = set()
for path in IMPLEMENTATIONS:
    tree = ast.parse(path.read_text(), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            methods.add(node.name)

missing = []
for path in WRITERS:
    tree = ast.parse(path.read_text(), filename=str(path))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        func = node.func
        value = func.value
        if (isinstance(value, ast.Attribute) and isinstance(value.value, ast.Name)
                and value.value.id == 'self' and value.attr == 'im'):
            if func.attr not in methods:
                missing.append((path.name, node.lineno, func.attr))

if missing:
    print('Writer API check: FAIL')
    for path, line, name in missing:
        print(f'  {path}:{line}: self.im.{name}() is not implemented')
    raise SystemExit(1)

print(f'Writer API check: PASS ({len(WRITERS)} writer modules, {len(methods)} available methods)')
