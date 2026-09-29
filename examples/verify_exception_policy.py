#!/usr/bin/python
"""Reject bare ``except:`` handlers in project Python source files.

Bare handlers also catch ``KeyboardInterrupt`` and ``SystemExit``.  This static
check keeps converter error handling explicit without requiring third-party
linters on lightweight CI runners.
"""

from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
SKIP_PARTS = {".git", ".venv", "build", "__pycache__"}
PATTERN = re.compile(r"^\s*except\s*:\s*(?:#.*)?$")

problems = []
for path in sorted(ROOT.rglob("*.py")):
    if any(part in SKIP_PARTS for part in path.parts):
        continue
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if PATTERN.match(line):
            problems.append((path.relative_to(ROOT), lineno))

if problems:
    for path, lineno in problems:
        print(f"Bare except forbidden: {path}:{lineno}")
    sys.exit(1)

print("Python exception policy check: PASS (no bare except handlers)")
