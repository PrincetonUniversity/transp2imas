#!/usr/bin/python
"""Fail when Python source contains duplicate literal dictionary keys.

Python silently keeps only the last occurrence of a duplicate key in a dict
literal.  That behavior is particularly dangerous for atomic/mapping tables,
so this source-only verifier catches the problem before conversion tests run.
"""

from __future__ import annotations

import ast
from pathlib import Path
import sys


def duplicate_keys(path: Path):
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Dict):
            continue
        seen = {}
        for key in node.keys:
            if not isinstance(key, ast.Constant):
                continue
            try:
                hash(key.value)
            except TypeError:
                continue
            if key.value in seen:
                yield key.value, seen[key.value], key.lineno
            else:
                seen[key.value] = key.lineno


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    failures = []
    for path in sorted(root.rglob("*.py")):
        if any(part in {".venv", "build", "__pycache__"} for part in path.parts):
            continue
        try:
            failures.extend((path.relative_to(root), *item) for item in duplicate_keys(path))
        except SyntaxError as exc:
            print(f"Syntax error while checking {path}: {exc}", file=sys.stderr)
            return 2

    if failures:
        print("Duplicate literal dictionary keys found:", file=sys.stderr)
        for path, key, first_line, duplicate_line in failures:
            print(
                f"  {path}:{duplicate_line}: {key!r} duplicates line {first_line}",
                file=sys.stderr,
            )
        return 1

    print("Python literal dictionary key check: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
