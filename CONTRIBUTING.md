# Contributing to transp2imas

Thank you for improving the TRANSP-to-IMAS converter. Physics correctness,
traceability, and reproducibility are more important than maximizing the number
of populated fields.

## Before changing a mapping

Read:

- `docs/mapping-workflow.md`
- `docs/technical-reference.md`
- `docs/testing.md`

A new mapping should identify the TRANSP meaning, units, grid, sign convention,
and exact IMAS destination and should include an independent validation where
possible.

## Repository conventions

- Every `*.py` file begins with exactly `#!/usr/bin/python`.
- Every Python module has a module docstring.
- New public behavior is documented in the same pull request.
- Pull-request titles should follow Conventional Commits (`fix:`, `feat:`,
  `docs:`, `test:`, etc.) because releases are automated with Release Please.

## Before opening a pull request

```bash
make check
make verify-all
```

Run `make verify-dd` in an IMAS environment for DD-sensitive changes.

The pull-request template contains the complete checklist.
