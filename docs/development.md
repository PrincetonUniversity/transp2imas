# Development guide

## Supported Python

The converter is Python 3 code. Repository policy nevertheless uses the exact
first line

```text
#!/usr/bin/python
```

for every tracked `*.py` file, including examples and verification scripts. The
header is enforced by `make check` and CI. On systems where `/usr/bin/python`
does not exist, invoke scripts explicitly with the active Python 3 interpreter
(for example `python3 convert.py`) or provide the site-appropriate `python`
alternative. Do not change individual files to a different shebang without
changing the repository-wide policy.

Every Python file should also contain a module-level docstring describing its
purpose.

## Local environment

A convenient local setup is:

```bash
make venv
make install
source .venv/bin/activate
```

For actual IDS creation, the environment must also provide a compatible
IMAS-Python installation and the requested Data Dictionary. Source-side checks
do not require IMAS.

## Common commands

```bash
make help
make check
make verify-all
make verify-dd
make test
```

`make check` is intended to be fast and should be run before every commit.
`make verify-all` exercises the bundled source regression cases. `make
verify-dd` requires an IMAS environment and should be run before merging changes
that touch schema-dependent code.

## Coding conventions

- Keep one physical IDS group per writer module.
- Put shared conversion logic in `nc2imas.py` rather than duplicating formulas.
- Preserve source units in comments/docstrings when the conversion is not
  obvious.
- Explicitly document sign conventions and radial grids.
- Prefer named helpers over repeated attribute-probing code.
- Treat missing optional source quantities as a supported condition.
- Do not place a TRANSP quantity into an IMAS node solely because the names look
  similar; establish physical semantics first.
- Add a source-side numerical check whenever TRANSP contains an independent
  identity, integral, component sum, or redundant representation that can
  validate the mapping.
- Avoid wildcard imports in production entry points; import writer classes and shared helpers explicitly.
- Do not change NumPy floating-point behavior process-wide with `np.seterr`; use a narrowly scoped `np.errstate` or safe `np.divide(..., where=...)` around converter-owned arithmetic.
- Do not use bare `except:` in converter code. Catch the expected exception types, keep the guarded block small, and log enough context to diagnose malformed optional data. A broader `except Exception` is acceptable only for deliberately isolated compatibility/optional-data fallbacks where the failure is logged or converted to a documented absence.
- `make check` enforces the bare-exception rule with `examples/verify_exception_policy.py`, so a new unqualified `except:` will fail CI.
- Resource-owning classes should support explicit close/context-manager usage. `netcdf_connection.__del__` is only a best-effort fallback; normal conversion closes the IMAS entry and input CDF deterministically.
- Dictionary literals used for physics/atomic metadata must not contain duplicate keys; Python silently retains only the last value.

## Git workflow

Development should normally occur on a branch and enter `main` through a pull
request. Pull-request titles/merge commits should follow Conventional Commits
because Release Please derives the next semantic version from them.

Examples:

```text
fix: preserve populated core_sources entries
feat: add full fast-ion distributions mapping
docs: document passive-structure requirements
test: add RF deposition regression case
```

Use `feat!:` or a `BREAKING CHANGE:` footer only for an intentional incompatible
change.

## Required change checklist

Before opening a pull request:

1. run `make check`;
2. run `make verify-all`;
3. run `make verify-dd` when an IMAS environment is available and the change
   touches DD-sensitive fields;
4. update the relevant mapping documentation;
5. update/add a verifier for new physics mappings;
6. confirm no generated database output or local virtual environment is staged.
