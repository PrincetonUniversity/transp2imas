# Testing and validation

## Validation layers

`transp2imas` uses several complementary levels of testing because a successful
Python run alone is not sufficient to establish a correct physics mapping.

### 1. Header and syntax checks

`make check` first verifies that every tracked Python source file starts with
`#!/usr/bin/python`, then compiles the Python tree. It also runs the quick
source-verification suite under `examples/`.

It also runs `examples/verify_python_literals.py` and `examples/verify_exception_policy.py`. The latter rejects bare `except:` handlers so `KeyboardInterrupt` and `SystemExit` are never swallowed by accidental catch-all blocks.

This catches repository-policy drift, syntax errors, and calls to helper methods
that do not exist.

### 2. Source-side physics checks

`make verify-all` runs independent checks against each bundled TRANSP regression
case. These tests operate on the CDF/TR.DAT source and do not require IMAS.

Examples of useful validation patterns include:

- profile/component closure;
- volume or area integral reconstruction;
- redundant TRANSP representations of the same quantity;
- expected finiteness and dimensionality;
- radial-grid compatibility;
- sign relationships implied by the TRANSP input controls.


### Mapping-coverage inventory

`make audit CASE=<runid>` runs `examples/audit_mapping_coverage.py` for one CDF
and writes complete Markdown, CSV, and JSON reports to `build/audit/`. This is a
development inventory rather than a pass/fail regression: a large unclassified
set is expected while coverage is being expanded.

The audit uses the same SciPy NetCDF reader as the portable source verifiers and
does not require IMAS. It performs static source-reference discovery, so its
`mapped/consumed` count is deliberately conservative for dynamically generated
TRANSP variable names. Do not make CI fail on the unclassified count unless the
project later adopts a maintained complete mapping registry.

### 3. Data Dictionary/API checks

`make verify-dd` inspects the installed IMAS object structure. It is intended to
catch field renames, array-shape differences, and DD3/DD4 incompatibilities.
This check requires the target IMAS Data Dictionary installation and therefore
is not part of the generic GitHub-hosted CI job.

### 4. End-to-end conversion

`make test` performs the source checks and converts the active regression cases
using the selected backend/environment. Run it in a configured IMAS environment
before a physics release when possible.

## GitHub CI

`.github/workflows/ci.yml` runs source-only checks on GitHub-hosted Ubuntu
runners for the supported Python matrix. It intentionally avoids tests that
require PPPL/ITER-specific IMAS installations.

A future self-hosted runner can add `make verify-dd` and end-to-end conversion
without weakening the portable source CI.

## Adding regression data

A new regression case should be chosen because it unlocks validation of a
specific previously untestable mapping. Useful future cases include:

- a run with full NUBEAM FBM phase-space output for `distributions`;
- a run with passive-structure currents and authoritative geometry for
  `pf_passive`;
- a run with spatial RF absorption/deposition for `waves`;
- a sawtoothing run with explicit crash-event diagnostics for `sawteeth`.

Do not add large binary fixtures without confirming repository size and data
redistribution constraints.
