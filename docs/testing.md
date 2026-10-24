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

`make verify-all` automatically discovers **all** bundled `examples/data/**/*.CDF`
files (10 in the complete local fixture set), verifies their grids, EC/IC wave
closures, PF channel mappings and XML input provenance, and runs detailed
ISolver PF regressions for the fixtures actually available. In GitHub CI the
large NSTX `142202U01.CDF` is absent, so `verify-all` reports it as `SKIP`
while still validating ITER and all 9 checked-in CDFs. Use
`make verify-isolver` to **require** both real ISolver CDFs locally. These
source tests do not require IMAS. A per-run Markdown,
JSON and log report is written to `build/verification/`. Missing source input
files and known numerical residuals are WARN, not silently ignored. Use
`make verify-deep` for the much slower legacy numerical diagnostic scripts
across the entire dataset, and `--strict-inputs` / `--strict-current` for
strict fail-on-warning versions. See `examples-verification.md`.

The quick suite also tests the ISolver fixture policy: CI may skip a missing
optional NSTX CDF without discarding ITER coverage, while strict local runs,
present CDFs missing their matching TR.DAT, and zero-fixture runs still fail.

`make verify-pf-passive` runs synthetic PS_n, machine-geometry and DD3/DD4
outline regression tests, and checks that the current real fixtures (none of
which has PS_n) never generate fictitious currents. It requires SciPy/NumPy
but no IMAS installation. It is also included in the fast `make check` path.
**The supplied ITER/NSTX ISolver cases now validate PS_n source channels; physical conductor geometry/signs still need authoritative machine data.**

`make verify-core-profiles` tests the electron and ion temperature-peaking
calculations against every example CDF. In particular, the TORAY 147634T13
case has 4 `TIME` points for `TE0`, `TI0` and `PVOL`, but 7 `TIME3` points
for profile data. The `core_profiles` writer must interpolate `TE0` before
forming `t_e_peaking`; the fast source check enforces that contract and the
full CDF check is also part of `make verify-all`.

`make verify-neoclassical` validates every CDF's `GFLNC_*`, `QFLNC_*`,
`QFLNCC_*` source divergence, unit/time-grid checks and finite-volume
reconstruction. It exercises the actual `core_transport` writer with IMAS-shaped
mocks; results are stored under `build/neoclassical-transport/`.

`make verify-neutrals` independently checks wall, volume and recombination
neutrals, their density SI conversion, already-eV temperature, and signed
rotation. Results are stored under `build/neutral-profiles/`. None of these
neutral quantities is treated as a turbulence-model coefficient.

`make verify-theory` tests unified MMM07, GLF23 and TGLF registration, native
`TIME3/XB` coefficient mapping, MMM07 model decomposition and reduced-model
`TIME3/X` outputs, and the actual source-side core_transport writer. GLF23
nonzero coefficients are tested on ITER 38530P94; DIII-D T16 validates MMM07;
other GLF23/MMM07 placeholder channels do not create false models. TGLF is
reported as awaiting a validated fixture. The per-case report is under
`build/theory-transport/`. `make verify-mmm` is an alias for compatibility,
but there is no separate MMM mapping module. The real IMAS DD4.1.1 write and
readback require `make verify-end-to-end` in an IMAS-equipped environment.
See [theory transport mapping notes](theory-transport.md).

`make verify-waves` runs the RF-specific source regressions against every bundled CDF, including the
DIII-D TORAY/TORBEAM/GENRAY EC, ITER EC+IC, and C-Mod IC cases. It validates component sums, finite-volume power
integrals, absorbed-power traces, and EC driven-current closure without
requiring an IMAS installation.

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

### 4. End-to-end conversion and real readback

`make verify-end-to-end DD_VERSION=4.1.1 BACKEND=HDF5` runs all available
examples in isolated real IMAS HDF5 databases, checks accepted writes and
reopens each stored IDS with `validate()`. It checks the native `TIME3` arrays
and associated profile/slice shape constraints. This requires IMAS-Python and
is deliberately separate from the CI source-only suite. The output includes
per-case JSON and Markdown reports; failures are never downgraded to skips.
See [end-to-end instructions](end-to-end.md).

### Tested matrix

What each check actually runs against (last verified on IMAS-Python 2.1.0.post1, Python 3.12.3, netCDF4 1.7.4, numpy 2.4.1, HDF5 backend):

| check | needs | DD 3.39.0 | DD 4.1.1 | backend |
|---|---|---|---|---|
| `make check` (CI) | numpy, scipy, matplotlib, PyYAML | source side only | source side only | none |
| `make verify-dd` | IMAS-Python | metadata | metadata | none |
| `make verify-readback` (includes `verify_registry_values.py`: 541 tests = 135 per case and DD on ITER 38530P94 and DIII-D 147634T16, plus the registry/test cross-check; about 15 s) | IMAS-Python, netCDF4 | yes | yes | HDF5 |
| `make verify-end-to-end` (every available bundled case, 12 with the local NSTX samples) | IMAS-Python, netCDF4 | `DD_VERSION=3.39.0` | default | HDF5 |
| `make verify-nubeam-fi` | IMAS-Python, netCDF4 (+ the local NSTX sample for the full run) | yes (synthetic and ITER cases) | yes | HDF5 |
| `python examples/verify_mapping_registry.py` | IMAS-Python for the Data Dictionary path check | paths checked | paths checked | none |

Tests that need a variable a case does not have skip with the reason instead of failing: ITER has no recombination
neutrals (`N0V0_*`, `S0RECO`), no lumped single impurity (`NIMPS_TOK`) and no crash events, so those entries are
verified on DIII-D 147634T16 or by the source-side scripts named in the registry.

MDSplus is not tested. The access layer is IMAS-Python's own (`imas.DBEntry`); the legacy access-layer
Python module is not used.

#### Which Data Dictionary a command uses

| command | `DD_VERSION` unset or `auto` | explicit `DD_VERSION=x` |
|---|---|---|
| `make verify-end-to-end` | 4.1.1 (pinned by the Makefile) | x |
| `make convert`, `make test`, `make verify` | `$IMAS_VERSION` or `$IMAS_DD_VERSION` if set, else the IMAS-Python default (3.39.0 here) | x |
| `python convert.py` without `--dd-version` | the same as the row above | x |

The end-to-end check once passed an empty version and silently used the IMAS-Python default; it now always passes an
explicit version, and the journals record it (`dd_version` in `build/end-to-end/*-conversion.json`). `make convert`
still follows the environment, so check `IMAS target DD:` in its output, or pass `DD_VERSION` explicitly, before
comparing results across DD versions.

#### Environment pitfalls

* An old Access Layer on `PYTHONPATH` (for example `AL-3.39.0/python/lib`) shadows IMAS-Python and fails with
  `No module named 'imp'`. Run `unset PYTHONPATH` before activating the IMAS-Python environment (see the README).
  `verify_mapping_registry.py` then prints `SKIP Data Dictionary check` instead of failing, so a registry
  destination typo passes unnoticed: look for that line in the output.
* `IMAS_VERSION` / `IMAS_DD_VERSION` in the environment select the DD for every command that does not pass one.

### 5. Legacy manual conversion

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
- the ITER/NSTX ISolver runs now supplied for `pf_active` and `pf_passive`
  currents; full geometry validation still requires authoritative device files;
- event-trigger mechanism and verified post-crash profiles for enhanced `sawteeth` mapping; C-Mod 89846 now covers native crash-history export.

Do not add large binary fixtures without confirming repository size and data
redistribution constraints.

## Active PF source-only regression

Run `make verify-pf-active` to exercise synthetic `KK_/VV_/KT_/VS_/CC_`
CDF channels, independent `TIME` and `TIME3` coordinates, measured-versus-target
current protection, DD3/DD4 geometry/circuit matrices and no-data skipping.
This regression is also run by `make check` through `examples/quick`. Real ITER/NSTX current-waveform regression is now provided by
`make verify-isolver`; an end-to-end IMAS write and physical-current sign/geometry
verification still require an authoritative device model. See
[`pf-active.md`](pf-active.md).

## Real ISolver PF source regression (optional large fixtures)

`make verify-isolver` checks the ITER 38530P94 (11 active coil circuits, 119
passive currents) and NSTX 142202U01 (9 active circuits, 82 passive currents)
fixtures against the mapping helpers and the native CDF time axes. It also
checks every `CC_* = N_effective * KK_*` closure and all `PFC_NAMES` from
TR.DAT. It writes CSV reports to `build/pf-isolver/`.

The NSTX CDF is roughly 166 MB and is intentionally not required by GitHub
CI: `make verify-all` runs the optional-fixture mode and reports a visible
`SKIP`, while direct `make verify-isolver` remains strict. See
[ISolver fixtures](isolver-fixtures.md) for storage policy and reproducible
commands.

## Equilibrium time-coordinate source regression

`make verify-equilibrium` validates all CDF files in `examples/data` and tests
resampling of `RMCB*/RMSB*/YMCB*/YMSB*` from source `TIME` into `TIME3`
boundary moments.  This test executes the actual boundary writer methods
through a SciPy-backed adapter, requiring no IMAS installation.  The
`equilibrium-time-grid` result is included in `make verify-all`.

See [TIME and TIME3](time-grids.md).  An IMAS DD4.1.1 runtime check must
still be performed on a workstation with IMAS-Python installed.


## Sawtooth history

`make verify-sawteeth` scans **every bundled CDF** for changes in the native
`TLASTSAW(TIME)` history, checks the `TIME` shape independently of `TIME3`,
retains only physically valid `XIQ1` values, and verifies the typed field
contract with mock IDS objects. C-Mod 89846 contains 38 observed new updates;
ITER has constant `TLASTSAW`, and the other cases have no crash history.
`make verify-all` includes a per-case `sawteeth-source` check. Real validation
requires `--ids sawteeth` with `verify_end_to_end.py` on an IMAS host.
See [sawteeth mapping](sawteeth.md).


PT_SOLVER multi-region inputs: `make verify-pt-solver` verifies all bundled
TR.DAT and PT_TEMPLATE inventories, per-equation X*BOUND ranges, scheduled
namelist updates and references to the single assembled core_transport IDS.
It is independent of IMAS-Python; use a real DD 4.1.1 conversion/readback for
final schema validation. See [pt-solver-regions.md](pt-solver-regions.md).
