# End-to-end TRANSP CDF to IMAS regression

The normal `make check` and `make verify-all` targets exercise **source-side**
CDF and mapping physics. They **do not** create or reopen IMAS database entries.
For DD validation of actual saved data, run the new **opt-in** end-to-end test on
an IMAS-equipped machine.

## Standard full validation

Activate the same IMAS-Python environment used for ordinary conversions, then
from the repository root run:

```bash
unset IMAS_AL_DISABLE_VALIDATE
make verify-dd PYTHON=python3
make verify-end-to-end PYTHON=python3 DD_VERSION=4.1.1 BACKEND=HDF5
```

This discovers **every existing** `examples/data/**/*.CDF` recursively. The
number of cases is whatever is installed (it is printed, and recorded as
`examples_discovered` in `summary.json`): large CDFs that are supplied locally,
such as NSTX 142202U01 and 120967E03, are included when present and absent from
a plain GitHub checkout, and companion CDFs such as `120967E03PH.CDF` are listed
but not converted. Tests do not silently invent or silently omit
cases. Conversion failures are reported as `FAIL`, even if some earlier IDSs
were already saved and can still be read back.

The end-to-end verifier uses its own temporary HDF5 database root, **not**
`~/public/imasdb`, and removes that test database after completion. It writes
persistent reports and converter logs under `build/end-to-end/`:

- `summary.json`: machine-readable aggregate results, individual writer
  outcomes and actual IDS readback outcomes for each CDF;
- `summary.md`: compact per-case status report;
- `<relative-case>-conversion.json`: conversion write journal;
- `<relative-case>-conversion.log`: stdout/stderr of the actual conversion.

`WROTE` means `DBEntry.put()` returned successfully, not merely that a writer
was called. `SKIP` means the writer made no `put()` (e.g. no supported EC/IC
waves). `FAIL` means conversion, schema validation or readback failed. A
successful first IDS is independently read back even when a later writer fails.
All time-grid failures remain failures; the runner refuses to proceed if the
environment sets `IMAS_AL_DISABLE_VALIDATE=1`.

## NumPy metadata in the conversion journal

TRANSP CDFs often expose `shot` and `run` as NumPy `int32` or `int64`, which
Python's default `json.dumps` cannot serialize.  The conversion journal
normalizes its database identity and accepted IDS occurrences to native Python
numbers.  A defensive JSON encoder also handles other NumPy scalar/array
values without converting them to strings.  Regression tests include the
`TypeError: Object of type int32 is not JSON serializable` failure previously
seen after otherwise successful IMAS writes in all 10 supplied examples.

If this error appears in an old checkout, update `conversion_plan.py` before
rerunning.  The original failure occurred **after** the IDS writers; no
schema/physics regression is implied by that traceback.  Nevertheless, an
actual IMAS readback is still required after applying the fix.

## Testing a single problematic run

```bash
python3 examples/verify_end_to_end.py \
    --case 147634T13 --dd-version 4.1.1 \
    --ids workflow,core_profiles,equilibrium \
    --report-dir build/e2e-t13
```

The reader opens each IDS recorded by the write journal and calls
`ids.validate()`. In particular, the `equilibrium.time` and
`core_profiles.time` arrays must match the **actual `TIME3` values** in the
CDF, and the number of equilibrium/core-profile slices must match `TIME3`.
The CDF's `TIME` scalar-history grid is **not** substituted for `TIME3`, even
when the two arrays happen to have the same length. Selected LCFS outline
and core-profile electron-temperature/radial-grid shapes are also checked.

This is a **schema, readback, and native time-coordinate regression**. It is
not yet a complete independent energy/current/particle conservation check; the
source-side power/current closures remain covered by `make verify-all` and
`make verify-waves`.

## Preserving test databases (optional)

By default, IMAS test entries are deleted after verification. For debugging:

```bash
python3 examples/verify_end_to_end.py \
    --case 38530P94 --dd-version 4.1.1 \
    --database-dir build/e2e-database --report-dir build/e2e-iter
```

**Use a new, dedicated directory.** The runner invokes the converter with
`-new` and must never target an existing production IMAS database. The
`--database-dir` argument disables automatic cleanup. Run with
`--fail-fast` to stop at the first error or `--timeout 3600` for longer cases.

To see which source cases *would* be processed without IMAS:

```bash
make verify-end-to-end-plan
```

This mode explicitly reports **`NOT_RUN`**; it does not count as validation.
In an environment with no IMAS-Python, normal end-to-end mode exits with a
configuration error rather than claiming a pass.

## Debugging with `--ids`

Ordinary `convert.py` still processes all supported IDSs by default. To focus
on a subset (executed in normal writer order):

```bash
python3 convert.py -i examples/data/d3d/147634T13/147634T13.CDF \
    -b HDF5 -new --dd-version 4.1.1 \
    --ids workflow,core_profiles,equilibrium \
    --report-json build/t13-conversion.json
```

Use the actual IMAS IDS names `lh_antennas` and `ic_antennas`; the older
`lh_antenna` and `ic_antenna` spellings are accepted as aliases.

The JSON distinguishes writers that did not output an IDS from successfully
stored ones. A JSON report is also saved if a writer raises an exception.
The option deliberately does **not** create other IDSs implicitly. Include
`workflow` explicitly when provenance should be written in a selective run.

## CI and Data Dictionary support

`examples/verify_end_to_end_contract.py` verifies selection, safe write
journaling, schema/time-grid validation logic, and mocked process/readback
failure handling in source-only CI. It runs without an installed IMAS package.
The real conversion/readback suite **requires** `imas-python` and the requested
Data Dictionary, and is opt-in; it must not be labeled PASS on GitHub-hosted
source-only runners. After DD 4.1.1 is green, repeat using
`DD_VERSION=3.39.0` (if installed), and audit DD-specific unsupported IDSs
separately. An absent optional IDS is SKIP only if no write was attempted.


## Sawtooth readback

The `sawteeth` writer (C-Mod `89846`) uses the native CDF `TIME` grid, not
`TIME3`. The end-to-end verifier compares stored `sawteeth.time`,
`diagnostics.previous_crash_time`, any inferred `previous_period`, and valid
`rho_tor_norm_q1` against the original source, and asserts that no trigger
mechanism was fabricated. To focus on the new IDS:

```bash
python3 examples/verify_end_to_end.py --case 89846 --ids sawteeth \
  --dd-version 4.1.1 --report-dir build/e2e-sawteeth
```
