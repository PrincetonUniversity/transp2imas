# Required versus optional mappings

A TRANSP CDF legitimately lacks many variables, and the converter then skips the matching IMAS
fields. A variable that **is** in the CDF but cannot be mapped (wrong radial size, unreadable, an
index or interpolation that does not fit, a schema mismatch) is a different thing: it is a defect in
the mapping or the data and must not look like "optional data absent". Earlier, both ended up as an
empty array or a log line (`get_ncvar` returned `[]` for any failure, and `fill_data` reported a size
mismatch only for time slice 0 and skipped the other slices silently).

## Policy

| Situation | Result |
|---|---|
| CDF variable not present | optional input: the field is skipped and the variable is recorded as `missing_optional` |
| Source data empty or invalid for a field | the field is skipped and recorded as `skipped` (`no source data`) |
| Variable present but cannot be read, indexed, interpolated, shaped or assigned | **broken mapping**: `MappingError` naming the writer, the CDF variable and the IDS target |
| IDS missing from the target DD (`IDSNameError`) | recorded as `unavailable_in_dd` |
| A field written after a documented adjustment of the source | recorded as `adjusted` with the reason (the LCFS contour level clamped into the PSIRZ range, see below) |
| A writer that writes nothing | always has a `skipped` entry with its reason |

`MappingError` aborts the conversion (the journal status is `FAIL` and records the message). With
`--lenient-mapping` the same failures are recorded as `skipped` fields with the exception text as
the reason instead, so a messy CDF can still be converted and nothing is hidden in either mode.
Best-effort provenance (file hashes, parameter files) never aborts a conversion but records what it
could not do.

## Where it is reported

`--report-json` writes a `mapping` object for every writer, and the journal has
`mapping_policy` (`strict` or `lenient`):

```json
{"ids": "core_profiles", "status": "WROTE",
 "mapping": {"missing_optional": {"count": 26, "sources": ["BDENS_3", "BDENS_4", "..."], "truncated": false},
             "skipped": {"count": 3, "fields": [{"target": "profiles_1d/ion/density_fast",
                                                 "reason": "no source data", "count": 12}]}}}
```

Counts are always complete; the listings are capped at 200 entries. The ledger is
`mapping_status.MappingLedger` (no IMAS or NetCDF dependency) and the converter exposes it as
`nci.mapping`.

## What changed in the code

* `netcdf_class.get_ncvar`: absent variable -> `missing_optional`; any other failure -> `MappingError`.
* `nc2imas.fill_data`: returns whether it wrote; absent data -> recorded skip; size mismatch (every
  time slice), scalar-for-array, shape or schema error -> `MappingError`.
* `nc2imas.check_profile_shape`: a `TIME3 x radial` profile of another radial size is a named broken mapping
  instead of a later numpy broadcast error or a silent drop.  `core_transport`, `summary`, `core_sources`, `waves`,
  `edge_sources`, `radiation` and `distribution_sources` apply it to every radial profile they read, and
  `core_profiles` and `equilibrium` to the radial profiles they combine before filling (`IdsWriter.profile`).
  `core_transport`, `core_sources` and `radiation` also pass `chk=` to `fill_data`, and so do the profile fills of
  `core_profiles` and `equilibrium` (the 2-D psi grid and the boundary scalars are not radial and are not checked
  this way); there the error names the field being filled (for example `zeff`, `elong`) and its IDS path.  A central
  check in `get_ncvar` was tried and removed: whether a variable is a radial profile is the reader's knowledge (a
  dimension named `XGRID_NPHI` is not radial), and in a real CDF a profile on `X` or `XB` cannot differ in length from `X`.
  Measured with one corrupted profile at a time under the strict policy on DIII-D 147634B13 (before -> after, of the
  radial profiles each writer read): `core_transport` 35 -> 58 of 59 fail naming the CDF variable (the other is `XB`
  itself); `summary` 0 -> 18 of 20 (the others are `X` and `XB`); `core_profiles` 12 -> 35 of 46, 10 still fail at fill
  time naming the field and IDS path, and `ZIMPS_TOK` (read, never written) is the one silent case; `equilibrium` 0 -> 29
  of 31 (the others are `X` and `XB`).  Before, 8 profiles (`UTHRM`, `UTOTL`, `CUROH`, `PTOWB`, `DRAV`, `DRAVFAC`, `GXI`,
  `TRIANG`) were skipped by a shape comparison without any record.
* 37 `except Exception` blocks were reduced to 8. Derivations that hid failures (parallel current,
  `rho_tor`, `rho_pol_norm`, `psi_norm`, beta, poloidal-field energy, thermal pressure, ...) now go
  through `mapping.broken`; metadata look-ups catch only what they can raise (`atom_base` raises
  `KeyError`); an absent vacuum field is an explicit recorded skip, never an exception. The 8 left are
  the top-level wrap-and-re-raise handlers (`convert`, `imas_class`) and best-effort provenance.
* Running all 12 bundled cases with the old code showed what the swallowed failures actually were:
  none of the CDF variables was present-but-unreadable and no size mismatch occurred, so strict mode
  changes no bundled conversion. The one real swallowed exception was `boundary_minor` on an empty
  boundary (397 time slices of CMOD 89846), which is now a recorded skip.

Tests: `examples/verify_mapping_status.py`.

## Measured reads and the mapping registry

The old static audit (`examples/audit_mapping_coverage.py`) guesses coverage from string constants in the
source, so a name that merely appears in a writer counts as covered. Two measured replacements sit next to it.

**Measured reads.** `netcdf_class` wraps the CDF dataset in `mapping_status.RecordingDataset`. Item access and
`get` on `variables` record the variable name against the running writer; membership tests, iteration and
`keys()` do not. The journal lists them per writer (`consumed`) and before any writer (`setup_reads`). Each writer's record also has `reads`: how often, not only which. `total` reads, `variables` (distinct), and `most_repeated`, the variables read more than ten times (up to five, most first). A variable read hundreds of times is a read inside a per-slice loop, which is what a caching decision needs to see.
A variable being read is not an export: for example `core_profiles` reads `VPOLD_NC`/`VTORD_NC` and writes
them only when they sit on the profile grid.

**Registry.** `mapping_registry.yaml` states each mapping: kind (`direct_typed`, `derived_typed`, `xml_only`,
`input_consumed`, `deferred`), writer, source variables or patterns, destination per DD major version, units,
transform, sign, verification level (`readback_verified`, `writer_verified`, `unverified`) and the tests that
pin it. `docs/mapping-registry.md` is generated from it. `examples/verify_mapping_registry.py` checks:

* schema, unique ids, that named tests and writers exist, and that the generated doc is current,
* every destination path exists in the Data Dictionary version it is registered for (DD 3.39.0 and 4.1.1),
* every typed source pattern matches a variable in a bundled CDF,
* with `--journals DIR`: a writer that wrote an IDS read at least one present source of each of its typed
  entries, and a deferred variable is read only by the writers its `read_by` allows.

`examples/report_coverage.py` (`make -C examples coverage`) classifies every CDF variable as registered typed,
registered input, deferred, read-but-unregistered or unknown, and reports three separate figures: all-source
reads, reads over the applicable set (numeric, time-dependent, not identically zero) and typed registry
entries over that set. The registry is a seed: it covers the mappings verified in this branch, so most read
variables are still "read, unregistered" and are not claimed as verified exports.

Finding from the first measured run: the audit listed `VPOLD_NC`/`VTORD_NC` as deferred while `core_profiles`
did read them; the registry now records both facts.

### Second case: DIII-D 147634T16

`verify_registry_values.py` runs the same tests on two cases. ITER 38530P94 is strict: a missing variable or IDS
is a failure, and every registered mapping is verified there. DIII-D 147634T16 is the bundled MMM07 reference
case; it runs the same tests where its data allow and skips, with the reason, what it does not contain (no NBI,
EC or IC data, no fast-ion density). On it, DD 4.1.1 verifies 87 tests and skips 19. It is also the only bundled
case that verifies `ct_mmm07_diffusivities` (all five MMM07 models: typed coefficients and auxiliary channels)
and `ct_theory_auxiliary_xml`; ITER verifies the same XML check for its GLF23 `DIFFIGLF`. Findings from adding it:
the aggregate impurity density comes from `NIMP` because `NIMPS_TOK` is an all-zero placeholder; impurity ions carry
`TX`, which differs from `TI` here; all-zero GLF23 channels give no `glf23` model; and the 2-D flux is clamped just
inside the transport-edge flux (`imas_equilibrium.trace_boundary` then traces the edge of the clamped region and records an `adjusted` entry), so the outline is known only to about one flux-map grid cell there, and the outline-versus-`boundary/psi`
check uses 1e-6 of the flux span on ITER and 1e-2 on this case.

### Verification levels

`readback_verified` means the value was written to HDF5, read back and compared with an independent
computation from the CDF; `writer_verified` and `unverified` are weaker or absent. The validator requires
that a verified entry lists test scripts that mention its id. `examples/verify_registry_values.py` is the
bulk readback test (ITER 38530P94, DD 4.1.1 and 3.39.0); a registry entry and a `@covers` test must agree
both ways. Signs (`ip_sign`, `q_sign`, `toroidal_rotation_sign`) are tested as a global sign with exact
magnitude; the sign choice itself belongs to the COCOS verifiers. The only entry still `unverified` is the EC frequency taken from the namelist (`wv_ec_frequency`). `pf_active`/`pf_passive` are not registered: they need a user-supplied geometry manifest and their tests use synthetic fixtures.
