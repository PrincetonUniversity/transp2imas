# Changelog

Notable changes to transp2imas are recorded here. Future entries and GitHub
releases are managed automatically by Release Please from Conventional Commit
messages.

- Verify the registry on a second case, DIII-D 147634T16. `examples/verify_registry_values.py` now runs on ITER
  (strict) and on 147634T16 (runs what its data allow, skips the rest with the reason): 425 tests at DD 4.1.1 and
  3.39.0. 147634T16 is the MMM07 reference, so `ct_mmm07_diffusivities` and `ct_theory_auxiliary_xml` are now
  readback-verified (all five MMM07 models); only `wv_ec_frequency` remains unverified. The new test of the 2-D flux
  checks that the boundary outline is a contour of `profiles_2d/psi` at `boundary/psi` (ITER: to 1e-7 of the span).
- Fix the LCFS outline of equilibria whose PSIRZ does not reach the contour level. The rectangular PSIRZ of
  DIII-D 147634T16 is clamped 0.1-0.2% below the transport-edge flux `PLFLXA`, so no contour exists at that level and
  whether a curve was found depended on the sign convention and on the contouring library: DD 4.1.1 got an outline,
  DD 3.39 got a 5-point fragment, rejected, and no outline. `trace_boundary` now uses the nearest level inside the data
  range, which traces the edge of the clamped region (identical outline in both DDs, within 0.14 mm of the earlier
  DD 4.1.1 result), and records it as an `adjusted` entry of the mapping journal (a new status kind next to
  `skipped`). Cases whose level is in range are unchanged.
- Start changing the writer interface so that constructing a writer does not do its work. `writer_base.IdsWriter` has
  `build()`, which returns a list of the IDS objects to write (or `(ids, occurrence)` pairs), and `run()`, which
  submits them; `run_writer` runs migrated writers through `run()` and leaves unmigrated ones alone, so
  `convert.make_writers` wraps every writer in it. All 22 writers are migrated, with unchanged stored output
  (`tf`, `wall`, `lh_antennas`, `ic_antennas`, `nbi`, `ec_launchers`, `sawteeth`, `dataset_fair`, `workflow`,
  `summary`, `transport_solver_numerics`, `pf_active`, `pf_passive`, `radiation`, `distribution_sources`,
  `edge_sources`, `waves`, `core_sources`, `core_profiles`, `equilibrium`, `distributions`, `core_transport`).
  Writers that read grids when constructed now read them at the start of `build()`; `core_transport` returns
  occurrence 0 and the region view occurrences as `(ids, occurrence)` pairs; `equilibrium` keeps its timing around
  `_build()`. A build that fails submits nothing. Tested by `examples/verify_writer_interface.py`, which also
  requires that every writer `convert.py` imports is migrated; the mock tests of the writers now check that
  construction does not write and call `run()`. The `core_profiles`/`equilibrium`/`core_transport` structure tests
  look for `build`/`_build`. The sawteeth, core-profiles and equilibrium time-grid verifiers skip companion CDFs
  without a `TIME` coordinate.
- Extend the radial-shape check to the readers of `waves`, `edge_sources`, `radiation` and `distribution_sources`.
  A profile present in the CDF with a different radial size than the writer's grid used to be dropped silently by
  those readers (indistinguishable from absent data); it is now a broken mapping that names the variable (strict:
  `MappingError`, lenient: a skipped field with the reason), as in `core_sources`. The profile fills of
  `core_profiles` (40) and `equilibrium` (33) now pass the radial size to `fill_data` too. Tested by corrupting one
  variable per writer in a bundled case (`EveryReaderShapeTests` in `examples/verify_mapping_status.py`, six writers).
- Fix documentation drift named in the review: the `pt_solver_regions` header said region settings are not separate
  `core_transport` occurrences although region view occurrences are written (it now describes occurrence 0 as the
  assembled solution and the region occurrences as views with no invented profiles); hard-coded example counts in
  the Makefile help and the docs (10, 9) are replaced by "every bundled example"; and `-new` is documented, in
  `--help` and the README, as creating a fresh entry for the selected database, pulse and run without allocating a
  run number, an existing entry being opened when it is omitted.
- Move the validity checks, array arithmetic and IDS-node assignment primitives out of `nc2imas` into
  `ids_assignment.py` (`is_valid`, `sum_arrays`, `operate_arrays`, `np_multiply_sum`, `check_1dim`, `fill_value`,
  `fill_data`), as functions that take the logger and the `MappingLedger` explicitly. The `nc2imas` methods of the
  same names keep their signatures and delegate. Strict and lenient mapping behaviour is unchanged and now has direct
  unit tests (`examples/verify_ids_assignment.py`, 19 tests, no IMAS needed).
- Move conversion provenance out of `nc2imas` into `code_provenance.py`: the portable source-file records (name with
  SHA-256, no absolute path, capped at 10), writing them into either DD layout, the `ids_properties/comment` text and
  the `code.parameters` XML, all as functions of plain values. About 300 lines leave `nc2imas.py` (1384 to 1089),
  including five unused legacy helpers (`_xml_safe_text`, `_clean_namelist_path`, `_resolve_input_file`,
  `_read_parameter_file`, `_append_namelist_xml`) that `input_provenance.py` had already replaced. The `comment`,
  `code.parameters` and provenance records stored for DD 4.1.1 and 3.39.0 are byte-identical to the previous code.
- Move the COCOS/direction sign policy out of `nc2imas` into `cocos_signs.py` (`resolve_signs` returns an immutable
  `CocosSigns`; `nc2imas` publishes the same attributes and prints the same line). Verified against a hand-written
  table of all direction and DD3/DD4 combinations (`examples/verify_cocos_signs.py`).
- Fix the default Data Dictionary of `make verify-end-to-end` run from `examples/`: with `DD_VERSION` unset it passed an
  empty version and used the IMAS-Python default (DD 3.39), so the 12-case check was not exercising DD 4.1.1. It now
  resolves `auto` to 4.1.1 like the root Makefile; `docs/testing.md` has a table of what each check covers.
- Read the FI_TAR archives of a NUBEAM run in parallel. Each archive is a gzip stream of hundreds of MB that must be
  decompressed to reach the equilibrium members, and they were read one after another. They are now read by up to
  four `spawn` worker processes (a bounded window of at most twice the worker count is held in memory, results are
  consumed in the original order and an error surfaces at its archive's turn); `TRANSP2IMAS_FI_WORKERS=1` restores
  the sequential read. On NSTX 120967E03 the `distributions` writer fell from 40 s to 21 s with identical snapshots.
- Record where conversion time goes and remove the largest cost in our own code. The journal now carries
  `write_seconds` per writer (time in the backend: IMAS validation, put and compression) next to `seconds`; on CMOD
  89846 the backend is about 60% of the run. Radial interpolation (`nc2imas.interp_radial_profiles`) now uses a
  vectorized operator cached by the content of the two coordinate arrays (`radial_interpolation.py`), because the
  `radiation` writer repeated the same X to XB mapping 2385 times with a per-time-slice `np.interp` loop. Grids that
  are not strictly increasing or contain NaN use the original loop. All 32,898 interpolations of every radial
  variable in the bundled CDFs equal the original to 1e-12, and the `radiation` build on CMOD fell from 5.4 s to 1.9 s.
- Start the composition refactor behind the unchanged `nc2imas` interface: Data Dictionary version logic moves to
  `dd_adapter.DDAdapter`, IDS puts become `WriteRequest` records queued by `ImasStore`, and the writer loop moves from
  `convert.py` to `ConversionRunner`, which flushes each writer's IDSs after the writer returned. Constructing a writer
  no longer writes to the database, and a writer that fails writes none of its own output (previously earlier puts of
  the same writer stayed). Outputs of the 12 bundled cases are unchanged.
- Fix mixed poloidal-flux references in `equilibrium`. `profiles_1d/psi` and `psi_magnetic_axis` are absolute
  (axis value `PSI0_TR`), but `boundary/psi` was written from `PBPMHD` and `profiles_2d/psi` from `PSIRZ`, both of
  which are relative to the magnetic axis (ITER: `boundary/psi` +75.4 Wb against -12.9 Wb at the profile edge).
  `boundary/psi` is now the flux of the traced outline (the profile edge, equal to `global_quantities/psi_boundary`),
  `profiles_2d/psi` is shifted by the axis value (the boundary is still traced on the relative field), and the DD3
  `boundary_separatrix/psi` is the absolute bounding-point flux; DD4 no longer reads `PBPMHD`.
- Register and read back the `equilibrium`, `summary`, `distribution_sources`, `waves` and the species-resolved,
  global and remaining `core_profiles` mappings (209 readback tests on ITER 38530P94 for DD 4.1.1 and 3.39.0).
- Grow the mapping registry to the hardware, heating, `core_profiles`, `core_transport` and `edge_sources` mappings and
  read them back: `examples/verify_registry_values.py` converts ITER 38530P94 for DD 4.1.1 and 3.39.0 and compares
  each registered value with an independent numpy computation from the CDF. A registry entry may only claim a
  test that mentions its id, and the test file and registry must list the same ids. Writing it corrected the
  registry's own over-claims (`PRAD_LI` is `emissivity_ion_total`, not electron emissivity; `VNDNC_E` is not an
  ion velocity; fast-ion temperature is the mean energy `EBEAM`/`EFUSN`) and showed the species-resolved
  `SIGF_*_*`/`SIRC_*_*` channels are read but not exported.
- Add measured mapping coverage and a declarative mapping registry. The CDF dataset records which variables each
  writer actually read (`consumed` and `setup_reads` in the `--report-json` journal);
  `mapping_registry.yaml` lists typed, XML-only, input and deliberately deferred mappings with units, transform,
  sign, verification level and tests, is validated against the Data Dictionary, the bundled CDFs and the
  measured reads (`examples/verify_mapping_registry.py`), and renders `docs/mapping-registry.md`.
  `examples/report_coverage.py` reports all-source reads, reads over the applicable set and typed registry
  entries separately, never counting a read as an export.
- Skip `edge_sources` GGD momentum components that the target DD ties to the radial component (DD 3.39
  `toroidal`) instead of failing validation; the old trim-by-resize had been discarding these sources.
- Separate absent optional data from broken required mappings. `get_ncvar` and `fill_data` no longer turn
  read, shape, index, interpolation or schema failures into empty arrays or slice-0-only log lines: they raise
  `MappingError` naming the writer, CDF variable and IDS target (or, with `--lenient-mapping`, record a skipped
  field with the reason). Absent optional variables, skipped fields and IDSs missing from the DD are recorded
  per writer in the `--report-json` journal (`mapping`), every skipped writer carries a reason, and the broad
  `except Exception` handlers were reduced from 37 to 8.
- Attribute regional PT_SOLVER models over time: keep every model selected in a region with its
  `initial`..`end` intervals (a GLF23 -> TGLF switch no longer becomes TGLF for the whole run, and a model
  disabled before the end is no longer lost), write one model slot per selected model without
  placeholders, add `effective_until` to the radial validity states, and stop relabelling the
  assembled neoclassical result with a selected model name (it is named `neoclassical`, with the
  selected models and their intervals as metadata).
- Add a `transport_solver_numerics` writer: PT_SOLVER `NPTRITER` and the `PT_*` timers as named
  control parameters (`newton_iterations`, `cumulative_cpu_time_hours`, `cumulative_wall_time_hours`)
  on the native `TIME` grid, for DD >= 3.42 (DD 3.39.0 has unnamed control parameters, so nothing is
  written); inactive channels, non-hour timers and non-integer counts are skipped.
- Keep the native `TIME` coordinate values with each PT_SOLVER runtime diagnostic channel (`NPTRITER`,
  `PT_CPTIME`, `PT_WTIME`): `TIME` is not the `TIME3` grid of `core_transport` and differs in length in
  some CDFs. Discover the channels independently of the PT_SOLVER input provenance (without a PT plan,
  only channels with nonzero data are recorded), decode their text units instead of storing bytes reprs,
  and run the PT verification scripts in `make quick`.
- Normalize `j_parallel` (`core_sources`, `core_profiles`, `equilibrium`) and `rho_tor` (`core_profiles`,
  `equilibrium`, shared grid helper) by the IDS's own `vacuum_toroidal_field/b0`, as the Data Dictionary
  defines them, instead of TRANSP's `BZ`, which differs from the exported `b0` by 5-38% in the bundled
  cases. Fill `vacuum_toroidal_field` in `equilibrium` (it was empty although its `j_parallel` is
  defined relative to it) and use the same `r0`/`b0` convention in `summary`.
- Fix IDS arrays of structures being erased on write: `radiation`, `distribution_sources` and
  `edge_sources` shrank a filled array with `resize(used)`, which IMAS-Python resets
  (`keep=False`), so they wrote empty entries while reporting success. Entries are now built in
  scratch IDSs and appended (`utilities.append_successful_entries`); regression tests read the
  entries back from HDF5.
- Make `vacuum_toroidal_field` consistent with the TF IDS: `b0` is now `BZXR / r0` with one
  constant reference radius, so `r0 * b0` equals the TF vacuum field times radius (it differed by
  up to 17% when the magnetic axis moved).
- Write the complete NBI toroidal momentum source to `core_sources`
  (`TQBE + TQBI + TQBTH + TQJXB`, plus the fast-ion anomalous-diffusion torque `TQABDF` when
  present) with `TQJXB` as the included JxB subset, instead of `TQBTH` alone, and fill the DD 3.39
  JxB node. The sum reproduces TRANSP's `TQIN` in every bundled case with beams.
- Write the predictive angular-momentum flux (`AMTR_MOD`) to the `transport_solver` model.
- End-to-end verification no longer treats `<run>PH.CDF` companion files as cases or regional
  settings-only `core_transport` models as missing profiles.

- Add TR.DAT-driven discovery of NUBEAM `DATA<n>` and `FI_TAR.GZ<n>` outputs,
  keeping their independent `outtim`/`fi_outtim` counts and excluding CPU
  particle-state files from conversion inputs. Keep the large NSTX case out
  of git. Write the `distributions` IDS (GGD on NUBEAM's MC zone grid x
  energy x pitch, pitch along B) from them, with tests and plotting script.

- Trace the equilibrium LCFS from `PSIRZ` and the enclosed boundary flux by
  default; reserve Fourier moment reconstruction for explicit `--moments yes`
  requests. Populate equilibrium geometric-axis coordinates from boundary and
  flux-surface zeroth moments, with an outline-centroid fallback. Use matching
  `TIME3` zone volumes for core-profile volume averages and preserve each
  impurity charge state's own charge label.

- Preserve PT_SOLVER axial, confinement and edge model selection as
  input-provenance metadata within a single core_transport occurrence.
  Include TR.DAT XIMIN_CONF/XIMAX_CONF, equation-specific XIBOUND/XNBOUND/
  XPHIBOUND, LPREDICT flags and ~update_time history; preserve referenced
  PT_TEMPLATE %active/%xanom/ExB options and model update blocks. Treat
  missing T16 TR.DAT and adjacent templates as unverified, not authoritative.
  Link configured model records without copying radial grids or fabricating
  region-specific profiles. Add verify-pt-solver on all bundled inputs and
  docs/pt-solver-regions.md.

- Refactor theory transport mapping into `theory_transport_mapping.py`:
  MMM07, GLF23 and neoclassical sources share one model-family registry;
  future TGLF requires validated source examples. Add source-verified GLF23
  electron/ion heat and toroidal momentum diffusivities from ITER 38530P94,
  preserving aggregate `DIFFIGLF` without inventing a species assignment.
  Separate neutral profiles into `neutral_profiles_mapping.py` with their own
  validation target, and split neoclassical and neutral verification workflows.
- Explicitly configure Release Please for conventional `feat:` minor and
  `fix:` patch bumps during 0.x releases, add release-config tests and explain
  squash-merge subjects and limits of source-only release verification.


- Implement source-validated `sawteeth` IDS with native `TIME` and typed
  `diagnostics.previous_crash_time`, `previous_period`, and
  `rho_tor_norm_q1`. Detect newly observed changes in C-Mod 89846
  `TLASTSAW` without inventing the initial event, crash-trigger mechanism,
  post-crash profiles or ad-hoc mixing-radius physics. Register `--ids sawteeth`,
  add `make verify-sawteeth`, all-case regressions and real-IMAS readback checks.


- Keep MMM07 `TIME3`, `X`, and `XB` coordinates in native `core_transport`
  fields instead of serializing them repeatedly in model XML. For MMM reduced
  DRIBM/ETG/MTM component outputs, store `X` in the otherwise-unused
  `profiles_1d.grid_flux.rho_tor_norm` (while keeping `XB` in `grid_d`), and
  reference native IDS coordinates from version 2 supplementary-output XML.
  Expand all-example writer and schema checks to verify native coordinate grids.
- Fix doubled paths when an absent top-level `TR.DAT` is supplied relative to
  the invocation working directory, notably DIII-D `147634T16`; preserve the
  genuine missing-file status and the separate unreferenced MMM/PT namelists.

- Map the reduced MMM07 anomalous-transport output to `core_transport` as
  a combined model plus physically distinct Weiland (CDF `W19`), DRIBM/DRBM,
  ETG, and reduced MTM contributions. Preserve source `TIME3/XB` coefficients
  in SI, model-specific `TIME3/X` reduced growth/frequency outputs as XML,
  and audit the component sum against total diffusivities on DIII-D T16.
  Do not create MMM models from zero-filled inactive placeholders. Add
  `make verify-mmm`, all-example source checks, writer regression and DD
  node checks. Clarify that MTM is an MMM **model output**, not a measurement.

- Fix real IMAS-Python DD4 `core_profiles.code.parameters` handling when
  appending neutral rotation XML: read the `IDSString0D.value` text before
  parsing, and always return a Python string even if no rotation channels
  are present. Preserve existing XML input-provenance nodes and add a mock
  `IDSString0D` regression reproducing the runtime `TypeError`.

- Extend `core_transport` with dedicated, species-resolved neoclassical
  particle/energy fluxes from `GFLNC_*` and `QFLNC_*`, plus a distinct
  `QFLNCC_*` model. Reconstruct boundary flux from source divergence using
  actual `DVOL`, `SURF`, and unit conversion. Do not invent aggregate species
  or derive undocumented effective diffusivities.
- Expand `core_profiles` neutral populations (wall, volume, recombination),
  correct neutral temperature eV scaling, and preserve `OM0*` toroidal neutral
  rotation with source/time-grid provenance in DD4 code XML (and DD3 typed
  velocity where available).
- Verify source units, native TIME3/X dimensions and flux closure on all
  bundled cases. Add `make verify-neoclassical-neutrals` and a source regression
  to each case in `make verify-all`.

- Fix end-to-end conversion journals failing *after* otherwise successful IMAS
  writes when NetCDF supplies NumPy `int32`/`int64` machine metadata. Normalize
  shot/run/occurrence IDs to native Python integers, serialize other NumPy
  scalars/lists safely, and test successful and partially failed write journals.
  This prevents all-example false failures from `json.dumps` and restores real
  IMAS readback to the end-to-end verification sequence.

## [0.12.1](https://github.com/PrincetonUniversity/transp2imas/compare/v0.12.0...v0.12.1) (2026-10-23)


### Bug Fixes

* Multiple verification scripts ([#116](https://github.com/PrincetonUniversity/transp2imas/issues/116)) ([15b40d1](https://github.com/PrincetonUniversity/transp2imas/commit/15b40d1bc54c6db1253ea2816bf82a6393b6813e))

## [0.12.0](https://github.com/PrincetonUniversity/transp2imas/compare/v0.11.1...v0.12.0) (2026-10-22)


### Features

* Nubeam related quantities ([#114](https://github.com/PrincetonUniversity/transp2imas/issues/114)) ([f1df6d0](https://github.com/PrincetonUniversity/transp2imas/commit/f1df6d0d91140923c2a8b03a726866872f24ab2e))

## [0.11.1](https://github.com/PrincetonUniversity/transp2imas/compare/v0.11.0...v0.11.1) (2026-10-19)


### Bug Fixes

* misc fixes for equilibrium and core profile IDSs ([#112](https://github.com/PrincetonUniversity/transp2imas/issues/112)) ([0e8ae1b](https://github.com/PrincetonUniversity/transp2imas/commit/0e8ae1b43141054defdf2f181d53a72b2d7abc33))

## [0.11.0](https://github.com/PrincetonUniversity/transp2imas/compare/v0.10.0...v0.11.0) (2026-10-19)


### Features

* Pt multi region ([#110](https://github.com/PrincetonUniversity/transp2imas/issues/110)) ([88f61d4](https://github.com/PrincetonUniversity/transp2imas/commit/88f61d420ce12fe711ce91a132d145125e1d9f19))

## 0.9.0 - 2026-09-28

Initial public pre-1.0 baseline of the IMAS-Python converter.

- Validate real ISolver `pf_active`/`pf_passive` source waveforms using ITER
  38530P94 (11 active / 119 passive) and NSTX 142202U01 (9 active / 82
  passive). Add `make verify-isolver`, per-case coil/loop CSV inventories and
  source-side `CC_* = N_effective * KK_*` closures. Require explicit,
  externally verified `cc_effective_turns` before using any `CC_*` ampere-turns
  as physical coil current; recognize and defer KD_/KS_ controller channels.
- Rebuild `pf_active` writer with distinct coil/supply names from the CDF,
  per-channel time and unit conversions, source-aware current and voltage
  assignments, and optional authoritative coil/circuit machine geometry.
  Correct the KT_ target-current bug without incorrectly treating targets
  as measured currents, remove fake `-999` times, and support DD3/DD4
  geometry/circuit conventions. Add `make verify-pf-active`, a synthetic
  device example, and `docs/pf-active.md`. Actual ISolver current waveforms are now tested; device geometry and winding
  conventions still require an authoritative machine model.

Code-quality fixes included in the baseline:

- fix duplicate atomic-species dictionary keys for lithium and Hassium/Holmium;
- correct the Dysprosium key and Technetium atomic metadata discovered in the same audit;
- remove process-global `numpy.seterr(all="raise")` from the conversion entry point;
- use explicit imports in `convert.py`;
- replace broad numerical exception handling in the temperature-peaking/shared-array paths with local NumPy error handling and specific expected exceptions;
- add a source-only duplicate-literal-key check to `make check`.

- Support DD 3.39 and DD 4.1.1 conversion paths with centralized COCOS/sign handling.
- Expand high-confidence TRANSP mappings across core profiles, equilibrium, summary,
  core sources, core transport, radiation, edge sources, actuator IDSs, and NBI
  distribution sources.
- Populate NBI parallel current from `CURB * CONPLJB / BZ`.
- Preserve populated `core_sources` entries with IMAS-Python arrays of structures.
- Add workflow, dataset FAIR metadata, source hashes, and per-IDS provenance.
- Add source-side regression checks for the bundled NSTX cases and static writer/API
  validation.
- Validate and extend the `waves` writer with bundled DIII-D TORBEAM/GENRAY EC
  and C-Mod multi-frequency ICRF cases. IC coherent waves now use resolved
  `PIC_Fj` absorbed-power profiles with `RFPWRABj` closure; `QICHA` is only a
  single-wave fallback. `PEICH`, `PIICH`, and `QRFBTOT` remain plasma-source
  quantities rather than being incorrectly summed as wave absorption.
- Add `examples/verify_waves_mapping.py` and `make verify-waves` for portable
  RF source-side regression checks.
- Extend `pf_passive`: preserve PS_n-native time bases, read A/kA/MA currents,
  allow audited external machine geometry and explicit source-channel sign
  correspondences, support all six DD geometries, and handle the DD3/DD4
  outline-closure difference. Keep geometry out of current-only mode and
  enable device-geometry-only export without PS_n.
- Add `examples/verify_pf_passive_mapping.py`, a synthetic device JSON example,
  `make verify-pf-passive`, and `docs/pf-passive.md`. Real PS_n regression is now
  available through `make verify-isolver`; physical geometry remains pending.
- Add `examples/audit_mapping_coverage.py` and `make audit` to classify CDF variables
  as converter-referenced, derived/support, validation-only, deliberately deferred,
  or still unclassified, with Markdown/CSV/JSON reports.

### Unreleased input provenance extension

- Preserve TR.DAT and recursively referenced external namelist contents in
  workflow XML with hashes, raw text, parsed assignments and missing-file flags.
- Add model-specific input XML to workflow components and waves/core_transport/NBI
  code libraries (TORIC, CQL3D, GENRAY, TORBEAM, FPP, LSC, PT_SOLVER, etc.).
- Store workflow provenance before physics writers; add source-only input XML
  verification and DD checks for wave-code library parameters.

### Unreleased all-example verification and plotting extension

- Dynamically discover every CDF under `examples/data/`; validate CDF grids,
  the available PF and waves mappings, and lossless XML model inputs.
  Write `build/verification/summary.md`, JSON and per-run diagnostic logs.
- Reorganize `make verify-all`, `make verify-waves` and `make verify-inputs` to
  include all 10 current examples; retain `make verify-deep` as an optional slow
  legacy diagnostic suite and run dedicated real-ISolver checks.
- Archive adjacent model input files with an explicit unreferenced marker;
  preserve explicit `PT_TEMPLATE` file authority and flag a missing T16 TR.DAT.
- Correct legacy source diagnostics where scalar `TIME` and profile `TIME3`
  differ (notably DIII-D TORAY `147634T13`).
- Replace fragile `plot_ids.py` positional/time-array logic with an argparse
  CLI that supports DD selection, equilibrium fields, correct R-Z plotting,
  zero-based time selection, and headless PNG/PDF output; add mock plot tests.
- Record the existing ITER 38530P94 EC current closure of 0.543% as a WARN
  (0.5% strict threshold; fail above 0.6%), not as an exact closure.

### Unreleased core-profile time-grid hotfix

- Fix DIII-D 147634T13 DD4.1.1 `core_profiles.global_quantities.t_e_peaking`
  broadcasting failure by interpolating scalar `TE0(TIME)` onto `TIME3`
  before dividing by the volume-averaged electron temperature.
- Add source-only verification for electron and ion peaking on every bundled
  CDF, a synthetic unequal-grid test, and a writer-contract check. This runs
  with `make verify-all` and the new `make verify-core-profiles` target.

### Unreleased equilibrium TIME/TIME3 coordinate fix

- Fix DD 4.1.1 `equilibrium` coordinate validation when CDF scalar TIME has
  fewer samples than profile TIME3 (DIII-D TORAY 147634T13: 4 versus 7).
  `equilibrium.time` and `time_slice` now both use native TIME3.
- Interpolate LCFS Fourier boundary coefficients and magnetic-axis helpers
  from their declared source time grids onto TIME3.  Keep profile-derived
  boundary shaping directly on TIME3; correct psi-axis time validation.
- Introduce dimension-aware `get_ncvar_on_time` and fix the shared time
  interpolation helpers to compare actual coordinates, not just lengths.
  Apply this distinction to core-profile scalar/profile mixing.
- Validate the boundary methods and time coordinates against all 10 bundled
  CDFs, including a synthetic case with equal-length but different TIME and
  TIME3 values.  See `make verify-equilibrium` and `docs/time-grids.md`.

## Unreleased

- Map the reduced MMM07 anomalous-transport output to `core_transport` as
  a combined model plus physically distinct Weiland (CDF `W19`), DRIBM/DRBM,
  ETG, and reduced MTM contributions. Preserve source `TIME3/XB` coefficients
  in SI, model-specific `TIME3/X` reduced growth/frequency outputs as XML,
  and audit the component sum against total diffusivities on DIII-D T16.
  Do not create MMM models from zero-filled inactive placeholders. Add
  `make verify-mmm`, all-example source checks, writer regression and DD
  node checks. Clarify that MTM is an MMM **model output**, not a measurement.

- End-to-end conversion/readback harness, selective IDS CLI and write journal.

## Historical development notes migrated from README

The following development descriptions previously appeared in the project README. They are preserved for traceability; some refer to earlier versions and should not be interpreted as current implementation guarantees.

## Conversion corrections

The following TRANSP-to-IMAS corrections are included:

- `CUR`, `CUROH`, and `CURBS` are converted from A/cm^2 to A/m^2 using `1.0e4`.
  The previous code accidentally used the cm^2-to-m^2 area factor (`1.0e-4`),
  making current densities smaller by `1.0e8`.
- TRANSP `DVOL` and `DAREA` are per-zone quantities. IMAS `volume` and `area` grid
  fields are enclosed/cumulative quantities, so cumulative sums are written. Raw
  zone values are still used internally for volume and area integrations.
- The existing ITER direction convention (`sgn=-1`) is applied consistently to
  plasma current, toroidal magnetic field, current-density profiles, and signed
  magnetic flux quantities in `core_profiles` and `equilibrium`.

 
## IMAS DD4 / COCOS17 note

When writing IMAS Data Dictionary major version 4, `transp2imas` targets COCOS17.
It reads `NLJCCW` and `NLBCCW` from the TRANSP TR.DAT file to determine the
physical plasma-current and toroidal-field directions.  See
the COCOS 17 sign-policy section for the full sign table.  For old runs without those
controls, use `--ip-direction {ccw,cw}` and `--bt-direction {ccw,cw}`.

Run `python examples/verify_cocos17.py <run>.CDF` for a source-side sign sanity check.

## Selecting the IMAS Data Dictionary (v0.3)

TRANSP2IMAS can now populate either a legacy DD3.39 schema or a DD4 schema
explicitly.  Use, for example:

```bash
python convert.py -i examples/data/nstx/141711Y01/141711Y01.CDF -new -b HDF5 --dd-version 3.39.0
python convert.py -i examples/data/nstx/141711Y01/141711Y01.CDF -new -b HDF5 --dd-version 4.1.1
```

With `--dd-version auto` (default), the converter follows the IMAS environment
or the default Data Dictionary bundled with IMAS-Python.  The exact selected
version is printed at startup.  IDSs are created from the DBEntry factory for
that same version rather than creating one DD and relying on cross-major
conversion at write time.

The EC launcher mapping is explicitly schema-aware: DD3 uses
`beam.identifier` + descriptive `beam.name` + `beam.mode`; DD4 uses the short
identifier in `beam.name`, descriptive text in `beam.description`, and (DD4.1+)
`beam.polarization.o_mode_fraction`.  See the Data Dictionary compatibility section for
the complete mapping table and run
`python examples/verify_dd_compat.py 3.39.0 4.1.1` to inspect the schemas available in
your IMAS-Python installation.

## Extended CDF mapping and provenance (v0.4)

Version v0.4 expands `core_profiles`, `equilibrium`, and `summary` with
high-confidence TRANSP CDF quantities and embeds the original TRANSP/model input
files in IMAS `code.parameters` XML when the IDS has a `code` structure.

Run the source-side checks for both active regression cases before conversion:

```bash
make verify-all
make verify-waves
make verify-dd
```

The source-only verifiers do not require IMAS. `examples/verify_dd_compat.py` exercises
the target DD object structure and therefore must be run in an environment with
the requested IMAS Data Dictionary versions installed.  Detailed source
mappings are documented in the core sources and transport mapping section,
the radiation mapping section, the edge sources mapping section, and the waves mapping section.

### Validated RF `waves` mapping

`make verify-waves` exercises three source-side RF regression CDFs without an
IMAS installation:

- DIII-D `147634B13`: six-gyrotron TORBEAM EC deposition/current drive;
- DIII-D `147634G14`: six-gyrotron GENRAY EC deposition/current drive;
- C-Mod `89846`: two-frequency ICRF/TORIC absorbed-power output.

For EC, `PEECHj(TIME3,X)` is the coherent-wave absorbed power density and
thermal-electron deposition. Its finite-volume integral provides wave power.
`ECCURj(TIME3,X)` is integrated over `DAREA` to obtain cumulative toroidal
driven current.

For IC, `PIC_Fj(TIME3,X)` is the authoritative frequency/antenna-resolved
absorbed wave-power profile and `RFPWRABj(TIME)` is the corresponding absorbed
power trace. `QICHA` is used only as a single-wave fallback if a resolved
`PIC_Fj` profile is unavailable. `PEICH`, `PIICH`, and `QRFBTOT` remain
plasma-heating/source quantities in `core_sources`; they are not summed to
construct wave absorption.

The C-Mod case has two IC frequencies but only an aggregate `XGRID_NPHI` /
`PIC_NPHI` spectrum, so `n_phi` is intentionally left empty rather than
duplicating an unresolved spectrum onto both coherent waves.

Referenced namelist templates (for example `PT_TEMPLATE`, `MMM_TEMPLATE`,
`TGLF_TEMPLATE`, etc.) are resolved relative to the TRANSP input directory.  If
present, they are embedded verbatim and also represented as parsed namelist XML;
if missing, their reference and missing status remain in the provenance payload.



### Audit CDF mapping coverage

To inventory which variables in a TRANSP CDF are referenced, validation-only,
intentionally deferred, or still unclassified, run:

```bash
make audit CASE=141711Y01
```

The command writes full Markdown, CSV, and JSON reports under `build/audit/`.
For a CDF outside the bundled regression-data tree, call the tool directly:

```bash
python examples/audit_mapping_coverage.py /path/to/run.CDF \
    --markdown build/audit/run.md --csv build/audit/run.csv
```

The audit is static and intentionally conservative: it is a planning tool for
choosing the next mapping group, not a claim that every referenced variable has
a complete one-to-one IMAS mapping. See
[the mapping workflow](docs/mapping-workflow.md) for interpretation and policy.

## Passive-structure currents and geometry (`pf_passive`)

The converter can export TRANSP `PS_n` passive-loop current histories, using
**each loop's own time axis** (not just top-level `pf_passive.time`). Without
an authoritative device description, it retains a clearly logged **current-only
legacy mode** and does not invent geometry, turns, resistance or connectivity.

With an explicitly sourced geometry manifest, the converter additionally writes
all `pf_passive/loop/element/geometry` shapes supported by the DD and maps a
current to a physical loop only after its identity and sign have been specified:

```bash
python3 convert.py -i /path/to/RUN.CDF -b HDF5 -new \
    --pf-passive-geometry /path/to/machine-pf-passive.json
make verify-pf-passive
```

The only manifest included in the repository is **synthetic test data**, not a
machine file. All current bundled CDF fixtures lack `PS_n`, so actual device
current/geometry correlation is not yet physically validated. See
[the complete manifest specification, DD3/DD4 outline rules and validation
status](docs/pf-passive.md).

## Pre-0.9.0 development milestones

The following 8.x labels were internal development snapshots created before the
first public pre-1.0 baseline. They are retained here only to document the
evolution of the converter.

### v0.8.3: NBI parallel current

v0.8.3 populates `core_sources/source[nbi]/profiles_1d/j_parallel` from the
TRANSP beam-driven toroidal current `CURB`. The converter changes A/cm² to
A/m², applies the physical plasma-current sign, and converts to the IMAS
`<J.B>/B0` representation using `CONPLJB/BZ`. The same shared conversion is
used for the bootstrap-current source.

This release also appends completed `core_sources` entries instead of shrinking
a preallocated IMAS array of structures. IMAS-Python recreates elements during
`resize()`, which previously left the correct number of source slots but erased
their populated contents before writing.

For `174783Q32`, `CURB` agrees with the sum of the fourteen `BDCnn` beam-current
profiles with correlation 0.9999993, and the resulting parallel-current profile
is finite.

### v0.8.2: workflow and provenance API fix

v0.8.2 fixes the metadata/provenance workflow path for DD 4.1.1.  The
workflow writer now uses `build_code_parameters_xml()` (the actual nc2imas
API) and the new `get_ncattr_text()` helper instead of the nonexistent
`build_transp_parameters_xml()` and `_decode_attr()` names.  The latter bug
had silently suppressed TRANSP/TGLF/GLF23 component metadata because `_attr()`
caught the exception; the former aborted conversion when `workflow` was
reached.  `dataset_fair` also goes through the common code-property helper when
the target DD supplies a code node.

Run `python examples/verify_writer_api.py` as a no-IMAS regression check; it statically
verifies every `self.im.<method>()` call made by writer modules against the
methods implemented by `nc2imas` and its local base classes.

### v0.8.1: distribution_sources namelist accessor fix

v0.8.1 fixes the per-beam species lookup in `imas_distribution_sources.py`.
The v0.8.0 writer accidentally called a nonexistent `get_namelist_value()`
method. It now uses the established `get_namelist_var()` API with the correct
zero-based index (`ib-1`) and treats missing species namelist metadata as
optional rather than aborting conversion. `distribution_sources` remains
independent of the full fast-ion `distributions` IDS: NBI birth/deposition
source moments can be exported from `BDEPnn_*` even when FBM phase-space output
was not saved.

### v0.8.0: NBI distribution sources and metadata/provenance

v0.8.0 adds `distribution_sources` for NUBEAM birth/deposition moments when
`BDEPnn_*` is present, including per-beam particle, energy, toroidal-momentum,
global power/torque, and shinethrough bookkeeping.  It does not synthesize a
full `distributions` IDS if FBM phase-space arrays are absent.

The converter also writes `workflow` and `dataset_fair` where supported, and
attaches whole-IDS source provenance to `ids_properties/provenance`.  Full
TR.DAT/template text and parsed XML remain in `code.parameters`.  See
the NBI kinetic-source mapping section and the metadata and provenance section.



## Documentation

Detailed mapping audits, compatibility notes, provenance design, and migration history are maintained in [the technical reference](docs/technical-reference.md).

## Versioning and releases

The public baseline version is `0.9.0`. Runtime provenance uses
`nc2imas.__version__`; `VERSION` and `TRANSP2IMAS_VERSION` are aliases.

Release Please uses the `simple` strategy because this repository is a
script-based Python project rather than an installable package with
`pyproject.toml`/`setup.py`. The canonical release file is `version.txt`; the
annotated version in `nc2imas.py` is updated in the same release PR.

- Conventional Commits such as `fix: ...` and `feat: ...` determine the next
  semantic version.
- Pushes to `main` update a release pull request and `CHANGELOG.md`.
- Merging that release pull request updates `version.txt`, `nc2imas.py`, and the
  manifest, then creates the corresponding `vX.Y.Z` GitHub release/tag.
- `0.9.0` is the bootstrap baseline. After the one-time `v0.9.0` GitHub release
  exists, routine development should not edit version files or the manifest
  manually.

### Model input XML in workflow and physics IDSs

The complete source TR.DAT and recursively referenced external model namelists
are embedded as XML under `workflow.time_loop.component[0].parameters`.
Individual configured wave codes (TORIC, GENRAY, TORBEAM, CQL3D, FPP, etc.)
are additionally recorded in separate workflow components and in
`waves.code.library[].parameters`. Referenced transport/NUBEAM model input
subsets are associated with their corresponding IDS code libraries. Missing
files are explicitly flagged; comments do not count as active model inputs.
Input provenance is written before physics IDS conversion starts.

See [docs/input-provenance.md](docs/input-provenance.md). Test using
`make verify-inputs` (source-only) and `make verify-dd` (IMAS-Python).


### Unified anomalous theory transport: MMM07 and GLF23

`core_transport.model[]` now preserves the combined MMM07 transport coefficients
and separate Weiland (source `W19`), DRIBM/DRBM, ETG and reduced-MTM
contributions as **model outputs**, not experimental diagnostics. Diffusivities
are mapped in SI units on the original `TIME3/XB` grid, while model growth
rates, frequencies, selected `ky*rho_s`, and untyped aggregate coefficients
remain in model-specific XML, but their `TIME3` and `X` coordinates are now
stored only in the native `core_transport.time` and
`model[].profiles_1d[].grid_flux.rho_tor_norm` fields. The transport
coefficient grid remains `model[].profiles_1d[].grid_d.rho_tor_norm` (`XB`). The combined result and its
components must not be added together. Only physically nonzero model families are instantiated; most examples
contain inactive, zero-filled placeholders. ITER `38530P94` has nonzero
GLF23 `XKEGLF23`, `XKIGLF23` and `ETPHIGLF` on `TIME3/XB` and is now
converted into its own `core_transport.model[glf23]`. The aggregate
`DIFFIGLF` remains supplementary rather than assigned to an invented species.
TGLF awaits a validated output example.

Verify with `make verify-theory`, followed by `make verify-all` and (where IMAS
is installed) `make verify-dd`. See [unified theory mapping details](docs/theory-transport.md).

### Separate neoclassical theory transport and neutral-state profiles

`core_transport.model[]` preserves neoclassical particle and heat fluxes from
`GFLNC_*`, `QFLNC_*`, and `QFLNCC_*` divergences with finite-volume
reconstruction. The same generic theory interface also exports nonzero
GLF23 heat and momentum transport from ITER 38530P94, and the MMM07 models
from DIII-D T16. TGLF requires a validated example. Neoclassical presence
alone does not prove that PT_SOLVER enabled it.

Independently, `core_profiles.neutral[]` contains wall, volume, and
recombination neutral densities and eV temperatures; DD4 toroidal neutral
rotation remains supplemental XML. These are not theory transport models.

Verify separately using `make verify-theory`, `make verify-neoclassical`,
`make verify-neutrals` and `make verify-all`. See the
[theory transport](docs/theory-transport.md),
[neoclassical](docs/neoclassical-transport.md), and
[neutral profiles](docs/neutral-profiles.md) guides.

### Sawtooth history

The `sawteeth` writer exports `TLASTSAW(TIME)` into the native IMAS
`diagnostics.previous_crash_time` history. C-Mod 89846 has 38 observed
updates; valid `XIQ1(TIME)` samples are mapped to the native `q=1` radius.
Trigger mechanisms, mixing/inversion radii, and post-crash profiles are **not**
inferred. The first recorded `TLASTSAW` value is not counted as a newly
observed crash. Run `make verify-sawteeth` or convert C-Mod with
`--ids sawteeth`. See [details](docs/sawteeth.md).
