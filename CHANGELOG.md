# Changelog

Notable changes to transp2imas are recorded here. Future entries and GitHub
releases are managed automatically by Release Please from Conventional Commit
messages.

## Unreleased

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
