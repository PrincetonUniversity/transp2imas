# Changelog

Notable changes to transp2imas are recorded here. Future entries and GitHub
releases are managed automatically by Release Please from Conventional Commit
messages.

## Unreleased

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
