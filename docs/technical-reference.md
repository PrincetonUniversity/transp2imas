# transp2imas technical reference

This document contains detailed mapping audits, compatibility notes, provenance design, and migration history.

The DIII-D `174783Q32` data referenced in historical audits below are no longer
distributed with the repository. The active, distributable regression fixtures
are the two NSTX cases under `examples/data/nstx/`, two DIII-D EC cases under
`examples/data/d3d/`, and the C-Mod ICRF case under `examples/data/cmod/`.



## TRANSP CDF -> IMAS coverage audit: 174783Q32

### Source inventory

The historical `174783Q32.CDF` contains **2465 variables**.  Most are detailed
TRANSP diagnostics, per-model intermediate quantities, balance terms, or
module-specific outputs that do not have a one-to-one IMAS node and therefore
should not be exported blindly.

The v8.4 extension focuses on quantities whose CDF meaning, units, radial grid,
and IMAS semantics are sufficiently clear to map without model-dependent
assumptions.

### Source checks used in the extension

- `PPLAS` is consistent with `(2/3)*UTHRM*1e6` to a median relative difference
  of about `1.36e-7`, confirming that `PPLAS` is the thermal pressure for this
  CDF while `UTHRM` is already a thermal energy density.
- `sum(DVOL)` agrees with `PVOL` within about `5.27 cm^3` over a plasma volume
  of order cubic metres.
- `sum(DAREA)` agrees with `PAREA` within about `7.2e-3 cm^2`.
- `VRPOT` is on the TRANSP `XB` grid; `SHAT`, `ETA_USE`, `CUR`, `PLJB`, and
  `PMASS` are on the `X` grid.
- `TRFLX` is stored in Webers; `PLFLX` and `PSI0_TR` are Wb/rad.

### Added to core_profiles

| IMAS quantity | TRANSP CDF source | Conversion / note |
|---|---|---|
| electrons/pressure_thermal, pressure | `UE` | `(2/3)*UE`, J/cm3 -> Pa |
| pressure_ion_total | `UI` | `(2/3)*UI`, J/cm3 -> Pa |
| pressure_thermal | `PPLAS` | direct Pa |
| pressure_perpendicular | `PPLAS`, `UFASTPP` | `PPLAS + UFASTPP*1e6` |
| pressure_parallel | `PPLAS`, `UFASTPA` | `PPLAS + 2*UFASTPA*1e6` |
| n_i_thermal_total | `NI` | N/cm3 -> m^-3 |
| t_i_average | `TI` | direct eV |
| ion/pressure_thermal | ion thermal density + `TI` | `n*T*e` |
| magnetic_shear | `SHAT` | direct |
| conductivity_parallel | `ETA_USE` | ohm*cm -> ohm*m, then reciprocal |
| e_field/parallel | `ETA_USE`, `PLJBOH` | `eta * j_ohmic_parallel` |
| phi_potential | `VRPOT` | XB -> X, volts |

The previous ad-hoc radial electric field reconstruction is removed.  The
potential is preserved, but `E_r` is intentionally left empty until a metric
mapping is selected explicitly.

### Added to equilibrium/profiles_1d

| IMAS quantity | TRANSP CDF source |
|---|---|
| psi_norm | signed `PLFLX` + `PSI0_TR` |
| rho_tor | signed `TRFLX` and physical `BZ` |
| dpsi_drho_tor | derivative of the above |
| dpressure_dpsi | `PTOWB` mapped X -> XB, derivative vs psi |
| f_df_dpsi | `GFUN*BZXR`, derivative vs psi |
| j_phi / legacy j_tor | `CUR`, A/cm2 -> A/m2, X -> XB |
| j_parallel | `PLJB/BZ`, X -> XB |
| magnetic_shear | `SHAT`, X -> XB |
| mass_density | `PMASS`, g/cm3 -> kg/m3, X -> XB |

### Rebuilt summary IDS

The previously disabled summary writer is enabled after correcting its energy
and sign handling.  It now includes global current/field/beta/volume/q95,
confinement quantities, direct integrated thermal/electron/ion/total energies,
fast anisotropic energies, diamagnetic and poloidal-field energies, ohmic and
radiated power, current components, and total/DD neutron rates where present.

Important correction: `UTOTL`, `UTHRM`, `UE`, and `UI` are **energy densities**
in the CDF, so they are integrated directly.  They are not multiplied by 1.5.

### Input/provenance conversion

For every active IDS that has a `code` structure, `code.parameters` now contains
XML with:

1. target DD and COCOS/sign metadata;
2. source CDF metadata (TRANSP version, DOI, build date, hashes);
3. complete source TRANSP `TR.DAT` text;
4. a structured XML representation of its namelist assignments;
5. every referenced `*TEMPLATE` file, if present;
6. adjacent run-specific `.nml` files, if present;
7. SHA-256 hashes for embedded files and explicit status for missing files.

For the historical audit input, the following files are referenced but are not in
the supplied repository snapshot:

- `PT_TEMPLATE = ./174783Q32_input_pt.nml`
- `MMM_TEMPLATE = ./174783Q32_input_mmm.nml`

They will be embedded automatically when present in the run directory.

### Subsequent IDS work and remaining gaps

- `core_sources` and `core_transport` were subsequently audited and enabled in
  v8.5/v8.5.1. In v8.5.2, predictive `*_MOD` transport is separated into a
  `transport_solver` model. See the core sources and transport mapping section.
- The dedicated `radiation` IDS was rebuilt and enabled in v8.5.2; see
  the radiation mapping section. Inspection of this CDF also corrected an important
  semantic error: `P0RECO` is **recombination-neutral kinetic source power**,
  not radiative recombination emission, and is therefore handled by
  `edge_sources`, not `radiation`.  v8.5.5 additionally uses the namelist
  `AIMP=12`, `XZIMP=6` to resolve generic `PRLS_TOK` / `PRLX_SINGL` line
  radiation to carbon; both aliases are identical to `PRAD_LI` in this run and
  only one is attached, preventing double counting.
- `edge_sources` was rebuilt in v8.5.3 around direct source-rate profiles and
  refined in v8.5.5 after validating the individual neutral-process terms on
  this CDF.  `SIGF_D_D = SEGF_D + SFCXGF_D` (~1e-12 relative) and
  `SIRC_D_D = SERC_D + SFCXRC_D` (~1e-6 relative), so ionization and charge
  exchange are now distinct IMAS source entries.  Their combined ion-source
  integrals agree with `GASD` and `RCYD` to better than 5e-7 relative.  The
  `PI-PCX` and `TQI-TQCX` sign conventions are independently reproduced by the
  TRANSP neutral power/momentum balance arrays.  Recombination recapture obeys
  the same decomposition: `SIV0D = SERECO + SFCXRECO`.  Neutral state profiles
  (`N0*`, `T0*`, `OM0*`) remain deliberately outside `edge_sources`. See
  the edge sources mapping section.
- `waves` was rebuilt and enabled in v8.5.4; see the waves mapping section.  This run
  contains no supported EC/IC absorbed/deposited wave payload: no `PEECH*`,
  `ECCUR*`, `PICHA*`, `FREQA*`, `RFPWRAB*`, `PEICH`, `PIICH`, `QRFBTOT`, or
  `XGRID_NPHI`.  `EC_MODEL='NONE'` and `NLLH=.false.`.  `PRFB_D`/`PRFFI` are
  summary fast-ion RF powers, not a spatial coherent-wave solution, so no
  `waves` IDS is written for 174783Q32.
- `NLFBM=.true.` in TR.DAT and the CDF contains `BMVOL` and `BDENS2_D`, but it
  contains no `FBM`, `EBA2PL`, or `EBA2PP`.  The phase-space distribution is
  therefore incomplete and must not be reconstructed from reduced moments.
  A full distribution writer should be validated on a run saved with the full
  FBM payload.
- `NLSAW=.false.` and the only explicitly sawtooth-named CDF quantity is scalar
  `TESAW`; `XIQ1` gives the q=1-surface location but is not a crash event.  No
  crash times or mixing diagnostics are present, so there is no sawteeth event
  to export for this run.
- 101-point `RMAJM` local poloidal/toroidal velocities remain excluded from the
  50-point flux-profile velocity nodes.
- A radial electric field is not inferred from `VRPOT` without an explicitly
  selected metric convention.
- `pf_passive` remains the final large audit item.  This attached CDF contains
  no `PS_1...PS_n` passive-structure current channels, so it cannot validate
  even the current-only legacy path, much less geometry, resistance,
  inductance, and connectivity.  A CDF/run with passive-structure output is
  required before this writer can be completed responsibly.

### v8.6.0 additional quantities and metadata

The real CDF exposes several high-confidence quantities that can be saved even
though its full FBM is absent:

- `PPHI(TIME3,X)` is now written to `core_profiles/profiles_1d/momentum_phi`
  (legacy `momentum_tor`).  The CDF independently verifies
  `PPHI = AMOI * OMEGA` to about 2.2e-7 relative RMS.  The scalar `AMOM` entry in
  this CDF is only a placeholder and is deliberately not used.
- Per-beam NBI `BDEPnn_*`, `EINJnn_*`, `BPSHnn_*`, and `TQTOTnn` are now written
  to `distribution_sources`, including radial particle/energy/momentum sources,
  global particles/power/torque, and shinethrough.  See
  the NBI kinetic-source mapping section.
- DD4.1+ `summary` metadata now marks the data as a simulation and records
  machine, pulse, relative simulation limits, and uniform time step when valid.
- `workflow` and `dataset_fair` IDSs are now written where supported.
- every output IDS receives whole-IDS source-file provenance in
  `ids_properties/provenance`, in addition to the full input text/XML already
  kept in `code.parameters`.  See the metadata and provenance section.

Potentially useful quantities that remain intentionally outside the current
mapping include neutral source-state temperatures/rotation (`T0*`, `OM0*`),
boundary energy transport (`PEEDG`, `PIEDG`), and reduced fast-ion moments when
the full FBM phase-space grids are absent.  They should be mapped only after a
semantically correct destination IDS/node is established rather than stored in
an unrelated structure for completeness alone.


## TRANSP2IMAS DD4 / COCOS 17 sign policy

IMAS Data Dictionary major version 4 uses **COCOS 17**.  In IMAS cylindrical
coordinates `(R, phi, Z)`, positive `phi` is counter-clockwise (CCW) when viewed
from above.  TRANSP stores the physical field-orientation controls in the input
namelist:

- `NLJCCW`: plasma current is CCW viewed from above.
- `NLBCCW`: toroidal magnetic field is CCW viewed from above.

`transp2imas` reads those controls (case-insensitively).  They can be overridden
with `--ip-direction {ccw,cw}` and `--bt-direction {ccw,cw}` when converting old
runs whose TR.DAT is incomplete.

### Sign table used by the converter

Let

- `s_I = +1` for `NLJCCW=.true.`, otherwise `-1`;
- `s_B = +1` for `NLBCCW=.true.`, otherwise `-1`.

For DD4/COCOS17 the converter applies:

| IMAS quantity | multiplier applied to TRANSP magnitude/profile |
|---|---:|
| plasma current `Ip` | `s_I` |
| toroidal current densities (`CUR`, `CUROH`, `CURBS`, ...) | `s_I` |
| loop voltage | `s_I` |
| toroidal field `B_phi`, vacuum `B0` | `s_B` |
| `F = R B_phi` | `s_B` |
| toroidal flux `Phi` | `s_B` |
| safety factor `q` | `s_I * s_B` |
| full poloidal flux `psi` | `-s_I` |
| poloidal field components `B_R`, `B_Z` | `s_I` |
| toroidal rotation referenced by TRANSP to the plasma-current direction | `s_I` |

COCOS11 and COCOS17 use the same cylindrical orientation and the same sign rule
for `q`, but reverse the reference sign of poloidal flux.  Thus the DD3 fallback
uses `psi multiplier = +s_I`, while DD4 uses `psi multiplier = -s_I`.

### Examples

#### MAST-U with clockwise Ip and clockwise Bt

`NLJCCW=.false.`, `NLBCCW=.false.` gives

- `Ip < 0`
- `B_phi < 0`
- `q > 0`
- COCOS17 `psi` increases numerically from an axis value of zero toward positive
  values at the edge (because the positive-Ip COCOS17 reference is reversed).

The same physical discharge written in COCOS11 would have the opposite sign of
`psi`.

#### Historical audit case 174783Q32 (DIII-D)

The historical `174783Q32TR.DAT` contains

```
NLJCCW = .true.
NLBCCW = .false.
```

so the expected DD4 signs are

- `Ip > 0`
- `B_phi < 0`
- `q < 0`
- `psi < 0` away from an axis offset of zero.

This makes the case useful for testing the sign conversion because Ip and Bt are
opposite.

### Flux units

TRANSP `TRFLX` is already **toroidal flux in Wb**.  It must not be multiplied by
`2*pi`.  `PLFLX` and `PSI0_TR` are in Wb/rad and are converted to full Wb with
`2*pi`; `PLFLX2PI` is already in Wb.

The `PSI0_TR` magnetic-axis offset is added to the signed enclosed poloidal flux
when writing absolute IMAS `psi` values.

### PF coil and passive-loop currents

PF circuit-current signs are **not** forced to follow `s_I`.  Their sign depends
on the independent circuit/winding orientation.  The converter therefore
preserves the sign stored in the TRANSP CDF.  A device-specific PF wiring audit
is required if those source signals are not already expressed in the IMAS
positive-current convention.

### Polar angle for inverse equilibrium

DD4 polar angle is `theta = -atan2(Z-Zaxis, R-Raxis)` (clockwise on the right-hand
poloidal cross-section).  TRANSP Fourier moments are evaluated using the
corresponding opposite geometrical angle, and the IMAS coordinate system is
identified as `inverse_rhotornorm_polar` (index 33) because TRANSP `XB` is
`sqrt(Phi/Phi_boundary)`.


### Current-density mapping in core_profiles

TRANSP `CUR`, `CUROH`, and `CURBS` are toroidal current densities in A/cm^2.
They are converted with `1e4` and the physical `Ip` sign.  In DD4,
`core_profiles/profiles_1d/j_phi` is the toroidal current-density node, while
`j_total`, `j_ohmic`, `j_non_inductive`, and `j_bootstrap` are parallel-current
definitions based on `<J.B>/B0`.  Therefore v8 writes `CUR` to `j_phi`, derives
`j_total` and `j_ohmic` from TRANSP `PLJB` and `PLJBOH`, obtains
`j_non_inductive = j_total - j_ohmic`, and maps the selected `CURBS` bootstrap
profile to the parallel convention using `CONPLJB/BZ`.


## TRANSP -> IMAS core_sources/core_transport mapping (v8.5.2)

### core_sources

| IMAS process | TRANSP CDF fields | Notes |
|---|---|---|
| nbi | PBE, PBI, PBTH, SBE, SBTH, TQBTH, TQJXB, CURB, CONPLJB, BZ, PBTH_*, SBTH_* | PBI + PBTH is total thermal-ion energy source; torques are converted from TRANSP co-current sign to physical +phi; CURB is converted to `<J.B>/B0` for `j_parallel`. |
| ec | PEECH | Electron energy source; ECCUR is not used as j_parallel because it is toroidal current density. |
| lh | PELH, PILH | Electron/ion energy sources when present. |
| ic | PEICH, PIICH | Electron/ion energy sources when present. |
| fusion | PFE_*, PFI_*, PFTH_*, SFTH_* | Species-resolved thermalization when present. |
| ohmic | POH, PLJBOH, BZ | POH -> electron energy; PLJBOH/BZ -> <J.B>/B0. |
| bremsstrahlung | PRAD_BR | Negative electron energy source. |
| synchrotron_radiation | PRAD_CY | Negative electron energy source. |
| line_radiation | PRAD_LI | Negative electron energy source. |
| collisional_equipartition | QIE | Electron = -QIE, ions = +QIE. |
| cold_neutrals | SCEW, SCEV, PION | Electron particles from wall+volume neutrals; ionization work is an electron energy loss. |
| bootstrap_current | CURBS, CONPLJB, BZ | Selected TRANSP bootstrap profile converted to approximate <J.B>/B0. |

All W/cm3 and particles/cm3/s profiles are converted to SI. `*_inside` and
global quantities are generated from cumulative/full volume integrals using
TRANSP `DVOL` converted to m3.

### core_transport

#### `model[combined]` (identifier index 1)

| IMAS field | TRANSP | Conversion |
|---|---|---|
| electrons.particles.d | DIFFE | cm2/s -> m2/s |
| electrons.particles.v | VELE | cm/s -> m/s |
| electrons.particles.flux | EPTR_OBS | cumulative integral(div*DVOL)/SURF |
| electrons.energy.d | CONDE | cm2/s -> m2/s |
| electrons.energy.v | VEL_TE | cm/s -> m/s |
| electrons.energy.flux | EETR_OBS | cumulative integral(div*DVOL)/SURF |
| total_ion_energy.d | CONDI | cm2/s -> m2/s |
| total_ion_energy.v | VEL_TI | cm/s -> m/s |
| total_ion_energy.flux | IETR_OBS | cumulative integral(div*DVOL)/SURF |
| momentum_phi/toroidal.d | CHPHI | cm2/s -> m2/s |
| momentum_phi/toroidal.v | VMO_NET | cm/s -> m/s |
| momentum_phi/toroidal.flux | AMTR_OBS | cumulative integral(div*DVOL)/SURF; physical toroidal sign applied |
| conductivity_parallel | ETA_USE | ohm*cm -> ohm*m; sigma=1/eta; X -> XB interpolation |
| ion[*].particles | DIFF*, VEL*, PTR*_OBS | H/D/T/He3/He4 when available |

#### `model[transport_solver]` (identifier index 2)

Predictive/model divergences are deliberately kept separate from the observed
balance/effective-coefficient model.

| IMAS field | TRANSP | Conversion |
|---|---|---|
| electrons.particles.flux | EPTR_MOD | cumulative integral(div*DVOL)/SURF |
| electrons.energy.flux | EETR_MOD | cumulative integral(div*DVOL)/SURF |
| total_ion_energy.flux | IETR_MOD | cumulative integral(div*DVOL)/SURF |
| ion[*].particles.flux | PTRH_MOD, PTRD_MOD, PTRT_MOD, PTRHE3_MOD/PTR3_MOD, PTRHE4_MOD/PTR4_MOD | species-resolved cumulative integral(div*DVOL)/SURF when present |

`IPTR_MOD` is the TRANSP divergence of the aggregate total-ion particle flux.
The IMAS `core_transport` schema has no aggregate total-ion particle-flux node,
so it is **not** assigned to a fictitious ion species. Its presence is logged;
species-resolved `PTR*_MOD` channels are exported instead when available. An
empty `transport_solver` model is not created when `IPTR_MOD` is the only model
quantity available.

The three IMAS radial grids (`grid_d`, `grid_v`, `grid_flux`) use the TRANSP
`XB` normalized toroidal-flux boundary grid. Enclosed volume/area and
flux-surface area are also populated on XB.

#### Intentionally not mapped

- `CURB`/`ECCUR` are toroidal current densities and are not written to
  `core_sources.j_parallel`.
- Global `current_parallel` is left empty because deriving it from `<J.B>/B0`
  requires the precise IMAS surface-integral geometry, not a simple zone-area
  sum.
- `VPOLD_NC`/`VTORD_NC` remain omitted because they are on the 101-point
  `RMAJM` grid rather than the 50-point flux grid.
- `E_r` is not inferred from `VRPOT` until the radial metric convention is
  explicitly defined.


## TRANSP2IMAS Data Dictionary compatibility

### Scope

TRANSP2IMAS v8.3 can create IDS objects directly in either the legacy IMAS
Data Dictionary 3.39 family or the IMAS DD4 family.  The target schema is
selected with `--dd-version` and is passed to both `imas.DBEntry` and its
`IDSFactory`, so the object that is populated has the same schema as the
entry that is written.

Examples:

```bash
## Legacy output
python convert.py -i examples/data/nstx/141711Y01/141711Y01.CDF -new -b HDF5 --dd-version 3.39.0

## Current output
python convert.py -i examples/data/nstx/141711Y01/141711Y01.CDF -new -b HDF5 --dd-version 4.1.1
```

`--dd-version auto` (the default) uses `IMAS_VERSION`, then
`IMAS_DD_VERSION`, and finally the default DD bundled with IMAS-Python.
The actual version selected by the DBEntry factory is printed at startup.

### Schema mappings used by the converter

| Quantity / IDS | DD 3.39.x | DD 4.x |
| --- | --- | --- |
| COCOS target | 11 | 17 |
| EC beam short identifier | `ec_launchers/beam/identifier` | `ec_launchers/beam/name` |
| EC descriptive text | `ec_launchers/beam/name` | `ec_launchers/beam/description` |
| EC dominant mode | `ec_launchers/beam/mode` | removed |
| EC O-mode fraction | beam-level field when present; DD3.39.0 itself may only have `mode` | DD4.1+: `beam/polarization/o_mode_fraction` with `beam/polarization/time` |
| NBI unit identifier | `nbi/unit/identifier` | `nbi/unit/name` |
| NBI description | `nbi/unit/name` | `nbi/unit/description` |
| Species/state text | DD3.39 `label` | DD4 `name` |
| Nuclear charge `z_n` | floating point | integer |
| Core total toroidal current density | `profiles_1d/j_tor` | `profiles_1d/j_phi` |
| TF vacuum field-times-radius | `b_field_tor_vacuum_r` | `b_field_phi_vacuum_r` |
| Equilibrium toroidal B | `b_field_tor` | `b_field_phi` |
| Equilibrium axis psi | `global_quantities/psi_axis` | DD4.1+: `global_quantities/psi_magnetic_axis` |
| Equilibrium normalized beta | `beta_normal` | `beta_tor_norm` |
| Selected boundary psi | `boundary_separatrix/psi` / legacy nodes | `boundary/psi` |
| Wall limiter closure | `unit/closed=1`, without duplicate endpoint | first `(R,Z)` point repeated as the last point |
| Waves toroidal mode number | `n_tor` | `n_phi` |

Some `_tor -> _phi` and `label -> name` renames occurred in late DD3.  The
converter contains version-tuple checks for those transitions, while the two
primary supported targets are explicitly 3.39.x and 4.x.

### EC launcher semantics

The EC change is semantic, not only syntactic.  In DD3.39, `beam/name` is the
descriptive beam name and `beam/identifier` is the short identifier.  In DD4,
`beam/name` is the short unique identifier and descriptive text belongs in
`beam/description`.  TRANSP2IMAS therefore does not use a generic `name`
assignment for EC beams.

TRANSP `RFMODECH` is interpreted as the fraction of launched EC power in
O-mode.  In DD3 the converter writes the legacy `mode` value and, where the
selected DD provides it, the beam-level `o_mode_fraction`.  In DD4.1+ it
writes the full fraction to `beam/polarization/o_mode_fraction` and fills the
polarization time base.  It does not write the removed DD3 `mode` node into a
DD4 object.

### Version selection and writing

Every IDS is created with:

```python
ids = dbentry.factory.new(ids_name)
```

The converter never creates a DD3 IDS and then relies on cross-major automatic
conversion to save it in DD4 (or vice versa).  This avoids ambiguity for
non-backwards-compatible schema changes and guarantees that hard-coded field
assignments are validated against the requested target schema immediately.

### Schema-only verification

Run the bundled verifier in an IMAS-Python environment:

```bash
python examples/verify_dd_compat.py 3.39.0 4.1.1
```

This creates in-memory IDS objects from both factories and prints the fields
that distinguish the two schemas.  It does not create or modify a database.
A non-zero exit status means one of the requested DD versions is not available
or the IDS layout could not be inspected.

### Physics conventions

DD3 output follows COCOS11 and DD4 output follows COCOS17.  The physical
plasma-current and toroidal-field directions are taken from TRANSP `NLJCCW`
and `NLBCCW`, with `--ip-direction` / `--bt-direction` overrides available.
See the COCOS 17 sign-policy section for the sign treatment and the source-side sanity
checks.


## NBI kinetic-source mapping (v8.6.1)

The supplied `174783Q32.CDF` does not contain a complete FBM phase-space
payload, but it does contain high-quality per-beam NUBEAM **birth-source
moments**.  These belong in `distribution_sources` and can be exported without
inventing a full distribution.

### Availability versus full fast-ion distributions

`distribution_sources` is intentionally independent of the full `distributions`
IDS. The former contains source terms for kinetic equations and can be written
from NUBEAM `BDEPnn_*`, `TQTOTnn`, and shinethrough outputs. The latter requires
the saved FBM phase-space distribution (for example `FBM`, `EBA2PL`, `EBA2PP`
and their grids). Therefore 174783Q32 can legitimately contain
`distribution_sources` while having no complete fast-ion `distributions` IDS.

Per-unit beam species metadata are read from `ABEAMA` and `XZBEAMA` with
TRANSP's 1-based beam number converted to the zero-based index expected by
`nc2imas.get_namelist_var()`. Missing namelist metadata does not suppress the
source itself.


### Per NBI unit

For each active `BDEPnn_TOT` beam, one `distribution_sources/source` is written.
Configured all-zero BDEP beams remain represented by the `nbi` hardware IDS but
are omitted from `distribution_sources`.

| IMAS quantity | TRANSP source | Conversion |
|---|---|---|
| profiles_1d/particles | `BDEPnn_TOT` | particles cm^-3 s^-1 -> m^-3 s^-1 |
| profiles_1d/energy | `sum_k BDEPnn_Ek * EINJnn_Ek` | eV -> J, giving W/m^3 |
| profiles_1d/momentum_phi | `TQTOTnn` | N m/cm^3 -> N/m^2; apply IMAS toroidal sign |
| global/particles | volume integral of particle profile | s^-1 |
| global/power | volume integral of energy profile | W |
| global/torque_phi | volume integral of `TQTOTnn` | N m |
| shinethrough/power | `sum_k BPSHnn_Ek` | W |
| shinethrough/particles | `sum_k BPSHnn_Ek/(EINJnn_Ek*e)` | s^-1 |

The source process is identified as NBI and records the NBI unit number.  The
source species is taken from `ABEAMA`/`XZBEAMA` when available.

### Validation on 174783Q32

Ten of the fourteen configured beam arrays are active: 01-06 and 11-14.

- `BDEPnn_TOT` agrees with the sum of the three energy components with typical
  normalized RMS residuals of about 2e-4 to 5e-4.
- `(birth power + shinethrough power) / PINJ` is close to unity for active beam
  periods, providing an independent energy bookkeeping check.
- The sum of all `TQTOTnn` profiles reproduces `TQIN` with about 2.45e-3
  relative RMS residual and correlation 0.9999968.
- The integrated sum of per-beam torques spans approximately -2.55 to
  +10.40 N m over this run.

The full `distributions` IDS remains intentionally unwritten because `FBM`,
`EBA2PL`, and `EBA2PP` are absent.  Reduced moments such as `BDENS2_D` are not
used to synthesize missing phase space.


## TRANSP -> IMAS `edge_sources` mapping (v8.5.5)

The v8.5.5 writer keeps TRANSP neutral **state** quantities separate from true
source-rate quantities and decomposes the neutral interaction terms into their
physical IMAS processes.  IMAS `edge_sources` particle nodes have units of
`m^-3 s^-1`; therefore `N0*` arrays (`N/cm^3`) are never written as particle
sources.

The IDS uses an explicit one-dimensional GGD on TRANSP `X`, interpreted as
`rho_tor_norm`.  The grid is `linear` (identifier 1), its primary space is
`primary_standard` (1), the coordinate is `rho_tor_norm` (12), and source
quantities use the implicit `nodes` subset (1).

### Gas-flow and recycling neutrals

For every resolved neutral isotope/group `G`, gas-flow (`GF`) and recycling
(`RC`) are split into two source entries instead of being merged under the
`ionization` identifier:

| Physical process | TRANSP profiles | IMAS representation |
| --- | --- | --- |
| Atomic ionization | `SEGF_G` / `SERC_G`; `PIGF_G` / `PIRC_G`; `TQIGF_G` / `TQIRC_G` | identifier `atomic_ionization` (707); +electron, +singly-charged ion, -neutral particle source; ionization power/torque to the ion and the opposite neutral term |
| Fast-ion charge exchange | `SFCXGF_G` / `SFCXRC_G`; `PCXGF_G` / `PCXRC_G`; `TQCXGF_G` / `TQCXRC_G` | identifier `charge_exchange` (305); +ion / -target-neutral conversion term; `PCX*` and `TQCX*` are explicitly transfers **to neutrals**, hence neutral + and aggregate-ion - energy/momentum |

TRANSP also supplies `SIGF_I_G` / `SIRC_I_G` (and aggregate `SISGF_G` /
`SISRC_G`).  They are retained as consistency diagnostics rather than being
written wholesale into the ionization entry because they contain both pieces.
For the supplied `174783Q32.CDF`:

- `SIGF_D_D = SEGF_D + SFCXGF_D` to about `1e-12` relative error;
- `SIRC_D_D = SERC_D + SFCXRC_D` to about `1e-6` relative error;
- volume-integrated `SIGF_D_D` agrees with `GASD` to `4.4e-7` relative;
- volume-integrated `SIRC_D_D` agrees with `RCYD` to `4.9e-7` relative.

The CDF neutral-balance diagnostics independently confirm the energy/momentum
sign convention: `BALE0_* ~= PI* - PCX* + PFLX0*` and
`TQBA0_* ~= TQI* - TQCX* + TQ0FL*` to floating-point accuracy/scaling.

`cm^-3` source densities and `W/cm^3` / `N m/cm^3` source densities are
multiplied by `1e6` before writing SI values.

### Recombination and recapture

`S0RECO`, `P0RECO`, and `TQ0RECO` are a recombination-neutral particle,
kinetic-energy, and toroidal-momentum source.  `P0RECO` is **not** photon
radiation.  When the CDF resolves exactly one recombining isotope through
`S0V0G` / `SIV0G`, v8.5.5 also writes the electron and ion particle sinks
implied by atomic recombination.  It does not invent the plasma energy split
corresponding to `P0RECO`/`TQ0RECO`.

The subsequent recapture of those recombination neutrals is also decomposed:

- `SERECO`, `PIRECO`, `TQIRECO` -> `atomic_ionization` (707);
- `SFCXRECO`, `PCXRECO`, `TQCXRECO` -> `charge_exchange` (305).

For `174783Q32`, `S0RECO == S0V0D` and
`SIV0D = SERECO + SFCXRECO` to about `9e-7` relative error, so the isotope and
decomposition are directly supported by the CDF.

### Quantities deliberately not forced into `edge_sources`

`N0SGF_*`, `N0SRC_*`, `N0GF_*_*`, `N0RC_*_*`, `N0V0_*`, `T0*`, and `OM0*`
are neutral density/temperature/rotation **state profiles**.  They belong in a
profile IDS (the TRANSP arrays are on the core `X` flux grid) rather than a
volumetric source node.

`T0GASFL`, `T0RECYC`, `OM0GASFL`, `OM0RECYC`, and `E0INR` are source/boundary
state metadata.  They remain preserved in TRANSP input provenance until a
schema-consistent neutral-state representation is selected.

`PEEDG` and `PIEDG` are integrated energy flow through the plasma boundary,
not volumetric source densities.  They should be handled with a boundary/edge
transport mapping rather than written into `edge_sources`.

The donor fast-ion species of the `SFCX*`/`PCX*`/`TQCX*` terms is not resolved
by these neutral-transport CDF arrays.  The born-ion particle conversion is
therefore identified by the target neutral isotope, while the CX energy and
momentum sink is kept at aggregate-ion level; no fictitious donor species is
created.

### Disposition of legacy edge comments (174783Q32 audit)

Several quantities listed as TODOs in the legacy writer do not belong in a
volumetric `edge_sources` mapping. `T0RECYC`, `T0GASFL`, `OM0RECYC`, and
`OM0GASFL` describe the neutral source/boundary state; `PEEDG` and `PIEDG` are
energy terms through the boundary. `GASIZ`/`GASIN` are absent from the supplied
CDF, while the available `GASD`/`RCYD` integrated rates are used to validate the
spatial source profiles. Neutral densities (`N0*`) are state/profile quantities,
not particle-source densities. These distinctions are recorded in
the next-IDS status section so that future work can place them in the
appropriate state/transport IDS rather than overloading `edge_sources`.


## Metadata and provenance mapping (v8.6.2)

v8.6.2 treats metadata/provenance as first-class converter output rather than
leaving provenance only in free-form comments.

### Per-IDS `ids_properties`

Every IDS passed through `fill_properties()` receives:

- `ids_properties/name = "TRANSP <runid>"` where supported;
- a UTC ISO-8601 `ids_properties/creation_date` for the conversion product;
- the existing run/build/COCOS description in `ids_properties/comment`;
- whole-IDS `ids_properties/provenance` with a single blank `path` node.

The provenance references are the direct source files used by the converter:

1. the TRANSP CDF;
2. the original `TR.DAT`;
3. resolved model/template `.nml` files, when present.

Each reference uses a portable file name plus SHA-256 digest.  On DD versions
with the `reference` AoS the source-file modification time is stored as its UTC
ISO-8601 timestamp.  Older DD layouts use the legacy `sources` string list.
The digest inventory is cached so the large CDF is hashed only once per
conversion, not once per IDS.

### `code.parameters`

For every output IDS that has a `code` structure, `code.parameters` continues
to contain the full provenance XML already introduced earlier:

- exact original `TR.DAT` text;
- parsed namelist XML;
- recursively referenced `*TEMPLATE` files when available;
- adjacent run-specific `.nml` files when available;
- explicit missing-file records;
- source hashes and CDF/TRANSP metadata;
- target DD and COCOS/sign convention metadata.

The converter version in this XML and in `code.version` is taken from the
single `TRANSP2IMAS_VERSION` constant, so it cannot drift between writers.

### `workflow` IDS

The converter now writes `workflow` when the target DD provides it.  It records:

- TRANSP as the workflow component;
- the TRANSP version when exposed by the CDF;
- the complete TRANSP parameter XML on the component;
- TGLF and GLF23 commit/hash information as workflow libraries when available;
- the end of the TRANSP simulation time loop;
- normal converter `code` provenance for the workflow IDS itself.

`workflow_cycle` is deliberately left empty: a post-processed TRANSP CDF does
not retain the actual component execution/freeze schedule at every internal
workflow cycle.

### `dataset_fair` IDS

The converter now writes `dataset_fair` when available.  It intentionally does
**not** invent a persistent dataset identifier, license, rights holder,
replacement chain, or validity interval.  Those values require an authoritative
catalog/publication source and are not present in the CDF/TR.DAT inputs.

The TRANSP software DOI is likewise not used as the dataset PID: it identifies
software/reference material, not this particular converted IMAS dataset.

### `summary` metadata

For DD 4.1+ the `summary` IDS is additionally identified as simulation data and
stores:

- `type = simulation`;
- machine/device name;
- pulse/shot number when numeric;
- simulation description;
- TRANSP relative simulation `time_begin` and `time_end`;
- `time_step` only when the exported time base is demonstrably uniform;
- a workflow note pointing to the detailed `workflow` IDS and `code.parameters`.

Experimental wall-clock pulse timestamps are not inferred from CDF file times.

### v8.6.2 workflow/provenance regression fix

The workflow writer now uses the public converter helpers
`get_ncattr_text()` and `build_code_parameters_xml()`.  This fixes two v8.6.0/
v8.6.1 API mismatches: the nonexistent private `_decode_attr()` call (which
silently suppressed TRANSP/TGLF/GLF metadata) and the nonexistent
`build_transp_parameters_xml()` call (which aborted conversion).

A static writer-API verifier is included so future `self.im.<method>()` calls
that are not implemented by `nc2imas` or its base classes are caught before
runtime.


## TRANSP2IMAS IMAS-Python migration and conversion fixes

### IMAS-Python API migration

Updated the converter from the legacy Access-Layer Python calling pattern to the current IMAS-Python API:

- `imas.imasdef.*_BACKEND` -> `imas.ids_defs.*_BACKEND`
- empty IDS construction via `DBEntry.factory.new(ids_name)`
- writes via `DBEntry.put(ids, occurrence=...)`
- `DBEntry.create()` / `open()` handled by exceptions rather than `(status, idx)` return values
- Data Dictionary version is optional and read from `IMAS_VERSION` or `IMAS_DD_VERSION`; otherwise IMAS-Python's default is used

The pulse/run `DBEntry` constructor is intentionally retained. It is supported by current IMAS-Python and preserves the existing transp2imas CLI and database layout.

### Dmitriy Yadykin report

#### 1. Ip / Bt direction and COCOS 17

The old converter used a machine-specific `sgn=-1` shortcut for ITER and applied
it inconsistently.  That is not sufficient for DD4.  IMAS Data Dictionary major
version 4 uses COCOS 17, so v8 now derives the physical directions from the
TRANSP namelist controls `NLJCCW` (plasma current) and `NLBCCW` (toroidal
field), with optional `--ip-direction` / `--bt-direction` overrides.

For DD4 the converter writes:

- `Ip`, toroidal current-density profiles and loop voltage with the `NLJCCW` sign;
- `B_phi`, `B0`, `F=RB_phi` and toroidal flux with the `NLBCCW` sign;
- `q` with `sign(Ip)*sign(Bphi)`;
- COCOS17 poloidal flux with `-sign(Ip)` (COCOS11 uses the opposite psi sign);
- poloidal magnetic-field components with the plasma-current sign;
- toroidal rotation quantities stored by TRANSP relative to plasma current with the plasma-current sign.

The sign policy is documented in the COCOS 17 sign-policy section and can be checked
against a CDF/TR.DAT pair with `examples/verify_cocos17.py`.  PF active/passive circuit
currents are no longer multiplied by the plasma-current sign because their sign
depends on independent circuit/winding orientation.

The historical 174783Q32 audit input contains `NLJCCW=.true.` and
`NLBCCW=.false.`, so expected DD4 signs are `Ip>0`, `Bphi<0`, `q<0` and
negative outward poloidal flux relative to a zero axis offset.

#### 2. Current-density normalization

TRANSP `CUR`, `CUROH`, and `CURBS` are in A/cm^2. The old `core_profiles` code used `1.0e-4`, which is the cm^2 -> m^2 conversion for an area, not the inverse-area conversion for a current density.

Correct conversion:

`1 A/cm^2 = 1.0e4 A/m^2`

Changing `1.0e-4` to `1.0e4` raises the output by exactly `1.0e8`, matching the reported discrepancy.

#### 3. Equilibrium/core-profile volume

TRANSP `DVOL` is a per-zone volume. IMAS `profiles_1d.volume` / `grid.volume` is the volume enclosed inside the magnetic surface. The converter now writes `cumsum(DVOL)` in SI units.

The same semantic mismatch existed for TRANSP `DAREA` (per-zone poloidal cross-sectional area) versus IMAS `area` (enclosed cross-sectional area), so `cumsum(DAREA)` is also written. Differential `DVOL` and `DAREA` are still retained internally for integrations.

### Verification performed

Using the historical `174783Q32.CDF`:

- `cumsum(DVOL)[:,-1]` agrees with `PVOL` to about `5.3e-6 m^3` maximum absolute error.
- `cumsum(DAREA)[:,-1]` agrees with `PAREA` to about `7.2e-7 m^2` maximum absolute error.
- old/new current-density scale ratio is exactly `1.0e8`.
- modified Python converter modules pass `py_compile`.
- the IMAS-Python DBEntry/factory/put path was exercised with an API-shaped mock.

A full end-to-end IDS write was not run in this environment because the IMAS-Python package and IMAS-Core backend are not installed here.

### DD 4.1.1 follow-up: TF rename and float precision

- `tf.b_field_tor_vacuum_r` was renamed to `tf.b_field_phi_vacuum_r` in newer Data Dictionaries. The converter now detects the available node and supports both names.
- The TF vacuum-field `data` and its node-local `time` coordinate are both filled.
- The physical toroidal-field sign (`nc2imas.bt_sign`) is applied to `BZXR` in the TF IDS, consistent with equilibrium/core_profiles.
- Floating-point values read from TRANSP CDF are promoted from `float32` to `float64` in `get_ncvar()`. This removes IMAS-Python automatic conversion messages without changing integer/string CDF data.


### DD 4.1.1 follow-up: NBI species and generic species names

- `nbi/unit/species/label` was renamed to `.../name` in DD 3.42.  The
  converter now uses the current `name` field (with a compatibility helper for
  older DDs) and strips the internal `_BEAM` suffix from names such as
  `D_BEAM`.
- `z_n` is an integer node in DD4; nuclear charge assignments are now written
  as integers.
- Applied the same `label` -> `name` compatibility handling to ion, ion-state,
  neutral, core-source/radiation helper paths, and edge-source neutral species
  so these do not fail later in the same conversion.
- Fixed the pre-existing `identifierindex` typo in `imas_edge_sources.py` to
  `identifier.index`.

### DD 4.1.1 follow-up: core_profiles velocity coordinates

IMAS-Python DD4 validation enforces the coordinate of
`core_profiles.profiles_1d[].ion[].velocity.{radial,poloidal,toroidal}` to be the
same radial grid as `profiles_1d[].grid.rho_tor_norm` (or an allowed alternative).
For the historical 174783Q32 CDF, TRANSP stores `X` on 50 points, but
`VPOL*_NC` and `VTOR*_NC` on the 101-point `RMAJM` equatorial-midplane grid.
Those quantities cannot be copied directly into the IMAS 1-D flux profile.

The converter now shape-checks these assignments against the IMAS radial grid.
The 101-point `RMAJM` poloidal/toroidal velocities are therefore intentionally
left empty instead of writing an invalid or physically ambiguous profile.
`rotation_frequency_tor` remains populated from TRANSP `OMEGA`/related X-grid
quantities, and the radial velocity is retained when it already matches the
radial grid. A future enhancement may map the midplane velocities to a clearly
defined IMAS representation, but that requires an explicit physics convention
rather than a blind interpolation.

`pf_active.time` is also assigned through the numeric-array API with an explicit
`float64` dtype to avoid DD4 automatic-conversion diagnostics.


### DD 4.1.1 v5
- Explicitly skip VPOLD_NC/VTORD_NC when their RMAJM grid length does not match the core_profiles radial grid; this no longer depends on nc2imas.fill_data(chk=...).
- Assign pf_active.time as an explicit float64 ndarray only when a real time coordinate exists; skip pf_active entirely when no KK/VV/VS source and no authoritative geometry are available (no sentinel).
- Added a startup patch-level print (`DD411-v5`) so the active checkout can be verified immediately.


### DD411-v6 equilibrium field-name migration

- Updated `equilibrium/time_slice/global_quantities/magnetic_axis/b_field_tor` to the DD 3.42+/DD4 name `b_field_phi`, with fallback for older DDs.
- Updated `equilibrium/time_slice/profiles_2d/b_field_tor` to `b_field_phi`, with fallback for older DDs.
- Prepared the currently-disabled `waves` converter for `n_tor` -> `n_phi`.

- Updated `equilibrium/time_slice/global_quantities/psi_axis` to the DD4 preferred `psi_magnetic_axis`, with fallback for older DDs.
- Mapped the removed `boundary_separatrix/psi` representation to DD4 `boundary/psi`, while retaining an old-DD fallback.

### DD 4.1.1 v7

- Migrated `imas_edge_sources.py` IDS construction from the legacy `self.im.ids.edge_sources` access pattern to `self.im.new_ids('edge_sources')`, which uses `DBEntry.factory.new`.
- Added an explicit `edge_sources` stage marker and patchlevel `DD411-v7` to `convert.py`.

### DD 4.1.1 v8: explicit COCOS17 audit

- Added central sign policy based on `NLJCCW` and `NLBCCW`, with CLI overrides.
- DD4 target is COCOS17; DD3 fallback is COCOS11.
- Separated Ip-like, Bphi-like, q-like, phi-like and psi-like signs instead of using one generic `sgn`.
- Fixed equilibrium toroidal flux: TRANSP `TRFLX` is already in Wb and must not be multiplied by `2*pi`.
- Made core/equilibrium poloidal flux absolute by adding the signed `PSI0_TR` axis offset after converting Wb/rad to Wb.
- Applied q sign consistently to q-axis/q95/profile values.
- Applied independent signs to 2-D `B_R/B_Z` (Ip-like) and `B_phi` (Bt-like).
- Applied the toroidal-rotation sign to thermal, impurity, beam and fusion-ion rotation frequencies.
- Updated inverse-equilibrium moments for the DD4 clockwise polar angle and identified the grid as `inverse_rhotornorm_polar` (index 33).
- Removed the plasma-current sign multiplier from PF active/passive circuit currents.
- Corrected `core_profiles.beta_tor_norm` scaling to use Ip in MA (`1e8` factor for beta fraction, Ip in A).
- Added the COCOS 17 sign-policy section and `examples/verify_cocos17.py`.


#### DD411-v8 COCOS17 sign audit

- DD4 output is explicitly treated as COCOS 17. Physical plasma-current and
  toroidal-field directions are read from `NLJCCW` and `NLBCCW`; command-line
  overrides are available for incomplete historical namelists.
- Converted signs are applied consistently to Ip/current density, Bphi/F/Phi,
  q, poloidal flux, loop voltage, and TRANSP toroidal rotation (which is stored
  positive co-current). PF circuit currents are intentionally not tied to Ip.
- `TRFLX` is written directly in Wb (the previous extra `2*pi` was removed).
- `beta_tor_norm` now includes the correct A-to-MA factor.
- X/XB radial-grid semantics were tightened: core/source profiles remain on X
  zone centres while boundary flux/geometric quantities are interpolated from
  XB; equilibrium profiles use the XB flux-surface grid and X quantities are
  mapped to it.
- The disabled `summary` path has the same Ip/B0/q sign treatment where obvious;
  `core_sources` remains disabled/unverified and is intentionally outside this v8 audit.

- DD4 core-profile current semantics are corrected: `CUR` is written to `j_phi`, while parallel `j_total`/`j_ohmic` are derived from `PLJB`/`PLJBOH`; bootstrap uses the selected `CURBS` with `CONPLJB/BZ`.

### DD411-v8.2 compatibility fix

- Restored `nc2imas.set_name_or_label()`, used by NBI, core_profiles, and edge_sources species metadata.
- The helper writes DD4 `name` and falls back to legacy `label` when running with older data dictionaries.
- v8.2 drop-in includes all modules that depend on this helper to avoid mixed-version installations.

### v8.3: explicit dual-DD schema support

- Added `--dd-version` (`auto`, `3.39.0`, `4.1.1`, etc.).  The selected value is
  passed to `imas.DBEntry`, and every IDS is created with that DBEntry's
  factory.  Startup output reports the exact factory version.
- Added version-tuple helpers (`dd_at_least`, `dd_before`, `is_dd3`, `is_dd4`)
  and centralized DD-dependent text/type mappings in `nc2imas`.
- Reworked EC launcher identity and polarization semantics:
  - DD3.39: short ID -> `beam.identifier`, description -> `beam.name`, legacy
    `beam.mode` retained (dominant mode inferred at O-mode fraction 0.5); beam-level
    `o_mode_fraction` is used only if present.
  - DD4: short ID -> `beam.name`, description -> `beam.description`; removed
    `beam.mode` is never accessed.  DD4.1+ writes `RFMODECH` to
    `beam.polarization.o_mode_fraction` with `beam.polarization.time`.
- NBI unit identity, species/state `label`/`name`, and `z_n` scalar type are
  selected according to the target DD.
- Added explicit DD-aware mappings for core-profile `j_tor`/`j_phi`, TF
  `b_field_tor_vacuum_r`/`b_field_phi_vacuum_r`, equilibrium toroidal-field
  names, axis/boundary psi and normalized beta, wall limiter closure, and the
  disabled waves `n_tor`/`n_phi` path.
- DD major version also selects the target coordinate convention: DD3 ->
  COCOS11, DD4 -> COCOS17.  The existing `NLJCCW`/`NLBCCW` physics sign policy
  is preserved.
- Added the Data Dictionary compatibility section and `examples/verify_dd_compat.py`.  The verifier
  creates schema objects only and can check 3.39.0 and 4.1.1 without writing a
  database.

### v8.4: CDF coverage and code-parameter provenance

- Audited the historical `174783Q32.CDF` (2465 variables) and added high-confidence
  mappings that were missing from the active IDS writers.
- `core_profiles` now adds:
  - electron thermal pressure and aggregate ion/thermal/perpendicular/parallel pressure,
  - total thermal ion density and `t_i_average`,
  - per-species thermal ion pressure when thermal density/temperature are available,
  - magnetic shear (`SHAT`), parallel conductivity (`1/ETA_USE`),
    electrostatic potential (`VRPOT`, XB -> X), and parallel electric field
    (`ETA_USE * j_ohmic`).
- Corrected `NI` from N/cm^3 to m^-3 before using it in IMAS profiles/global averages.
- Removed the old ad-hoc radial-electric-field reconstruction; `phi_potential` is
  retained, while `E_r` is left empty until an explicit physical metric is chosen.
- `equilibrium/profiles_1d` now also writes `psi_norm`, metric `rho_tor`,
  `dpsi_drho_tor`, `dpressure_dpsi`, `f_df_dpsi`, toroidal/parallel current
  density, magnetic shear, and mass density (`PMASS`, g/cm^3 -> kg/m^3).
- Rebuilt and enabled the `summary` IDS.  Energy densities `UTOTL`, `UTHRM`,
  `UE`, and `UI` are volume-integrated directly (the old extra factor 1.5 was
  removed).  Fast parallel/perpendicular energies follow the IMAS pressure
  definitions.  Current components, ohmic/radiated power, confinement metrics,
  and neutron-rate decomposition are populated when present.
- `code.parameters` is now an XML provenance payload containing:
  - target DD/COCOS/sign metadata,
  - TRANSP CDF build/version/DOI/model hashes,
  - the complete TRANSP input (`TR.DAT`) as escaped source text and parsed XML,
  - every directly referenced `*TEMPLATE` file when present,
  - adjacent run-specific `.nml` files when present,
  - SHA-256 hashes and explicit `missing` status for unresolved references.
- Added `examples/verify_cdf_extension.py` and expanded `examples/verify_dd_compat.py` for the new
  schema nodes.

### v8.5.1 hotfix

- Fixed `netcdf_connection.check_var_name()` so it returns a boolean instead of comparing a list with an integer. This unblocks `core_sources` and `core_transport` variable-presence checks under Python 3.
- Removed a duplicated `SBTH` read in the NBI source builder.

### v8.5.2 transport/radiation

- Added a separate `core_transport` `transport_solver` model (identifier 2) for
  predictive `EPTR_MOD`, `EETR_MOD`, `IETR_MOD`, and species-resolved
  `PTR*_MOD` fluxes. The existing `combined` observed/effective model is left
  unchanged.
- `IPTR_MOD` is deliberately not forced into a pseudo-species because
  `core_transport` has no aggregate total-ion particle-flux node.
- Rebuilt and enabled the dedicated `radiation` IDS for `PRAD_BR`, `PRAD_CY`,
  `PRAD_LI`, `P0RECO`, and recognized `PRLS_*`/`PRLX_*` impurity detail.
  Emissivities are positive emitted power and `power_inside` is now a true
  cumulative `emissivity*DVOL` integral.
- Fixed the generic `nc2imas.operate_arrays(..., '+', ...)` branch, which
  incorrectly performed subtraction.

### v8.5.3 edge sources

- Rebuilt `edge_sources` around TRANSP **source-rate** profiles rather than
  neutral state/density profiles.  Gas-flow and recycling ionization now use
  `SESGF`/`SISGF_*` and `SESRC`/`SISRC_*` (with species-resolved aliases when
  available), with matched electron, ion, and consumed-neutral particle terms.
- Added ionization power and toroidal-momentum transfer from `PISGF`/`PISRC`
  and `TQISGF`/`TQISRC`, including opposite neutral terms where the same
  interaction is represented.
- Added recombination-neutral production from `S0RECO`, `P0RECO`, and
  `TQ0RECO`.  Real-CDF inspection established that `P0RECO` is kinetic power
  carried by the recombination-neutral source, **not** recombination-radiation
  emissivity.  The earlier v8.5.2 `radiation` interpretation of `P0RECO` was
  removed.
- `N0*`, `T0*`, and `OM0*` neutral quantities are retained as state/profile
  information and are no longer mislabeled as `edge_sources` source rates.
- Added direct CDF numerical checks for source integrals and alias consistency.
  On 174783Q32, integrated `SISGF_D` agrees with `GASD` to better than 5e-7
  relative error.

### v8.5.4 waves

- Rebuilt and enabled `waves` for supported EC and IC absorbed/deposited RF
  outputs.  EC maps `PEECHj`, `ECCURj`, and `FREQECH`; IC maps `RFPWRABj`,
  `FREQAj`, and, for a single coherent wave, aggregate `PEICH`, `PIICH`, and
  `QRFBTOT` deposition.
- Wave power-inside is calculated as a cumulative finite-volume integral using
  `DVOL`.  EC toroidal current-inside is calculated from `ECCURj*DAREA` with
  the configured COCOS current sign.
- `ECCUR` is not written to `current_parallel_density`; that IMAS quantity
  requires `<j.B>/B0`, not TRANSP toroidal current density.
- The writer does not create a `waves` IDS for launcher-only or RF-summary-only
  runs.  The supplied 174783Q32 validation run has no supported EC/IC
  deposition payload, so this behavior is exercised explicitly.
- IC antenna names and time bases were aligned with coherent-wave references.
- Fixed the `nc2imas.get_namelist_var()` compatibility wrapper so it calls the
  instance method rather than an undefined global function.
- Added the waves mapping section and `examples/verify_next_ids_inventory.py` to report whether
  a CDF contains enough wave, full-FBM, or sawtooth-event data to exercise the
  next IDS mappings.

### v8.5.5 neutral-process decomposition and attached-run validation

- Re-audited `edge_sources` against the supplied `174783Q32.CDF` and split
  neutral interactions by physical process.  Gas-flow and recycling
  electron-impact ionization are now `atomic_ionization` (707), while
  `SFCX*`/`PCX*`/`TQCX*` are emitted as `charge_exchange` (305).
- Corrected the previous mixed-process interpretation of `SIGF_*_*` and
  `SIRC_*_*`.  On the supplied CDF these arrays satisfy
  `SI = SE + SFCX`; they are now used as consistency diagnostics rather than
  being written wholesale as an ionization-only source.
- Added charge-exchange energy and toroidal-momentum transfer with the signs
  stated by the CDF: `PCX*` and `TQCX*` are positive **to neutrals**, with the
  opposite aggregate-ion sink.
- Made recombination particle stoichiometry species-resolved when the CDF has a
  unique `S0V0G`/`SIV0G` isotope.  For 174783Q32 this is D; recapture is split
  into `SERECO` ionization and `SFCXRECO` charge exchange.  No unresolved plasma
  energy sink is invented for `P0RECO`.
- Extended `examples/verify_sources_transport.py` with the TRANSP identities and neutral
  energy/momentum-balance equations.  On 174783Q32, `SIGF_D_D` equals
  `SEGF_D+SFCXGF_D` to ~1e-12 relative, `SIRC_D_D` equals
  `SERC_D+SFCXRC_D` to ~1e-6, and their volume integrals agree with `GASD` and
  `RCYD` to better than 5e-7 relative.
- Improved single-impurity line-radiation metadata.  Generic `PRLS_TOK` /
  `PRLX_SINGL` profiles are no longer discarded when `AIMP`/`XZIMP` uniquely
  identify the species.  The supplied namelist has `AIMP=12`, `XZIMP=6`, so
  the identical `PRLS_TOK=PRLX_SINGL=PRAD_LI` profile is attached to carbon,
  without treating `TOK` or `SINGL` as element names or double counting the
  aliases.
- Re-ran the next-IDS inventory on the supplied inputs.  The run has no
  supported coherent-wave deposition payload, no full FBM phase-space payload,
  no sawtooth crash-event payload, and no `PS_*` passive-structure currents.
  Those IDSs are therefore skipped/conservatively deferred for this validation
  case rather than populated from incomplete diagnostics.


## Next-IDS status from 174783Q32.CDF / 174783Q32TR.DAT

This note records what can and cannot be implemented responsibly from the supplied
TRANSP run after the v8.5.5 `edge_sources`/radiation audit.

### edge_sources: implemented and numerically checked

The old `edge_sources` comments mixed three different kinds of quantities.  Only
actual volumetric source terms are written to `edge_sources`.

- Gas-flow atomic ionization uses `SEGF_G`, `PIGF_G`, `TQIGF_G`.
- Recycling atomic ionization uses `SERC_G`, `PIRC_G`, `TQIRC_G`.
- Fast-ion charge exchange is kept separate using `SFCXGF_G`/`PCXGF_G`/`TQCXGF_G`
  and `SFCXRC_G`/`PCXRC_G`/`TQCXRC_G`.
- Recombination neutral production uses `S0RECO`, `P0RECO`, `TQ0RECO`; subsequent
  recapture is split into `SERECO` ionization and `SFCXRECO` charge exchange.

For the supplied D-only neutral population, the CDF directly verifies

- `SIGF_D_D = SEGF_D + SFCXGF_D` to about 1e-12 relative;
- `SIRC_D_D = SERC_D + SFCXRC_D` to about 1e-6 relative;
- the volume integrals of those combined ion sources reproduce `GASD` and `RCYD`
  to better than 5e-7 relative;
- `SIV0D = SERECO + SFCXRECO` to about 9e-7 relative.

The following legacy-comment quantities are *not* forced into `edge_sources`:

- `N0SGF_*`, `N0SRC_*`, `N0GF_*_*`, `N0RC_*_*` are neutral **densities**, not
  source rates.  The existing total wall/volume-neutral representation in
  `core_profiles` remains the proper starting point for neutral-state output;
  source-origin-resolved neutral populations need a separate state/profile audit.
- `T0RECYC`, `T0GASFL`, `OM0RECYC`, `OM0GASFL` are source/boundary neutral-state
  descriptors.  `edge_sources` has particle/energy/momentum source terms rather
  than a generic source-temperature/source-rotation field.  The spatial `P*` and
  `TQ*` source profiles are therefore exported instead of fabricating a state node.
- `PEEDG` and `PIEDG` are electron/ion energy terms "via BDY".  They are boundary
  transport terms, not volumetric source profiles, and are not relabelled as an
  `edge_sources` source.
- `GASIZ`/`GASIN` are not present in this CDF.  `GASD` and `RCYD` are present and
  are used as independent integral checks of the spatial source profiles rather
  than duplicated as a second source representation.

### waves: writer enabled, this run correctly produces no waves IDS

`EC_MODEL='NONE'`, `NLLH=.false.`, and the CDF has no supported EC/IC absorbed
wave payload (`PEECH*`, `ECCUR*`, `RFPWRAB*`, `FREQA*`, `PEICH`, `PIICH`,
`QRFBTOT`, `XGRID_NPHI`).  `PRFB_D` and `PRFFI` are summary RF-to-fast-ion powers,
not a coherent-wave spatial solution.  The v8.5.4+ writer therefore skips this run.

### fast-ion distributions: input is incomplete

`NLFBM=.true.` and the main CDF contains `BMVOL` and `BDENS2_D`, plus reduced
midplane moments, but it contains no `FBM`, `EBA2PL`, or `EBA2PP`.  A phase-space
distribution must not be reconstructed from these reduced moments.

For implementation/validation of the distributions IDS, provide a TRANSP fast-ion
ACFILE/get_fbm output containing the full distribution and its pitch/energy/spatial
grid description (typically the generated `*_fi_*.cdf` product).

### sawteeth: no event payload in this run

`NLSAW=.false.`.  `TESAW` is a scalar sawtooth parameter and `XIQ1(TIME)` is q=1
surface context, but there are no crash times or mixing-radius/pre/post-crash arrays.
The sawteeth IDS remains disabled rather than writing an empty or inferred event.
A sawtooth-enabled run with actual crash output (and, if prescribed, the PRESAW/
EXTSAW event input) is required for validation.

### pf_passive: implemented with explicit machine-data provenance

None of the bundled CDFs contains `PS_1...PS_n` passive-current channels, so
**real passive-current physics remains unvalidated**. The current writer detects
`PS_n` using strictly matched numeric indices, retains its physical-current
sign independent of plasma `Ip`/`Bphi`, converts supported A/kA/MA units, and
uses each time series' native timebase in `pf_passive/loop[i]/time`.

In the legacy mode (`--pf-passive-geometry` absent), it writes any available
`PS_n` histories as current-only loops and warns that the resulting IDS lacks
authoritative mechanical/electrical properties. With an external, audited
machine-specific JSON description it fills geometries, signed element turns,
optional area/resistance/resistivity and trace-to-device mappings. Current sign
conversion in this mode is **explicit**, not inferred from the CDF label.
All six geometry types are supported; outline vertices are normalized open
for DD3 and explicitly closed for DD4. The JSON source, SHA256 and reference
are recorded in IDS provenance. A full geometry-only IDS can also be written
when `PS_n` output is absent.

The portable `make verify-pf-passive` suite uses synthetic CDFs and manifest
values only; no NSTX, DIII-D or C-Mod geometry is invented. Refer to the
[full passive-conductor mapping guide](pf-passive.md) for schema, sign and time
conventions, and the requirements for real-machine validation.

### Remaining inputs that would unlock the next implementation steps

1. A full `*_fi_*.cdf` (or equivalent ACFILE/get_fbm product) for distributions.
2. A run with `NLSAW=.true.` and saved sawtooth crash/mixing output for sawteeth.
3. A run containing `PS_n` passive currents plus the matching authoritative machine
   description / `pf_passive` geometry if a complete passive-loop IDS is desired.
4. The bundled DIII-D and C-Mod EC/IC cases now validate `waves`; no further
   RF fixture is required for the already supported coherent-wave quantities.

### v8.6.0 metadata/provenance and additional saved moments

This run does contain enough information for two additional physics mappings
that do not require a full FBM:

- `PPHI` -> total `core_profiles` toroidal momentum density.  Its identity with
  `AMOI*OMEGA` is verified at ~2.2e-7 relative RMS on this CDF.
- per-NBI-unit BDEP birth particle/energy source and `TQTOTnn` torque ->
  `distribution_sources`, with `BPSHnn_*` shinethrough bookkeeping.

Dataset provenance is now explicit at three levels: common per-IDS
`ids_properties/provenance`, detailed original-input XML/text in
`code.parameters`, and dataset-level `workflow` + `dataset_fair` IDSs.  Missing
FAIR catalog values (dataset PID, license, rights holder, replacement chain)
remain blank rather than being guessed.


## TRANSP -> IMAS `radiation` mapping (v8.5.5)

The dedicated `radiation` IDS writes **positive emitted power/emissivity**. The
corresponding radiation losses in `core_sources` retain the opposite sign.

| IMAS radiation process | TRANSP CDF profile | Spatial mapping |
| --- | --- | --- |
| bremsstrahlung (8) | `PRAD_BR` | `electrons.emissivity`; cumulative volume integral -> `electrons.power_inside`; final integral -> global emitted power |
| synchrotron_radiation (9) | `PRAD_CY` | same electron mapping as bremsstrahlung |
| line_radiation (10) | `PRAD_LI` | aggregate line emissivity in `emissivity_ion_total`; cumulative volume integral -> `power_inside_ion_total`; final integral -> global emitted power |

TRANSP stores these profiles in `W/cm^3`; the converter multiplies by `1e6` to
obtain `W/m^3`, and integrates once over `DVOL` converted from `cm^3` to `m^3`.

### Species-resolved line radiation

Recognized `PRLS_<element>` profiles are attached to the corresponding ion
species. Recognized `PRLX_<element>_<charge>` profiles are attached to
charge-state entries.

Older/single-impurity output can instead use generic model labels
`PRLS_TOK`/`PRLX_SINGL`.  These strings are **not** chemical species.  v8.5.5
now resolves them only when the TRANSP namelist contains a unique single
impurity `AIMP`/`XZIMP` pair.  The pair is matched against the converter's
periodic-table data; only one generic alias is used, preferring `PRLS_TOK`, so
identical aliases are never double counted.

For `174783Q32TR.DAT`, `AIMP=12` and `XZIMP=6`, which resolves to carbon.  In
the attached CDF, `PRLS_TOK`, `PRLX_SINGL`, and `PRAD_LI` are numerically
identical.  The line-radiation process can therefore carry a carbon-resolved
ion entry without pretending that `TOK` or `SINGL` is an element name.

If `PRAD_LI` is absent, recognized species-resolved line profiles can be summed
to form the aggregate line emissivity.

`PRAD_LI` itself does not encode whether every line photon originates from an
ion or a neutral state. The aggregate mapping to `emissivity_ion_total` is
retained as the least-lossy representation of TRANSP impurity line radiation;
richer charge-state CDFs should still be used to validate state-level mapping.

### Recombination-neutral power is not radiation

`P0RECO` is **not** treated as recombination radiation. In `174783Q32.CDF` its
long name is `RECOMB NEUTRAL SCE POWER`; paired arrays `S0RECO` and `TQ0RECO`
show that it is kinetic power carried by a recombination-neutral source. It is
handled by `edge_sources`. A radiation process 11 should only be added when an
explicit radiative recombination-emission quantity is identified.


## TRANSP -> IMAS `waves` mapping (validated for 0.9.0)

The `waves` writer is enabled in v8.5.4, but it is intentionally conservative:
it creates a coherent wave only when the TRANSP CDF contains absorbed/deposited
wave output.  Launcher-only quantities remain in `ec_launchers` or
`ic_antennas` and are not relabelled as absorbed power.

### Electron-cyclotron waves

| TRANSP | IMAS `waves/coherent_wave` | Conversion / rule |
| --- | --- | --- |
| `PEECHj(TIME3,X)` | `profiles_1d/power_density` | W/cm3 -> W/m3 (`x1e6`) |
| `PEECHj(TIME3,X)` | `profiles_1d/electrons/power_density_thermal` | same profile; EC deposition is assigned to thermal electrons |
| volume integral of `PEECHj` | `profiles_1d/power_inside`, `global_quantities/power` | cumulative finite-volume integral using `DVOL` |
| `ECCURj(TIME3,X)` | `profiles_1d/current_phi_inside` (legacy `current_tor_inside`) | cumulative surface integral using `DAREA`; COCOS toroidal-current sign applied |
| outer value of integrated `ECCURj` | `global_quantities/current_phi` (legacy `current_tor`) | A |
| `FREQECH(j)` from TR.DAT | `global_quantities/frequency` | Hz |
| EC solver selected by `EC_MODEL` | `wave_solver_type` | beam/ray-tracing category (index 1) |

For a run with exactly one EC wave, aggregate `PEECH` and `ECCUR` are accepted
as fallbacks when per-launcher profiles are absent.  Aggregate profiles are
never duplicated across multiple gyrotrons.

`ECCUR` is TRANSP toroidal current density in A/cm2.  It is **not** written to
IMAS `current_parallel_density`, whose definition is `<j.B>/B0`.  Only the
cumulative toroidal current fields are populated.

### Ion-cyclotron waves

| TRANSP | IMAS `waves/coherent_wave` | Conversion / rule |
| --- | --- | --- |
| `FREQAj(TIME)` | `global_quantities/frequency` | interpolated onto `TIME3` |
| `RFPWRABj(TIME)` | `global_quantities/power` | absorbed wave power, interpolated onto `TIME3` |
| `PIC_Fj(TIME3,X)` | `profiles_1d/power_density` | frequency/antenna-resolved absorbed RF power; W/cm3 -> W/m3 |
| cumulative integral of `PIC_Fj` | `profiles_1d/power_inside` | using finite-zone `DVOL` |
| `QICHA(TIME3,X)` | fallback `profiles_1d/power_density` | only when exactly one IC coherent wave exists and `PIC_Fj` is unavailable |
| `ICCUR(TIME3,X)` | cumulative `current_phi_inside` / legacy `current_tor_inside` | only when exactly one IC coherent wave exists; aggregate current is not duplicated across multiple waves |
| `XGRID_NPHI` | `n_phi` (legacy `n_tor`) | only for a single IC coherent wave because the CDF grid is not antenna-resolved |

IC uses wave identifier `IC` index **3** and the full-wave solver category
index **2**.  The antenna reference is `antenna{j}`, matching the IC antenna
writer.

`PEICH`, `PIICH`, and `QRFBTOT` are not interpreted as additive wave-absorption
partitions.  In the C-Mod validation case their sum greatly exceeds the
absorbed antenna power because these profiles describe plasma heating/source
channels after RF/Fokker-Planck processing.  They remain mapped through
`core_sources`, where their source semantics belong.

For multi-frequency IC runs, a CDF may contain an aggregate `XGRID_NPHI` and
`PIC_NPHI` spectrum without a frequency-resolved spectral decomposition.  The
writer does not copy that same aggregate spectrum onto every coherent wave.  A
single-wave case may still use the aggregate `XGRID_NPHI` representation.

### Time base and grid

The writer uses `TIME3` for all coherent-wave globals and profiles and marks the
IDS as homogeneous in time.  Scalar `TIME` traces are interpolated to `TIME3`.
The radial grid is the TRANSP `X` finite-volume grid and reuses the standard
`fill_profiles1d_grid()` geometry/flux mapping.

### RF regression results

`examples/verify_waves_mapping.py` provides the source-side regression contract
for the mapping and is run by `make verify-waves`.

For DIII-D `147634B13` (TORBEAM) and `147634G14` (GENRAY):

- six `PEECHj` profiles sum to aggregate `PEECH` at numerical precision;
- the finite-volume integral of aggregate `PEECH` reproduces `PECHT` with
  relative error of order `1e-7`;
- integrated `ECCURj` contributions reproduce aggregate `ECCUR` within about
  `2.3e-3` relative error or better.

For C-Mod `89846`:

- `integral(PIC_F1)` reproduces `RFPWRAB1` within `2.3e-7` relative error;
- `integral(PIC_F2)` reproduces `RFPWRAB2` within `3.8e-7` relative error;
- `PIC_F1 + PIC_F2` reproduces aggregate `QICHA` with normalized maximum error
  below `5e-7`;
- the CDF contains two frequencies (80 and 78 MHz) and an aggregate
  `XGRID_NPHI/PIC_NPHI` spectrum, so per-frequency `n_phi` is intentionally not
  fabricated.

### Not yet mapped

- Lower-hybrid coherent-wave output.  No suitable LH payload is present in the
  supplied validation cases, so no unvalidated mapping is enabled.
- Parallel current density from EC current drive, because it requires the
  proper `<j.B>/B0` quantity rather than direct reuse of `ECCUR`.
- Frequency-resolved IC `n_phi` spectra when only aggregate `PIC_NPHI` output
  is present.

## Active PF circuits and DD3/DD4 geometry (0.9.0 extension)

`imas_pf_active.py` now delegates to `pf_active_mapping.py`. The writer
collects KK_ and VV_ by union of coil names, independently collects VS_
supplies, and records each waveform's native 1D time and source units.
The old bug writing a KT_ *target current* to the supply-voltage array
is removed; KT_ is not a measured current and is reserved for later
`pulse_schedule` mapping. CC_ equilibrium-current output is used only
when explicitly bound by an authoritative machine manifest. The manifest
provides geometry, resistance, turns/sign, and optionally supply current
relations and DD3/DD4 circuit connectivity. Without it the IDS is
electrical-output-only and explicitly incomplete for field reconstruction.
The previously bundled six general CDFs lack active PF coil currents. Two
additional real ISolver CDFs now provide active currents and passive loops:
ITER 38530P94 (11 coil circuits, 119 passive channels) and NSTX 142202U01
(9 coil circuits, 82 passive channels). The new `verify-isolver` source regression
checks source/time identity and CC/KK integer-like turns, while authoritative
machine geometry remains a separate requirement.
See [pf_active details](pf-active.md).

### Per-model input provenance (0.9.0 extension)

In addition to the full original TR.DAT and referenced model input XML in
`workflow` and `code.parameters`, `workflow.time_loop.component[].parameters`
and `waves.code.library[].parameters` now separately retain configured
TORIC/GENRAY/TORBEAM/CQL3D/FPP/LSC model input sets, including hashes and
explicit missing-file records. Transport and NUBEAM model input sets are
associated with their respective IDS code libraries. See
`docs/input-provenance.md` for resolution rules and safety limits.

## Neoclassical and neutral state mapping (2026-10)

The new CDF-to-IMAS mapping explicitly treats neoclassical particle/energy
quantities as **divergences** on TRANSP `TIME3, X`, and volume-integrates
them to boundary fluxes in SI. `QFLNCC_*` is kept separate from `QFLNC_*`
so no neoclassical/classical double counting is introduced. Aggregate `_X`
and `_TOK` channels are not treated as physical ions; `_I` thermal-ion heat
is represented as `total_ion_energy.flux`.

Neutral `DN0W*`, `DN0V*` and `N0V0_*` represent distinct origins rather than
a single density. `T0*` is already eV and is no longer incorrectly scaled by
1e6. `OM0*` is in rad/s; DD4.1.1 lacks a typed neutral rotation field, so
that information is preserved in supplementary `core_profiles.code.parameters`
XML with source channel, grid and physical sign. See
[the detailed mapping specification](neoclassical-neutrals.md).



## MMM07 anomalous-transport decomposition (October 2026)

TRANSP's MMM07 `core_transport` outputs are *reduced turbulence-model
predictions*: the combined `XKEMMM07`/`XKIMMM07` coefficients and individual
Weiland `XKEW19`/`XKIW19`, DRIBM `XKEDRBM`/`XKIDRBM`, ETG `XKEETG`, and MTM
`XKEMTM` components. They are not experimental diagnostics. Their CDF
`TIME3/XB` data are converted from cm^2/s to m^2/s without flux integration.
`GAMMAETG`, `GAMMAMTM`, `OMEGA*`, `KYRS*`, `XDBMTM` and untyped aggregate
momentum/impurity coefficients are preserved in the corresponding model's
`code.parameters` XML at native coordinates. See [the precise mapping and
closure tests](mmm-transport.md). No MMM model is fabricated from all-zero
placeholder arrays.
