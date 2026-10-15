# Unified theory transport: MMM07, GLF23, neoclassical, and TGLF

The common `theory_transport_mapping.py` source contract supplies a registry
for four theory families. `imas_core_transport.py` creates separate model
records and never double-counts the total MMM07 and its submodels.

| Family | Source-supported status | Output |
|---|---|---|
| MMM07 | Active in DIII-D 147634T16 | total, Weiland, DRIBM, ETG, MTM |
| GLF23 | Active in ITER 38530P94 | electron/ion heat diffusivity and toroidal momentum diffusivity |
| Neoclassical | Source fields present in all ten CDFs | separate particle and heat flux models from divergences |
| TGLF | No validated TGLF coefficient fixture | reserved family only; no false outputs |

GLF23 mapping uses `XKEGLF23 -> electrons.energy.d`,
`XKIGLF23 -> total_ion_energy.d`, and `ETPHIGLF -> momentum_phi.d`
(`momentum_tor` in older DD3) with cm²/s -> m²/s conversion on `TIME3/XB`.
`DIFFIGLF` is explicitly an **aggregate ion particle diffusivity** and remains
supplementary GLF23 model XML until species-specific interpretation is known.
It is not copied into an arbitrary ion `particles.d` field. The XML references
native IDS time and `grid_d` coordinates; no duplicated time/grid arrays.

GLF23 channels are zero-filled placeholders in other fixtures. A nonzero
`TGLF_DR` scalar is not sufficient to infer that a TGLF transport model was
run, because its source meaning is interchange-instability location. Future
TGLF coefficient mappings will be added once representative files and input
namelists are available.

See [neoclassical theory transport](neoclassical-transport.md) for the
volume-integrated flux mapping. Neutrals belong independently in
[core_profiles](neutral-profiles.md).

## MMM07 reduced anomalous components

The TRANSP Multi-Mode Model (MMM) combines **reduced anomalous transport
models**. The `GAMMAMTM`, `OMEGAMTM`, `GAMMAETG`, etc. variables are predictions
from those reduced models, **not measured diagnostics**. The primary validated
example is DIII-D `147634T16.CDF`; the adjacent `mmm.nml` identifies switches
for Weiland, DRIBM, ETGM and MTM and the selected model parameters. Because
`147634T16TR.DAT` is absent, the adjacent `mmm.nml` is retained as a nearby
input file but must not be represented as conclusively proven to have been
used in the TRANSP execution.

## IMAS representation and units

`core_transport.model[]` holds separate model entries, using **negative private
identifier indices** rather than assigning arbitrary public enumeration indices.
The individual contributions are **components of MMM**, and must not be added
to the combined MMM result a second time. All transport profiles are on the
source `TIME3, XB` grid. `1 cm^2/s = 1e-4 m^2/s`, `1 cm/s = 1e-2 m/s`.

| Model identifier | Physical model | TRANSP CDF | Typed IMAS fields |
|---|---|---|---|
| `mmm07_total` | Combined MMM07 output | `XKDMMM07`, `XKEMMM07`, `XKIMMM07` | `electrons.particles.d`, `electrons.energy.d`, `total_ion_energy.d` |
| `mmm07_total` | Combined MMM07 momentum | `XPTMMM07`, `VPHIMMM07` | `momentum_phi.d`, `momentum_phi.v` (DD4; `momentum_tor` in early DD3) |
| `mmm07_weiland` | Weiland model (CDF label Weiland19, namelist options `W20`) | `XKEW19`, `XKIW19` | `electrons.energy.d`, `total_ion_energy.d` |
| `mmm07_dribm` | Drift-resistive ballooning/interchange model (`DRIBM` in namelist, `DRBM` in CDF) | `XKEDRBM`, `XKIDRBM` | `electrons.energy.d`, `total_ion_energy.d` |
| `mmm07_etg` | Reduced electron-temperature-gradient model | `XKEETG` | `electrons.energy.d` |
| `mmm07_mtm` | Reduced microtearing model | `XKEMTM` | `electrons.energy.d` |

Do **not** reconstruct energy *fluxes* from these coefficients without the
corresponding gradients, pinch contributions, and conventions. The new models
therefore populate `d` (and the measured MMM toroidal convection `v`) rather
than inventing heat/particle fluxes. In contrast, the separate neoclassical
`GFLNC_*`/`QFLNC_*` inputs are *flux divergences* and require finite-volume
integration; they are not handled by this MMM map.

## Reduced-model output retained as model-specific XML

IMAS DD4.1.1 has no dedicated typed fields for reduced-model growth rates,
frequencies or the selected `k_y rho_s`. To avoid losing scientifically useful
MMM output, the converter preserves these arrays in
`core_transport.model[].code.parameters`, tagged as **reduced-model output**.
These are not observation arrays and not native IMAS typed transport fields.
The XML records the original variable name, physical meaning, SI units,
array shape, and **references** to native IMAS coordinates. Coordinate arrays
are no longer duplicated inside the XML.

### Native IMAS grids (version 2 XML)

- `core_transport.time` and each `model[].profiles_1d[it].time` contain the
  native TRANSP `TIME3` values (s).
- `model[].profiles_1d[it].grid_d.rho_tor_norm` contains native `XB` for MMM
  transport coefficients (`d`, and their SI conversion).
- For DRIBM, ETG and MTM components with additional reduced-mode results on
  `X`, `model[].profiles_1d[it].grid_flux.rho_tor_norm` contains **native `X`**.
  This is a secondary grid only: no physical fluxes are fabricated or written
  to the `*.flux` nodes of these models. The grid is provided solely as a
  native coordinate reference for the supplementary reduced-mode arrays.
- The combined MMM model, with only `XB`-based supplementary coefficients,
  retains `XB` for its `grid_flux` as before.

IMAS's `grid_flux` node is designated as the **grid for flux quantities**;
using it as a second radial grid for supplementary reduced-model outputs is a
practical interim convention, explicitly documented here, not a claim that
these model growth rates are transport fluxes. An eventual dedicated typed
reduced-model structure would be preferable if a suitable standard becomes
available. No regridding from `X` to `XB` is performed and the difference in
radial grids is not lost.

The XML `supplementary_model_outputs/coordinate_reference` elements point to
these IDS arrays. The numerical arrays in XML are the **model outputs only**
(`GAMMAMTM`, `OMEGAMTM`, `KYRSMTM`, etc.), not coordinate copies.


| MMM component | Retained extra output | Description |
|---|---|---|
| ETG | `GAMMAETG`, `OMEGAETG`, `KYRSETG` | Reduced ETG growth rate (s^-1), frequency (rad/s), normalized wavenumber |
| MTM | `GAMMAMTM`, `OMEGAMTM`, `KYRSMTM`, `XDBMTM` | Reduced MTM growth rate, real frequency, normalized wavenumber, estimated delta-B/B |
| DRIBM/DRBM | `KYRSDRBM` | Selected normalized wavenumber |
| MMM combined | `XKZMMM07`, `XPPMMM07` | Unresolved aggregate MMM chi(Z) and poloidal-momentum chi(Theta) coefficients (converted to m^2/s) |

The aggregate `XKZMMM07` channel is described in the CDF only as
`MMM07 MODEL CHI(Z)`; its precise transport-equation interpretation is **not
assumed**, and it is **not** assigned to a fictitious ion species. `XPPMMM07` represents aggregate poloidal momentum
transport and has no unambiguous aggregate poloidal-momentum field in
`core_transport`; both are preserved as supplementary XML profiles. Standalone
`XKECDBM`, `XKEGLF23`, `XKEPALEO`, `XKINC`, etc. are **not** claimed as MMM
components. Negative source coefficients (including small Weiland negative
values) are preserved rather than clipped.

## Source closure from the DIII-D T16 fixture

Across the available nine `TIME3` profiles, the component decomposition is:

```
XKEMMM07 ~= XKEW19 + XKEDRBM + XKEETG + XKEMTM
XKIMMM07 ~= XKIW19 + XKIDRBM
```

The maximal residual normalized by the maximal absolute *total* coefficient is
5.42e-7 for electron heat and 3.41e-7 for ion heat (source CDF, cm^2/s).
These are numerical closure tests, **not** a claim that the total modeled
energy flux is simply additive after any subsequent solver corrections. The
converter retains all source profiles unchanged, and its source verifier
checks closure without silently expanding tolerances.

Other bundled examples declare zero-valued MMM placeholder arrays even when
other turbulence models are used; the converter **does not create** MMM models
for those runs. It validates the source dimensions/units of all declared MMM
variables but only exports MMM when data are actually present.

## Verification

```
make check
make verify-theory
make verify-all
make verify-dd
python3 convert.py -i examples/data/d3d/147634T16/147634T16.CDF \
    -b HDF5 -new --dd-version 4.1.1 --ids core_transport
```

`make verify-theory` checks every supplied CDF and writes
`build/theory-transport/summary.json` and `summary.md`. It tests the **actual
writer methods** with an IMAS-shaped test double. `make verify-dd` requires
IMAS-Python. Actual DD4.1.1 write/readback must still be confirmed on a
workstation using `make verify-end-to-end ... E2E_IDS=core_transport`.

## T16 input-file path warning

For `-i examples/data/d3d/147634T16/147634T16.CDF`, the top-level
`147634T16TR.DAT` filename is a path relative to the *current directory*.
`InputInventory` now resolves it from that directory exactly once, instead of
erroneously appending it to the containing case directory again. The original
`TR.DAT` is genuinely absent from the supplied T16 case, so the inventory
retains a **missing** record with the correct absolute path. Adjacent
`mmm.nml` and `t16_pt.nml` remain labelled
`ADJACENT_UNREFERENCED`, since their use by the original run is not proven.

## End-to-end coordinate readback

In an IMAS-enabled environment, `make verify-end-to-end` (or the selective
`examples/verify_end_to_end.py --case 147634T16 --ids core_transport` test)
reads `core_transport` back and compares all stored MMM `grid_d` and
`grid_flux` values with the CDF's native `XB` and `X` arrays.  This includes
all time slices, not just array shape checks.  The source-only tests use
schema-shaped mocks and cannot replace this live DD4.1.1 validation.
