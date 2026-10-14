# MMM07 reduced anomalous transport in `core_transport`

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
The XML records the original variable name, physical meaning, SI units, native
`TIME3/X` or `TIME3/XB` grids and the array shape.

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
make verify-mmm
make verify-all
make verify-dd
python3 convert.py -i examples/data/d3d/147634T16/147634T16.CDF \
    -b HDF5 -new --dd-version 4.1.1 --ids core_transport
```

`make verify-mmm` checks every supplied CDF and writes
`build/mmm-transport/summary.json` and `summary.md`. It tests the **actual
writer methods** with an IMAS-shaped test double. `make verify-dd` requires
IMAS-Python. Actual DD4.1.1 write/readback must still be confirmed on a
workstation using `make verify-end-to-end ... E2E_IDS=core_transport`.
