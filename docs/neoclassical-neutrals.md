# Neoclassical transport and neutral-state CDF mappings

This extension is source-validated on **every bundled TRANSP CDF**. It does
not imply that all per-species physics is available on every machine, and
the DD4.1.1 write/readback must still be run with installed IMAS-Python.

## Neoclassical transport

| TRANSP CDF, dimensions | Physical quantity | IMAS field |
|---|---|---|
| `GFLNC_E(TIME3,X)` | electron neoclassical particle flux divergence | `core_transport.model[neoclassical].profiles_1d.electrons.particles.flux` |
| `GFLNC_<ion>` | identified ion particle flux divergence | `.ion[<ion>].particles.flux` |
| `QFLNC_E` | electron neoclassical heat flux divergence | `.electrons.energy.flux` |
| `QFLNC_<ion>` | identified ion neoclassical heat flux divergence | `.ion[<ion>].energy.flux` |
| `QFLNC_I` | summed thermal-ion neoclassical heat flux divergence | `.total_ion_energy.flux` |
| `QFLNCC_E`, `_I`, or `<ion>` | independent TRANSP NC-classical heat diagnostic | analogous heat flux in **separate** `neoclassical_classical_heat` model |

TRANSP metadata defines the source quantities as **divergences**, in
`N/CM3/SEC` and `WATTS/CM3`. With `DVOL` in cm3 per zone and `SURF` in cm2
at its outer face, the SI flux at radial boundary `j` is

```text
flux[j] = 10000 * sum(divergence[k] * DVOL[k], k=0..j) / SURF[j]
```

Particle flux units are particles/m2/s, energy flux W/m2. Missing/zero-area
boundaries yield zero rather than infinities. The source's `TIME3, X`
dimensions and physical units are checked before integration; output is on
the `XB` boundary grid. Flux values are independent of the aggregate
`combined` (observed effective coefficients) or `transport_solver` model
contributions, and should not be added together without a physics rationale.

`GFLNC_I` is an aggregate particle divergence without a faithful
`total_ion_particles.flux` destination, so it remains out of the typed IDS.
`GFLNC_X`, `QFLNC_X`, `*_TOK` are aggregates and are not invented as
individual ion species. Unknown ion tokens are skipped with a warning.
`QFLNCC_*` is named using TRANSP's own **NC-classical** designation; no
claim of additivity to `QFLNC_*` is implied. No neoclassical effective
diffusivity or convection is fabricated from a divergence.

## Neutral state profiles

| TRANSP CDF (TIME3,X) | Population | IMAS destination | Units |
|---|---|---|---|
| `DN0W<species>`, `T0W<species>` | wall-related | `core_profiles.profiles_1d.neutral[]` (distinct named population) | density cm^-3 to m^-3; T already eV |
| `DN0V<species>`, `T0V<species>` | volume-related | separate `neutral[]` | same |
| `N0V0_<species>`, `T0V0_<species>` | recombination | separate `neutral[]` | same |
| `OM0W<species>`, `OM0V<species>`, `OM0V0_<species>` | corresponding population toroidal angular rotation | see below | rad/s |

The wall, volume and recombination populations are exported separately.
They are **not** added, because their overlap/normalization should first
be established from the relevant TRANSP neutral-source model. Names have
origin suffixes (`D_WALL`, `D_VOLUME`, `D_RECOMBINATION`). Nuclear masses,
charges and number of atoms refer to the associated neutral chemical species.

**Neutral temperature bug fixed:** earlier `get_neutral_data_type` used
`get_ncvar(T0*, 1e6)` even though NetCDF units are `EV`. The new code
reads it directly as eV; the ×1e6 scaling applies only to cm^-3 densities.

**Neutral rotation on IMAS DD4.1.1:** `core_profiles.profiles_1d.neutral[]`
contains density and temperature, but no neutral toroidal velocity or
angular-frequency field. We preserve each **signed** `OM0*` profile in a
clearly labeled supplementary XML `<neutral_toroidal_rotation>` group in
`core_profiles.code.parameters`. The group documents the source channel,
axis names (`TIME3`, `X`), units (rad/s), toroidal sign, shape and numeric
values. This is an extension for provenance/reconstruction, **not** a typed
IMAS neutral-rotation field. In early DD3, a typed toroidal neutral velocity
may be available; when so, the code uses outboard `RMJMP` to compute
`v_tor = omega * R` after XB -> X interpolation. XML is preserved in
both versions. The converter reads the existing XML through
`IDSString0D.value` before appending the neutral data and assigns the
resulting Python XML string back to `core_profiles.code.parameters`. Input
provenance is preserved, and the no-rotation path also yields plain text.

## Tests and remaining validation

```bash
make verify-neoclassical-neutrals
make check
make verify-all
make verify-dd
make verify-end-to-end DD_VERSION=4.1.1 BACKEND=HDF5
```

`make check` executes synthetic/unit tests and real-writer IMAS-shaped
mock tests. The `verify-neoclassical-neutrals` target checks all available
CDFs; `verify-all` adds the same per-run source check to the standard CDF
inventory. The final two commands require a workstation with IMAS-Python.
