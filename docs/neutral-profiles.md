# Neutral density, temperature and rotation profiles

The neutral-state mapping is independent of transport-theory models.
`neutral_profiles_mapping.py` handles the CDF contracts and
`imas_core_profiles.py` owns IMAS storage. None of these quantities is
misidentified as an MMM07/GLF23/TGLF transport coefficient.

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


## Verification

```bash
make verify-neutrals
make verify-all
```

Neutral state is checked independently from neoclassical flux-divergence
reconstruction. Real schema validation requires IMAS-Python.
