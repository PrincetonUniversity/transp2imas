# Neoclassical theory transport in core_transport

Neoclassical source quantities are represented in the unified
`theory_transport_mapping.py` family registry, while actual IDS assignments
remain in `imas_core_transport.py`. Presence of neoclassical CDF channels
**does not prove PT_SOLVER selected them**; only reliable solver input
provenance can establish that association. Neoclassical fluxes remain
physically separate from predicted MMM07/GLF23 transport coefficients.

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


## Verification

```bash
make verify-neoclassical
make verify-all
```

No neutral profiles are audited by this target. IMAS write/readback requires
`make verify-end-to-end E2E_IDS=core_transport DD_VERSION=4.1.1` on an
IMAS-equipped machine.
