# TRANSP to IMAS mapping registry

Generated from `mapping_registry.yaml` by `python examples/verify_mapping_registry.py --write-docs`;
do not edit by hand.  Recognising a variable name is not an export: only typed entries are exports.

## direct_typed (2)

| id | writer | source | destination (DD4 / DD3) | units | verification | tests |
|---|---|---|---|---|---|---|
| tf_vacuum_field | tf | `BZXR` | `b_field_phi_vacuum_r` / `b_field_tor_vacuum_r` | T*cm -> T*m | unverified | - |
| pt_runtime_diagnostics_typed | transport_solver_numerics | `NPTRITER`, `PT_CPTIME`, `PT_WTIME` | `solver_1d/control_parameters` | count, hour | readback_verified | `verify_transport_solver_numerics.py` |

## derived_typed (8)

| id | writer | source | destination (DD4 / DD3) | units | verification | tests |
|---|---|---|---|---|---|---|
| vacuum_field_reference | core_profiles | `RAXIS`, `BZXR` | `vacuum_toroidal_field/b0` | cm, T*cm -> m, T | readback_verified | `verify_vacuum_field.py` |
| core_profiles_rho_tor | core_profiles | `TRFLX`, `BZXR` | `profiles_1d/grid/rho_tor` | Wb -> m | readback_verified | `verify_parallel_current.py` |
| equilibrium_rho_tor | equilibrium | `TRFLX`, `BZXR` | `time_slice/profiles_1d/rho_tor` | Wb -> m | readback_verified | `verify_parallel_current.py` |
| j_parallel_core_profiles | core_profiles | `PLJB`, `PLJBOH`, `CURBS`, `CONPLJB` | `profiles_1d/j_total` | A*T/cm^2 -> A/m^2 | readback_verified | `verify_parallel_current.py` |
| j_parallel_equilibrium | equilibrium | `PLJB` | `time_slice/profiles_1d/j_parallel` | A*T/cm^2 -> A/m^2 | readback_verified | `verify_parallel_current.py` |
| nbi_momentum_source | core_sources | `TQBE`, `TQBI`, `TQBTH`, `TQJXB`, `TQABDF` | `source/profiles_1d/momentum_phi` / `source/profiles_1d/momentum_tor` | N*m/cm^3 -> N*m/m^3 | readback_verified | `verify_nbi_momentum.py` |
| momentum_flux_solver_model | core_transport | `AMTR_MOD`, `SURF`, `DVOL` | `model/profiles_1d/momentum_phi/flux` / `model/profiles_1d/momentum_tor/flux` | N*m -> flux through the surface | readback_verified | `verify_momentum_transport.py` |
| momentum_flux_observed | core_transport | `AMTR_OBS` | `model/profiles_1d/momentum_phi/flux` / `model/profiles_1d/momentum_tor/flux` | N*m -> flux through the surface | unverified | - |

## xml_only (2)

| id | writer | source | destination (DD4 / DD3) | units | verification | tests |
|---|---|---|---|---|---|---|
| pt_runtime_diagnostics_xml | core_transport | `NPTRITER`, `PT_CPTIME`, `PT_WTIME` | `model/code/parameters` | - | readback_verified | `verify_runtime_diagnostics_readback.py` |
| pt_regional_models | core_transport | TR.DAT, PT_SOLVER template | `model/code/parameters` | - | readback_verified | `verify_regional_model_timeline.py` |

## input_consumed (2)

| id | source | reason |
|---|---|---|
| input_time_grids | `TIME`, `TIME3` | native time bases for every writer |
| input_geometry | `X`, `XB`, `DVOL`, `DAREA`, `SURF`, `BZ` | coordinates and zone geometry used to normalise and integrate |

## deferred (6)

| id | source | reason |
|---|---|---|
| deferred_boundary_sources | `T0RECYC`, `T0GASFL`, `OM0RECYC`, `OM0GASFL`, `E0INR` | neutral/source-state quantities with no audited IMAS destination |
| deferred_boundary_energy | `PEEDG`, `PIEDG` | energy crossing the plasma boundary, not a volumetric source |
| deferred_neoclassical_velocity | `VPOL?_NC`, `VTOR?_NC` | RMAJM-grid (101 midplane points) velocities; core_profiles reads them and exports only if a CDF provides them on the profile grid |
| deferred_aggregate_fluxes | `GFLNC_I`, `GFLNC_X`, `GFLNC_TOK`, `QFLNC_X`, `QFLNC_TOK`, `QFLNCC_X`, `QFLNCC_TOK` | aggregate divergences; no typed total destination and no extra ion species is invented |
| deferred_ic_spectrum | `PIC_NPHI`, `PICF??N??` | frequency/n_phi-resolved IC families need a dedicated spectral mapping audit |
| deferred_ripple_torque | `TQRPL`, `FTQRIP` | ripple torque weight FTQRIP is not in the CDF; TQRPL is read only to warn that it is excluded from the NBI momentum sum |

