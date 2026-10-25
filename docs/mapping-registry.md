# TRANSP to IMAS mapping registry

Generated from `mapping_registry.yaml` by `python examples/verify_mapping_registry.py --write-docs`;
do not edit by hand.  Recognising a variable name is not an export: only typed entries are exports.

## direct_typed (77)

| id | writer | source | destination (DD4 / DD3) | units | verification | tests |
|---|---|---|---|---|---|---|
| tf_vacuum_field | tf | `BZXR` | `b_field_phi_vacuum_r` / `b_field_tor_vacuum_r` | T*cm -> T*m | readback_verified | `verify_registry_values.py` |
| pt_runtime_diagnostics_typed | transport_solver_numerics | `NPTRITER`, `PT_CPTIME`, `PT_WTIME` | `solver_1d/control_parameters` | count, hour | readback_verified | `verify_transport_solver_numerics.py` |
| wall_limiter_outline | wall | `RLIM`, `YLIM` | `description_2d/limiter/unit/outline/r` | cm -> m | readback_verified | `verify_registry_values.py` |
| ic_antenna_frequency | ic_antennas | `FREQA?` | `antenna/frequency/data` | Hz -> Hz | readback_verified | `verify_registry_values.py` |
| ic_antenna_power | ic_antennas | `PICHA?` | `antenna/power_launched/data` | W -> W | readback_verified | `verify_registry_values.py` |
| ic_total_power | ic_antennas | `PICHTOT` | `power_launched/data` | W -> W | readback_verified | `verify_registry_values.py` |
| ec_beam_power | ec_launchers | `PECIN?` | `beam/power_launched/data` | W -> W | readback_verified | `verify_registry_values.py` |
| nbi_unit_power | nbi | `PINJ??` | `unit/power_launched/data` | W -> W | readback_verified | `verify_registry_values.py` |
| nbi_unit_energy | nbi | `EINJ??_E1` | `unit/energy/data` | eV -> eV | readback_verified | `verify_registry_values.py` |
| nbi_native_time | nbi | `TIME` | `unit/power_launched/time` | s -> s | readback_verified | `verify_registry_values.py`, `verify_nbi_mapping.py` |
| radiation_electron_emissivity | radiation | `PRAD_BR`, `PRAD_CY` | `process/profiles_1d/electrons/emissivity` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| radiation_line_emissivity | radiation | `PRAD_LI` | `process/profiles_1d/emissivity_ion_total` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| core_sources_nbi_electron_energy | core_sources | `PBE` | `source/profiles_1d/electrons/energy` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| core_sources_nbi_electron_particles | core_sources | `SBE` | `source/profiles_1d/electrons/particles` | 1/(cm^3 s) -> 1/(m^3 s) | readback_verified | `verify_registry_values.py` |
| core_sources_ohmic_energy | core_sources | `POH` | `source/profiles_1d/electrons/energy` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| cp_electron_temperature | core_profiles | `TE` | `profiles_1d/electrons/temperature` | eV -> eV | readback_verified | `verify_registry_values.py` |
| cp_electron_density | core_profiles | `NE` | `profiles_1d/electrons/density` | cm^-3 -> m^-3 | readback_verified | `verify_registry_values.py` |
| cp_ion_average_temperature | core_profiles | `TI` | `profiles_1d/t_i_average` | eV -> eV | readback_verified | `verify_registry_values.py` |
| cp_ion_temperature | core_profiles | `TI` | `profiles_1d/ion/temperature` | eV -> eV | readback_verified | `verify_registry_values.py` |
| cp_fast_ion_temperature | core_profiles | `EBEAM_D`, `EBEAM_T`, `EBEAM_H`, `EFUSN_T`, `EFUSN_3`, `EFUSN_4`, `EFUSN_P`, `TMINI_D`, `TMINI_T`, `TMINI_H` | `profiles_1d/ion/temperature` | eV -> eV | readback_verified | `verify_registry_values.py` |
| cp_ion_density_thermal | core_profiles | `ND`, `NT`, `NHE?` | `profiles_1d/ion/density_thermal` | cm^-3 -> m^-3 | readback_verified | `verify_registry_values.py` |
| cp_ion_radial_velocity | core_profiles | `VNDNC_D`, `VNDNC_T`, `VNDNC_H`, `VNDNC_HE3`, `VNDNC_HE4` | `profiles_1d/ion/velocity/radial` | cm/s -> m/s | readback_verified | `verify_registry_values.py` |
| cp_ion_rotation_frequency | core_profiles | `OMEGA` | `profiles_1d/ion/rotation_frequency_tor` | rad/s -> rad/s | readback_verified | `verify_registry_values.py` |
| cp_zeff | core_profiles | `ZEFFP` | `profiles_1d/zeff` | 1 | readback_verified | `verify_registry_values.py` |
| cp_magnetic_shear | core_profiles | `SHAT` | `profiles_1d/magnetic_shear` | 1 | readback_verified | `verify_registry_values.py` |
| ct_combined_electron_particle_dv | core_transport | `DIFFE`, `VELE` | `model/profiles_1d/electrons/particles/d` | cm^2/s, cm/s -> m^2/s, m/s | readback_verified | `verify_registry_values.py` |
| ct_combined_electron_heat_dv | core_transport | `CONDE`, `VEL_TE` | `model/profiles_1d/electrons/energy/d` | cm^2/s, cm/s -> m^2/s, m/s | readback_verified | `verify_registry_values.py` |
| ct_combined_ion_heat_dv | core_transport | `CONDI`, `VEL_TI` | `model/profiles_1d/total_ion_energy/d` | cm^2/s, cm/s -> m^2/s, m/s | readback_verified | `verify_registry_values.py` |
| ct_combined_momentum_dv | core_transport | `CHPHI`, `VMO_NET` | `model/profiles_1d/momentum_phi/d` / `model/profiles_1d/momentum_tor/d` | cm^2/s, cm/s -> m^2/s, m/s | readback_verified | `verify_registry_values.py` |
| ct_combined_species_dv | core_transport | `DIFFH`, `DIFFD`, `DIFFT`, `DIFFHE3`, `DIFFHE4`, `VELH`, `VELD`, `VELT`, `VELHE3`, `VELHE4` | `model/profiles_1d/ion/particles/d` | cm^2/s, cm/s -> m^2/s, m/s | readback_verified | `verify_registry_values.py` |
| ct_glf23_diffusivities | core_transport | `XKEGLF23`, `XKIGLF23`, `ETPHIGLF` | `model/profiles_1d/electrons/energy/d` | cm^2/s -> m^2/s | readback_verified | `verify_registry_values.py` |
| ct_mmm07_diffusivities | core_transport | `XKDMMM07`, `XKEMMM07`, `XKIMMM07`, `XPTMMM07`, `VPHIMMM07`, `XKEW19`, `XKIW19`, `XKEDRBM`, `XKIDRBM`, `XKEETG`, `XKEMTM` | `model/profiles_1d/electrons/energy/d` | cm^2/s, cm/s -> m^2/s, m/s | readback_verified | `verify_registry_values.py` |
| es_ionization_particles | edge_sources | `SEGF_?`, `SERC_?`, `SERECO` | `source/ggd/electrons/particles/values` | 1/(cm^3 s) -> 1/(m^3 s) | readback_verified | `verify_registry_values.py` |
| es_ionization_energy | edge_sources | `PIGF_?`, `PIRC_?`, `PIRECO` | `source/ggd/ion/energy/values` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| es_ionization_momentum | edge_sources | `TQIGF_?`, `TQIRC_?`, `TQIRECO` | `source/ggd/ion/momentum/phi` / `-` | N*m/cm^3 -> N*m/m^3 | readback_verified | `verify_registry_values.py` |
| es_cx_particles | edge_sources | `SFCXGF_?`, `SFCXRC_?`, `SFCXRECO` | `source/ggd/ion/particles/values` | 1/(cm^3 s) -> 1/(m^3 s) | readback_verified | `verify_registry_values.py` |
| es_cx_energy | edge_sources | `PCXGF_?`, `PCXRC_?`, `PCXRECO` | `source/ggd/total_ion_energy/values` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| es_cx_momentum | edge_sources | `TQCXGF_?`, `TQCXRC_?`, `TQCXRECO` | `source/ggd/momentum/phi` / `-` | N*m/cm^3 -> N*m/m^3 | readback_verified | `verify_registry_values.py` |
| es_recombination_energy | edge_sources | `P0RECO` | `source/ggd/neutral/energy/values` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| es_recombination_momentum | edge_sources | `TQ0RECO` | `source/ggd/neutral/momentum/phi` / `-` | N m/cm^3 -> N m/m^3 | readback_verified | `verify_registry_values.py` |
| eq_plasma_current | equilibrium | `PCURC` | `time_slice/global_quantities/ip` | A -> A | readback_verified | `verify_registry_values.py` |
| eq_magnetic_axis | equilibrium | `RAXIS`, `YAXIS` | `time_slice/global_quantities/magnetic_axis/r` | cm -> m | readback_verified | `verify_registry_values.py` |
| eq_psi_magnetic_axis | equilibrium | `PSI0_TR` | `time_slice/global_quantities/psi_magnetic_axis` / `time_slice/global_quantities/psi_axis` | Wb/rad -> Wb | readback_verified | `verify_registry_values.py` |
| eq_q_scalars | equilibrium | `Q95`, `Q0` | `time_slice/global_quantities/q_95` | 1 | readback_verified | `verify_registry_values.py` |
| eq_li3 | equilibrium | `LI_3` | `time_slice/global_quantities/li_3` | 1 | readback_verified | `verify_registry_values.py` |
| eq_plasma_volume | equilibrium | `PVOL` | `time_slice/global_quantities/volume` | cm^3 -> m^3 | readback_verified | `verify_registry_values.py` |
| eq_boundary_shape | equilibrium | `ELONG`, `TRIANG`, `TRIANGU`, `TRIANGL` | `time_slice/boundary/elongation` | 1 | readback_verified | `verify_registry_values.py` |
| eq_profile_shaping | equilibrium | `ELONG`, `TRIANGU`, `TRIANGL`, `SQUARE_UO`, `SQUARE_LO` | `time_slice/profiles_1d/elongation` | 1 | readback_verified | `verify_registry_values.py` |
| eq_profile_phi | equilibrium | `TRFLX` | `time_slice/profiles_1d/phi` | Wb -> Wb | readback_verified | `verify_registry_values.py` |
| eq_profile_q | equilibrium | `Q` | `time_slice/profiles_1d/q` | 1 | readback_verified | `verify_registry_values.py` |
| sm_plasma_current | summary | `PCURC` | `global_quantities/ip/value` | A -> A | readback_verified | `verify_registry_values.py` |
| sm_direct_scalars | summary | `LI_3`, `VSURC`, `PVOL`, `Q95`, `TAUA1`, `H98Y2`, `TAUE98Y2`, `LHMODE` | `global_quantities/li_3/value` / `global_quantities/li/value` | V, 1, m^3, s | readback_verified | `verify_registry_values.py` |
| wv_ec_frequency | waves | TR.DAT (FREQECH) | `coherent_wave/global_quantities/frequency` | Hz -> Hz | unverified | - |
| cp_fast_ion_density | core_profiles | `BDENS_D`, `BDENS_T`, `BDENS_H`, `FDENS_T`, `FDENS_3`, `FDENS_4`, `FDENS_P`, `NMINI_D`, `NMINI_T`, `NMINI_H` | `profiles_1d/ion/density_fast` | cm^-3 -> m^-3 | readback_verified | `verify_registry_values.py` |
| cp_fast_ion_rotation | core_profiles | `OMEGB_D`, `OMEGB_T`, `OMEGB_H`, `OMEGF_T`, `OMEGF_3`, `OMEGF_4`, `OMEGF_P` | `profiles_1d/ion/rotation_frequency_tor` | rad/s -> rad/s | readback_verified | `verify_registry_values.py` |
| cp_impurity_density | core_profiles | `NIMPS_BE`, `NIMPS_NE`, `NIMPS_W`, `NIMPS_AR`, `NIMPS_B`, `NIMPS_GA`, `NIMP` | `profiles_1d/ion/density` | cm^-3 -> m^-3 | readback_verified | `verify_registry_values.py` |
| cp_impurity_temperature_rotation | core_profiles | `TX`, `OMGX` | `profiles_1d/ion/temperature` | eV, rad/s | readback_verified | `verify_registry_values.py` |
| cp_impurity_charge_states | core_profiles | `NIMP_BE_?`, `NIMP_NE_?`, `NIMP_NE_??`, `NIMP_W_??`, `NIMP_AR_??`, `NIMP_B_?`, `NIMP_GA_??` | `profiles_1d/ion/state/density` | cm^-3 -> m^-3 | readback_verified | `verify_registry_values.py` |
| cp_neutral_profiles | core_profiles | `DN0V?`, `DN0W?`, `DN0V??`, `DN0W??`, `DN0V???`, `DN0W???`, `T0V?`, `T0W?`, `T0V??`, `T0W??`, `T0V???`, `T0W???` | `profiles_1d/neutral/density` | cm^-3 -> m^-3; eV | readback_verified | `verify_registry_values.py` |
| cp_momentum_phi | core_profiles | `PPHI` | `profiles_1d/momentum_phi` / `profiles_1d/momentum_tor` | N*m*s/cm^3 -> kg/(m*s) | readback_verified | `verify_registry_values.py` |
| cp_global_scalars | core_profiles | `PCURC`, `VSURC`, `LI_3` | `global_quantities/ip` | A, V, 1 | readback_verified | `verify_registry_values.py` |
| cp_profile_j_tor | core_profiles | `CUR` | `profiles_1d/j_phi` / `profiles_1d/j_tor` | A/cm^2 -> A/m^2 | readback_verified | `verify_registry_values.py` |
| pf_active_coil_current | pf_active | `KK_*` | `coil/current/data` | A (CDF units attribute; kA/MA scaled) -> A | writer_verified | `verify_pf_active_mapping.py` |
| pf_active_coil_voltage | pf_active | `VV_*` | `coil/voltage/data` | V (CDF units attribute; mV/kV/MV scaled) -> V | writer_verified | `verify_pf_active_mapping.py` |
| pf_active_supply_voltage | pf_active | `VS_*` | `supply/voltage/data` | V (CDF units attribute; mV/kV/MV scaled) -> V | writer_verified | `verify_pf_active_mapping.py` |
| pf_passive_loop_current | pf_passive | `PS_[1-9]*` | `loop/current` | A (CDF units attribute) -> A | writer_verified | `verify_pf_passive_mapping.py` |
| core_sources_ec_electron_energy | core_sources | `PEECH` | `source/profiles_1d/electrons/energy` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| core_sources_ic_electron_energy | core_sources | `PEICH` | `source/profiles_1d/electrons/energy` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| core_sources_ic_ion_energy | core_sources | `PIICH` | `source/profiles_1d/total_ion_energy` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| core_sources_nbi_ion_particles | core_sources | `SBTH` | `source/profiles_1d/total_ion_particles` / `-` | 1/cm^3/s -> 1/m^3/s | readback_verified | `verify_registry_values.py` |
| cp_neutral_recombination | core_profiles | `N0V0_*`, `T0V0_*` | `profiles_1d/neutral/density` | cm^-3 -> m^-3; eV | readback_verified | `verify_registry_values.py` |
| cp_single_impurity_mass | core_profiles | `AIMP` | `profiles_1d/ion/element/a` | amu -> amu | readback_verified | `verify_registry_values.py` |
| cp_single_impurity_charge | core_profiles | `XZIMP` | `profiles_1d/ion/element/z_n` | e -> e | readback_verified | `verify_registry_values.py` |
| eq_profile_geometric_axis | equilibrium | `RMC00`, `YMC00`, `RMJMP` | `time_slice/profiles_1d/geometric_axis/r` | cm -> m | readback_verified | `verify_registry_values.py` |
| eq_profile_2d_b_field | equilibrium | `B_FIELD` | `time_slice/profiles_2d/b_field_r` | T -> T | readback_verified | `verify_registry_values.py` |
| sawteeth_previous_crash_time | sawteeth | `TLASTSAW` | `diagnostics/previous_crash_time` | s -> s | writer_verified | `verify_sawteeth_mapping.py` |
| sawteeth_q1_radius | sawteeth | `XIQ1` | `diagnostics/rho_tor_norm_q1` | 1 -> 1 | writer_verified | `verify_sawteeth_mapping.py` |

## derived_typed (89)

| id | writer | source | destination (DD4 / DD3) | units | verification | tests |
|---|---|---|---|---|---|---|
| vacuum_field_reference | core_profiles | `RLIM`, `RAXIS`, `BZXR` | `vacuum_toroidal_field/b0` | cm, T*cm -> m, T | readback_verified | `verify_vacuum_field.py` |
| core_profiles_rho_tor | core_profiles | `TRFLX`, `BZXR` | `profiles_1d/grid/rho_tor` | Wb -> m | readback_verified | `verify_parallel_current.py` |
| equilibrium_rho_tor | equilibrium | `TRFLX`, `BZXR` | `time_slice/profiles_1d/rho_tor` | Wb -> m | readback_verified | `verify_parallel_current.py` |
| j_parallel_core_profiles | core_profiles | `PLJB`, `PLJBOH`, `CURBS`, `CONPLJB` | `profiles_1d/j_total` | A*T/cm^2 -> A/m^2 | readback_verified | `verify_parallel_current.py` |
| j_parallel_equilibrium | equilibrium | `PLJB` | `time_slice/profiles_1d/j_parallel` | A*T/cm^2 -> A/m^2 | readback_verified | `verify_parallel_current.py` |
| nbi_momentum_source | core_sources | `TQBE`, `TQBI`, `TQBTH`, `TQJXB`, `TQABDF` | `source/profiles_1d/momentum_phi` / `source/profiles_1d/momentum_tor` | N*m/cm^3 -> N*m/m^3 | readback_verified | `verify_nbi_momentum.py` |
| nbi_momentum_ripple | core_sources | `TQRPL`, `FTQRIP` | `source/profiles_1d/momentum_phi` / `source/profiles_1d/momentum_tor` | N*m/cm^3 -> N*m/m^3 | writer_verified | `verify_nbi_momentum.py` |
| nbi_unit_species | nbi | `PINJ??`, `PINJ_*` | `unit/species/a` | 1 -> 1 | writer_verified | `verify_nbi_mapping.py` |
| momentum_flux_solver_model | core_transport | `AMTR_MOD`, `SURF`, `DVOL` | `model/profiles_1d/momentum_phi/flux` / `model/profiles_1d/momentum_tor/flux` | N m cm^-3 (a torque density) -> N m m^-2, i.e. kg s^-2 | readback_verified | `verify_momentum_transport.py` |
| momentum_flux_observed | core_transport | `AMTR_OBS` | `model/profiles_1d/momentum_phi/flux` / `model/profiles_1d/momentum_tor/flux` | N m cm^-3 (a torque density) -> N m m^-2, i.e. kg s^-2 | readback_verified | `verify_registry_values.py` |
| ec_steering_angles | ec_launchers | `PHAIECH?`, `THETECH?` | `beam/steering_angle_tor` | deg -> rad | readback_verified | `verify_registry_values.py` |
| core_sources_nbi_ion_energy | core_sources | `PBI`, `PBTH` | `source/profiles_1d/total_ion_energy` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| core_sources_current_drive_ohmic_nbi | core_sources | `PLJBOH`, `CURB`, `CONPLJB` | `source/profiles_1d/j_parallel` | A*T/cm^2 or A/cm^2 -> A/m^2 | readback_verified | `verify_parallel_current.py` |
| core_sources_bootstrap_current | core_sources | `CURBS`, `CONPLJB` | `source/profiles_1d/j_parallel` | A/cm^2 -> A/m^2 | readback_verified | `verify_registry_values.py` |
| core_sources_radiation_sink | core_sources | `PRAD_BR`, `PRAD_CY`, `PRAD_LI` | `source/profiles_1d/electrons/energy` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| cp_ion_to_electron_density | core_profiles | `NI`, `NE` | `profiles_1d/n_i_total_over_n_e` | 1 | readback_verified | `verify_registry_values.py` |
| cp_pressure_thermal | core_profiles | `PPLAS`, `UTHRM` | `profiles_1d/pressure_thermal` | Pa -> Pa; J/cm^3 -> Pa | readback_verified | `verify_registry_values.py` |
| cp_pressure_electron | core_profiles | `UE` | `profiles_1d/electrons/pressure_thermal` | J/cm^3 -> Pa | readback_verified | `verify_registry_values.py` |
| cp_pressure_ion_total | core_profiles | `UI` | `profiles_1d/pressure_ion_total` | J/cm^3 -> Pa | readback_verified | `verify_registry_values.py` |
| cp_pressure_anisotropy | core_profiles | `PPLAS`, `UTHRM`, `UFASTPP`, `UFASTPA` | `profiles_1d/pressure_perpendicular` | J/cm^3 -> Pa | readback_verified | `verify_registry_values.py` |
| cp_q_profile | core_profiles | `Q`, `Q0`, `X`, `XB` | `profiles_1d/q` | 1 | readback_verified | `verify_registry_values.py` |
| cp_parallel_conductivity | core_profiles | `ETA_USE` | `profiles_1d/conductivity_parallel` | ohm*cm -> 1/(ohm*m) | readback_verified | `verify_registry_values.py` |
| cp_phi_potential | core_profiles | `VRPOT`, `X`, `XB` | `profiles_1d/phi_potential` | V -> V | readback_verified | `verify_registry_values.py` |
| ct_combined_conductivity | core_transport | `ETA_USE`, `X`, `XB` | `model/profiles_1d/conductivity_parallel` | ohm*cm -> 1/(ohm*m) | readback_verified | `verify_registry_values.py` |
| ct_observed_fluxes | core_transport | `EPTR_OBS`, `EETR_OBS`, `IETR_OBS`, `DVOL`, `SURF` | `model/profiles_1d/electrons/particles/flux` | per cm^3 -> per m^2 | readback_verified | `verify_registry_values.py` |
| ct_observed_species_flux | core_transport | `PTRH_OBS`, `PTRD_OBS`, `PTRT_OBS`, `PTRHE3_OBS`, `PTRHE4_OBS`, `DVOL`, `SURF` | `model/profiles_1d/ion/particles/flux` | per cm^3 -> per m^2 | readback_verified | `verify_registry_values.py` |
| ct_solver_fluxes | core_transport | `EPTR_MOD`, `EETR_MOD`, `IETR_MOD`, `DVOL`, `SURF` | `model/profiles_1d/electrons/particles/flux` | per cm^3 -> per m^2 | readback_verified | `verify_registry_values.py` |
| ct_solver_species_flux | core_transport | `PTRH_MOD`, `PTRD_MOD`, `PTRT_MOD`, `PTRHE3_MOD`, `PTRHE4_MOD`, `DVOL`, `SURF` | `model/profiles_1d/ion/particles/flux` | per cm^3 -> per m^2 | readback_verified | `verify_registry_values.py` |
| ct_neoclassical_fluxes | core_transport | `GFLNC_E`, `GFLNC_H`, `GFLNC_D`, `GFLNC_T`, `GFLNC_3`, `GFLNC_4`, `GFLNC_B`, `GFLNC_BE`, `GFLNC_NE`, `GFLNC_AR`, `GFLNC_GA`, `GFLNC_W`, `QFLNC_E`, `QFLNC_H`, `QFLNC_D`, `QFLNC_T`, `QFLNC_3`, `QFLNC_4`, `QFLNC_B`, `QFLNC_BE`, `QFLNC_NE`, `QFLNC_AR`, `QFLNC_GA`, `QFLNC_W`, `QFLNC_I`, `DVOL`, `SURF` | `model/profiles_1d/electrons/particles/flux` | per cm^3 -> per m^2 | readback_verified | `verify_registry_values.py` |
| ct_neoclassical_classical_heat | core_transport | `QFLNCC_E`, `QFLNCC_H`, `QFLNCC_D`, `QFLNCC_T`, `QFLNCC_3`, `QFLNCC_4`, `QFLNCC_B`, `QFLNCC_BE`, `QFLNCC_NE`, `QFLNCC_AR`, `QFLNCC_GA`, `QFLNCC_W`, `QFLNCC_I`, `DVOL`, `SURF` | `model/profiles_1d/electrons/energy/flux` | per cm^3 -> per m^2 | readback_verified | `verify_registry_values.py` |
| es_recombination_particles | edge_sources | `S0RECO`, `S0V0?`, `S0V0??` | `source/ggd/neutral/particles/values` | 1/(cm^3 s) -> 1/(m^3 s) | readback_verified | `verify_registry_values.py` |
| eq_betas | equilibrium | `PTOWB`, `DVOL`, `PCURC`, `BZXR`, `RLIM`, `RAXIS`, `RMCB?`, `RMCB??`, `RMSB?`, `RMSB??` | `time_slice/global_quantities/beta_tor` | Pa, m^3, A, T -> 1 | readback_verified | `verify_registry_values.py` |
| eq_profile_psi | equilibrium | `PLFLX2PI`, `PSI0_TR` | `time_slice/profiles_1d/psi` | Wb | readback_verified | `verify_registry_values.py` |
| eq_profile_f | equilibrium | `GFUN`, `BZXR` | `time_slice/profiles_1d/f` | T*m | readback_verified | `verify_registry_values.py` |
| eq_profile_pressure | equilibrium | `PTOWB`, `X`, `XB` | `time_slice/profiles_1d/pressure` | Pa | readback_verified | `verify_registry_values.py` |
| eq_profile_geometry | equilibrium | `DVOL`, `DAREA`, `SURF` | `time_slice/profiles_1d/volume` | cm^3, cm^2 -> m^3, m^2 | readback_verified | `verify_registry_values.py` |
| eq_profile_field_averages | equilibrium | `GB1`, `BMAX`, `BMIN`, `GB2`, `GR2I`, `GRI`, `GXI`, `DRAV`, `DRAVFAC`, `X`, `XB` | `time_slice/profiles_1d/b_field_average` | T, 1/m, 1/m^2 | readback_verified | `verify_registry_values.py` |
| eq_profile_j_tor | equilibrium | `CUR`, `X`, `XB` | `time_slice/profiles_1d/j_phi` / `time_slice/profiles_1d/j_tor` | A/cm^2 -> A/m^2 | readback_verified | `verify_registry_values.py` |
| eq_profile_shear_mass | equilibrium | `SHAT`, `PMASS`, `X`, `XB` | `time_slice/profiles_1d/magnetic_shear` | 1, g/cm^3 -> kg/m^3 | readback_verified | `verify_registry_values.py` |
| eq_boundary_psi | equilibrium | `PLFLX2PI`, `PSI0_TR` | `time_slice/boundary/psi` | Wb | readback_verified | `verify_registry_values.py` |
| eq_boundary_separatrix_psi | equilibrium | `PBPMHD`, `PSI0_TR` | `-` / `time_slice/boundary_separatrix/psi` | Wb/rad -> Wb | readback_verified | `verify_registry_values.py` |
| eq_profile_2d_psi | equilibrium | `PSIRZ`, `PSI0_TR`, `RGRID`, `ZGRID` | `time_slice/profiles_2d/psi` | Wb/rad -> Wb | readback_verified | `verify_registry_values.py` |
| sm_vacuum_field | summary | `RLIM`, `RAXIS`, `BZXR` | `global_quantities/b0/value` | cm, T*cm -> m, T | readback_verified | `verify_registry_values.py` |
| sm_betas | summary | `PTOWB`, `DVOL`, `PCURC`, `BZXR`, `RLIM`, `RAXIS`, `RMCB?`, `RMCB??`, `RMSB?`, `RMSB??` | `global_quantities/beta_tor_norm/value` | Pa, m^3, A, T -> 1 | readback_verified | `verify_registry_values.py` |
| sm_energies | summary | `UTOTL`, `UTHRM`, `UE`, `UI`, `PTOWB`, `PPLAS`, `UFASTPP`, `UFASTPA`, `DVOL` | `global_quantities/energy_total/value` | J/cm^3 -> J | readback_verified | `verify_registry_values.py` |
| sm_poloidal_field_energy | summary | `BPOL`, `X`, `XB`, `DVOL` | `global_quantities/energy_b_field_pol/value` | T -> J | readback_verified | `verify_registry_values.py` |
| sm_currents | summary | `CUROH`, `CUR`, `CURBS`, `DAREA` | `global_quantities/current_ohm/value` | A/cm^2 * cm^2 -> A | readback_verified | `verify_registry_values.py` |
| sm_powers | summary | `POHT`, `PRAD`, `PRAD_BR`, `PRAD_LI`, `PRAD_CY`, `PINJ`, `DVOL` | `global_quantities/power_ohm/value` | W, W/cm^3 -> W | readback_verified | `verify_registry_values.py` |
| sm_neutron_rates | summary | `NEUTT`, `NEUTX`, `NEUTX_DD`, `BTNTS_DD`, `BBNTS_DD` | `fusion/neutron_rates/total/value` | N/s | readback_verified | `verify_registry_values.py` |
| ds_birth_particles | distribution_sources | `BDEP??_TOT` | `source/profiles_1d/particles` | 1/(cm^3 s) -> 1/(m^3 s) | readback_verified | `verify_registry_values.py` |
| ds_birth_energy | distribution_sources | `BDEP??_E1`, `BDEP??_E2`, `BDEP??_E3`, `EINJ??_E1`, `EINJ??_E2`, `EINJ??_E3` | `source/profiles_1d/energy` | cm^-3 s^-1, eV -> W/m^3 | readback_verified | `verify_registry_values.py` |
| ds_birth_momentum | distribution_sources | `TQTOT??` | `source/profiles_1d/momentum_phi` / `source/profiles_1d/momentum_tor` | N*m/cm^3 -> N/m^2 | readback_verified | `verify_registry_values.py` |
| ds_global_quantities | distribution_sources | `BDEP??_TOT`, `BDEP??_E1`, `BDEP??_E2`, `BDEP??_E3`, `TQTOT??`, `DVOL` | `source/global_quantities/particles` | -> 1/s, W, N*m | readback_verified | `verify_registry_values.py` |
| ds_shinethrough | distribution_sources | `BPSH??_E1`, `BPSH??_E2`, `BPSH??_E3`, `EINJ??_E1`, `EINJ??_E2`, `EINJ??_E3` | `source/global_quantities/shinethrough/power` | W -> W, 1/s | readback_verified | `verify_registry_values.py` |
| wv_ec_power | waves | `PEECH?` | `coherent_wave/profiles_1d/power_density` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| wv_ec_current | waves | `ECCUR?`, `DAREA` | `coherent_wave/profiles_1d/current_phi_inside` / `coherent_wave/profiles_1d/current_tor_inside` | A/cm^2 * cm^2 -> A | readback_verified | `verify_registry_values.py` |
| wv_ic_power | waves | `RFPWRAB?`, `PIC_F?`, `DVOL` | `coherent_wave/global_quantities/power` | W, W/cm^3 -> W, W/m^3 | readback_verified | `verify_registry_values.py` |
| wv_ic_frequency | waves | `FREQA?` | `coherent_wave/global_quantities/frequency` | Hz -> Hz | readback_verified | `verify_registry_values.py` |
| cp_ion_density_total | core_profiles | `ND`, `NT`, `NH`, `NHE3`, `NHE4`, `BDENS_D`, `BDENS_T`, `BDENS_H`, `FDENS_T`, `FDENS_3`, `FDENS_4`, `FDENS_P`, `NMINI_D`, `NMINI_T`, `NMINI_H` | `profiles_1d/ion/density` | cm^-3 -> m^-3 | readback_verified | `verify_registry_values.py` |
| cp_betas | core_profiles | `PTOWB`, `DVOL`, `PCURC`, `BZXR`, `RLIM`, `RAXIS`, `RMCB?`, `RMCB??`, `RMSB?`, `RMSB??` | `global_quantities/beta_tor` | Pa, m^3, A, T -> 1 | readback_verified | `verify_registry_values.py` |
| cp_global_currents | core_profiles | `CUR`, `CUROH`, `CURBS`, `DAREA` | `global_quantities/current_non_inductive` | A/cm^2 * cm^2 -> A | readback_verified | `verify_registry_values.py` |
| cp_energy_diamagnetic | core_profiles | `PTOWB`, `DVOL` | `global_quantities/energy_diamagnetic` | Pa*cm^3 -> J | readback_verified | `verify_registry_values.py` |
| cp_volume_averages | core_profiles | `TE`, `TI`, `NE`, `NI`, `DVOL` | `global_quantities/t_e_volume_average` | eV, cm^-3 -> m^-3 | readback_verified | `verify_registry_values.py` |
| cp_peaking | core_profiles | `TE0`, `TI0`, `TE`, `TI`, `DVOL` | `global_quantities/t_e_peaking` | 1 | readback_verified | `verify_registry_values.py` |
| cp_j_bootstrap_and_non_inductive | core_profiles | `CURBS`, `CONPLJB`, `PLJB`, `PLJBOH` | `profiles_1d/j_bootstrap` | A/cm^2, A*T/cm^2 -> A/m^2 | readback_verified | `verify_registry_values.py` |
| radiation_line_ion_emissivity | radiation | `PRLS_*` | `process/profiles_1d/ion/emissivity` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| radiation_line_state_emissivity | radiation | `PRLX_*` | `process/profiles_1d/ion/state/emissivity` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| core_sources_fusion_electron_energy | core_sources | `PFE_*` | `source/profiles_1d/electrons/energy` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| core_sources_fusion_ion_energy | core_sources | `PFI_*`, `PFTH_*` | `source/profiles_1d/total_ion_energy` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| core_sources_fusion_ion_particles | core_sources | `SFTH_*` | `source/profiles_1d/total_ion_particles` / `-` | 1/cm^3/s -> 1/m^3/s | readback_verified | `verify_registry_values.py` |
| core_sources_fusion_ion_species | core_sources | `SFTH_*`, `PFTH_*` | `source/profiles_1d/ion/particles` | 1/cm^3/s -> 1/m^3/s ; W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| core_sources_nbi_ion_species | core_sources | `SBTH_?`, `PBTH_?` | `source/profiles_1d/ion/particles` | 1/cm^3/s -> 1/m^3/s ; W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| core_sources_equipartition | core_sources | `QIE` | `source/profiles_1d/electrons/energy` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| core_sources_cold_neutral_particles | core_sources | `SCEW`, `SCEV` | `source/profiles_1d/electrons/particles` | 1/cm^3/s -> 1/m^3/s | readback_verified | `verify_registry_values.py` |
| core_sources_cold_neutral_energy | core_sources | `PION` | `source/profiles_1d/electrons/energy` | W/cm^3 -> W/m^3 | readback_verified | `verify_registry_values.py` |
| dist_orbit_cx_losses | distributions | `BPLIM_*`, `BPCXI_*`, `BPCXX_*` | `distribution/global_quantities/power_first_orbit` / `-` | W -> W | writer_verified | `verify_nubeam_fi.py` |
| dist_collisions | distributions | `BPTE_*`, `BPTI_*`, `PBE_*`, `PBI_*` | `distribution/global_quantities/collisions/electrons/power_thermal` | W -> W; W/cm^3 -> W/m^3 | writer_verified | `verify_nubeam_fi.py` |
| dist_thermalization | distributions | `BPTH_*`, `BSTH_*`, `PBTH_*` | `distribution/global_quantities/thermalization/power` / `distribution/global_quantities/thermalisation/power` | W -> W; 1/s -> 1/s; W/cm^3 -> W/m^3 | writer_verified | `verify_nubeam_fi.py` |
| dist_nbi_source | distributions | `SBDEPMC_*`, `PBDEPMC_*`, `BPHDEP0_*`, `BDEP_*` | `distribution/global_quantities/source/particles` | 1/s, W, N m; 1/cm^3/s -> 1/m^3/s | writer_verified | `verify_nubeam_fi.py` |
| dist_trapped_density | distributions | `BTRAP_*`, `BDENS_*` | `distribution/profiles_1d/trapped/density_fast` / `-` | - -> m^-3 | writer_verified | `verify_nubeam_fi.py` |
| dist_torque_current | distributions | `TQJXB_*`, `UCURB`, `BPHXB*` | `distribution/global_quantities/torque_tor_j_radial` | N m cm^-3 -> N m (volume integral), N m^-2 (profile) | writer_verified | `verify_nubeam_fi.py` |
| eq_energy_mhd | equilibrium | `PMHD`, `DVOL` | `time_slice/global_quantities/energy_mhd` | Pa, cm^3 -> J | readback_verified | `verify_registry_values.py` |
| eq_boundary_geometric_axis | equilibrium | `RMCB0`, `YMCB0` | `time_slice/boundary/geometric_axis/r` | cm -> m | readback_verified | `verify_registry_values.py` |
| eq_boundary_outline | equilibrium | `PLFLXA`, `PBPMHD`, `PSIRZ`, `RGRID`, `ZGRID`, `RAXIS`, `YAXIS` | `time_slice/boundary/outline/r` | Wb/rad, cm -> m | writer_verified | `verify_equilibrium_boundary_trace.py` |
| grid_poloidal_flux_core_sources | core_sources | `PLFLX`, `PSI0_TR`, `X`, `XB` | `source/profiles_1d/grid/rho_pol_norm` | Wb/rad -> Wb, 1 | readback_verified | `verify_registry_values.py` |
| grid_poloidal_flux_radiation | radiation | `PLFLX`, `PSI0_TR`, `X`, `XB` | `process/profiles_1d/grid/rho_pol_norm` | Wb/rad -> Wb, 1 | readback_verified | `verify_registry_values.py` |
| grid_poloidal_flux_waves | waves | `PLFLX`, `PSI0_TR`, `X`, `XB` | `coherent_wave/profiles_1d/grid/rho_pol_norm` | Wb/rad -> Wb, 1 | readback_verified | `verify_registry_values.py` |
| grid_poloidal_flux_distribution_sources | distribution_sources | `PLFLX`, `PSI0_TR`, `X`, `XB` | `source/profiles_1d/grid/rho_pol_norm` | Wb/rad -> Wb, 1 | readback_verified | `verify_registry_values.py` |
| sawteeth_previous_period | sawteeth | `TLASTSAW` | `diagnostics/previous_period` | s -> s | writer_verified | `verify_sawteeth_mapping.py` |

## xml_only (4)

| id | writer | source | destination (DD4 / DD3) | units | verification | tests |
|---|---|---|---|---|---|---|
| pt_runtime_diagnostics_xml | core_transport | `NPTRITER`, `PT_CPTIME`, `PT_WTIME` | `model/code/parameters` | - | readback_verified | `verify_runtime_diagnostics_readback.py` |
| pt_regional_models | core_transport | TR.DAT, PT_SOLVER template | `model/code/parameters` | - | readback_verified | `verify_regional_model_timeline.py` |
| ct_theory_auxiliary_xml | core_transport | `KYRSDRBM`, `GAMMAETG`, `OMEGAETG`, `KYRSETG`, `GAMMAMTM`, `OMEGAMTM`, `KYRSMTM`, `XDBMTM`, `XKZMMM07`, `XPPMMM07`, `DIFFIGLF` | `model/code/parameters` | - | readback_verified | `verify_registry_values.py` |
| cp_neutral_rotation_xml | core_profiles | `OM0V*`, `OM0W*` | `code/parameters` | rad/s -> rad/s | readback_verified | `verify_registry_values.py` |

## input_consumed (2)

| id | source | reason |
|---|---|---|
| input_time_grids | `TIME`, `TIME3` | native time bases for every writer |
| input_geometry | `X`, `XB`, `DVOL`, `DAREA`, `SURF`, `BZ` | coordinates and zone geometry used to normalise and integrate |

## deferred (10)

| id | source | reason |
|---|---|---|
| deferred_edge_species_ion_sources | `SIGF_?_?`, `SIRC_?_?` | species-resolved ion sources are read to build a neutral-loss term that is not exported; the electron source SE* sets the particle stoichiometry because SI* = SE* + SFCX* |
| deferred_boundary_sources | `T0RECYC`, `T0GASFL`, `OM0RECYC`, `OM0GASFL`, `E0INR` | neutral/source-state quantities with no audited IMAS destination |
| deferred_boundary_energy | `PEEDG`, `PIEDG` | energy crossing the plasma boundary, not a volumetric source |
| deferred_neoclassical_velocity | `VPOL?_NC`, `VTOR?_NC` | RMAJM-grid (101 midplane points) velocities; core_profiles reads them and exports only if a CDF provides them on the profile grid |
| deferred_aggregate_fluxes | `GFLNC_I`, `GFLNC_X`, `GFLNC_TOK`, `QFLNC_X`, `QFLNC_TOK`, `QFLNCC_X`, `QFLNCC_TOK` | aggregate divergences; no typed total destination and no extra ion species is invented |
| deferred_ic_spectrum | `PIC_NPHI`, `PICF??N??` | frequency/n_phi-resolved IC families need a dedicated spectral mapping audit |
| deferred_ripple_torque | `TQRPL`, `FTQRIP` | ripple torque weight FTQRIP is not in the CDF; TQRPL is read only to warn that it is excluded from the NBI momentum sum |
| deferred_impurity_charge_profile | `ZIMPS_*` | mean impurity charge profile is read for bookkeeping, but ion/z_ion_1d is deliberately not written (the per-charge-state NIMP_* densities carry the charge resolution) |
| deferred_symmetric_major_radius | `RMJSYM` | read only for an angular-velocity to velocity conversion that is not used; rotation is exported as frequency |
| deferred_sawtooth_adhoc_boundary | `XSAWBOUND` | TRANSP labels it an ad-hoc sawtooth boundary; nothing proves it is a mixing or inversion radius, so it is only counted and a warning is logged |

