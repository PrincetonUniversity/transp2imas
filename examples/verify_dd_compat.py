#!/usr/bin/python
"""Verify IMAS schema layouts used by TRANSP2IMAS.

This is a schema-only diagnostic. It creates IDS objects from the requested
IDSFactory versions but does not open, create, or modify an IMAS database.

Examples
--------
    python examples/verify_dd_compat.py 3.39.0 4.1.1
    python examples/verify_dd_compat.py 4.1.1
"""

import sys
import imas


def version_tuple(version):
    out = []
    for token in str(version).split(".")[:3]:
        digits = "".join(ch for ch in token if ch.isdigit())
        out.append(int(digits) if digits else 0)
    while len(out) < 3:
        out.append(0)
    return tuple(out)


def yn(value):
    return "yes" if value else "no"


def check(condition, message, errors):
    if not condition:
        errors.append(message)


def inspect_version(version):
    factory = imas.IDSFactory(version)
    actual = str(factory.version)
    vt = version_tuple(actual)
    errors = []

    print(f"\nDD {actual}")
    print("-" * (3 + len(actual)))

    ec = factory.new("ec_launchers")
    ec.beam.resize(1)
    beam = ec.beam[0]
    ec_fields = {
        "identifier": hasattr(beam, "identifier"),
        "name": hasattr(beam, "name"),
        "description": hasattr(beam, "description"),
        "mode": hasattr(beam, "mode"),
        "o_mode_fraction": hasattr(beam, "o_mode_fraction"),
        "polarization": hasattr(beam, "polarization"),
    }
    print("ec_launchers/beam:")
    for key, value in ec_fields.items():
        print(f"  {key:24s} {yn(value)}")
    if ec_fields["polarization"]:
        pfrac = hasattr(beam.polarization, "o_mode_fraction")
        ptime = hasattr(beam.polarization, "time")
        print(f"  {'polarization/o_mode_fraction':24s} {yn(pfrac)}")
        print(f"  {'polarization/time':24s} {yn(ptime)}")
    else:
        pfrac = ptime = False

    if vt[0] == 3:
        check(ec_fields["identifier"], "DD3 EC beam lacks identifier", errors)
        check(ec_fields["name"], "DD3 EC beam lacks name", errors)
        check(ec_fields["mode"], "DD3 EC beam lacks mode", errors)
    elif vt[0] >= 4:
        check(ec_fields["name"], "DD4 EC beam lacks name", errors)
        check(ec_fields["description"], "DD4 EC beam lacks description", errors)
        check(not ec_fields["mode"], "DD4 EC beam unexpectedly still has mode", errors)
        if vt >= (4, 1, 0):
            check(pfrac and ptime, "DD4.1+ EC polarization fields are incomplete", errors)

    nbi = factory.new("nbi")
    nbi.unit.resize(1)
    unit = nbi.unit[0]
    nbi_fields = {
        "identifier": hasattr(unit, "identifier"),
        "name": hasattr(unit, "name"),
        "description": hasattr(unit, "description"),
        "species_label": hasattr(unit.species, "label"),
        "species_name": hasattr(unit.species, "name"),
    }
    print("nbi/unit:")
    for key, value in nbi_fields.items():
        print(f"  {key:24s} {yn(value)}")
    if vt[0] == 3:
        check(nbi_fields["identifier"], "DD3 NBI unit lacks identifier", errors)
        if vt < (3, 42, 0):
            check(nbi_fields["species_label"], "early DD3 NBI species lacks label", errors)
    elif vt[0] >= 4:
        check(nbi_fields["name"] and nbi_fields["description"],
              "DD4 NBI unit lacks name/description", errors)
        check(nbi_fields["species_name"], "DD4 NBI species lacks name", errors)

    cp = factory.new("core_profiles")
    cp.profiles_1d.resize(1)
    prof = cp.profiles_1d[0]
    prof.ion.resize(1)
    prof.neutral.resize(1)
    cp_fields = {
        "ion_label": hasattr(prof.ion[0], "label"),
        "ion_name": hasattr(prof.ion[0], "name"),
        "neutral_label": hasattr(prof.neutral[0], "label"),
        "neutral_name": hasattr(prof.neutral[0], "name"),
        "neutral_density": hasattr(prof.neutral[0], "density"),
        "neutral_temperature": hasattr(prof.neutral[0], "temperature"),
        "neutral_toroidal_velocity": (hasattr(prof.neutral[0], "velocity")
                                      and hasattr(prof.neutral[0].velocity, "toroidal")),
        "j_tor": hasattr(prof, "j_tor"),
        "j_phi": hasattr(prof, "j_phi"),
        "pressure_thermal": hasattr(prof, "pressure_thermal"),
        "pressure_perpendicular": hasattr(prof, "pressure_perpendicular"),
        "pressure_parallel": hasattr(prof, "pressure_parallel"),
        "n_i_thermal_total": hasattr(prof, "n_i_thermal_total"),
        "t_i_average": hasattr(prof, "t_i_average"),
        "magnetic_shear": hasattr(prof, "magnetic_shear"),
        "conductivity_parallel": hasattr(prof, "conductivity_parallel"),
        "phi_potential": hasattr(prof, "phi_potential"),
        "e_field_parallel": hasattr(prof, "e_field") and hasattr(prof.e_field, "parallel"),
        "momentum_phi": hasattr(prof, "momentum_phi"),
        "momentum_tor": hasattr(prof, "momentum_tor"),
        "code_parameters": hasattr(cp, "code") and hasattr(cp.code, "parameters"),
    }
    print("core_profiles:")
    for key, value in cp_fields.items():
        print(f"  {key:24s} {yn(value)}")
    if vt < (3, 42, 0):
        check(cp_fields["ion_label"] and cp_fields["neutral_label"],
              "early DD3 core_profiles species labels missing", errors)
        check(cp_fields["j_tor"], "early DD3 core_profiles lacks j_tor", errors)
    elif vt >= (3, 42, 0):
        check(cp_fields["ion_name"] and cp_fields["neutral_name"],
              "late-DD3/DD4 core_profiles species names missing", errors)
        check(cp_fields["j_phi"], "late-DD3/DD4 core_profiles lacks j_phi", errors)
    for field in ("neutral_density", "neutral_temperature"):
        check(cp_fields[field], f"core_profiles lacks {field}", errors)
    if vt[0] >= 4:
        print("  DD4 neutral rotation: supplementary core_profiles.code.parameters XML")
    elif not cp_fields["neutral_toroidal_velocity"]:
        print("  WARNING: DD3 neutral velocity.toroidal unavailable; rotation retained as XML")
    for field in ("pressure_thermal", "pressure_perpendicular", "pressure_parallel",
                  "n_i_thermal_total", "t_i_average", "magnetic_shear",
                  "conductivity_parallel", "phi_potential", "e_field_parallel",
                  "code_parameters"):
        check(cp_fields[field], f"core_profiles lacks {field}", errors)
    if vt < (3, 42, 0):
        check(cp_fields["momentum_tor"], "early-DD3 core_profiles lacks momentum_tor", errors)
    else:
        check(cp_fields["momentum_phi"], "late-DD3/DD4 core_profiles lacks momentum_phi", errors)

    # All populated IDSs receive whole-IDS source provenance.  The reference
    # AoS replaced the legacy string sources array in the DD evolution.
    prov = cp.ids_properties.provenance
    prov.node.resize(1)
    pnode = prov.node[0]
    provenance_fields = {
        "node_path": hasattr(pnode, "path"),
        "reference": hasattr(pnode, "reference"),
        "sources": hasattr(pnode, "sources"),
    }
    if provenance_fields["reference"]:
        pnode.reference.resize(1)
        pref = pnode.reference[0]
        provenance_fields["reference_name"] = hasattr(pref, "name")
        provenance_fields["reference_timestamp"] = hasattr(pref, "timestamp")
    else:
        provenance_fields["reference_name"] = False
        provenance_fields["reference_timestamp"] = False
    print("ids_properties/provenance:")
    for key, value in provenance_fields.items():
        print(f"  {key:24s} {yn(value)}")
    check(provenance_fields["node_path"], "ids_properties provenance node lacks path", errors)
    if provenance_fields["reference"]:
        check(provenance_fields["reference_name"], "provenance reference lacks name", errors)
    else:
        check(provenance_fields["sources"], "provenance lacks both reference and legacy sources", errors)

    eq = factory.new("equilibrium")
    eq.time_slice.resize(1)
    ts = eq.time_slice[0]
    gq = ts.global_quantities
    ma = gq.magnetic_axis
    ts.profiles_2d.resize(1)
    p2 = ts.profiles_2d[0]
    ep1 = ts.profiles_1d
    eq_fields = {
        "axis_b_field_tor": hasattr(ma, "b_field_tor"),
        "axis_b_field_phi": hasattr(ma, "b_field_phi"),
        "psi_axis": hasattr(gq, "psi_axis"),
        "psi_magnetic_axis": hasattr(gq, "psi_magnetic_axis"),
        "beta_normal": hasattr(gq, "beta_normal"),
        "beta_tor_norm": hasattr(gq, "beta_tor_norm"),
        "boundary_separatrix": hasattr(ts, "boundary_separatrix"),
        "p2_b_field_tor": hasattr(p2, "b_field_tor"),
        "p2_b_field_phi": hasattr(p2, "b_field_phi"),
        "psi_norm": hasattr(ep1, "psi_norm"),
        "rho_tor": hasattr(ep1, "rho_tor"),
        "dpsi_drho_tor": hasattr(ep1, "dpsi_drho_tor"),
        "dpressure_dpsi": hasattr(ep1, "dpressure_dpsi"),
        "f_df_dpsi": hasattr(ep1, "f_df_dpsi"),
        "j_parallel": hasattr(ep1, "j_parallel"),
        "magnetic_shear": hasattr(ep1, "magnetic_shear"),
        "mass_density": hasattr(ep1, "mass_density"),
    }
    print("equilibrium:")
    for key, value in eq_fields.items():
        print(f"  {key:24s} {yn(value)}")
    if vt < (3, 42, 0):
        check(eq_fields["axis_b_field_tor"] and eq_fields["p2_b_field_tor"],
              "early DD3 equilibrium lacks b_field_tor", errors)
    else:
        check(eq_fields["axis_b_field_phi"] and eq_fields["p2_b_field_phi"],
              "late-DD3/DD4 equilibrium lacks b_field_phi", errors)
    if vt >= (4, 1, 0):
        check(eq_fields["psi_magnetic_axis"], "DD4.1+ lacks psi_magnetic_axis", errors)
    if vt[0] >= 4:
        check(eq_fields["beta_tor_norm"], "DD4 equilibrium lacks beta_tor_norm", errors)
    else:
        check(eq_fields["beta_normal"], "DD3 equilibrium lacks beta_normal", errors)
    for field in ("rho_tor", "dpsi_drho_tor", "dpressure_dpsi",
                  "f_df_dpsi", "j_parallel", "magnetic_shear", "mass_density"):
        check(eq_fields[field], f"equilibrium profiles_1d lacks {field}", errors)
    if vt[0] >= 4:
        check(eq_fields["psi_norm"], "DD4 equilibrium profiles_1d lacks psi_norm", errors)

    summary = factory.new("summary")
    sgq = summary.global_quantities
    summary_fields = {
        "energy_total": hasattr(sgq, "energy_total"),
        "energy_thermal": hasattr(sgq, "energy_thermal"),
        "energy_ion_total_thermal": hasattr(sgq, "energy_ion_total_thermal"),
        "energy_electrons_thermal": hasattr(sgq, "energy_electrons_thermal"),
        "energy_fast_perpendicular": hasattr(sgq, "energy_fast_perpendicular"),
        "energy_fast_parallel": hasattr(sgq, "energy_fast_parallel"),
        "power_ohm": hasattr(sgq, "power_ohm"),
        "power_radiated": hasattr(sgq, "power_radiated"),
        "power_radiated_inside_lcfs": hasattr(sgq, "power_radiated_inside_lcfs"),
        "power_bremsstrahlung": hasattr(sgq, "power_bremsstrahlung"),
        "power_line": hasattr(sgq, "power_line"),
        "power_synchrotron": hasattr(sgq, "power_synchrotron"),
        "power_launched_nbi": hasattr(summary, "heating_current_drive") and hasattr(summary.heating_current_drive, "power_launched_nbi"),
        "current_ohm": hasattr(sgq, "current_ohm"),
        "current_non_inductive": hasattr(sgq, "current_non_inductive"),
        "current_bootstrap": hasattr(sgq, "current_bootstrap"),
        "fusion_neutron_rates": hasattr(summary, "fusion") and hasattr(summary.fusion, "neutron_rates"),
        "code_parameters": hasattr(summary, "code") and hasattr(summary.code, "parameters"),
        "type": hasattr(summary, "type"),
        "description": hasattr(summary, "description"),
        "simulation": hasattr(summary, "simulation"),
        "machine": hasattr(summary, "machine"),
        "pulse": hasattr(summary, "pulse"),
    }
    print("summary:")
    for key, value in summary_fields.items():
        print(f"  {key:24s} {yn(value)}")
    for field, present in summary_fields.items():
        if field in ("type", "description", "simulation", "machine", "pulse") and vt < (4, 1, 0):
            continue
        check(present, f"summary lacks {field}", errors)
    if vt >= (4, 1, 0):
        check(summary_fields["type"] and hasattr(summary.type, "name") and hasattr(summary.type, "index"),
              "DD4.1+ summary data type identifier is incomplete", errors)
        check(summary_fields["simulation"] and hasattr(summary.simulation, "time_begin")
              and hasattr(summary.simulation, "time_end") and hasattr(summary.simulation, "time_step")
              and hasattr(summary.simulation, "workflow"),
              "DD4.1+ summary simulation metadata is incomplete", errors)

    # NUBEAM fast-ion birth/deposition moments are exported to
    # distribution_sources even when a full phase-space FBM is unavailable.
    ds = factory.new("distribution_sources")
    ds.source.resize(1)
    dsrc = ds.source[0]
    dsrc.process.resize(1)
    dsrc.profiles_1d.resize(1)
    dsrc.global_quantities.resize(1)
    dsp = dsrc.profiles_1d[0]
    dsg = dsrc.global_quantities[0]
    ds_fields = {
        "process_type": hasattr(dsrc.process[0], "type"),
        "process_nbi_energy": hasattr(dsrc.process[0], "nbi_energy"),
        "process_nbi_unit": hasattr(dsrc.process[0], "nbi_unit"),
        "species": hasattr(dsrc, "species"),
        "profile_particles": hasattr(dsp, "particles"),
        "profile_energy": hasattr(dsp, "energy"),
        "profile_momentum_phi": hasattr(dsp, "momentum_phi"),
        "profile_momentum_tor": hasattr(dsp, "momentum_tor"),
        "global_particles": hasattr(dsg, "particles"),
        "global_power": hasattr(dsg, "power"),
        "global_torque_phi": hasattr(dsg, "torque_phi"),
        "global_torque_tor": hasattr(dsg, "torque_tor"),
        "shinethrough_power": hasattr(dsg, "shinethrough") and hasattr(dsg.shinethrough, "power"),
        "shinethrough_particles": hasattr(dsg, "shinethrough") and hasattr(dsg.shinethrough, "particles"),
        "code_parameters": hasattr(ds, "code") and hasattr(ds.code, "parameters"),
    }
    print("distribution_sources:")
    for key, value in ds_fields.items():
        print(f"  {key:32s} {yn(value)}")
    for field in ("process_type", "process_nbi_energy", "process_nbi_unit", "species",
                  "profile_particles", "profile_energy", "global_particles", "global_power",
                  "shinethrough_power", "shinethrough_particles", "code_parameters"):
        check(ds_fields[field], f"distribution_sources lacks {field}", errors)
    if vt < (3, 42, 0):
        check(ds_fields["profile_momentum_tor"] and ds_fields["global_torque_tor"],
              "early-DD3 distribution_sources lacks toroidal torque fields", errors)
    else:
        check(ds_fields["profile_momentum_phi"] and ds_fields["global_torque_phi"],
              "late-DD3/DD4 distribution_sources lacks phi torque fields", errors)

    ct = factory.new("core_transport")
    ct.model.resize(2)
    ctm = ct.model[1]
    ctm.profiles_1d.resize(1)
    ctp = ctm.profiles_1d[0]
    ctp.ion.resize(1)
    ct_fields = {
        "model_identifier": hasattr(ctm, "identifier"),
        "flux_multiplier": hasattr(ctm, "flux_multiplier"),
        "electron_particle_flux": hasattr(ctp.electrons.particles, "flux"),
        "electron_energy_flux": hasattr(ctp.electrons.energy, "flux"),
        "total_ion_energy_flux": hasattr(ctp.total_ion_energy, "flux"),
        "ion_particle_flux": hasattr(ctp.ion[0].particles, "flux"),
        "ion_energy_flux": hasattr(ctp.ion[0].energy, "flux"),
        "model_comment": hasattr(ctm, "comment"),
        "aggregate_total_ion_particles": hasattr(ctp, "total_ion_particles"),
    }
    print("core_transport:")
    for key, value in ct_fields.items():
        print(f"  {key:32s} {yn(value)}")
    for field in ("model_identifier", "flux_multiplier", "model_comment", "electron_particle_flux",
                  "electron_energy_flux", "total_ion_energy_flux", "ion_particle_flux",
                  "ion_energy_flux"):
        check(ct_fields[field], f"core_transport lacks {field}", errors)
    # The v8.5.2 mapper intentionally relies on there being no aggregate ion
    # particle-flux node; IPTR_MOD must not be forced into ion[0].
    check(not ct_fields["aggregate_total_ion_particles"],
          "core_transport unexpectedly has aggregate total_ion_particles; revisit IPTR_MOD mapping", errors)

    radiation = factory.new("radiation")
    radiation.process.resize(1)
    rp = radiation.process[0]
    rp.global_quantities.resize(1)
    rp.profiles_1d.resize(1)
    r1 = rp.profiles_1d[0]
    r1.ion.resize(1)
    r1.ion[0].state.resize(1)
    r1.neutral.resize(1)
    ri = r1.ion[0]
    rs = ri.state[0]
    rn = r1.neutral[0]
    rad_fields = {
        "identifier": hasattr(rp, "identifier"),
        "global_power": hasattr(rp.global_quantities[0].inside_lcfs, "power"),
        "electron_emissivity": hasattr(r1.electrons, "emissivity"),
        "electron_power_inside": hasattr(r1.electrons, "power_inside"),
        "emissivity_ion_total": hasattr(r1, "emissivity_ion_total"),
        "power_inside_ion_total": hasattr(r1, "power_inside_ion_total"),
        "emissivity_neutral_total": hasattr(r1, "emissivity_neutral_total"),
        "power_inside_neutral_total": hasattr(r1, "power_inside_neutral_total"),
        "ion_emissivity": hasattr(ri, "emissivity"),
        "ion_power_inside": hasattr(ri, "power_inside"),
        "ion_state_z_min": hasattr(rs, "z_min"),
        "ion_state_z_max": hasattr(rs, "z_max"),
        "ion_state_emissivity": hasattr(rs, "emissivity"),
        "ion_state_power_inside": hasattr(rs, "power_inside"),
        "neutral_emissivity": hasattr(rn, "emissivity"),
        "neutral_power_inside": hasattr(rn, "power_inside"),
    }
    print("radiation:")
    for key, value in rad_fields.items():
        print(f"  {key:32s} {yn(value)}")
    for field, present in rad_fields.items():
        check(present, f"radiation lacks {field}", errors)

    edge = factory.new("edge_sources")
    edge.grid_ggd.resize(1)
    eg = edge.grid_ggd[0]
    eg.space.resize(1)
    esp = eg.space[0]
    esp.objects_per_dimension.resize(1)
    eobj = esp.objects_per_dimension[0]
    eobj.object.resize(1)
    eg.grid_subset.resize(1)
    edge.source.resize(1)
    esrc = edge.source[0]
    esrc.ggd.resize(1)
    eggd = esrc.ggd[0]
    eggd.electrons.particles.resize(1)
    eggd.ion.resize(1)
    eggd.ion[0].particles.resize(1)
    eggd.ion[0].energy.resize(1)
    eggd.ion[0].momentum.resize(1)
    eggd.neutral.resize(1)
    eggd.neutral[0].particles.resize(1)
    eggd.neutral[0].energy.resize(1)
    eggd.neutral[0].momentum.resize(1)
    eggd.total_ion_energy.resize(1)
    eggd.momentum.resize(1)
    emom = eggd.ion[0].momentum[0]
    nmom = eggd.neutral[0].momentum[0]
    tmom = eggd.momentum[0]
    edge_fields = {
        "grid_identifier": hasattr(eg, "identifier"),
        "space_identifier": hasattr(esp, "identifier"),
        "geometry_type": hasattr(esp, "geometry_type"),
        "coordinates_type": hasattr(esp, "coordinates_type"),
        "geometry_content": hasattr(eobj, "geometry_content"),
        "object_geometry": hasattr(eobj.object[0], "geometry"),
        "grid_subset_identifier": hasattr(eg.grid_subset[0], "identifier"),
        "source_identifier": hasattr(esrc, "identifier"),
        "source_species": hasattr(esrc, "species"),
        "electron_particles_values": hasattr(eggd.electrons.particles[0], "values"),
        "ion_particles_values": hasattr(eggd.ion[0].particles[0], "values"),
        "neutral_particles_values": hasattr(eggd.neutral[0].particles[0], "values"),
        "ion_energy_values": hasattr(eggd.ion[0].energy[0], "values"),
        "neutral_energy_values": hasattr(eggd.neutral[0].energy[0], "values"),
        "total_ion_energy_values": hasattr(eggd.total_ion_energy[0], "values"),
        "ion_momentum_phi": hasattr(emom, "phi"),
        "ion_momentum_toroidal": hasattr(emom, "toroidal"),
        "neutral_momentum_phi": hasattr(nmom, "phi"),
        "neutral_momentum_toroidal": hasattr(nmom, "toroidal"),
        "total_momentum_phi": hasattr(tmom, "phi"),
        "total_momentum_toroidal": hasattr(tmom, "toroidal"),
    }
    print("edge_sources:")
    for key, value in edge_fields.items():
        print(f"  {key:32s} {yn(value)}")
    for field in ("grid_identifier", "space_identifier", "geometry_type", "coordinates_type",
                  "geometry_content", "object_geometry", "grid_subset_identifier",
                  "source_identifier", "source_species", "electron_particles_values",
                  "ion_particles_values", "neutral_particles_values", "ion_energy_values",
                  "neutral_energy_values", "total_ion_energy_values"):
        check(edge_fields[field], f"edge_sources lacks {field}", errors)
    if vt < (3, 42, 0):
        check(edge_fields["ion_momentum_toroidal"] and edge_fields["neutral_momentum_toroidal"]
              and edge_fields["total_momentum_toroidal"],
              "early-DD3 edge_sources lacks toroidal momentum fields", errors)
    else:
        check(edge_fields["ion_momentum_phi"] and edge_fields["neutral_momentum_phi"]
              and edge_fields["total_momentum_phi"],
              "late-DD3/DD4 edge_sources lacks phi momentum fields", errors)

    # coordinates_type changes from INT_1D in DD3 to identifier AoS in DD4.
    try:
        esp.coordinates_type.resize(1)
        if vt[0] >= 4:
            check(hasattr(esp.coordinates_type[0], "index") and hasattr(esp.coordinates_type[0], "name"),
                  "DD4 edge_sources coordinates_type is not identifier AoS", errors)
    except Exception as exc:
        errors.append(f"edge_sources coordinates_type resize failed: {exc}")

    waves = factory.new("waves")
    waves.coherent_wave.resize(1)
    cw = waves.coherent_wave[0]
    cw.global_quantities.resize(1)
    wg = cw.global_quantities[0]
    cw.profiles_1d.resize(1)
    wp = cw.profiles_1d[0]
    waves_fields = {
        "identifier_type": hasattr(cw, "identifier") and hasattr(cw.identifier, "type"),
        "antenna_name": hasattr(cw.identifier, "antenna_name"),
        "index_in_antenna": hasattr(cw.identifier, "index_in_antenna"),
        "wave_solver_type": hasattr(cw, "wave_solver_type"),
        "global_frequency": hasattr(wg, "frequency"),
        "global_power": hasattr(wg, "power"),
        "global_n_phi": hasattr(wg, "n_phi"),
        "global_n_tor": hasattr(wg, "n_tor"),
        "global_current_phi": hasattr(wg, "current_phi"),
        "global_current_tor": hasattr(wg, "current_tor"),
        "global_electron_power": hasattr(wg.electrons, "power_thermal"),
        "profile_power_density": hasattr(wp, "power_density"),
        "profile_power_inside": hasattr(wp, "power_inside"),
        "profile_n_phi": hasattr(wp, "n_phi"),
        "profile_n_tor": hasattr(wp, "n_tor"),
        "profile_current_phi_inside": hasattr(wp, "current_phi_inside"),
        "profile_current_tor_inside": hasattr(wp, "current_tor_inside"),
        "profile_current_parallel_density": hasattr(wp, "current_parallel_density"),
        "electron_power_density_thermal": hasattr(wp.electrons, "power_density_thermal"),
        "electron_power_inside_thermal": hasattr(wp.electrons, "power_inside_thermal"),
        "profile_time": hasattr(wp, "time"),
        "vacuum_toroidal_field": hasattr(waves, "vacuum_toroidal_field"),
        "generic_time": hasattr(waves, "time"),
        "code_parameters": hasattr(waves, "code") and hasattr(waves.code, "parameters"),
        "model_libraries": hasattr(waves, "code") and hasattr(waves.code, "library"),
    }
    print("waves:")
    for key, value in waves_fields.items():
        print(f"  {key:32s} {yn(value)}")
    for field in ("identifier_type", "antenna_name", "index_in_antenna",
                  "wave_solver_type", "global_frequency", "global_power",
                  "global_electron_power", "profile_power_density",
                  "profile_power_inside", "profile_current_parallel_density",
                  "electron_power_density_thermal", "electron_power_inside_thermal",
                  "profile_time", "vacuum_toroidal_field", "generic_time",
                  "code_parameters"):
        check(waves_fields[field], f"waves lacks {field}", errors)
    check(waves_fields["model_libraries"],
          "waves.code.library unavailable for model-specific input XML", errors)
    if waves_fields["model_libraries"]:
        waves.code.library.resize(1)
        model_parameters = hasattr(waves.code.library[0], "parameters")
        print(f"  model_library_parameters         {yn(model_parameters)}")
        check(model_parameters, "waves.code.library[].parameters unavailable", errors)
    if vt < (3, 42, 0):
        check(waves_fields["global_n_tor"] and waves_fields["profile_n_tor"],
              "early-DD3 waves lacks n_tor fields", errors)
        check(waves_fields["global_current_tor"] and waves_fields["profile_current_tor_inside"],
              "early-DD3 waves lacks toroidal-current fields", errors)
    else:
        check(waves_fields["global_n_phi"] and waves_fields["profile_n_phi"],
              "late-DD3/DD4 waves lacks n_phi fields", errors)
        check(waves_fields["global_current_phi"] and waves_fields["profile_current_phi_inside"],
              "late-DD3/DD4 waves lacks phi-current fields", errors)

    # Dataset-level metadata/provenance IDSs. They are optional only for DDs
    # predating the IDS itself; supported DDs should expose the structures used
    # by the converter without requiring fabricated PID/license information.
    try:
        fair = factory.new("dataset_fair")
        fair_fields = {
            "identifier": hasattr(fair, "identifier"),
            "replaces": hasattr(fair, "replaces"),
            "is_replaced_by": hasattr(fair, "is_replaced_by"),
            "rights_holder": hasattr(fair, "rights_holder"),
            "license": hasattr(fair, "license"),
            "provenance": hasattr(fair.ids_properties, "provenance"),
        }
        print("dataset_fair:")
        for key, value in fair_fields.items():
            print(f"  {key:32s} {yn(value)}")
        for field, present in fair_fields.items():
            check(present, f"dataset_fair lacks {field}", errors)
    except Exception as exc:
        if vt >= (3, 30, 0):
            errors.append(f"dataset_fair unavailable in supported DD: {exc}")
        else:
            print("dataset_fair: unavailable (predates IDS)")

    try:
        workflow = factory.new("workflow")
        workflow.time_loop.component.resize(1)
        wcomp = workflow.time_loop.component[0]
        workflow_fields = {
            "component_name": hasattr(wcomp, "name"),
            "component_version": hasattr(wcomp, "version"),
            "component_parameters": hasattr(wcomp, "parameters"),
            "component_library": hasattr(wcomp, "library"),
            "time_end": hasattr(workflow.time_loop, "time_end"),
            "code_parameters": hasattr(workflow, "code") and hasattr(workflow.code, "parameters"),
        }
        print("workflow:")
        for key, value in workflow_fields.items():
            print(f"  {key:32s} {yn(value)}")
        for field, present in workflow_fields.items():
            if field == "component_library" and vt[0] < 4:
                continue
            check(present, f"workflow lacks {field}", errors)
    except Exception as exc:
        if vt >= (3, 34, 0):
            errors.append(f"workflow unavailable in supported DD: {exc}")
        else:
            print("workflow: unavailable (predates IDS)")

    # Active PF hardware: each waveform has its own time axis, coil geometry is
    # explicitly element-based and circuit connections differ by DD version.
    active = factory.new("pf_active")
    active.coil.resize(1)
    active.supply.resize(1)
    active.circuit.resize(1)
    ac = active.coil[0]
    ac.element.resize(1)
    asupply = active.supply[0]
    acircuit = active.circuit[0]
    active_fields = {
        "coil_current_data": hasattr(ac.current, "data"),
        "coil_current_time": hasattr(ac.current, "time"),
        "coil_voltage_time": hasattr(ac.voltage, "time"),
        "coil_elements": hasattr(ac.element[0], "turns_with_sign"),
        "coil_geometry": hasattr(ac.element[0], "geometry") and
                         hasattr(ac.element[0].geometry, "rectangle"),
        "supply_current_time": hasattr(asupply.current, "time"),
        "supply_voltage_time": hasattr(asupply.voltage, "time"),
        "circuit_connections": hasattr(acircuit, "connections"),
        "circuit_current_time": hasattr(acircuit.current, "time"),
        "generic_time": hasattr(active, "time"),
        "code_parameters": hasattr(active, "code") and hasattr(active.code, "parameters"),
    }
    print("pf_active:")
    for key, present in active_fields.items():
        print(f"  {key:32s} {yn(present)}")
        check(present, f"pf_active lacks {key}", errors)

    tf = factory.new("tf")
    tf_tor = hasattr(tf, "b_field_tor_vacuum_r")
    tf_phi = hasattr(tf, "b_field_phi_vacuum_r")
    print("tf:")
    print(f"  {'b_field_tor_vacuum_r':24s} {yn(tf_tor)}")
    print(f"  {'b_field_phi_vacuum_r':24s} {yn(tf_phi)}")
    if vt < (3, 42, 0):
        check(tf_tor, "early DD3 TF lacks b_field_tor_vacuum_r", errors)
    else:
        check(tf_phi, "late-DD3/DD4 TF lacks b_field_phi_vacuum_r", errors)

    wall = factory.new("wall")
    wall.description_2d.resize(1)
    wall.description_2d[0].limiter.unit.resize(1)
    unit = wall.description_2d[0].limiter.unit[0]
    closed = hasattr(unit, "closed")
    print("wall:")
    print(f"  {'limiter/unit/closed':24s} {yn(closed)}")
    if vt[0] == 3:
        check(closed, "DD3 wall limiter lacks closed flag", errors)
    elif vt[0] >= 4:
        check(not closed, "DD4 wall unexpectedly has legacy limiter/unit/closed", errors)

    if errors:
        print("  RESULT: FAIL")
        for error in errors:
            print(f"    - {error}")
        return False

    print("  RESULT: PASS")
    return True


if __name__ == "__main__":
    versions = sys.argv[1:] or ["3.39.0", "4.1.1"]
    ok = True
    for version in versions:
        try:
            ok = inspect_version(version) and ok
        except Exception as exc:
            ok = False
            print(f"\nDD {version}: ERROR: {type(exc).__name__}: {exc}")
    raise SystemExit(0 if ok else 1)
