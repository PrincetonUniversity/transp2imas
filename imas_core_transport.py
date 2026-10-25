#!/usr/bin/python
"""Populate the IMAS core_transport IDS from observed and modeled TRANSP transport quantities."""

import xml.etree.ElementTree as ET

import numpy as np
from nc2imas import *
from writer_base import IdsWriter
from pt_solver_diagnostics import append_runtime_diagnostics_xml
from pt_solver_regions import (
    build_pt_region_plan, append_pt_regions_xml, annotate_theory_model_xml,
    append_pt_regions_reference, interval_attributes, occurrence_region_models,
    region_occurrence_xml, region_validity_xml, regional_model_slots,
)
from theory_transport_mapping import (
    COEFFICIENT_MODEL_SPECS, THEORY_AUXILIARY_SPECS, build_theory_plans,
    theory_model_xml, validate_profile, divergence_to_surface_flux,
)


class imas_core_transport(IdsWriter):
    """Map TRANSP effective/observed transport and predictive model fluxes."""

    def _load_grid(self):
        """Read the grids the models are written on (done in ``build``, not at construction)."""
        self.time=np.asarray(self.cdf.get_ncvar('TIME3'),dtype=np.float64)
        self.x=np.asarray(self.cdf.get_ncvar('X'),dtype=np.float64)
        self.xb=np.asarray(self.cdf.get_ncvar('XB'),dtype=np.float64)
        self.nt=len(self.time)
        self.nx=self.xb.shape[1]
        self.dvol=self._get('DVOL',self.grid.vconv)
        self.surf=self._get('SURF',self.grid.aconv)

    def _get(self,name,mult=None):
        """A CDF profile as float64; one whose radial size is not the writer's is a named broken mapping."""
        a=self.cdf.get_ncvar(name)
        if not self.fill.is_valid(a): return np.array([])
        a=np.asarray(self.check_profile(name,np.asarray(a,dtype=np.float64),self.nt,self.nx))
        return a if mult is None else a*mult

    def _first_available(self,names,mult=None):
        for name in names:
            if self.cdf.check_var_name(name):
                return name,self._get(name,mult)
        return '',np.array([])

    def _flux_from_divergence(self,name,mult,physical_sign=1.0):
        """Convert a TRANSP flux divergence on X to average flux on XB.

        TRANSP stores EPTR_*, EETR_*, IETR_* and related quantities as
        divergences per unit volume.  The outward surface-averaged flux is the
        enclosed volume integral divided by the flux-surface area SURF.
        """
        div=self._get(name,mult)
        if not self.fill.is_valid(div,self.dvol,self.surf): return np.array([])
        integrated=np.cumsum(div*self.dvol,axis=1)
        return physical_sign*np.divide(integrated,self.surf,out=np.zeros_like(integrated),where=np.abs(self.surf)>0)

    def _solver_flux(self,name,physical_sign=1.0):
        """Flux of a predictive ``*_MOD`` channel from its divergence, or empty when the channel is all zero.

        A run that did not solve a channel saves zeros for it (every NSTX case has all-zero EPTR/EETR/AMTR_MOD,
        and NSTX 141711Y05 solves only the ion heat).  Writing those would state that the solver predicted no
        flux, so an inactive channel is left out.  The flux is the volume integral of the divergence over the
        surface area, in N m m^-2 for the angular momentum and W m^-2, m^-2 s^-1 for the heat and particles.
        """
        if not self._has_nonzero(name):
            return np.array([])
        return self._flux_from_divergence(name,self.grid.vol_inv,physical_sign=physical_sign)

    def _has_nonzero(self,name):
        """True when CDF variable ``name`` exists and holds a nonzero value somewhere."""
        values=self._get(name,None)
        return bool(values.size and np.any(values))

    def _enclosed_geometry(self):
        """``(volume, area)`` enclosed by each XB surface for all times, computed once per grid.

        ``_fill_grid`` runs for every model and time slice and used to re-read DAREA and re-sum the whole
        time-by-radius arrays each time (5800 reads on NSTX 142202U01).  The cache is keyed on the identity of
        ``self.dvol``, so a writer given another grid recomputes.
        """
        cache=self.__dict__.get('_enclosed_cache')
        if cache is None or cache[0] is not self.dvol:
            area=self._get('DAREA',self.grid.aconv)
            cache=self._enclosed_cache=(self.dvol,self.grid.cumulative_zone_quantity(self.dvol),
                                        self.grid.cumulative_zone_quantity(area))
        return cache[1],cache[2]

    def _fill_grid(self,grid,it):
        self.fill.fill_data(grid.rho_tor_norm,self.xb,ind=it,chk=self.nx)
        # Add geometry where readily available on the same XB grid.
        vol_b,area_b=self._enclosed_geometry()
        self.fill.fill_data(grid.volume,vol_b,ind=it,chk=self.nx)
        self.fill.fill_data(grid.area,area_b,ind=it,chk=self.nx)
        self.fill.fill_data(grid.surface,self.surf,ind=it,chk=self.nx)

    def _fill_ion_identity(self,node,sp):
        atom=atom_base(sp)
        node.element.resize(1)
        node.element[0].a=float(atom.a)
        self.fill.set_nuclear_charge(node.element[0],atom.z_n)
        node.element[0].atoms_n=int(atom.atoms_n)
        node.z_ion=float(atom.z_n)
        self.fill.set_species_name(node,atom.label)

    def _species_observed(self):
        specs=[]
        candidates=[
            ('H','DIFFH','VELH','PTRH_OBS'),
            ('D','DIFFD','VELD','PTRD_OBS'),
            ('T','DIFFT','VELT','PTRT_OBS'),
            ('HE3','DIFFHE3','VELHE3','PTRHE3_OBS'),
            ('HE4','DIFFHE4','VELHE4','PTRHE4_OBS'),
        ]
        for sp,dn,vn,fn in candidates:
            if self.cdf.check_var_name(dn) or self.cdf.check_var_name(vn) or self.cdf.check_var_name(fn):
                specs.append((sp,self._get(dn,self.grid.aconv),self._get(vn,self.grid.cm2m),
                              self._flux_from_divergence(fn,self.grid.vol_inv),fn))
        return specs

    def _species_model(self):
        """Return species-resolved predictive particle fluxes when available.

        TRANSP's aggregate IPTR_MOD is the divergence of the *total* ion
        particle flux.  core_transport has no aggregate total-ion-particle
        field, so IPTR_MOD must not be represented as a fictitious species.
        Instead, export the physically meaningful PTR*_MOD channels whenever
        TRANSP saved them.  Multiple historical spellings for He3/He4 are
        accepted.
        """
        candidates=[
            ('H',['PTRH_MOD']),
            ('D',['PTRD_MOD']),
            ('T',['PTRT_MOD']),
            ('HE3',['PTRHE3_MOD','PTR3_MOD']),
            ('HE4',['PTRHE4_MOD','PTR4_MOD']),
        ]
        specs=[]
        for sp,names in candidates:
            name,arr=self._first_available(names)
            if name and arr.size and np.any(arr):          # an all-zero channel was not solved: see _solver_flux
                specs.append((sp,self._flux_from_divergence(name,self.grid.vol_inv),name))
        return specs

    def _prepare_model(self,model,name,index,description,comment,flux_grid='XB'):
        model.identifier.name=name
        model.identifier.index=index
        model.identifier.description=description
        model.comment=comment
        model.flux_multiplier=0.0
        if hasattr(model,'code'):
            model.code.name='TRANSP'
            if hasattr(model.code,'description'):
                model.code.description=description
            if hasattr(model.code, 'parameters'):
                # A valid XML root also makes regional configuration references
                # safe when no earlier model-specific parameters were set.
                model.code.parameters = '<transp2imas_model_output/>'
        model.profiles_1d.resize(self.nt)
        for it in range(self.nt):
            p=model.profiles_1d[it]
            p.time=float(self.time[it])
            self._fill_grid(p.grid_d,it)
            self._fill_grid(p.grid_v,it)
            if flux_grid == 'X':
                # MMM reduced-mode arrays (gamma, omega, ky*rhos) reside on
                # TIME3/X, unlike transport diffusivities on TIME3/XB.
                # No typed fluxes are assigned for these MMM component models.
                # Store X once as the secondary native radial coordinate;
                # do not fabricate geometry/flux values on this grid.
                self.fill.fill_data(p.grid_flux.rho_tor_norm, self.x,
                                  ind=it, chk=self.x.shape[1])
            else:
                self._fill_grid(p.grid_flux,it)

    def _fill_combined_model(self,model):
        # Effective coefficients used in the TRANSP balance analysis.
        de=self._get('DIFFE',self.grid.aconv)
        ve=self._get('VELE',self.grid.cm2m)
        che=self._get('CONDE',self.grid.aconv)
        vhe=self._get('VEL_TE',self.grid.cm2m)
        chi=self._get('CONDI',self.grid.aconv)
        vhi=self._get('VEL_TI',self.grid.cm2m)
        chphi=self._get('CHPHI',self.grid.aconv)
        vmom=self._get('VMO_NET',self.grid.cm2m)

        epart_flux=self._flux_from_divergence('EPTR_OBS',self.grid.vol_inv)
        een_flux=self._flux_from_divergence('EETR_OBS',self.grid.vol_inv)
        ien_flux=self._flux_from_divergence('IETR_OBS',self.grid.vol_inv)
        mom_flux=self._flux_from_divergence('AMTR_OBS',self.grid.vol_inv,physical_sign=self.signs.ip_sign)

        # ETA_USE is ohm*cm on X; conductivity in core_transport is on grid_d (XB).
        eta=self._get('ETA_USE',0.01)
        sigma=np.array([])
        if eta.size:
            eta_b=self.grid.interp_radial_profiles(eta,self.x,self.xb)
            sigma=np.divide(1.0,eta_b,out=np.zeros_like(eta_b),where=np.abs(eta_b)>0)

        species_specs=self._species_observed()
        self._prepare_model(
            model,'combined',1,
            'Combination of available TRANSP effective transport and observed balance fluxes',
            'D and v are effective coefficients reported by TRANSP. Fluxes are reconstructed '
            'from observed divergence profiles by cumulative volume integration divided by SURF. '
            'Energy fluxes therefore already include the particle-convective contribution.'
        )

        for it in range(self.nt):
            p=model.profiles_1d[it]
            self.fill.fill_data(p.conductivity_parallel,sigma,ind=it,chk=self.nx)

            self.fill.fill_data(p.electrons.particles.d,de,ind=it,chk=self.nx)
            self.fill.fill_data(p.electrons.particles.v,ve,ind=it,chk=self.nx)
            self.fill.fill_data(p.electrons.particles.flux,epart_flux,ind=it,chk=self.nx)
            self.fill.fill_data(p.electrons.energy.d,che,ind=it,chk=self.nx)
            self.fill.fill_data(p.electrons.energy.v,vhe,ind=it,chk=self.nx)
            self.fill.fill_data(p.electrons.energy.flux,een_flux,ind=it,chk=self.nx)
            self.fill.fill_data(p.total_ion_energy.d,chi,ind=it,chk=self.nx)
            self.fill.fill_data(p.total_ion_energy.v,vhi,ind=it,chk=self.nx)
            self.fill.fill_data(p.total_ion_energy.flux,ien_flux,ind=it,chk=self.nx)

            mom=p.momentum_phi if self.dd.dd_at_least(3,42,0) else p.momentum_tor
            self.fill.fill_data(mom.d,chphi,ind=it,chk=self.nx)
            self.fill.fill_data(mom.v,vmom,ind=it,chk=self.nx)
            self.fill.fill_data(mom.flux,mom_flux,ind=it,chk=self.nx)

            p.ion.resize(len(species_specs))
            for ji,(sp,d,v,flux,_fn) in enumerate(species_specs):
                inode=p.ion[ji]
                self._fill_ion_identity(inode,sp)
                self.fill.fill_data(inode.particles.d,d,ind=it,chk=self.nx)
                self.fill.fill_data(inode.particles.v,v,ind=it,chk=self.nx)
                self.fill.fill_data(inode.particles.flux,flux,ind=it,chk=self.nx)
        return len(species_specs)

    def _fill_transport_solver_model(self,model):
        """Write predictive/model fluxes separately from observed balance data."""
        epart_flux=self._solver_flux('EPTR_MOD')
        een_flux=self._solver_flux('EETR_MOD')
        ien_flux=self._solver_flux('IETR_MOD')
        # Angular-momentum counterpart of the three channels above, converted exactly like the
        # observed AMTR_OBS (cumulative volume integral over SURF, IMAS toroidal sign).  It is a
        # divergence, so it is never copied into the flux node.  An all-zero AMTR_MOD only says
        # that the momentum model was inactive, so it is not written.
        mom_flux=self._solver_flux('AMTR_MOD',physical_sign=self.signs.ip_sign)
        species_specs=self._species_model()

        present=[]
        for name,arr in [('EPTR_MOD',epart_flux),('EETR_MOD',een_flux),('IETR_MOD',ien_flux),('AMTR_MOD',mom_flux)]:
            if arr.size: present.append(name)
        present.extend(name for _sp,_flux,name in species_specs)
        has_iptr=self.cdf.check_var_name('IPTR_MOD')
        # Do not create an empty solver model merely because IPTR_MOD exists:
        # that aggregate quantity has no faithful core_transport destination.
        if not present:
            return False,[]

        iptr_note=(
            ' IPTR_MOD is present in the CDF but is not written directly because IMAS core_transport '
            'has no aggregate total-ion particle flux node; species-resolved PTR*_MOD channels are '
            'written when available.' if has_iptr else ''
        )
        self._prepare_model(
            model,'transport_solver',2,
            'Predictive transport-solver fluxes used by TRANSP',
            'Fluxes (including the angular-momentum flux from AMTR_MOD) are reconstructed from TRANSP '
            '*_MOD divergence profiles by cumulative volume integration divided by SURF. No effective '
            'D/v coefficients are inferred from these modeled divergences.'+iptr_note
        )
        # This is the single assembled predictive-flux model, not three
        # independent radial-region model outputs or IDS occurrences.
        if self._pt_region_plan is not None:
            model.code.parameters = append_pt_regions_reference(
                model.code.parameters, self._pt_region_plan)

        for it in range(self.nt):
            p=model.profiles_1d[it]
            self.fill.fill_data(p.electrons.particles.flux,epart_flux,ind=it,chk=self.nx)
            self.fill.fill_data(p.electrons.energy.flux,een_flux,ind=it,chk=self.nx)
            self.fill.fill_data(p.total_ion_energy.flux,ien_flux,ind=it,chk=self.nx)
            if mom_flux.size:
                mom=p.momentum_phi if self.dd.dd_at_least(3,42,0) else p.momentum_tor
                self.fill.fill_data(mom.flux,mom_flux,ind=it,chk=self.nx)
            p.ion.resize(len(species_specs))
            for ji,(sp,flux,_name) in enumerate(species_specs):
                inode=p.ion[ji]
                self._fill_ion_identity(inode,sp)
                self.fill.fill_data(inode.particles.flux,flux,ind=it,chk=self.nx)

        return True,present

    def _neoclassical_flux(self, source, kind):
        """Integrate the explicitly labelled NC divergence on its native grid.

        Never treat ``GFLNC_*`` or ``QFLNC*`` as a flux density directly:
        the CDF metadata describes them as *divergences* on X.  Their
        volume integrals are fluxes on the XB grid in IMAS SI units.
        """
        raw = validate_profile(self.cdf.nc, source, kind, required_shape=self.dvol.shape)
        volume_cm3 = np.asarray(self.dvol, dtype=np.float64) / self.grid.vconv
        surf_cm2 = np.asarray(self.surf, dtype=np.float64) / self.grid.aconv
        return divergence_to_surface_flux(raw, volume_cm3, surf_cm2, kind)

    def _fill_neoclassical_model(self, model, channels, classical=False):
        """Dedicated neoclassical vs NC+classical *heat* flux contributions.

        The source contains separate aggregate ion ('I'), impurity ('X'),
        and species channels.  'I' is used for total_ion_energy only; it is
        never invented as an actual species.  'X'/'TOK' are aggregate and
        therefore never counted as additional ions alongside real elements.
        """
        particle_names = {} if classical else channels['particles']
        heat_names = channels['heat_nc_classical' if classical else 'heat']
        particle_names = {key: name for key, name in particle_names.items()
                          if key not in ('I', 'X', 'TOK')}
        heat_names = {key: name for key, name in heat_names.items()
                      if key not in ('X', 'TOK')}
        if not particle_names and not heat_names:
            return False, []
        species = sorted((set(particle_names) | set(heat_names)) - {'E', 'I', 'X', 'TOK'})
        known_species = []
        for token in species:
            atom = atom_base(token)
            if float(atom.a) <= 0 or float(atom.z_n) <= 0:
                self.logger.warning('core_transport: unidentified neoclassical species %s: skipped', token)
            else:
                known_species.append(token)
        used = list(particle_names.values()) + list(heat_names.values())
        notes = ('TRANSP NC-classical heat-flux divergence QFLNCC_* (kept separate from QFLNC_*); '
                 'not summed with QFLNC_* to avoid double counting.' if classical else
                 'Neoclassical particle (GFLNC_*) and heat (QFLNC_*) flux divergences.')
        notes += (' Source channels: ' + ', '.join(sorted(used)) + '. '
                  'Each is integrated with TRANSP DVOL and divided by SURF. '
                  'The *_I heat channel represents the aggregate thermal-ion energy, '
                  'whereas *_X, *_TOK, and the aggregate *_I particle channel '
                  'are not represented as fictitious ion species.')
        self._prepare_model(
            model, 'neoclassical_classical_heat' if classical else 'neoclassical',
            4 if classical else 3,
            'TRANSP NC-classical heat flux (QFLNCC)' if classical else
            'TRANSP neoclassical particle and heat transport', notes)

        if not classical and self._pt_region_plan is not None:
            model.code.parameters = annotate_theory_model_xml(
                model.code.parameters, self._pt_region_plan, 'neoclassical')

        particle_flux = {name: self._neoclassical_flux(name, 'particles')
                         for name in particle_names.values()}
        heat_kind = 'heat_nc_classical' if classical else 'heat'
        energy_flux = {name: self._neoclassical_flux(name, heat_kind)
                       for name in heat_names.values()}
        for it in range(self.nt):
            p = model.profiles_1d[it]
            if 'E' in particle_names:
                self.fill.fill_data(p.electrons.particles.flux,
                                  particle_flux[particle_names['E']], ind=it, chk=self.nx)
            if 'E' in heat_names:
                self.fill.fill_data(p.electrons.energy.flux,
                                  energy_flux[heat_names['E']], ind=it, chk=self.nx)
            if 'I' in heat_names:
                self.fill.fill_data(p.total_ion_energy.flux,
                                  energy_flux[heat_names['I']], ind=it, chk=self.nx)
            p.ion.resize(len(known_species))
            for j, sp in enumerate(known_species):
                ion = p.ion[j]
                self._fill_ion_identity(ion, sp)
                if sp in particle_names:
                    self.fill.fill_data(ion.particles.flux,
                                      particle_flux[particle_names[sp]], ind=it, chk=self.nx)
                if sp in heat_names:
                    self.fill.fill_data(ion.energy.flux,
                                      energy_flux[heat_names[sp]], ind=it, chk=self.nx)
        return True, used

    def _fill_theory_coefficient_model(self, model, family, spec, plan):
        """Populate one physical theory-model result, never an observation.

        Model coefficients are already flux-surface-boundary (XB) profiles;
        unlike TRANSP flux *divergences*, they must NOT be volume-integrated.
        Preserve source separation and convert cm^2/s and cm/s to SI.
        """
        fields = plan['fields']
        sources = [name for name in (*spec['coefficients'], *spec['auxiliary'])
                   if name in fields]
        if not sources:
            return False
        is_total = spec['name'] == 'mmm07_total'
        comment = (
            spec['label'] + '. TRANSP theory transport MODEL OUTPUT. '
            + ('MMM07 total is decomposed by child models and must NOT be '
               'summed with its components. ' if is_total else
               'MMM07 component of the total; do NOT double count. '
               if family == 'mmm07' else
               'Independent GLF23 model: do NOT add to MMM07 unless '
               'PT_SOLVER input explicitly specifies that combination. ')
            + 'CDF coefficients on TIME3/XB are in cm^2/s or cm/s; '
              'typed IMAS coefficients are SI. Supplemental reduced-model '
              'outputs reference native IMAS time/radial coordinates. '
            + 'Sources: ' + ', '.join(sources) + '.'
        )
        use_native_x = any(
            source in fields and THEORY_AUXILIARY_SPECS[source][1] == 'X'
            for source in spec['auxiliary']
        )
        self._prepare_model(model, spec['name'], spec['index'],
                            spec['description'], comment,
                            flux_grid='X' if use_native_x else 'XB')
        model.code.name = 'TRANSP/' + family.upper()
        model.code.description = spec['label']
        model.code.parameters = annotate_theory_model_xml(
            theory_model_xml(plan, spec, family), self._pt_region_plan,
            spec['name'])
        converted = {}
        for source, (destination, kind) in spec['coefficients'].items():
            if source not in fields:
                continue
            multiplier = 1.0e-4 if kind == 'diffusivity' else 1.0e-2
            converted[source] = (destination, fields[source] * multiplier)
        for it in range(self.nt):
            prof = model.profiles_1d[it]
            for destination, values in converted.values():
                target = prof
                for section in destination.split('.'):
                    if section == 'momentum_phi' and not self.dd.dd_at_least(3, 42, 0):
                        section = 'momentum_tor'
                    target = getattr(target, section)
                self.fill.fill_data(target, values, ind=it, chk=self.nx)
        return True

    def build(self):
        self._load_grid()
        ct=self.dd.new_ids('core_transport')
        ct=self.fill.fill_properties(ct)
        self.fill.fill_code_properties(ct)
        self.fill.fill_external_model_inputs(ct, 'core_transport')
        # PT_SOLVER region settings and TR.DAT boundaries are input metadata.
        # Store them once in the assembled core_transport code record, while
        # leaving the native TIME3/X/XB output grids in the typed model nodes.
        self._pt_region_plan = (build_pt_region_plan(self.namelist.input_inventory())
                                if hasattr(self.namelist, 'input_inventory') else None)
        if self._pt_region_plan is not None:
            ct.code.parameters = append_pt_regions_xml(ct.code.parameters,
                                                       self._pt_region_plan)
        # Runtime diagnostics are CDF output: discover them whether or not a PT_SOLVER model
        # namelist is available.  Channels of a run without a PT plan are kept only when they hold
        # nonzero data (all-zero channels exist in every CDF).
        ct.code.parameters = append_runtime_diagnostics_xml(
            ct.code.parameters, getattr(self.cdf, 'nc', None),
            record_inactive=self._pt_region_plan is not None)
        ct.ids_properties.homogeneous_time=1
        ct.time=np.asarray(self.time,dtype=np.float64)
        self.fill.fill_vac_tor_field(ct)

        # Important: resizing an IMAS array of structures after populating it
        # can recreate *all* children, including the combined model. Allocate
        # the FINAL number of slots once, BEFORE filling any model.
        plans = build_theory_plans(self.cdf.nc)
        channels = plans['neoclassical']
        solver_present = any(self._has_nonzero(n) for n in
                             ('EPTR_MOD', 'EETR_MOD', 'IETR_MOD', 'AMTR_MOD', 'PTRH_MOD',
                              'PTRD_MOD', 'PTRT_MOD', 'PTRHE3_MOD', 'PTR3_MOD',
                              'PTRHE4_MOD', 'PTR4_MOD'))
        nc_present = bool(channels['particles'] or channels['heat'])
        classical_present = bool(channels['heat_nc_classical'])
        theory_slots = []
        for family, specs in COEFFICIENT_MODEL_SPECS.items():
            plan = plans[family]
            if not plan['active']:
                continue
            for quantity, check in plan['closure'].items():
                if check['status'] == 'FAIL':
                    self.logger.warning(
                        'core_transport: %s %s component closure mismatch '
                        '(normalized discrepancy %.4g); preserving source data',
                        family, quantity, check['max_normalized_difference'])
            for spec in specs:
                if any(name in plan['fields'] for name in
                       (*spec['coefficients'], *spec['auxiliary'])):
                    theory_slots.append((family, spec, plan))
        expected_models = (1 + int(solver_present) + int(nc_present) +
                           int(classical_present) + len(theory_slots))
        ct.model.resize(expected_models)
        nobs = self._fill_combined_model(ct.model[0])
        next_model = 1
        has_solver, present = False, []
        if solver_present:
            has_solver, present = self._fill_transport_solver_model(ct.model[next_model])
            if not has_solver:
                raise RuntimeError('Predictive flux model was expected but empty')
            next_model += 1
        has_nc, nc_outputs = False, []
        if nc_present:
            has_nc, nc_outputs = self._fill_neoclassical_model(
                ct.model[next_model], channels, classical=False)
            if not has_nc:
                raise RuntimeError('Neoclassical model was expected but empty')
            next_model += 1
        has_classical, cl_present = False, []
        if classical_present:
            has_classical, cl_present = self._fill_neoclassical_model(
                ct.model[next_model], channels, classical=True)
            if not has_classical:
                raise RuntimeError('Classical heat model was expected but empty')
            next_model += 1
        for family, spec, plan in theory_slots:
            if not self._fill_theory_coefficient_model(
                    ct.model[next_model], family, spec, plan):
                raise RuntimeError('Registered transport model is empty: ' + spec['name'])
            next_model += 1
            self.logger.info('core_transport: %s model: %s', family, spec['name'])
        if next_model != expected_models or len(ct.model) != expected_models:
            raise RuntimeError('core_transport model allocation mismatch')
        if len(ct.model[0].profiles_1d) != self.nt:
            raise RuntimeError('Combined core_transport model lost its time profiles')
        if ct.model[0].identifier.name != 'combined':
            raise RuntimeError('Combined core_transport model identifier was lost')
        if self.cdf.check_var_name('IPTR_MOD'):
            self.logger.info('core_transport: IPTR_MOD present but aggregate total-ion particle flux has no IMAS core_transport node')
        self.logger.info(f'core_transport: combined model, {nobs} observed ion particle species')
        if has_solver:
            self.logger.info('core_transport: transport_solver model channels: '+', '.join(present))
        if has_nc:
            self.logger.info('core_transport: neoclassical channels: '+', '.join(nc_outputs))
        if has_classical:
            self.logger.info('core_transport: NC+classical heat channels: '+', '.join(cl_present))
        # Occurrence 0 is the complete assembled/effective TRANSP result.
        # Region occurrences are additional views with two *model slots*,
        # populated only when matching native CDF output channels exist.
        if self.ctx.supports_occurrences:
            return [(ct, 0)] + self._region_occurrences(plans)
        return [ct]

    # Native TRANSP output channels exist for these anomalous model families only.
    NATIVE_ANOMALOUS_FAMILIES = {'glf23': 'glf23', 'mmm': 'mmm07'}

    def _region_occurrences(self, plans):
        """``[(ids, occurrence)]`` for the region view occurrences (1, 2, ...)."""
        occurrences = []
        for occurrence, item in enumerate(occurrence_region_models(self._pt_region_plan), 1):
            regional = self.dd.new_ids('core_transport')
            regional = self.fill.fill_properties(regional)
            self.fill.fill_code_properties(regional)
            regional.code.parameters = region_occurrence_xml(
                self._pt_region_plan, item, occurrence)
            regional.ids_properties.homogeneous_time = 1
            regional.time = np.asarray(self.time, dtype=np.float64)
            self.fill.fill_vac_tor_field(regional)
            channels = plans['neoclassical']
            slots = regional_model_slots(item, bool(channels['particles'] or channels['heat']))
            regional.model.resize(len(slots))          # final size, once, before filling
            for index, (node, spec) in enumerate(zip(regional.model, slots), 1):
                self._fill_regional_slot(node, index, item, spec, plans)
            occurrences.append((regional, occurrence))
        return occurrences

    def _fill_regional_slot(self, slot, index, item, spec, plans):
        category, selected = spec['category'], spec['selected']
        names = [s['model'] for s in selected]
        assembled = spec.get('assembled', False)
        # Identity of the OUTPUT written in this slot, which is not necessarily the selected model's:
        # the assembled neoclassical result is named by what it is.
        identity = 'neoclassical' if assembled else names[0].upper()
        listed = ', '.join(n.upper() for n in names)
        description = (f'Assembled TRANSP neoclassical output; selected neoclassical model(s) {listed} '
                       f'in pt_{item["region"]} are not independently attributed in the CDF'
                       if assembled else
                       f'{category} transport model {listed} in pt_{item["region"]}')
        self._prepare_model(slot, identity, index, description,
                            'PT_SOLVER regional selection; outputs present only where separately '
                            'available in TRANSP CDF')
        populated = False
        if assembled:
            populated, _ = self._fill_neoclassical_model(slot, plans['neoclassical'], classical=False)
        elif category == 'anomalous' and names[0] in self.NATIVE_ANOMALOUS_FAMILIES:
            family = self.NATIVE_ANOMALOUS_FAMILIES[names[0]]
            matching = [model for model in COEFFICIENT_MODEL_SPECS[family]
                        if names[0] != 'mmm' or model['name'] == 'mmm07_total']
            for model in matching[:1]:
                populated = self._fill_theory_coefficient_model(slot, family, model, plans[family])
        if not populated:
            # No physical CDF arrays: keep an explicit settings-only model, not a
            # misleading populated time/radius grid.
            slot.profiles_1d.resize(0)
        slot.identifier.name = identity
        slot.identifier.description = description

        original = getattr(slot.code.parameters, 'value', slot.code.parameters)
        root = ET.fromstring(original) if str(original).strip() else ET.Element('model')
        if self._pt_region_plan is not None:
            root.append(region_validity_xml(self._pt_region_plan, item['region']))
        scope = ET.SubElement(root, 'pt_solver_occurrence', {
            'region': item['region'], 'model': identity, 'category': category,
            'selected_models': ','.join(n.upper() for n in names),
            'cdf_model_output_available': str(bool(populated)).lower(),
            'output_scope': 'native_whole_CDF_profile_not_region_clipped' if populated else 'settings_only',
            'output_attribution': ('assembled_result_not_attributed_to_a_selected_model' if assembled and populated
                                   else 'native_model_channels' if populated else 'selection_only_no_output'),
        })
        for selection in selected:
            node = ET.SubElement(scope, 'selected_model', {'name': selection['model'].upper()})
            for start, end in selection['intervals']:
                ET.SubElement(node, 'interval', interval_attributes(start, end))
        slot.code.parameters = ET.tostring(root, encoding='unicode')


if __name__ == '__main__':
    print('Please use: python convert.py')
