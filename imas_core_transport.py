#!/usr/bin/python
"""Populate the IMAS core_transport IDS from observed and modeled TRANSP transport quantities."""

import numpy as np
from nc2imas import *


class imas_core_transport:
    """Map TRANSP effective/observed transport and predictive model fluxes."""

    def __init__(self, ncimas:nc2imas):
        self.im=ncimas
        self.logger=ncimas.logger
        self.time=np.asarray(self.im.get_ncvar('TIME3'),dtype=np.float64)
        self.x=np.asarray(self.im.get_ncvar('X'),dtype=np.float64)
        self.xb=np.asarray(self.im.get_ncvar('XB'),dtype=np.float64)
        self.dvol=self.im.get_ncvar('DVOL',self.im.vconv)
        self.surf=self.im.get_ncvar('SURF',self.im.aconv)
        self.nt=len(self.time)
        self.nx=self.xb.shape[1]
        self.create_ids()

    def _get(self,name,mult=None):
        a=self.im.get_ncvar(name)
        if not self.im.is_valid(a): return np.array([])
        a=np.asarray(a,dtype=np.float64)
        return a if mult is None else a*mult

    def _first_available(self,names,mult=None):
        for name in names:
            if self.im.check_var_name(name):
                return name,self._get(name,mult)
        return '',np.array([])

    def _flux_from_divergence(self,name,mult,physical_sign=1.0):
        """Convert a TRANSP flux divergence on X to average flux on XB.

        TRANSP stores EPTR_*, EETR_*, IETR_* and related quantities as
        divergences per unit volume.  The outward surface-averaged flux is the
        enclosed volume integral divided by the flux-surface area SURF.
        """
        div=self._get(name,mult)
        if not self.im.is_valid(div,self.dvol,self.surf): return np.array([])
        integrated=np.cumsum(div*self.dvol,axis=1)
        return physical_sign*np.divide(integrated,self.surf,out=np.zeros_like(integrated),where=np.abs(self.surf)>0)

    def _fill_grid(self,grid,it):
        self.im.fill_data(grid.rho_tor_norm,self.xb,ind=it,chk=self.nx)
        # Add geometry where readily available on the same XB grid.
        vol_b=self.im.cumulative_zone_quantity(self.dvol)
        area_b=self.im.cumulative_zone_quantity(self.im.get_ncvar('DAREA',self.im.aconv))
        self.im.fill_data(grid.volume,vol_b,ind=it,chk=self.nx)
        self.im.fill_data(grid.area,area_b,ind=it,chk=self.nx)
        self.im.fill_data(grid.surface,self.surf,ind=it,chk=self.nx)

    def _fill_ion_identity(self,node,sp):
        atom=atom_base(sp)
        node.element.resize(1)
        node.element[0].a=float(atom.a)
        self.im.set_nuclear_charge(node.element[0],atom.z_n)
        node.element[0].atoms_n=int(atom.atoms_n)
        node.z_ion=float(atom.z_n)
        self.im.set_species_name(node,atom.label)

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
            if self.im.check_var_name(dn) or self.im.check_var_name(vn) or self.im.check_var_name(fn):
                specs.append((sp,self._get(dn,self.im.aconv),self._get(vn,self.im.cm2m),
                              self._flux_from_divergence(fn,self.im.vol_inv),fn))
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
            if name:
                specs.append((sp,self._flux_from_divergence(name,self.im.vol_inv),name))
        return specs

    def _prepare_model(self,model,name,index,description,comment):
        model.identifier.name=name
        model.identifier.index=index
        model.identifier.description=description
        model.comment=comment
        model.flux_multiplier=0.0
        if hasattr(model,'code'):
            model.code.name='TRANSP'
            if hasattr(model.code,'description'):
                model.code.description=description
        model.profiles_1d.resize(self.nt)
        for it in range(self.nt):
            p=model.profiles_1d[it]
            p.time=float(self.time[it])
            self._fill_grid(p.grid_d,it)
            self._fill_grid(p.grid_v,it)
            self._fill_grid(p.grid_flux,it)

    def _fill_combined_model(self,model):
        # Effective coefficients used in the TRANSP balance analysis.
        de=self._get('DIFFE',self.im.aconv)
        ve=self._get('VELE',self.im.cm2m)
        che=self._get('CONDE',self.im.aconv)
        vhe=self._get('VEL_TE',self.im.cm2m)
        chi=self._get('CONDI',self.im.aconv)
        vhi=self._get('VEL_TI',self.im.cm2m)
        chphi=self._get('CHPHI',self.im.aconv)
        vmom=self._get('VMO_NET',self.im.cm2m)

        epart_flux=self._flux_from_divergence('EPTR_OBS',self.im.vol_inv)
        een_flux=self._flux_from_divergence('EETR_OBS',self.im.vol_inv)
        ien_flux=self._flux_from_divergence('IETR_OBS',self.im.vol_inv)
        mom_flux=self._flux_from_divergence('AMTR_OBS',self.im.vol_inv,physical_sign=self.im.ip_sign)

        # ETA_USE is ohm*cm on X; conductivity in core_transport is on grid_d (XB).
        eta=self._get('ETA_USE',0.01)
        sigma=np.array([])
        if eta.size:
            eta_b=self.im.interp_radial_profiles(eta,self.x,self.xb)
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
            self.im.fill_data(p.conductivity_parallel,sigma,ind=it,chk=self.nx)

            self.im.fill_data(p.electrons.particles.d,de,ind=it,chk=self.nx)
            self.im.fill_data(p.electrons.particles.v,ve,ind=it,chk=self.nx)
            self.im.fill_data(p.electrons.particles.flux,epart_flux,ind=it,chk=self.nx)
            self.im.fill_data(p.electrons.energy.d,che,ind=it,chk=self.nx)
            self.im.fill_data(p.electrons.energy.v,vhe,ind=it,chk=self.nx)
            self.im.fill_data(p.electrons.energy.flux,een_flux,ind=it,chk=self.nx)
            self.im.fill_data(p.total_ion_energy.d,chi,ind=it,chk=self.nx)
            self.im.fill_data(p.total_ion_energy.v,vhi,ind=it,chk=self.nx)
            self.im.fill_data(p.total_ion_energy.flux,ien_flux,ind=it,chk=self.nx)

            mom=p.momentum_phi if self.im.dd_at_least(3,42,0) else p.momentum_tor
            self.im.fill_data(mom.d,chphi,ind=it,chk=self.nx)
            self.im.fill_data(mom.v,vmom,ind=it,chk=self.nx)
            self.im.fill_data(mom.flux,mom_flux,ind=it,chk=self.nx)

            p.ion.resize(len(species_specs))
            for ji,(sp,d,v,flux,_fn) in enumerate(species_specs):
                inode=p.ion[ji]
                self._fill_ion_identity(inode,sp)
                self.im.fill_data(inode.particles.d,d,ind=it,chk=self.nx)
                self.im.fill_data(inode.particles.v,v,ind=it,chk=self.nx)
                self.im.fill_data(inode.particles.flux,flux,ind=it,chk=self.nx)
        return len(species_specs)

    def _fill_transport_solver_model(self,model):
        """Write predictive/model fluxes separately from observed balance data."""
        epart_flux=self._flux_from_divergence('EPTR_MOD',self.im.vol_inv)
        een_flux=self._flux_from_divergence('EETR_MOD',self.im.vol_inv)
        ien_flux=self._flux_from_divergence('IETR_MOD',self.im.vol_inv)
        species_specs=self._species_model()

        present=[]
        for name,arr in [('EPTR_MOD',epart_flux),('EETR_MOD',een_flux),('IETR_MOD',ien_flux)]:
            if arr.size: present.append(name)
        present.extend(name for _sp,_flux,name in species_specs)
        has_iptr=self.im.check_var_name('IPTR_MOD')
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
            'Fluxes are reconstructed from TRANSP *_MOD divergence profiles by cumulative volume '
            'integration divided by SURF. No effective D/v coefficients are inferred from these '
            'modeled divergences.'+iptr_note
        )

        for it in range(self.nt):
            p=model.profiles_1d[it]
            self.im.fill_data(p.electrons.particles.flux,epart_flux,ind=it,chk=self.nx)
            self.im.fill_data(p.electrons.energy.flux,een_flux,ind=it,chk=self.nx)
            self.im.fill_data(p.total_ion_energy.flux,ien_flux,ind=it,chk=self.nx)
            p.ion.resize(len(species_specs))
            for ji,(sp,flux,_name) in enumerate(species_specs):
                inode=p.ion[ji]
                self._fill_ion_identity(inode,sp)
                self.im.fill_data(inode.particles.flux,flux,ind=it,chk=self.nx)

        return True,present

    def create_ids(self):
        ct=self.im.new_ids('core_transport')
        ct=self.im.fill_properties(ct)
        self.im.fill_code_properties(ct)
        ct.ids_properties.homogeneous_time=1
        ct.time=np.asarray(self.time,dtype=np.float64)
        self.im.fill_vac_tor_field(ct)

        # Allocate both model slots, then shrink if no predictive/model channels exist.
        ct.model.resize(2)
        nobs=self._fill_combined_model(ct.model[0])
        has_solver,present=self._fill_transport_solver_model(ct.model[1])
        if not has_solver:
            ct.model.resize(1)

        if self.im.check_var_name('IPTR_MOD'):
            self.logger.info('core_transport: IPTR_MOD present but aggregate total-ion particle flux has no IMAS core_transport node')
        self.logger.info(f'core_transport: combined model, {nobs} observed ion particle species')
        if has_solver:
            self.logger.info('core_transport: transport_solver model channels: '+', '.join(present))
        self.im.putdb_imas(ct)


if __name__ == '__main__':
    print('Please use: python convert.py')
