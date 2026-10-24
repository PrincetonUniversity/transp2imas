#!/usr/bin/python
"""Populate the IMAS core_sources IDS from TRANSP heating, particle, current, and torque sources."""

import numpy as np
from nc2imas import *
from writer_base import IdsWriter


class imas_core_sources(IdsWriter):
    """Translate high-confidence TRANSP thermal source terms to core_sources."""

    SOURCE_IDS = {
        'nbi': (2, 'Source from Neutral Beam Injection'),
        'ec': (3, 'Sources from electron cyclotron heating and current drive'),
        'lh': (4, 'Sources from lower hybrid heating and current drive'),
        'ic': (5, 'Sources from heating at the ion cyclotron range of frequencies'),
        'fusion': (6, 'Sources from fusion reactions'),
        'ohmic': (7, 'Source from ohmic heating'),
        'bremsstrahlung': (8, 'Bremsstrahlung radiation loss'),
        'synchrotron_radiation': (9, 'Cyclotron/synchrotron radiation loss'),
        'line_radiation': (10, 'Line radiation loss'),
        'collisional_equipartition': (11, 'Collisional electron-ion equipartition'),
        'cold_neutrals': (12, 'Sources and ionization work from wall/volume neutrals'),
        'bootstrap_current': (13, 'Bootstrap current source'),
    }

    def _load_grid(self):
        """Read the profile grids the sources are written on (done in ``build``, not at construction)."""
        self.time = np.asarray(self.cdf.get_ncvar('TIME3'), dtype=np.float64)
        self.x = self.cdf.get_ncvar('X')
        self.dvol = self.cdf.get_ncvar('DVOL', self.grid.vconv)
        self.nt = len(self.time)
        self.nx = self.x.shape[1]

    def _has(self, *names):
        return any(self.cdf.check_var_name(n) for n in names)

    def _get(self, name, mult=None):
        a = self.cdf.get_ncvar(name)
        if not self.fill.is_valid(a):
            return np.array([])
        a = np.asarray(a, dtype=np.float64)
        a = self.fill.check_profile_shape(name, a, self.nt, self.nx)
        if not a.size:
            return a
        return a if mult is None else a * mult

    @staticmethod
    def _sum_valid(*arrays):
        good = [np.asarray(a, dtype=np.float64) for a in arrays if isinstance(a, np.ndarray) and a.size]
        if not good:
            return np.array([])
        out = np.zeros_like(good[0], dtype=np.float64)
        for a in good:
            out += a
        return out

    def _inside(self, profile):
        if not self.fill.is_valid(profile, self.dvol):
            return np.array([])
        return np.cumsum(np.asarray(profile) * self.dvol, axis=1)

    def _integral(self, profile):
        inside = self._inside(profile)
        return inside[:, -1] if inside.size else np.array([])

    def _parallel_from_toroidal_current(self, variable):
        """Convert a TRANSP toroidal current profile to ``<J.B>/B0``.

        ``variable`` is in A/cm^2. ``CONPLJB`` approximates ``<J.B>/Jtor``
        and IMAS normalizes by the IDS's own ``b0`` (``<J.B>/B0``, B0 =
        ``vacuum_toroidal_field/b0``), not TRANSP's ``BZ``. The plasma current sign
        maps the resulting profile to the physical IMAS +phi convention.
        """
        current = self._get(variable, self.grid.area_inv)
        conpljb = self._get('CONPLJB')
        b0 = self.signs.vacuum_b0_magnitude(self.time)
        if b0 is None:
            self.ctx.mapping.skipped('vacuum toroidal field (RAXIS/BZXR) is not in the CDF',
                                    source=variable, target='source/profiles_1d/j_parallel')
            return np.array([])
        if not self.fill.is_valid(current, conpljb):
            return np.array([])
        denominator = b0[:, None]
        geometry = np.divide(
            conpljb,
            denominator,
            out=np.zeros_like(conpljb, dtype=np.float64),
            where=np.isfinite(denominator) & (np.abs(denominator) > 0.0),
        )
        return self.signs.ip_sign * current * geometry

    def _set_identifier(self, source, name):
        idx, desc = self.SOURCE_IDS[name]
        source.identifier.name = name
        source.identifier.index = idx
        source.identifier.description = desc
        source.species.type.name = 'unspecified'
        source.species.type.index = 0
        source.species.type.description = 'N/A'

    def _set_source_species_neutral(self, source, sp):
        atom = atom_base(sp)
        source.species.type.name = 'neutral'
        source.species.type.index = 4
        source.species.type.description = 'Neutral species in a single/average state'
        n = source.species.neutral
        n.element.resize(1)
        n.element[0].a = float(atom.a)
        self.fill.set_nuclear_charge(n.element[0], atom.z_n)
        n.element[0].atoms_n = int(atom.atoms_n)
        self.fill.set_species_name(n, atom.label.replace('_BEAM',''))

    def _fill_ion_identity(self, node, sp, fusion=False):
        key = sp
        if fusion:
            key = {'P':'P_FUSION','P14':'P14_FUSION','T':'T_FUSION',
                   'HE3':'HE3_FUSION','HE4':'HE4_FUSION'}.get(sp, sp)
        atom = atom_base(key)
        node.element.resize(1)
        node.element[0].a = float(atom.a)
        self.fill.set_nuclear_charge(node.element[0], atom.z_n)
        node.element[0].atoms_n = int(atom.atoms_n)
        node.z_ion = float(atom.z_n)
        self.fill.set_species_name(node, atom.label.replace('_BEAM',''))

    def _prepare_source(self, source, name):
        self._set_identifier(source, name)
        source.profiles_1d.resize(self.nt)
        source.global_quantities.resize(self.nt)
        self.fill.fill_profiles1d_grid(source)
        for it in range(self.nt):
            source.profiles_1d[it].time = float(self.time[it])
            source.global_quantities[it].time = float(self.time[it])
        if hasattr(source, 'code'):
            source.code.name = 'TRANSP'
            if hasattr(source.code, 'description'):
                source.code.description = self.SOURCE_IDS[name][1]

    def _fill_channel(self, source, electron_energy=np.array([]), electron_particles=np.array([]),
                      ion_energy=np.array([]), ion_particles=np.array([]),
                      momentum=np.array([]), momentum_jxb=np.array([]),
                      j_parallel=np.array([]), ions=None):
        ions = ions or []
        e_pow_in = self._inside(electron_energy)
        e_part_in = self._inside(electron_particles)
        i_pow_in = self._inside(ion_energy)
        i_part_in = self._inside(ion_particles)
        torq_in = self._inside(momentum)

        e_pow_g = self._integral(electron_energy)
        e_part_g = self._integral(electron_particles)
        i_pow_g = self._integral(ion_energy)
        i_part_g = self._integral(ion_particles)
        torque_g = self._integral(momentum)
        total_pow_g = self._sum_valid(e_pow_g, i_pow_g)
        # Enclosed particles and power of each ion species do not depend on the time slice: computed once.
        ion_inside = [(self._inside(part), self._inside(en)) for _sp, part, en, _fusion in ions]

        for it in range(self.nt):
            p = source.profiles_1d[it]
            g = source.global_quantities[it]
            self.fill.fill_data(p.electrons.energy, electron_energy, ind=it, chk=self.nx)
            if hasattr(p.electrons, 'power_inside'):
                self.fill.fill_data(p.electrons.power_inside, e_pow_in, ind=it, chk=self.nx)
            self.fill.fill_data(p.electrons.particles, electron_particles, ind=it, chk=self.nx)
            if hasattr(p.electrons, 'particles_inside'):
                self.fill.fill_data(p.electrons.particles_inside, e_part_in, ind=it, chk=self.nx)
            self.fill.fill_data(p.total_ion_energy, ion_energy, ind=it, chk=self.nx)
            if hasattr(p, 'total_ion_power_inside'):
                self.fill.fill_data(p.total_ion_power_inside, i_pow_in, ind=it, chk=self.nx)
            if hasattr(p, 'total_ion_particles'):
                self.fill.fill_data(p.total_ion_particles, ion_particles, ind=it, chk=self.nx)
            if hasattr(p, 'total_ion_particles_inside'):
                self.fill.fill_data(p.total_ion_particles_inside, i_part_in, ind=it, chk=self.nx)

            mom_node = p.momentum_phi if self.dd.dd_at_least(3,42,0) else p.momentum_tor
            self.fill.fill_data(mom_node, momentum, ind=it, chk=self.nx)
            # DD4 / late DD3 spell it phi, DD 3.39 tor.
            jxb_node = next((getattr(p, name) for name in ('momentum_phi_j_cross_b_field',
                                                           'momentum_tor_j_cross_b_field')
                             if hasattr(p, name)), None)
            if jxb_node is not None:
                self.fill.fill_data(jxb_node, momentum_jxb, ind=it, chk=self.nx)
            torq_node = p.torque_phi_inside if self.dd.dd_at_least(3,42,0) else p.torque_tor_inside
            self.fill.fill_data(torq_node, torq_in, ind=it, chk=self.nx)
            self.fill.fill_data(p.j_parallel, j_parallel, ind=it, chk=self.nx)

            p.ion.resize(len(ions))
            for ji, iondat in enumerate(ions):
                sp, part, en, fusion = iondat
                inode = p.ion[ji]
                self._fill_ion_identity(inode, sp, fusion=fusion)
                self.fill.fill_data(inode.particles, part, ind=it, chk=self.nx)
                if hasattr(inode, 'particles_inside'):
                    self.fill.fill_data(inode.particles_inside, ion_inside[ji][0], ind=it, chk=self.nx)
                self.fill.fill_data(inode.energy, en, ind=it, chk=self.nx)
                if hasattr(inode, 'power_inside'):
                    self.fill.fill_data(inode.power_inside, ion_inside[ji][1], ind=it, chk=self.nx)

            if total_pow_g.size: g.power = float(total_pow_g[it])
            if i_pow_g.size: g.total_ion_power = float(i_pow_g[it])
            if i_part_g.size: g.total_ion_particles = float(i_part_g[it])
            if e_pow_g.size: g.electrons.power = float(e_pow_g[it])
            if e_part_g.size: g.electrons.particles = float(e_part_g[it])
            if torque_g.size:
                if self.dd.dd_at_least(3,42,0): g.torque_phi = float(torque_g[it])
                else: g.torque_tor = float(torque_g[it])

    # TRANSP's total input torque (xplbal.f90, TQINP -> CDF TQIN) is
    #   TQJXB + FTQRIP*TQRPL + TQBE + TQBI + TQBTH + TQDIFB + TQANOM + TQNTV,
    # of which TQBE/TQBI (collisional torque on electrons/ions), TQBTH (thermalization) and
    # TQJXB (return-current JxB) are the beam terms; TQDIFB is the fast-ion anomalous-diffusion
    # torque, written to the CDF as TQABDF.  The sum of the first five reproduces TQIN to
    # single precision in every bundled case with beams (ITER needs TQABDF); TQANOM and TQNTV
    # are not beam torques.  The IMAS JxB field is a subset already included in the total.
    NBI_REQUIRED_TORQUES = ('TQBE', 'TQBI', 'TQBTH', 'TQJXB')
    NBI_OPTIONAL_TORQUES = ('TQABDF',)

    def _nbi_momentum(self):
        """``(total, jxb)`` NBI toroidal momentum source on the thermal plasma (IMAS sign), SI.

        Both are empty unless every required beam torque is in the CDF: a partial sum
        must not be written as the total.
        """
        terms = {name: self._get(name, self.grid.vol_inv)
                 for name in self.NBI_REQUIRED_TORQUES + self.NBI_OPTIONAL_TORQUES}
        missing = [name for name in self.NBI_REQUIRED_TORQUES if not terms[name].size]
        if missing:
            self.logger.warning('core_sources: nbi momentum source not written; missing '
                                + ', '.join(missing) + ' (a partial sum would not be the total torque)')
            return np.array([]), np.array([])
        ripple = self._get('TQRPL', self.grid.vol_inv)
        if ripple.size and np.any(ripple):
            self.logger.warning('core_sources: TQRPL (beam ripple JxB torque) is nonzero but its '
                                'weight FTQRIP is not in the CDF; it is excluded from the nbi momentum')
        total = self._sum_valid(*terms.values())
        sign = self.signs.ip_sign
        return sign * total, sign * terms['TQJXB']

    def _build_nbi(self, source):
        if not self._has('PBE','PBI','PBTH','SBTH','SBE'):
            return False
        self._prepare_source(source, 'nbi')
        ee = self._get('PBE', self.grid.vol_inv)
        ep = self._get('SBE', self.grid.vol_inv)
        ie = self._sum_valid(self._get('PBI', self.grid.vol_inv), self._get('PBTH', self.grid.vol_inv))
        ip = self._get('SBTH', self.grid.vol_inv)
        mom, jxb = self._nbi_momentum()
        ions=[]
        species=self.cdf.get_species_names('SBTH_?', cut='SBTH_')
        for sp in species:
            part=self._get('SBTH_'+sp.replace('HE3','3').replace('HE4','4'), self.grid.vol_inv)
            en=self._get('PBTH_'+sp.replace('HE3','3').replace('HE4','4'), self.grid.vol_inv)
            if part.size or en.size: ions.append((sp,part,en,False))
        if len(species)==1:
            try: self._set_source_species_neutral(source, species[0])
            except KeyError:                    # atom_base: unknown species token
                self.ctx.mapping.skipped(f'unknown beam species {species[0]!r}', target='source/species')
        jp = self._parallel_from_toroidal_current('CURB')
        self._fill_channel(source, electron_energy=ee, electron_particles=ep,
                           ion_energy=ie, ion_particles=ip, momentum=mom,
                           momentum_jxb=jxb, j_parallel=jp, ions=ions)
        return True

    def _build_wave(self, source, name, evar, ivar=None):
        if not self._has(evar, ivar or '__none__'):
            return False
        self._prepare_source(source, name)
        ee=self._get(evar,self.grid.vol_inv)
        ie=self._get(ivar,self.grid.vol_inv) if ivar else np.array([])
        self._fill_channel(source,electron_energy=ee,ion_energy=ie)
        return True

    def _build_ohmic(self, source):
        if not self._has('POH','PLJBOH'):
            return False
        self._prepare_source(source,'ohmic')
        ee=self._get('POH',self.grid.vol_inv)
        jb=self._get('PLJBOH',self.grid.area_inv)
        b0=self.signs.vacuum_b0_magnitude(self.time)       # <J.B>/B0 with the IDS's own b0
        jp=np.array([])
        if jb.size and b0 is not None:
            jp=self.signs.jparallel_sign*np.divide(jb,b0[:,None],out=np.zeros_like(jb),where=np.abs(b0[:,None])>0)
        self._fill_channel(source,electron_energy=ee,j_parallel=jp)
        return True

    def _build_radiation(self, source, name, var):
        if not self._has(var): return False
        self._prepare_source(source,name)
        ee=-self._get(var,self.grid.vol_inv)
        self._fill_channel(source,electron_energy=ee)
        return True

    def _build_equip(self, source):
        if not self._has('QIE'): return False
        self._prepare_source(source,'collisional_equipartition')
        q=self._get('QIE',self.grid.vol_inv)
        self._fill_channel(source,electron_energy=-q,ion_energy=q)
        return True

    def _build_neutrals(self, source):
        if not self._has('SCEW','SCEV','PION'): return False
        self._prepare_source(source,'cold_neutrals')
        ep=self._sum_valid(self._get('SCEW',self.grid.vol_inv),self._get('SCEV',self.grid.vol_inv))
        ee=-self._get('PION',self.grid.vol_inv)
        self._fill_channel(source,electron_energy=ee,electron_particles=ep)
        return True

    def _build_bootstrap(self, source):
        if not self._has('CURBS','CONPLJB'): return False
        self._prepare_source(source,'bootstrap_current')
        jp = self._parallel_from_toroidal_current('CURBS')
        self._fill_channel(source,j_parallel=jp)
        return True

    def _build_fusion(self, source):
        species=self.cdf.get_species_names('SFTH_?', cut='SFTH_')
        if not species: return False
        self._prepare_source(source,'fusion')
        eprof=[]; iprof=[]; pprof=[]; ions=[]
        for sp in species:
            raw=sp.replace('HE3','3').replace('HE4','4')
            sf=self._get('SFTH_'+raw,self.grid.vol_inv)
            pfth=self._get('PFTH_'+raw,self.grid.vol_inv)
            pfe=self._get('PFE_'+raw,self.grid.vol_inv)
            pfi=self._get('PFI_'+raw,self.grid.vol_inv)
            if sf.size: pprof.append(sf)
            if pfe.size: eprof.append(pfe)
            if pfi.size: iprof.append(pfi)
            if pfth.size: iprof.append(pfth)
            if sf.size or pfth.size: ions.append((sp,sf,pfth,True))
        self._fill_channel(source,electron_energy=self._sum_valid(*eprof),
                           ion_energy=self._sum_valid(*iprof),
                           ion_particles=self._sum_valid(*pprof),ions=ions)
        return True

    def build(self):
        self._load_grid()
        cs=self.dd.new_ids('core_sources')
        cs=self.fill.fill_properties(cs)
        self.fill.fill_code_properties(cs)
        cs.ids_properties.homogeneous_time=1
        cs.time=np.asarray(self.time,dtype=np.float64)
        self.fill.fill_vac_tor_field(cs)

        builders=[
            self._build_nbi,
            lambda s:self._build_wave(s,'ec','PEECH'),
            lambda s:self._build_wave(s,'lh','PELH','PILH'),
            lambda s:self._build_wave(s,'ic','PEICH','PIICH'),
            self._build_fusion,
            self._build_ohmic,
            lambda s:self._build_radiation(s,'bremsstrahlung','PRAD_BR'),
            lambda s:self._build_radiation(s,'synchrotron_radiation','PRAD_CY'),
            lambda s:self._build_radiation(s,'line_radiation','PRAD_LI'),
            self._build_equip,
            self._build_neutrals,
            self._build_bootstrap,
        ]
        # IMAS-Python resize() recreates every AoS element, so resizing a
        # populated source array would erase it. Build each candidate in a
        # temporary IDS and append only successful sources to the output IDS.
        used=0
        for fill_source in builders:
            candidate = self.dd.new_ids('core_sources')
            candidate.source.resize(1)
            if fill_source(candidate.source[0]):
                cs.source.append(candidate.source[0])
                used += 1
        self.logger.info(f'core_sources: built {used} source processes')
        return [cs]


if __name__ == '__main__':
    print('Please use: python convert.py')
