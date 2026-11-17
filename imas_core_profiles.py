#!/usr/bin/python
"""Populate the IMAS core_profiles IDS from TRANSP profiles and metadata."""

from typing import List
import numpy as np
from nc2imas import *
from neutral_profiles_mapping import (
    discover_neutral_channels, validate_profile, attach_rotation_xml,
)
from writer_base import IdsWriter

class ion_states:
    def __init__(self, sdens, label=''):
        #arrays for all times
        self.sdens = sdens
        self.label = label

class ion_data(atom_base):
    def __init__(self, key):
        super().__init__(key)
        #arrays for times
        self.dens = []
        self.dens_fast = []
        self.dens_th = []
        self.temp = []
        self.vel_rad = []
        self.vel_tor = []
        self.vel_pol = []
        self.vel_par = []
        self.z_ion_1d = []
        self.rot_freq_tor = []
        self.aas = []
        self.zns = []
        self.zions = []
        #same as z_n
        self.z_ion = self.z_n
        self.states = []


    def is_valid(self, logger=None):
        check1 = self.a>0 and self.z_n>0 and self.z_ion>0
        check2 = len(self.aas)>0 and len(self.zns)>0 and len(self.zions)>0
        check3 = self.atoms_n>0 and len(self.dens)>0
        valid = check3 and (check1 or check2)
        if logger:
            if not valid:
               logger.info(f'\t invalid ion: dens {len(self.dens)} '
                           f' {self.label} : a {self.a} z_n {self.z_n} z_ion {self.z_ion}')
        return valid

class neutrals_data(atom_base):
    """Separate wall, volume or recombination neutral population on TIME3/X."""

    def __init__(self, symb, label, origin):
        super().__init__(symb)
        self.label = label
        self.origin = origin
        self.dens = np.array([])
        self.temp = np.array([])
        self.omega = np.array([])
        self.rotation_source = ''
        self.rotation_sign = 1.0
        self.density_source = ''
        self.temperature_source = ''

    def is_valid(self, logger=None):
        valid = self.a > 0 and self.z_n > 0 and self.atoms_n > 0 \
            and isinstance(self.dens, np.ndarray) and self.dens.ndim == 2 and self.dens.size > 0
        if not valid and logger is not None:
            logger.info('invalid neutral: %s (a=%g, z_n=%g)', self.label, self.a, self.z_n)
        return valid


class imas_core_profiles(IdsWriter):
    #total ions
    def get_thermal_ion_data(self, atomsDef) -> List:
        ionDict = atomsDef.ions
        ions = []
        #RMJSYM(TIME3, RMJSYM)
        rmaj_ = self.cdf.get_ncvar('RMJSYM')
        rmaj = np.average(rmaj_) if self.fill.is_valid(rmaj_) else 0
        self.logger.info(f' ions::')
        for short, key in ionDict.items():
            symb = 'HE' + short if short in ['3', '4'] else short
            # only thermal; symb needed
            dens = self.get_ion_density('thermal', ext=symb)
            if len(dens) == 0:
                continue
            ion = ion_data(key)
            self.logger.info(f'   ion: {symb}=> {ion.label}:')
            ion.dens_th = dens
            self.logger.info(f'\t{symb}: N{symb} {ion.dens_th.shape}')
            ion.dens =  self.get_total_ion_density(short, symb)
            ion.temp = self.cdf.get_ncvar(f'TI')
            self.logger.info(f'\t{symb}: TI {ion.temp.shape}')
            ion.vel_rad = self.cdf.get_ncvar(f'VNDNC_{symb}', 0.01)
            self.logger.info(f'\t{symb} vel_rad {ion.vel_rad.shape}')
            ion.vel_pol = self.cdf.get_ncvar(f'VPOL{short}_NC', 0.01)
            self.logger.info(f'\t{symb} vel_pol {ion.vel_pol.shape}')
            ion.vel_tor = self.signs.toroidal_rotation_sign * self.cdf.get_ncvar(f'VTOR{short}_NC', 0.01)
            self.logger.info(f'\t{symb} vel_tor {ion.vel_tor.shape}')
            #ion.vel_par = [] #UPLLHATI/(BZXR/RAXIS)*1.0D-2
            #ion.z_ion_1d = []
            # RAD/SEC
            ion.rot_freq_tor = self.signs.toroidal_rotation_sign * self.cdf.get_ncvar(f'OMEGA') #CMG_H CMG_D

            #if self.fill.is_valid(ang_vel) and rmaj>0:
            #    ion.rot_freq_tor = ang_vel / rmaj
            if ion.is_valid(self.logger):
               ions.append(ion)
        return ions

    #impurities
    def get_imp_ion_data(self) -> List:
        #NIMPS_ total for all charge states, per species
        self.logger.info(' impurities:')
        ions = []
        imp_atoms = self.cdf.get_species_names('NIMPS_?', cut='NIMPS_')
        #multiple species if not TOK
        for symb in imp_atoms:
            #total
            dens = 'NIMP' if symb=='TOK' else f'NIMPS_{symb}'
            # (TIME3, X) N/CM3 ;TOTAL IMPURITY DENSITY
            den = self.cdf.get_ncvar(dens, scale=1.0e6)
            if len(den)==0:
               continue
            ion = ion_data(symb)
            self.logger.info(f'  {symb} => {ion.label}:')
            ion.dens = den
            ion.dens_th = den
            #temperature
            ion.temp = self.profile('TX')
            #charge states. Only if NIMP_symb_n is present
            dstr = f'NIMP_{symb}_'
            dnames, ncs = self.cdf.get_species_data(f'{dstr}?', cut=dstr)

            #ZIMPS_TOK or ZIMPS_<symb>(TIME3, X) <symb> Avg. Mult. Impurity Z
            zprof = self.cdf.get_ncvar(f'ZIMPS_{symb}')
            ion.z_ion_1d = zprof
            #a,z_n from definition in atoms_class, except for TOK

            #OMGX(TIME3, X)
            ion.rot_freq_tor = self.signs.toroidal_rotation_sign * self.cdf.get_ncvar('OMGX')

            ai, zn = 0,0
            if symb=='TOK':
               #a   AIMP(TIME)
               ion.aas = self.cdf.get_ncvar(f'AIMP', interp=True)
               #z_n  XZIMP(TIME)
               ion.zns = self.cdf.get_ncvar(f'XZIMP', interp=True)
               if len(ion.aas)>0:
                  ai = ion.aas[0]
               if len(ion.zns)>0:
                  zn = ion.zns[0]
               #z_ion
               ion.zions = ion.zns #ion charge, fully ionized ?

            self.logger.info(f'\t{symb}: {ai}, {zn}, {dens} {ion.dens.shape}')
            self.logger.info(f'\t{symb}: charges: {ncs} NIMP_ {dnames} ZIMPS_{symb} {zprof.shape}')
            #charge states. eg:NIMP_BE_4
            for dname, ncharge in zip(dnames, ncs):
               self.logger.info(f'\t{symb}: {dname} ch: {ncharge}')
               sdens = self.cdf.get_ncvar(dname, scale=1.0e6)
               istate = ion_states(sdens, f'{symb}+{ncharge}')
               ion.states.append(istate)
            if ion.is_valid(self.logger):
               ions.append(ion)
        return ions

    def get_ion_density(self, cat, ext, name_only=False):
        if cat == 'thermal':
           dstr = 'N'
        elif cat == 'beam':
           dstr = 'BDENS_'
        elif cat == 'fusion':
           dstr = 'FDENS_'
        elif cat == 'mini':
           dstr = 'NMINI_'
        if name_only is True:
           return dstr
        return self.cdf.get_ncvar(f'{dstr}{ext}', 1.0e6)

    def get_total_ion_density(self, short, symb):
        tdens = self.get_ion_density('thermal', symb)
        bdens = self.get_ion_density('beam', short)
        fdens = self.get_ion_density('fusion', short)
        mdens = self.get_ion_density('mini', short)
        tot = self.fill.sum_arrays(tdens, bdens, fdens, mdens)
        return tot

    def get_ion_data_type(self, atomsDef, cat) -> List:
        #cat= 'beam', 'fusion', 'mini'
        #BDENS_{T,D,H,3,4}
        ionDict = ''
        dstr, tstr, freqs = '', '', ''
        if cat == 'beam': #_{T,D,H,3,4}
            ionDict = atomsDef.beams
            tstr = 'EBEAM'
            freqs = 'OMEGB'
        elif cat == 'fusion': #_{P,P14,T,3,4}
            ionDict = atomsDef.fusion_ions
            tstr = 'EFUSN'
            freqs = 'OMEGF'
        elif cat == 'mini': #_{T,D,H,3,4}
            ionDict = atomsDef.minority_ions
            tstr = 'TMINI'
            freqs = 'INVALID'
        ions = []
        self.logger.info(f' {cat}-ions::')
        for short, key in ionDict.items():
            symb = 'HE' + short if short in ['3', '4'] else short
            # no thermal called from here; short is needed
            dens = self.get_ion_density(cat, short)
            dstr = self.get_ion_density(cat, short, name_only=True)
            if len(dens) == 0:
               continue
            ion = ion_data(key)
            ion.dens_fast = dens
            ion.dens = self.get_total_ion_density(short, symb)
            ion.temp = self.cdf.get_ncvar(f'{tstr}_{short}')
            # (TIME3, X) RAD/SEC
            ion.rot_freq_tor = self.signs.toroidal_rotation_sign * self.cdf.get_ncvar(f'{freqs}_{short}')

            self.logger.info(f'\t {symb} => {ion.label} :')
            self.logger.info(f'\t{symb} : {dstr}{short} {ion.dens_fast.shape}')
            if ion.is_valid(self.logger):
               ions.append(ion)
        return ions


    def get_ion_data(self) -> List:
        atomsDef = atom_base('')
        iions = self.get_imp_ion_data()
        tions = self.get_thermal_ion_data(atomsDef)
        bions = self.get_ion_data_type(atomsDef, cat='beam')
        fions = self.get_ion_data_type(atomsDef, cat='fusion')
        mions = self.get_ion_data_type(atomsDef, cat='mini')
        for lst in [iions, bions, fions, mions]:
            if len(lst)>0:
               for ion in lst:
                  tions.append(ion)
        return tions

    def get_ion_data_occur(self, prof_id=0) -> List:
        atomsDef = atom_base('')
        res = self.get_thermal_ion_data(atomsDef)
        iions = self.get_imp_ion_data()
        bions = self.get_beam_ion_data(atomsDef)
        fions = []
        if prof_id==1:
            fions = self.get_fusion_ion_data(atomsDef)
        elif prof_id>1:
            raise Exception('Error: invalid core_profile occuurence')
        for ions in [iions, bions, fions]:
            if len(ions)>0:
               for ion in ions:
                  res.append(ion)
        return res


    def get_neutral_data_type(self, atomsDef, cat) -> List:
        """Load independent neutral populations with verified units and grids.

        Wall DN0W*, volume DN0V*, and recombination N0V0_* densities are
        distinct CDF populations.  Temperatures are ALREADY eV (the legacy
        1e6 temperature factor was incorrect).  OM0* are signed toroidal
        angular frequencies in rad/s, not temperatures or particle fluxes.
        """
        names = discover_neutral_channels(self.cdf.vars_)
        result = []
        shape = np.asarray(self.cdf.get_ncvar('X')).shape
        for spec in names:
            if spec['origin'] != cat:
                continue
            symb = spec['species']
            # atom_base knows real element/isotope masses and nuclear charges,
            # even when the neutral originated in the recombination solver.
            atom = atom_base(symb)
            if atom.a <= 0 or atom.z_n <= 0:
                self.logger.warning('skipping unknown TRANSP neutral %s', symb)
                continue
            label = f'{symb}_{spec["origin"].upper()}'
            neut = neutrals_data(symb, label, spec['origin'])
            neut.density_source = spec['density']
            neut.dens = validate_profile(self.cdf.nc, spec['density'], 'density', shape)*1.0e6
            if spec['temperature']:
                neut.temperature_source = spec['temperature']
                neut.temp = validate_profile(
                    self.cdf.nc, spec['temperature'], 'temperature', shape)
            if spec['rotation']:
                neut.rotation_source = spec['rotation']
                neut.rotation_sign = float(self.signs.toroidal_rotation_sign)
                neut.omega = neut.rotation_sign * validate_profile(
                    self.cdf.nc, spec['rotation'], 'rotation', shape)
            if neut.is_valid(self.logger):
                result.append(neut)
        return result

    def get_neutral_data(self):
        atomsDef = atom_base('')
        return [n for origin in ('wall', 'volume', 'recombination')
                for n in self.get_neutral_data_type(atomsDef, origin)]

    def _neutral_midplane_radius(self):
        """Outboard R (m) on neutral X grid, if provided by TRANSP.

        Needed only for the pre-DD4 neutral toroidal *velocity* field:
        v_phi = omega_tor * R_midplane.  If no consistent R exists, do not
        invent a neutral velocity profile.
        """
        rmj = self.cdf.get_ncvar('RMJMP', self.grid.cm2m)
        xb = self.cdf.get_ncvar('XB')
        x = self.cdf.get_ncvar('X')
        if not self.fill.is_valid(rmj, xb, x):
            return np.array([])
        if rmj.shape != xb.shape or x.ndim != 2 or x.shape[0] != rmj.shape[0]:
            return np.array([])
        radius = self.grid.interp_radial_profiles(rmj, xb, x)
        return radius if radius.shape == x.shape and np.isfinite(radius).all() else np.array([])


    def build(self):
        self.logger.info('\ncore_profiles')
        cp = self.dd.new_ids('core_profiles')
        cm2m = 0.01
        aconv = 1.0e-4
        vconv = 1.0e-6
        pi = constants.pi  #scipy
        ip_sign = self.signs.ip_sign
        bt_sign = self.signs.bt_sign
        psi_sign = self.signs.psi_sign
        q_sign = self.signs.q_sign

        #scalar time  len(TIME)
        ntimes = self.grid.ntimes
        #profile time  len(TIME3)
        ntimep = self.grid.ntimep
        coord_x = self.cdf.get_ncvar('X')
        coord_xb = self.cdf.get_ncvar('XB')
        if not self.fill.is_valid(coord_x, coord_xb, name=r'coord_x/xb'):
            print('\t Error: invalid coord_x/xb', coord_x, coord_xb)
            return []

        cp = self.fill.fill_properties(cp)
        self.fill.fill_code_properties(cp)
        cp.ids_properties.homogeneous_time = 1
        cp.time = self.grid.timep

        #NOTE: nc2imas.fill_data and fill_data_nc functions do interpolation. These
        #functions handle only 1D data. If 2D data is passed or read from NC variable,
        # a valid index is required for selecting data at that index.

        #If profile time(TIME3) and scalar time(TIME) are different, interp onto profile time
        #Since homogeneous time=1, all quantities should have the same time.
        interp = (ntimes != ntimep)

        # Geometry related quantities
        # grid area. DAREA (TIME3, X) CM**2 ZONE CROSS SECTIONAL AREA
        grid_area = self.profile('DAREA', aconv)
        area_boundary = self.grid.cumulative_zone_quantity(grid_area)
        area_enclosed = self.grid.interp_radial_profiles(area_boundary, coord_xb, coord_x, axis_value=0.0)
        #DVOL(TIME3, X) CM**3  ZONE VOLUME. Cumulative values live at XB.
        gvol = self.profile('DVOL', vconv)
        volume_boundary = self.grid.cumulative_zone_quantity(gvol)
        volume_enclosed = self.grid.interp_radial_profiles(volume_boundary, coord_xb, coord_x, axis_value=0.0)
        #SURF is defined on XB; map it to the X grid used by core_profiles.
        surface_boundary = self.profile('SURF', aconv)
        surface = self.grid.interp_radial_profiles(surface_boundary, coord_xb, coord_x, axis_value=0.0)
        # PAREA(TIME) - plasma CS area[cm^2]
        parea = self.cdf.get_ncvar('PAREA', aconv)

        #vacuum toroidal field
        self.fill.fill_vac_tor_field(cp)

        #ions for profiles_1d
        ions = self.get_ion_data()
        #ions = self.get_ion_data_occur(prof_id=occur)

        # neutrals for profile_id
        neutrals = self.get_neutral_data()
        if hasattr(cp, 'code') and hasattr(cp.code, 'parameters'):
            cp.code.parameters = attach_rotation_xml(
                cp.code.parameters, neutrals, np.asarray(self.grid.timep, dtype=np.float64))
        radius_for_neutral_velocity = (self._neutral_midplane_radius()
                                       if not self.dd.dd_at_least(4, 0, 0)
                                       else np.array([]))


        #global_quantities
        #------------------------------------------------------------------
        self.fill.fill_data_nc(cp.global_quantities.ip, 'PCURC', mult=ip_sign, interp=True)

        # TRANSP CUR/CUROH/CURBS are TOROIDAL current densities [A/cm^2].
        # Keep them for toroidal/global-current quantities.  DD4 core_profiles
        # j_total/j_ohmic/j_non_inductive/j_bootstrap, however, are PARALLEL
        # current densities defined as <J.B>/B0.  Use PLJB/PLJBOH for those
        # fields and write CUR separately to j_phi.
        cur_phi = self.profile('CUR', self.grid.area_inv)
        if self.fill.is_valid(cur_phi):
            cur_phi = ip_sign * cur_phi

        curoh_phi = self.profile('CUROH', self.grid.area_inv)
        if self.fill.is_valid(curoh_phi):
            curoh_phi = ip_sign * curoh_phi

        # Global non-inductive current is explicitly the toroidal component.
        j_nonind_phi = []
        if self.fill.is_valid(cur_phi, curoh_phi, grid_area):
            j_nonind_phi = self.fill.operate_arrays(cur_phi, '-', curoh_phi)
            curt = self.fill.np_multiply_sum(j_nonind_phi, grid_area, axis=1)
            self.fill.fill_data(cp.global_quantities.current_non_inductive, curt)

        # Bootstrap global current is also a toroidal component.
        curbs_phi = self.profile('CURBS', self.grid.area_inv)
        if self.fill.is_valid(curbs_phi):
            curbs_phi = ip_sign * curbs_phi
        if len(curbs_phi) > 0:
            curbsxa = self.fill.np_multiply_sum(curbs_phi, grid_area, axis=1)
            self.fill.fill_data(cp.global_quantities.current_bootstrap, curbsxa)

        # Parallel current profiles. TRANSP PLJB* are <J.B> in A*T/cm^2.
        # IMAS defines these as <J.B>/B0 with B0 the IDS's own vacuum_toroidal_field/b0
        # (not TRANSP's BZ); dividing by |b0| and multiplying by ip_sign gives that
        # definition in the physical +phi convention.
        bz_parallel = self.signs.vacuum_b0_magnitude(self.grid.timep)

        def _parallel_from_jdotb(varname):
            jb = self.cdf.get_ncvar(varname, self.grid.area_inv)
            if not self.fill.is_valid(jb) or bz_parallel is None:
                return []
            try:
                return self.signs.jparallel_sign * np.divide(jb, bz_parallel[:,None])
            except (ValueError, TypeError, FloatingPointError) as exc:
                self.ctx.mapping.broken(exc, source=varname, target='profiles_1d/j_*')
                return []

        j_total = _parallel_from_jdotb('PLJB')
        j_ohmic = _parallel_from_jdotb('PLJBOH')
        j_nonind = []
        if self.fill.is_valid(j_total, j_ohmic):
            j_nonind = self.fill.operate_arrays(j_total, '-', j_ohmic)

        # CURBS is the bootstrap model actually used by TRANSP. Convert that
        # selected toroidal profile to the approximate <J.B>/B0 representation
        # with CONPLJB = approximate <J.B>/Jtor.  This avoids accidentally using
        # a diagnostic PLJBSNC/PLJBSNEO profile from a bootstrap model that was
        # not selected for the run.
        conpljb = self.cdf.get_ncvar('CONPLJB')
        j_bootstrap = []
        if bz_parallel is not None and self.fill.is_valid(curbs_phi, conpljb, bz_parallel):
            try:
                # curbs_phi already has the physical Ip sign. CONPLJB and BZ
                # are internal magnitudes, so their ratio only supplies geometry.
                j_bootstrap = curbs_phi * np.divide(conpljb, bz_parallel[:,None])
            except (ValueError, TypeError, FloatingPointError) as exc:
                self.ctx.mapping.broken(exc, source='CONPLJB', target='profiles_1d/j_bootstrap')
                j_bootstrap = []

        #vloop VSURC(TIME)
        self.fill.fill_data_nc(cp.global_quantities.v_loop, 'VSURC', mult=self.signs.vloop_sign, interp=True)

        #li  LI_3(TIME)
        self.fill.fill_data_nc(cp.global_quantities.li_3, 'LI_3', interp=True)

        # The betas by the Data Dictionary definitions, against this IDS's own b0 and r0 and, for the normalised
        # one, the boundary minor radius (see ``plasma_betas``).  TRANSP's BTEQ and BETAT use other normalisations.
        betas = self.signs.plasma_betas()
        if not betas:
            self.ctx.mapping.skipped('betas need PTOWB, DVOL, PCURC, the vacuum field and the RMCB boundary moments',
                                     target='global_quantities/beta_tor')
        for name in ('beta_tor', 'beta_pol', 'beta_tor_norm'):
            if name in betas and hasattr(cp.global_quantities, name):
                self.fill.fill_data(getattr(cp.global_quantities, name), betas[name])

        #Plasma diamagnetic energy content = 3/2*integral over the plasma volume of
        #  the total perpendicular pressure [J].
        #PPLAS(TIME3, X) PASCALS PLASMA PRESSURE
        pplas = self.profile('PTOWB') #'PPLAS'
        if len(pplas)> 0:
            edia = self.fill.np_multiply_sum(pplas, gvol, axis=1, scale=1.5)
            self.fill.fill_data(cp.global_quantities.energy_diamagnetic, edia)

        # Use the sum of the TIME3 zone volumes as the averaging denominator.
        # PVOL is a separate TIME history and can be out of sync with DVOL, or
        # differ slightly from its sum; the profiles and denominator must use
        # the same time samples and geometry.
        volume_profile = np.sum(gvol, axis=1) if self.fill.is_valid(gvol) else []

        #t_e_peaking
        # TE0 is sampled on scalar TIME, while TE and DVOL use profile TIME3.
        # Interpolate the axis temperature before computing the peaking factor.
        te0 = self.cdf.get_ncvar('TE0', interp=True)  # magnetic axis, TIME3
        tel = self.profile('TE')  #TE(TIME3,X) eV
        tev = []
        # DVOL and TE are both sampled on TIME3; their summed zone volume is
        # the consistent denominator for the volume average.

        te0_tev = []
        if self.fill.is_valid(tel, gvol, volume_profile, te0):
            with np.errstate(divide='ignore', invalid='ignore', over='raise'):
                thermal_energy = np.sum(np.multiply(tel, gvol), axis=1)
                tev = np.divide(
                    thermal_energy,
                    volume_profile,
                    out=np.zeros_like(thermal_energy, dtype=float),
                    where=np.asarray(volume_profile) != 0,
                )
                te0_tev = np.divide(
                    te0,
                    tev,
                    out=np.zeros_like(np.asarray(te0), dtype=float),
                    where=np.asarray(tev) != 0,
                )

        self.fill.fill_data(cp.global_quantities.t_e_peaking, te0_tev)

        #t_i_average_peaking
        ti0 = self.cdf.get_ncvar('TI0') #magnetic axis
        tion = self.profile('TI')  #TE(TIME3,X)EV
        tiv = []

        if self.fill.is_valid(ti0, tion, volume_profile):
            ti_ = self.fill.np_multiply_sum(tion, gvol, axis=1)
            tiv = self.fill.operate_arrays(ti_, '/', volume_profile)
            if not np.array_equal(self.grid.timep, self.grid.timesc):
               ti0_pt = self.grid.interp1d_time(self.grid.timep, self.grid.timesc, ti0)
            else:
               ti0_pt = ti0
            ti0_ti = self.fill.operate_arrays(ti0_pt, '/', tiv)
            self.fill.fill_data(cp.global_quantities.t_i_average_peaking, ti0_ti,
               name='t_i_average_peak')

        #resistive_psi_losses
        #ejima

        #t_e_volume_average
        if len(tev) >0:
            self.fill.fill_data(cp.global_quantities.t_e_volume_average, tev)

        #n_e_volume_average
        #NE(TIME3, X) N/CM3
        nel = self.profile('NE', self.grid.vol_inv)
        if self.fill.is_valid(nel, gvol, volume_profile):
            ne_ = self.fill.np_multiply_sum(nel, gvol, axis=1)
            ne_ = self.fill.operate_arrays(ne_, '/', volume_profile)
            self.fill.fill_data(cp.global_quantities.n_e_volume_average, ne_)

        #IONS
        #NI(TIME3, X)
        # NI is stored in N/cm^3 in TRANSP; IMAS requires m^-3.
        nion = self.profile('NI', self.grid.vol_inv)
        cp.global_quantities.ion.resize(1)
        #t_i_volume_average
        if len(tiv) >0:
            self.fill.fill_data(cp.global_quantities.ion[0].t_i_volume_average, tiv, name='tiv')

        if self.fill.is_valid(nion, gvol, volume_profile):
            ni_ = self.fill.np_multiply_sum(nion, gvol, axis=1)
            ni_ = self.fill.operate_arrays(ni_, '/', volume_profile)
            self.fill.fill_data(cp.global_quantities.ion[0].n_i_volume_average, ni_, name='niv')

        #ion_time_slice
        # ----------end of global -----------------------------------------

        #profiles_1d
        # X(TIME3,X) is sqrt(norm. toroidal flux)
        rho_tor_norm = coord_x

        #TRFLX(TIME3, XB)
        bflux_tor = self.cdf.get_ncvar('TRFLX')
        rho_tor = []
        # rho_tor = sqrt(phi/(pi*B0)) with B0 the IDS's own vacuum_toroidal_field/b0
        b0_rho = self.signs.vacuum_b0_magnitude(self.grid.timep)

        if b0_rho is None:
            self.ctx.mapping.skipped('vacuum toroidal field (RAXIS/BZXR) is not in the CDF',
                                    source='BZXR', target='profiles_1d/grid/rho_tor')
        else:
            try:
                bflux_tor_x = self.grid.interp_radial_profiles(bflux_tor, coord_xb, coord_x, axis_value=0.0)
                rho_tor = np.sqrt(np.clip(bflux_tor_x/(pi*b0_rho[:,None]), 0.0, None))
            except (ValueError, TypeError, IndexError, FloatingPointError) as exc:
                self.ctx.mapping.broken(exc, source='TRFLX', target='profiles_1d/grid/rho_tor')

        #(TIME3, XB)
        psib = self.profile('PLFLX2PI')
        if not self.fill.is_valid(psib):
            psib = self.cdf.get_ncvar('PLFLX', 2.*pi) #(TIME3, XB)
        if self.fill.is_valid(psib):
            psib = psi_sign * psib
        psi_rel = []
        for i in range(ntimep):
            # PLFLX/PLFLX2PI is the enclosed poloidal-flux increment relative
            # to the magnetic axis. Interpolate that increment onto X centers.
            psi_rel.append(np.interp(coord_x[i],
                                      np.concatenate(([0.0], coord_xb[i])),
                                      np.concatenate(([0.0], psib[i]))))
        psi_rel = np.array(psi_rel)

        # PSI0_TR is the magnetic-axis offset in Wb/rad. Convert to full Wb
        # and apply the target COCOS sign before constructing absolute psi.
        psi_axis = psi_sign * self.cdf.get_ncvar('PSI0_TR', 2.*pi)
        if not np.array_equal(self.grid.timep, self.grid.timesc):
            psi_axis = self.grid.interp1d_time(self.grid.timep, self.grid.timesc, psi_axis)

        psi = psi_rel + psi_axis[:,None] if self.fill.is_valid(psi_rel, psi_axis) else psi_rel
        psi_boundary = psib[:,-1] + psi_axis if self.fill.is_valid(psib, psi_axis) else []

        # Normalised poloidal flux coordinate.  The ratio is positive for both
        # COCOS11 and COCOS17 because numerator and denominator transform alike.
        rho_pol_norm = []
        try:
            denom = psi_boundary - psi_axis
            rho_pol_norm = np.sqrt(np.clip((psi - psi_axis[:,None]) / denom[:,None], 0.0, None))
        except (ValueError, TypeError, IndexError, FloatingPointError) as exc:
            self.ctx.mapping.broken(exc, source='PLFLX', target='profiles_1d/grid/rho_pol_norm')

        nrho = rho_tor_norm.shape[0]

        self.logger.info(f'  ntimep {ntimep} nrho {nrho}')
        self.logger.info(f'\tb0 {len(b0_rho)} rho_tor {len(rho_tor)} bflux_tor {len(bflux_tor)}')

        n_i_total_over_n_e = []
        if self.fill.is_valid(nion, nel):
            n_i_total_over_n_e = nion/nel

        # High-confidence thermodynamic quantities directly available in the CDF.
        # UE/UI/UTHRM are energy densities in J/cm^3.  For an isotropic thermal
        # species p = 2/3 U.  PPLAS is the TRANSP thermal plasma pressure and
        # agrees numerically with 2/3*UTHRM for the bundled verification case.
        ue_m3 = self.profile('UE', self.grid.vol_inv)
        ui_m3 = self.profile('UI', self.grid.vol_inv)
        utherm_m3 = self.profile('UTHRM', self.grid.vol_inv)
        pressure_e_thermal = (2.0/3.0) * ue_m3 if self.fill.is_valid(ue_m3) else []
        pressure_i_total = (2.0/3.0) * ui_m3 if self.fill.is_valid(ui_m3) else []
        pressure_thermal = self.profile('PPLAS')
        if not self.fill.is_valid(pressure_thermal) and self.fill.is_valid(utherm_m3):
            pressure_thermal = (2.0/3.0) * utherm_m3

        # TRANSP UFASTPP is perpendicular fast-ion kinetic energy density; for
        # two perpendicular degrees of freedom this equals p_perp. UFASTPA is
        # the parallel kinetic energy density = p_parallel/2.
        ufast_perp_m3 = self.profile('UFASTPP', self.grid.vol_inv)
        ufast_parallel_m3 = self.profile('UFASTPA', self.grid.vol_inv)
        pressure_perpendicular = []
        pressure_parallel = []
        if self.fill.is_valid(pressure_thermal, ufast_perp_m3):
            pressure_perpendicular = pressure_thermal + ufast_perp_m3
        if self.fill.is_valid(pressure_thermal, ufast_parallel_m3):
            pressure_parallel = pressure_thermal + 2.0*ufast_parallel_m3

        # Additional profile quantities with direct IMAS counterparts.
        magnetic_shear = self.profile('SHAT')
        eta_parallel = self.cdf.get_ncvar('ETA_USE', 0.01)  # ohm*cm -> ohm*m
        conductivity_parallel = []
        e_field_parallel = []
        if self.fill.is_valid(eta_parallel):
            with np.errstate(divide='ignore', invalid='ignore'):
                conductivity_parallel = np.where(np.abs(eta_parallel) > 0.0, 1.0/eta_parallel, 0.0)
            if self.fill.is_valid(j_ohmic):
                try:
                    e_field_parallel = eta_parallel * j_ohmic
                except (ValueError, TypeError, FloatingPointError) as exc:
                    self.ctx.mapping.broken(exc, source='ETA', target='profiles_1d/e_field/parallel')
                    e_field_parallel = []

        # Total toroidal angular-momentum density.  PPHI is the time/radius
        # profile (AMOM is only a scalar placeholder in some CDFs, including
        # 174783Q32).  TRANSP Nt*m*s/cm^3 converts by 1e6 to
        # kg*m^-1*s^-1, exactly the IMAS momentum_phi unit.
        momentum_phi = self.signs.toroidal_rotation_sign * self.profile('PPHI', self.grid.vol_inv)

        #ZEFFP(TIME3, X) PLASMA COMPOSITION ZEFF PROFILE
        zeff = self.profile('ZEFFP')

        #VRPOT(TIME3, XB); V; RADIAL ELECTRICAL POTENTIAL. Map XB -> X.
        vrpot_b = self.cdf.get_ncvar('VRPOT')
        phi_potential = []
        if self.fill.is_valid(vrpot_b):
            try:
                phi_potential = self.grid.interp_radial_profiles(
                    vrpot_b, coord_xb, coord_x, axis_value=vrpot_b[:,0])
            except (ValueError, TypeError, IndexError, FloatingPointError) as exc:
                self.ctx.mapping.broken(exc, source='VRPOT', target='profiles_1d/phi_potential')
                phi_potential = []
        # A radial electric field is deliberately not reconstructed here.
        # VRPOT is preserved as phi_potential; converting -grad(phi) to a
        # physical radial component requires an explicit metric convention.
        # The older implementation mixed zone area/volume factors and could
        # generate inconsistent shapes, so it has been removed.

        # Q profile at XB ; Q(TIME3, XB)
        qtmp = q_sign * self.profile('Q')
        q0 = q_sign * self.cdf.get_ncvar('Q0')
        if not np.array_equal(self.grid.timep, self.grid.timesc):
            q0 = self.grid.interp1d_time(self.grid.timep, self.grid.timesc, q0)
        qprof = []
        if len(qtmp)>1:
           for i in range(ntimep):
              qprof.append(np.interp(coord_x[i],
                                     np.concatenate(([0.0], coord_xb[i])),
                                     np.concatenate(([q0[i]], qtmp[i]))))
        qprof = np.array(qprof)

        self.logger.info(f'  Filling profile1d time slices: {ntimep}')
        cp.profiles_1d.resize(ntimep)

        #fill profiles_1d[*].grid
        #self.fill.fill_profiles1d_grid(cp) #TODO use this

        for i in range(ntimep):
          nrad = len(rho_tor_norm[i])        # radial size every profile of this slice must have
          #grid
          self.fill.fill_data(cp.profiles_1d[i].grid.rho_tor_norm, rho_tor_norm, ind=i, name='tho_tor_norm', chk=nrad)
          self.fill.fill_data(cp.profiles_1d[i].grid.rho_tor, rho_tor, ind=i, name='rho_tor', chk=nrad)
          self.fill.fill_data(cp.profiles_1d[i].grid.rho_pol_norm, rho_pol_norm, ind=i, name='rho_pol_norm', chk=nrad)
          self.fill.fill_data(cp.profiles_1d[i].grid.volume, volume_enclosed, ind=i, name='volume_enclosed', chk=nrad)
          self.fill.fill_data(cp.profiles_1d[i].grid.psi, psi, ind=i, name='psi', chk=nrad)
          self.fill.fill_data(cp.profiles_1d[i].grid.area, area_enclosed, ind=i, name='area_enclosed', chk=nrad)
          self.fill.fill_data(cp.profiles_1d[i].grid.surface, surface, ind=i, name='surface', chk=nrad)

          if len(psi_axis) > i:
              cp.profiles_1d[i].grid.psi_magnetic_axis = psi_axis[i].astype(float)
          if len(psi) > i:
              cp.profiles_1d[i].grid.psi_boundary = psi_boundary[i].astype(float)

          #electrons
          self.fill.fill_data(cp.profiles_1d[i].electrons.temperature, tel, ind=i, name='tel', chk=nrad)
          self.fill.fill_data(cp.profiles_1d[i].electrons.density, nel, ind=i, name='nel', chk=nrad)
          eprof = cp.profiles_1d[i].electrons
          if hasattr(eprof, 'pressure_thermal'):
              self.fill.fill_data(eprof.pressure_thermal, pressure_e_thermal, ind=i, name='pe_thermal', chk=nrad)
          if hasattr(eprof, 'pressure'):
              # No separate non-thermal electron population is represented here.
              self.fill.fill_data(eprof.pressure, pressure_e_thermal, ind=i, name='pe', chk=nrad)

          # Aggregate thermal/anisotropic quantities.
          p1d = cp.profiles_1d[i]
          if hasattr(p1d, 't_i_average'):
              self.fill.fill_data(p1d.t_i_average, tion, ind=i, name='t_i_average', chk=nrad)
          if hasattr(p1d, 'n_i_thermal_total'):
              self.fill.fill_data(p1d.n_i_thermal_total, nion, ind=i, name='n_i_thermal_total', chk=nrad)
          if hasattr(p1d, 'pressure_ion_total'):
              self.fill.fill_data(p1d.pressure_ion_total, pressure_i_total, ind=i, name='pressure_ion_total', chk=nrad)
          if hasattr(p1d, 'pressure_thermal'):
              self.fill.fill_data(p1d.pressure_thermal, pressure_thermal, ind=i, name='pressure_thermal', chk=nrad)
          if hasattr(p1d, 'pressure_perpendicular'):
              self.fill.fill_data(p1d.pressure_perpendicular, pressure_perpendicular, ind=i, name='pressure_perpendicular', chk=nrad)
          if hasattr(p1d, 'pressure_parallel'):
              self.fill.fill_data(p1d.pressure_parallel, pressure_parallel, ind=i, name='pressure_parallel', chk=nrad)
          if hasattr(p1d, 'magnetic_shear'):
              self.fill.fill_data(p1d.magnetic_shear, magnetic_shear, ind=i, name='magnetic_shear', chk=nrad)
          if hasattr(p1d, 'conductivity_parallel'):
              self.fill.fill_data(p1d.conductivity_parallel, conductivity_parallel, ind=i, name='conductivity_parallel', chk=nrad)
          if hasattr(p1d, 'e_field') and hasattr(p1d.e_field, 'parallel'):
              self.fill.fill_data(p1d.e_field.parallel, e_field_parallel, ind=i, name='e_field_parallel', chk=nrad)
          if hasattr(p1d, 'phi_potential'):
              self.fill.fill_data(p1d.phi_potential, phi_potential, ind=i, name='phi_potential', chk=nrad)
          # momentum_tor was renamed momentum_phi in DD 3.42.
          if self.dd.dd_at_least(3, 42, 0):
              if hasattr(p1d, 'momentum_phi'):
                  self.fill.fill_data(p1d.momentum_phi, momentum_phi, ind=i, name='momentum_phi', chk=nrad)
          elif hasattr(p1d, 'momentum_tor'):
              self.fill.fill_data(p1d.momentum_tor, momentum_phi, ind=i, name='momentum_tor', chk=nrad)

          #ions
          cp.profiles_1d[i].ion.resize(len(ions))
          if i==0:
              self.logger.info(f'\tnumber of ions {len(ions)}')
          for ij in range(len(ions)):
              ion = ions[ij]
              if i==0:
                  self.logger.info(f'\tfilling ion {ion.label}')
              cp.profiles_1d[i].ion[ij].element.resize(1)
              cp.profiles_1d[i].ion[ij].element[0].atoms_n = ion.atoms_n
              self.fill.set_species_name(cp.profiles_1d[i].ion[ij], ion.label)
              arr = len(ion.aas)>i and len(ion.zns)>i and len(ion.zions)>i
              aa = ion.aas[i].astype(float) if arr else ion.a
              z_n = ion.zns[i].astype(float) if arr else ion.z_n
              zion = ion.zions[i].astype(float) if arr else ion.z_ion
              cp.profiles_1d[i].ion[ij].z_ion = zion
              cp.profiles_1d[i].ion[ij].element[0].a = aa
              self.fill.set_nuclear_charge(cp.profiles_1d[i].ion[ij].element[0], z_n)

              #self.fill.fill_data(cp.profiles_1d[i].ion[ij].z_ion_1d, ion.z_ion_1d, ind=i, name='z_ion_1d')
              self.fill.fill_data(cp.profiles_1d[i].ion[ij].temperature, ion.temp, ind=i, name='temp', chk=nrad)
              self.fill.fill_data(cp.profiles_1d[i].ion[ij].density, ion.dens, ind=i, name='dens', chk=nrad)
              self.fill.fill_data(cp.profiles_1d[i].ion[ij].density_fast, ion.dens_fast, ind=i, name='densf', chk=nrad)
              self.fill.fill_data(cp.profiles_1d[i].ion[ij].density_thermal, ion.dens_th, ind=i, name='densth', chk=nrad)
              # Thermal species pressure p=n*T, with T in eV and n in m^-3.
              ion_node = cp.profiles_1d[i].ion[ij]
              if hasattr(ion_node, 'pressure_thermal') and self.fill.is_valid(ion.dens_th, ion.temp):
                  try:
                      pth = ion.dens_th * ion.temp * constants.e
                  except (ValueError, TypeError, FloatingPointError) as exc:
                      self.ctx.mapping.broken(exc, source='ion density/temperature',
                                             target='profiles_1d/ion/pressure_thermal')
                  else:
                      self.fill.fill_data(ion_node.pressure_thermal, pth, ind=i, name='ion_pressure_thermal', chk=nrad)
              # Radial NCLASS velocity is stored by TRANSP on the XB grid and has
              # the same number of points as the core-profile X grid.  Keep it for
              # now; a later cleanup can interpolate XB -> X explicitly.
              self.fill.fill_data(cp.profiles_1d[i].ion[ij].velocity.radial,
                                ion.vel_rad, ind=i, name='vel_rad', chk=len(rho_tor_norm[i]))

              # VPOL*_NC and VTOR*_NC are not X-grid flux profiles in TRANSP.
              # They are written on RMAJM, the 101-point equatorial-midplane grid
              # (inboard edge -> magnetic axis -> outboard edge), whereas IMAS
              # core_profiles ion velocity components use the profiles_1d radial
              # flux coordinate (50 X points for this run).  Older access layers
              # accepted the inconsistent shapes, but IMAS-Python DD4 validation
              # correctly rejects them.  Do not invent a midplane-to-flux-surface
              # averaging convention here: only fill them if a future TRANSP CDF
              # provides them directly on the core-profile grid.
              # Be explicit here rather than relying on fill_data(chk=...).  This
              # keeps the converter safe even if nc2imas.py from an older checkout
              # is accidentally mixed with this file.
              _nr = len(rho_tor_norm[i])
              _vtor = ion.vel_tor[i] if self.fill.is_valid(ion.vel_tor) and len(ion.vel_tor) > i else None
              _vpol = ion.vel_pol[i] if self.fill.is_valid(ion.vel_pol) and len(ion.vel_pol) > i else None
              if _vtor is not None and not np.isscalar(_vtor) and len(_vtor) == _nr:
                  self.fill.fill_data(cp.profiles_1d[i].ion[ij].velocity.toroidal,
                                    ion.vel_tor, ind=i, name='vel_tor')
              elif i == 0 and _vtor is not None and not np.isscalar(_vtor):
                  self.logger.info(f'\t\tskipping vel_tor: source grid has {len(_vtor)} points, core_profiles grid has {_nr}')

              if _vpol is not None and not np.isscalar(_vpol) and len(_vpol) == _nr:
                  self.fill.fill_data(cp.profiles_1d[i].ion[ij].velocity.poloidal,
                                    ion.vel_pol, ind=i, name='vel_pol')
              elif i == 0 and _vpol is not None and not np.isscalar(_vpol):
                  self.logger.info(f'\t\tskipping vel_pol: source grid has {len(_vpol)} points, core_profiles grid has {_nr}')

              # OMEGA/OMEGF/OMEGB are TRANSP X-grid quantities and therefore map
              # directly to the IMAS radial coordinate.
              self.fill.fill_data(cp.profiles_1d[i].ion[ij].rotation_frequency_tor,
                                ion.rot_freq_tor, ind=i, name='rftor', chk=len(rho_tor_norm[i]))
              #self.fill.fill_data(cp.profiles_1d[i].ion[ij].velocity.parallel, ion.vel_par, ind=i, name='vel_par')

              #states
              nstates = len(ion.states)
              cp.profiles_1d[i].ion[ij].state.resize(nstates)
              for ist, istate in enumerate(ion.states):
                  self.fill.fill_data(cp.profiles_1d[i].ion[ij].state[ist].density, istate.sdens, ind=i, name='sdens', chk=nrad)
                  self.fill.set_species_name(cp.profiles_1d[i].ion[ij].state[ist], istate.label)

          # neutrals
          cp.profiles_1d[i].neutral.resize(len(neutrals))
          if i==0:
              self.logger.info(f'\tnumber of neutrals {len(neutrals)}')
          for j in range(len(neutrals)):
              neut = neutrals[j]
              if i==0:
                  self.logger.info(f'\t\t {neut.label} temp {len(neut.temp)} dens {len(neut.dens)}')
              cp.profiles_1d[i].neutral[j].element.resize(1)
              cp.profiles_1d[i].neutral[j].element[0].atoms_n = neut.atoms_n
              self.fill.set_species_name(cp.profiles_1d[i].neutral[j], neut.label)
              cp.profiles_1d[i].neutral[j].element[0].a = neut.a
              self.fill.set_nuclear_charge(cp.profiles_1d[i].neutral[j].element[0], neut.z_n)

              self.fill.fill_data(cp.profiles_1d[i].neutral[j].temperature, neut.temp, ind=i, name='n0temp', chk=nrad)
              self.fill.fill_data(cp.profiles_1d[i].neutral[j].density, neut.dens, ind=i, name='n0dens', chk=nrad)
              # DD 3.39 has neutral.velocity/toroidal; DD 4.1.1 does not.
              # The original (signed) angular frequency is preserved for
              # both DD families in core_profiles.code.parameters XML.
              node = cp.profiles_1d[i].neutral[j]
              if (not self.dd.dd_at_least(4, 0, 0) and
                  hasattr(node, 'velocity') and hasattr(node.velocity, 'toroidal') and
                  neut.omega.size and radius_for_neutral_velocity.size):
                  self.fill.fill_data(node.velocity.toroidal,
                                    neut.omega * radius_for_neutral_velocity,
                                    ind=i, name='neutral_v_tor')

          #
          self.fill.fill_data_nc(cp.profiles_1d[i].t_i_average, 'TI', ind=i)
          self.fill.fill_data(cp.profiles_1d[i].q, qprof, ind=i, name='qprof', chk=nrad)
          self.fill.fill_data(cp.profiles_1d[i].zeff, zeff, ind=i, name='zeff', chk=nrad)
          self.fill.fill_data(cp.profiles_1d[i].n_i_total_over_n_e, n_i_total_over_n_e, ind=i, name='ni_o_ne', chk=nrad)

          self.fill.fill_data(cp.profiles_1d[i].j_bootstrap, j_bootstrap, ind=i, name='j_bootstrap_parallel', chk=nrad)
          self.fill.fill_data(cp.profiles_1d[i].j_ohmic, j_ohmic, ind=i, name='j_ohmic_parallel', chk=nrad)
          self.fill.fill_data(cp.profiles_1d[i].j_total, j_total, ind=i, name='j_total_parallel', chk=nrad)
          self.fill.fill_data(cp.profiles_1d[i].j_non_inductive, j_nonind, ind=i, name='j_nonind_parallel', chk=nrad)
          # Total toroidal current density was renamed j_tor -> j_phi
          # in late DD3.  DD3.39 uses j_tor; DD4.1.1 uses j_phi.
          if self.dd.dd_at_least(3, 42, 0):
              self.fill.fill_data(cp.profiles_1d[i].j_phi, cur_phi, ind=i, name='j_phi', chk=nrad)
          else:
              self.fill.fill_data(cp.profiles_1d[i].j_tor, cur_phi, ind=i, name='j_tor', chk=nrad)

          # e_field.radial intentionally left empty: VRPOT is preserved as
          # phi_potential, but converting -grad(phi) requires an explicit
          # radial metric convention.

        self.logger.info(f'  Done core_profiles')
        return [cp]

if __name__ == '__main__':
    print("Please use: python convert.py ")
