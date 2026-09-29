#!/usr/bin/python
"""Populate the IMAS core_profiles IDS from TRANSP profiles and metadata."""

from typing import List
from nc2imas import *

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
     def __init__(self, symb):
        super().__init__(symb)
        #arrays for times
        self.dens = []
        self.temp = []
        self.pres = []
     def is_valid(self, logger=None):
        check1 = self.a>0 and self.z_n>0 and self.z_ion>0
        check2 = self.atoms_n>0 and len(self.dens)>0
        valid = check1 and check2
        if logger:
            if not valid:
               logger.info(f'\t invalid neutral: {self.label} : a {self.a}'
                           f' z_n {self.z_n} z_ion {self.z_ion}')
        return valid


class imas_core_profiles:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        #for i in range(2):
        self.create_ids() #occur=i)

    #total ions
    def get_thermal_ion_data(self, atomsDef) -> List:
        ionDict = atomsDef.ions
        ions = []
        #RMJSYM(TIME3, RMJSYM)
        rmaj_ = self.im.get_ncvar('RMJSYM')
        rmaj = np.average(rmaj_) if self.im.is_valid(rmaj_) else 0
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
            ion.temp = self.im.get_ncvar(f'TI')
            self.logger.info(f'\t{symb}: TI {ion.temp.shape}')
            ion.vel_rad = self.im.get_ncvar(f'VNDNC_{symb}', 0.01)
            self.logger.info(f'\t{symb} vel_rad {ion.vel_rad.shape}')
            ion.vel_pol = self.im.get_ncvar(f'VPOL{short}_NC', 0.01)
            self.logger.info(f'\t{symb} vel_pol {ion.vel_pol.shape}')
            ion.vel_tor = self.im.toroidal_rotation_sign * self.im.get_ncvar(f'VTOR{short}_NC', 0.01)
            self.logger.info(f'\t{symb} vel_tor {ion.vel_tor.shape}')
            #ion.vel_par = [] #UPLLHATI/(BZXR/RAXIS)*1.0D-2
            #ion.z_ion_1d = []
            # RAD/SEC
            ion.rot_freq_tor = self.im.toroidal_rotation_sign * self.im.get_ncvar(f'OMEGA') #CMG_H CMG_D

            #if self.im.is_valid(ang_vel) and rmaj>0:
            #    ion.rot_freq_tor = ang_vel / rmaj
            if ion.is_valid(self.logger):
               ions.append(ion)
        return ions

    #impurities
    def get_imp_ion_data(self) -> List:
        #NIMPS_ total for all charge states, per species
        self.logger.info(' impurities:')
        ions = []
        imp_atoms = self.im.get_species_names('NIMPS_?', cut='NIMPS_')
        #multiple species if not TOK
        for symb in imp_atoms:
            #total
            dens = 'NIMP' if symb=='TOK' else f'NIMPS_{symb}'
            # (TIME3, X) N/CM3 ;TOTAL IMPURITY DENSITY
            den = self.im.get_ncvar(dens, scale=1.0e6)
            if len(den)==0:
               continue
            ion = ion_data(symb)
            self.logger.info(f'  {symb} => {ion.label}:')
            ion.dens = den
            ion.dens_th = den
            #temperature
            ion.temp = self.im.get_ncvar('TX')
            #charge states. Only if NIMP_symb_n is present
            dstr = f'NIMP_{symb}_'
            dnames, ncs = self.im.get_species_data(f'{dstr}?', cut=dstr)

            #ZIMPS_TOK or ZIMPS_<symb>(TIME3, X) <symb> Avg. Mult. Impurity Z
            zprof = self.im.get_ncvar(f'ZIMPS_{symb}')
            ion.z_ion_1d = zprof
            #a,z_n from definition in atoms_class, except for TOK

            #OMGX(TIME3, X)
            ion.rot_freq_tor = self.im.toroidal_rotation_sign * self.im.get_ncvar('OMGX')

            ai, zn = 0,0
            if symb=='TOK':
               #a   AIMP(TIME)
               ion.aas = self.im.get_ncvar(f'AIMP', interp=True)
               #z_n  XZIMP(TIME)
               ion.zns = self.im.get_ncvar(f'XZIMP', interp=True)
               if len(ion.aas)>0:
                  ai = ion.aas[0]
               if len(ion.zns)>0:
                  zn = ion.zns[0]
               #z_ion
               ion.zions = ion.zns #ion charge, fully ionized ?

            self.logger.info(f'\t{symb}: {ai}, {zn}, {dens} {ion.dens.shape}')
            self.logger.info(f'\t{symb}: charges: {ncs} NIMP_ {dnames} ZIMPS_{symb} {zprof.shape}')
            #charge states. eg:NIMP_BE_4
            for dname in dnames:
               ncharge = ncs[0]
               self.logger.info(f'\t{symb}: {dname} ch: {ncharge}')
               sdens = self.im.get_ncvar(dname, scale=1.0e6)
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
        return self.im.get_ncvar(f'{dstr}{ext}', 1.0e6)

    def get_total_ion_density(self, short, symb):
        tdens = self.get_ion_density('thermal', symb)
        bdens = self.get_ion_density('beam', short)
        fdens = self.get_ion_density('fusion', short)
        mdens = self.get_ion_density('mini', short)
        tot = self.im.sum_arrays(tdens, bdens, fdens, mdens)
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
            ion.temp = self.im.get_ncvar(f'{tstr}_{short}')
            # (TIME3, X) RAD/SEC
            ion.rot_freq_tor = self.im.toroidal_rotation_sign * self.im.get_ncvar(f'{freqs}_{short}')

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

    def get_ion_data_global(ions:ion_data) ->ion_data:
        iong = []
        #for ion in ions:


    #volume/wall neutral density
    def get_neutral_data_type(self, atomsDef, cat) -> List:
        dstr, tstr = '',''
        if cat=='wall':
            dstr = 'DN0W'
            tstr = 'T0W'
            neutDict = atomsDef.neutrals_wall
        elif cat=='volume':
            dstr = 'DN0V'
            tstr = 'T0V'
            neutDict = atomsDef.neutrals_volume
        neuts = []
        self.logger.info(f' {cat}-neutrals::')
        for symb, label in neutDict.items():
            dens = self.im.get_ncvar(f'{dstr}{symb}', 1.0e6)
            if len(dens)==0:
                continue
            self.logger.info(f'\tneutral: {symb} => {label}')
            neut = neutrals_data(label)
            neut.dens = dens
            neut.temp = self.im.get_ncvar(f'{tstr}{symb}', 1.0e6)
            self.logger.info(f'\t{symb}, dens {len(neut.dens)} temp {len(neut.temp)}')
            if neut.is_valid(self.logger):
               neuts.append(neut)
        return neuts


    def get_neutral_data(self):
        atomsDef = atom_base('')
        neuts = self.get_neutral_data_type(atomsDef, cat='wall')
        vneuts = self.get_neutral_data_type(atomsDef, cat='volume')
        if len(vneuts)>0:
            for neut in vneuts:
               neuts.append(neut)
        return neuts


    def create_ids(self): #, occur=0):
        self.im.logger.info('\ncore_profiles')
        cp = self.im.new_ids('core_profiles')
        cm2m = 0.01
        aconv = 1.0e-4
        vconv = 1.0e-6
        pi = constants.pi  #scipy
        ip_sign = self.im.ip_sign
        bt_sign = self.im.bt_sign
        psi_sign = self.im.psi_sign
        q_sign = self.im.q_sign

        #scalar time  len(TIME)
        ntimes = self.im.ntimes
        #profile time  len(TIME3)
        ntimep = self.im.ntimep
        coord_x = self.im.get_ncvar('X')
        coord_xb = self.im.get_ncvar('XB')
        if not self.im.is_valid(coord_x, coord_xb, name=r'coord_x/xb'):
            print('\t Error: invalid coord_x/xb', coord_x, coord_xb)
            return

        cp = self.im.fill_properties(cp)
        self.im.fill_code_properties(cp)
        cp.ids_properties.homogeneous_time = 1
        cp.time.resize(ntimep)
        cp.time = self.im.timep

        #NOTE: nc2imas.fill_data and fill_data_nc functions do interpolation. These
        #functions handle only 1D data. If 2D data is passed or read from NC variable,
        # a valid index is required for selecting data at that index.

        #If profile time(TIME3) and scalar time(TIME) are different, interp onto profile time
        #Since homogeneous time=1, all quantities should have the same time.
        interp = (ntimes != ntimep)

        # Geometry related quantities
        # grid area. DAREA (TIME3, X) CM**2 ZONE CROSS SECTIONAL AREA
        grid_area = self.im.get_ncvar('DAREA', aconv)
        area_boundary = self.im.cumulative_zone_quantity(grid_area)
        area_enclosed = self.im.interp_radial_profiles(area_boundary, coord_xb, coord_x, axis_value=0.0)
        #DVOL(TIME3, X) CM**3  ZONE VOLUME. Cumulative values live at XB.
        gvol = self.im.get_ncvar('DVOL', vconv)
        volume_boundary = self.im.cumulative_zone_quantity(gvol)
        volume_enclosed = self.im.interp_radial_profiles(volume_boundary, coord_xb, coord_x, axis_value=0.0)
        #SURF is defined on XB; map it to the X grid used by core_profiles.
        surface_boundary = self.im.get_ncvar('SURF', aconv)
        surface = self.im.interp_radial_profiles(surface_boundary, coord_xb, coord_x, axis_value=0.0)
        # PAREA(TIME) - plasma CS area[cm^2]
        parea = self.im.get_ncvar('PAREA', aconv)

        #vacuum toroidal field
        self.im.fill_vac_tor_field(cp)

        #ions for profiles_1d
        ions = self.get_ion_data()
        #ions = self.get_ion_data_occur(prof_id=occur)
        #ion_glob = self.get_ion_data_global(ions)

        # neutrals for profile_id
        neutrals = self.get_neutral_data()


        #global_quantities
        #------------------------------------------------------------------
        self.im.fill_data_nc(cp.global_quantities.ip, 'PCURC', mult=ip_sign, interp=True)

        # TRANSP CUR/CUROH/CURBS are TOROIDAL current densities [A/cm^2].
        # Keep them for toroidal/global-current quantities.  DD4 core_profiles
        # j_total/j_ohmic/j_non_inductive/j_bootstrap, however, are PARALLEL
        # current densities defined as <J.B>/B0.  Use PLJB/PLJBOH for those
        # fields and write CUR separately to j_phi.
        cur_phi = self.im.get_ncvar('CUR', self.im.area_inv)
        if self.im.is_valid(cur_phi):
            cur_phi = ip_sign * cur_phi

        curoh_phi = self.im.get_ncvar('CUROH', self.im.area_inv)
        if self.im.is_valid(curoh_phi):
            curoh_phi = ip_sign * curoh_phi

        # Global non-inductive current is explicitly the toroidal component.
        j_nonind_phi = []
        if self.im.is_valid(cur_phi, curoh_phi, grid_area):
            j_nonind_phi = self.im.operate_arrays(cur_phi, '-', curoh_phi)
            curt = self.im.np_multiply_sum(j_nonind_phi, grid_area, axis=1)
            self.im.fill_data(cp.global_quantities.current_non_inductive, curt)

        # Bootstrap global current is also a toroidal component.
        curbs_phi = self.im.get_ncvar('CURBS', self.im.area_inv)
        if self.im.is_valid(curbs_phi):
            curbs_phi = ip_sign * curbs_phi
        if len(curbs_phi) > 0:
            curbsxa = self.im.np_multiply_sum(curbs_phi, grid_area, axis=1)
            self.im.fill_data(cp.global_quantities.current_bootstrap, curbsxa)

        # Parallel current profiles. TRANSP PLJB* are <J.B> in A*T/cm^2.
        # TRANSP's scalar BZ is an unsigned/internal toroidal-field magnitude;
        # dividing by it and multiplying by ip_sign gives the IMAS definition
        # <J.B>/B0 in the physical +phi convention.
        bz_parallel = self.im.get_ncvar('BZ')
        if ntimep != ntimes:
            bz_parallel = self.im.interp1d_time(self.im.timep, self.im.timesc, bz_parallel)

        def _parallel_from_jdotb(varname):
            jb = self.im.get_ncvar(varname, self.im.area_inv)
            if not self.im.is_valid(jb, bz_parallel):
                return []
            try:
                return self.im.jparallel_sign * np.divide(jb, bz_parallel[:,None])
            except Exception:
                return []

        j_total = _parallel_from_jdotb('PLJB')
        j_ohmic = _parallel_from_jdotb('PLJBOH')
        j_nonind = []
        if self.im.is_valid(j_total, j_ohmic):
            j_nonind = self.im.operate_arrays(j_total, '-', j_ohmic)

        # CURBS is the bootstrap model actually used by TRANSP. Convert that
        # selected toroidal profile to the approximate <J.B>/B0 representation
        # with CONPLJB = approximate <J.B>/Jtor.  This avoids accidentally using
        # a diagnostic PLJBSNC/PLJBSNEO profile from a bootstrap model that was
        # not selected for the run.
        conpljb = self.im.get_ncvar('CONPLJB')
        j_bootstrap = []
        if self.im.is_valid(curbs_phi, conpljb, bz_parallel):
            try:
                # curbs_phi already has the physical Ip sign. CONPLJB and BZ
                # are internal magnitudes, so their ratio only supplies geometry.
                j_bootstrap = curbs_phi * np.divide(conpljb, bz_parallel[:,None])
            except Exception:
                j_bootstrap = []

        #vloop VSURC(TIME)
        self.im.fill_data_nc(cp.global_quantities.v_loop, 'VSURC', mult=self.im.vloop_sign, interp=True)

        #li  LI_3(TIME)
        self.im.fill_data_nc(cp.global_quantities.li_3, 'LI_3', interp=True)

        #beta_tor BTEQ(TIME) EQUILIBRIUM BETA(TOROIDAL)
        self.im.fill_data_nc(cp.global_quantities.beta_tor, 'BTEQ', interp=True)

        #beta_tor_norm
        bz = self.im.get_ncvar('BZ') #(TIME) T. BZ outside plasma
        btor = self.im.get_ncvar('BTEQ') #toroidal beta
        beta_tor_norm = []
        try:
            # 100 * beta_tor * a[m] * B0 [T] / ip [MA]   ???
            #PCURC(TIME)
            currents = ip_sign * self.im.get_ncvar('PCURC')
            bz_signed = bt_sign * bz
            a = np.sqrt(parea/pi)
            beta_tor_norm = 1.0e8*btor*a*bz_signed/currents
        except:
            pass

        #beta_pol BETAT(TIME)
        self.im.fill_data_nc(cp.global_quantities.beta_pol, 'BETAT', interp=True)

        #Plasma diamagnetic energy content = 3/2*integral over the plasma volume of
        #  the total perpendicular pressure [J].
        #PPLAS(TIME3, X) PASCALS PLASMA PRESSURE
        pplas = self.im.get_ncvar('PTOWB') #'PPLAS'
        if len(pplas)> 0:
            edia = self.im.np_multiply_sum(pplas, gvol, axis=1, scale=1.5)
            self.im.fill_data(cp.global_quantities.energy_diamagnetic, edia)

        #t_e_peaking
        te0 = self.im.get_ncvar('TE0') #magnetic axis
        tel = self.im.get_ncvar('TE')  #TE(TIME3,X) eV
        tev = []
        #volume PVOL(TIME) CM^3 : PLASMA VOLUME
        pvol_pt = self.im.get_ncvar('PVOL', vconv)
        #interpolate pvol from scalar time onto profile times
        if ntimep != ntimes:
            pvol_pt = self.im.interp1d_time(self.im.timep, self.im.timesc, pvol_pt)

        te0_tev = []
        try:
            tev = np.sum(np.multiply(tel, gvol), axis=1) / pvol_pt
            te0_tev = te0/tev
        except:
            pass

        self.im.fill_data(cp.global_quantities.t_e_peaking, te0_tev)

        #t_i_average_peaking
        ti0 = self.im.get_ncvar('TI0') #magnetic axis
        tion = self.im.get_ncvar('TI')  #TE(TIME3,X)EV
        tiv = []

        if self.im.is_valid(ti0, tion, pvol_pt):
            ti_ = self.im.np_multiply_sum(tion, gvol, axis=1)
            tiv = self.im.operate_arrays(ti_, '/', pvol_pt)
            if ntimep != ntimes:
               ti0_pt = self.im.interp1d_time(self.im.timep, self.im.timesc, ti0)
            else:
               ti0_pt = ti0
            ti0_ti = self.im.operate_arrays(ti0_pt, '/', tiv)
            self.im.fill_data(cp.global_quantities.t_i_average_peaking, ti0_ti,
               name='t_i_average_peak')

        #resistive_psi_losses
        #ejima

        #t_e_volume_average
        if len(tev) >0:
            self.im.fill_data(cp.global_quantities.t_e_volume_average, tev)

        #n_e_volume_average
        #NE(TIME3, X) N/CM3
        nel = self.im.get_ncvar('NE', self.im.vol_inv)
        if self.im.is_valid(nel, gvol, pvol_pt):
            ne_ = self.im.np_multiply_sum(nel, gvol, axis=1)
            ne_ = self.im.operate_arrays(ne_, '/', pvol_pt)
            self.im.fill_data(cp.global_quantities.n_e_volume_average, ne_)

        #IONS
        #NI(TIME3, X)
        # NI is stored in N/cm^3 in TRANSP; IMAS requires m^-3.
        nion = self.im.get_ncvar('NI', self.im.vol_inv)
        cp.global_quantities.ion.resize(1)
        #t_i_volume_average
        if len(tiv) >0:
            self.im.fill_data(cp.global_quantities.ion[0].t_i_volume_average, tiv, name='tiv')

        if self.im.is_valid(nion, gvol, pvol_pt):
            ni_ = self.im.np_multiply_sum(nion, gvol, axis=1)
            ni_ = self.im.operate_arrays( ni_, '/', pvol_pt)
            self.im.fill_data(cp.global_quantities.ion[0].n_i_volume_average, ni_, name='niv')

        #ion_time_slice
        # ----------end of global -----------------------------------------

        #profiles_1d
        # X(TIME3,X) is sqrt(norm. toroidal flux)
        rho_tor_norm = coord_x

        #TRFLX(TIME3, XB)
        bflux_tor = self.im.get_ncvar('TRFLX')
        rho_tor = []
        #interpolate onto profile time
        if ntimep != ntimes:
            bz = self.im.interp1d_time(self.im.timep, self.im.timesc, bz)

        try:
            bflux_tor_x = self.im.interp_radial_profiles(bflux_tor, coord_xb, coord_x, axis_value=0.0)
            rho_tor = np.sqrt(np.clip(bflux_tor_x/(pi*bz[:,None]), 0.0, None))
        except Exception:
            self.im.logger.info('error calculating rho_tor')

        #(TIME3, XB)
        psib = self.im.get_ncvar('PLFLX2PI')
        if not self.im.is_valid(psib):
            psib = self.im.get_ncvar('PLFLX', 2.*pi) #(TIME3, XB)
        if self.im.is_valid(psib):
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
        psi_axis = psi_sign * self.im.get_ncvar('PSI0_TR', 2.*pi)
        if ntimep != ntimes:
            psi_axis = self.im.interp1d_time(self.im.timep, self.im.timesc, psi_axis)

        psi = psi_rel + psi_axis[:,None] if self.im.is_valid(psi_rel, psi_axis) else psi_rel
        psi_boundary = psib[:,-1] + psi_axis if self.im.is_valid(psib, psi_axis) else []

        # Normalised poloidal flux coordinate.  The ratio is positive for both
        # COCOS11 and COCOS17 because numerator and denominator transform alike.
        rho_pol_norm = []
        try:
            denom = psi_boundary - psi_axis
            rho_pol_norm = np.sqrt(np.clip((psi - psi_axis[:,None]) / denom[:,None], 0.0, None))
        except Exception:
            pass

        nrho = rho_tor_norm.shape[0]

        self.logger.info(f'  ntimep {ntimep} nrho {nrho}')
        self.logger.info(f'\tbz {len(bz)} rho_tor {len(rho_tor)} bflux_tor {len(bflux_tor)}')

        n_i_total_over_n_e = []
        if self.im.is_valid(nion, nel):
            n_i_total_over_n_e = nion/nel

        # High-confidence thermodynamic quantities directly available in the CDF.
        # UE/UI/UTHRM are energy densities in J/cm^3.  For an isotropic thermal
        # species p = 2/3 U.  PPLAS is the TRANSP thermal plasma pressure and
        # agrees numerically with 2/3*UTHRM for the bundled verification case.
        ue_m3 = self.im.get_ncvar('UE', self.im.vol_inv)
        ui_m3 = self.im.get_ncvar('UI', self.im.vol_inv)
        utherm_m3 = self.im.get_ncvar('UTHRM', self.im.vol_inv)
        pressure_e_thermal = (2.0/3.0) * ue_m3 if self.im.is_valid(ue_m3) else []
        pressure_i_total = (2.0/3.0) * ui_m3 if self.im.is_valid(ui_m3) else []
        pressure_thermal = self.im.get_ncvar('PPLAS')
        if not self.im.is_valid(pressure_thermal) and self.im.is_valid(utherm_m3):
            pressure_thermal = (2.0/3.0) * utherm_m3

        # TRANSP UFASTPP is perpendicular fast-ion kinetic energy density; for
        # two perpendicular degrees of freedom this equals p_perp. UFASTPA is
        # the parallel kinetic energy density = p_parallel/2.
        ufast_perp_m3 = self.im.get_ncvar('UFASTPP', self.im.vol_inv)
        ufast_parallel_m3 = self.im.get_ncvar('UFASTPA', self.im.vol_inv)
        pressure_perpendicular = []
        pressure_parallel = []
        if self.im.is_valid(pressure_thermal, ufast_perp_m3):
            pressure_perpendicular = pressure_thermal + ufast_perp_m3
        if self.im.is_valid(pressure_thermal, ufast_parallel_m3):
            pressure_parallel = pressure_thermal + 2.0*ufast_parallel_m3

        # Additional profile quantities with direct IMAS counterparts.
        magnetic_shear = self.im.get_ncvar('SHAT')
        eta_parallel = self.im.get_ncvar('ETA_USE', 0.01)  # ohm*cm -> ohm*m
        conductivity_parallel = []
        e_field_parallel = []
        if self.im.is_valid(eta_parallel):
            with np.errstate(divide='ignore', invalid='ignore'):
                conductivity_parallel = np.where(np.abs(eta_parallel) > 0.0, 1.0/eta_parallel, 0.0)
            if self.im.is_valid(j_ohmic):
                try:
                    e_field_parallel = eta_parallel * j_ohmic
                except Exception:
                    e_field_parallel = []

        # Total toroidal angular-momentum density.  PPHI is the time/radius
        # profile (AMOM is only a scalar placeholder in some CDFs, including
        # 174783Q32).  TRANSP Nt*m*s/cm^3 converts by 1e6 to
        # kg*m^-1*s^-1, exactly the IMAS momentum_phi unit.
        momentum_phi = self.im.toroidal_rotation_sign * self.im.get_ncvar('PPHI', self.im.vol_inv)

        #ZEFFP(TIME3, X) PLASMA COMPOSITION ZEFF PROFILE
        zeff = self.im.get_ncvar('ZEFFP')

        #VRPOT(TIME3, XB); V; RADIAL ELECTRICAL POTENTIAL. Map XB -> X.
        vrpot_b = self.im.get_ncvar('VRPOT')
        phi_potential = []
        if self.im.is_valid(vrpot_b):
            try:
                phi_potential = self.im.interp_radial_profiles(
                    vrpot_b, coord_xb, coord_x, axis_value=vrpot_b[:,0])
            except Exception:
                phi_potential = []
        # A radial electric field is deliberately not reconstructed here.
        # VRPOT is preserved as phi_potential; converting -grad(phi) to a
        # physical radial component requires an explicit metric convention.
        # The older implementation mixed zone area/volume factors and could
        # generate inconsistent shapes, so it has been removed.

        # Q profile at XB ; Q(TIME3, XB)
        qtmp = q_sign * self.im.get_ncvar('Q')
        q0 = q_sign * self.im.get_ncvar('Q0')
        if ntimep != ntimes:
            q0 = self.im.interp1d_time(self.im.timep, self.im.timesc, q0)
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
        #self.im.fill_profiles1d_grid(cp) #TODO use this

        for i in range(ntimep):
          #grid
          self.im.fill_data(cp.profiles_1d[i].grid.rho_tor_norm, rho_tor_norm, ind=i, name='tho_tor_norm')
          self.im.fill_data(cp.profiles_1d[i].grid.rho_tor, rho_tor, ind=i, name='rho_tor')
          self.im.fill_data(cp.profiles_1d[i].grid.rho_pol_norm, rho_pol_norm, ind=i, name='rho_pol_norm')
          self.im.fill_data(cp.profiles_1d[i].grid.volume, volume_enclosed, ind=i, name='volume_enclosed')
          self.im.fill_data(cp.profiles_1d[i].grid.psi, psi, ind=i, name='psi')
          self.im.fill_data(cp.profiles_1d[i].grid.area, area_enclosed, ind=i, name='area_enclosed')
          self.im.fill_data(cp.profiles_1d[i].grid.surface, surface, ind=i, name='surface')

          if len(psi_axis) > i:
              cp.profiles_1d[i].grid.psi_magnetic_axis = psi_axis[i].astype(float)
          if len(psi) > i:
              cp.profiles_1d[i].grid.psi_boundary = psi_boundary[i].astype(float)

          #electrons
          self.im.fill_data(cp.profiles_1d[i].electrons.temperature, tel, ind=i, name='tel')
          self.im.fill_data(cp.profiles_1d[i].electrons.density, nel, ind=i, name='nel')
          eprof = cp.profiles_1d[i].electrons
          if hasattr(eprof, 'pressure_thermal'):
              self.im.fill_data(eprof.pressure_thermal, pressure_e_thermal, ind=i, name='pe_thermal')
          if hasattr(eprof, 'pressure'):
              # No separate non-thermal electron population is represented here.
              self.im.fill_data(eprof.pressure, pressure_e_thermal, ind=i, name='pe')

          # Aggregate thermal/anisotropic quantities.
          p1d = cp.profiles_1d[i]
          if hasattr(p1d, 't_i_average'):
              self.im.fill_data(p1d.t_i_average, tion, ind=i, name='t_i_average')
          if hasattr(p1d, 'n_i_thermal_total'):
              self.im.fill_data(p1d.n_i_thermal_total, nion, ind=i, name='n_i_thermal_total')
          if hasattr(p1d, 'pressure_ion_total'):
              self.im.fill_data(p1d.pressure_ion_total, pressure_i_total, ind=i, name='pressure_ion_total')
          if hasattr(p1d, 'pressure_thermal'):
              self.im.fill_data(p1d.pressure_thermal, pressure_thermal, ind=i, name='pressure_thermal')
          if hasattr(p1d, 'pressure_perpendicular'):
              self.im.fill_data(p1d.pressure_perpendicular, pressure_perpendicular, ind=i, name='pressure_perpendicular')
          if hasattr(p1d, 'pressure_parallel'):
              self.im.fill_data(p1d.pressure_parallel, pressure_parallel, ind=i, name='pressure_parallel')
          if hasattr(p1d, 'magnetic_shear'):
              self.im.fill_data(p1d.magnetic_shear, magnetic_shear, ind=i, name='magnetic_shear')
          if hasattr(p1d, 'conductivity_parallel'):
              self.im.fill_data(p1d.conductivity_parallel, conductivity_parallel, ind=i, name='conductivity_parallel')
          if hasattr(p1d, 'e_field') and hasattr(p1d.e_field, 'parallel'):
              self.im.fill_data(p1d.e_field.parallel, e_field_parallel, ind=i, name='e_field_parallel')
          if hasattr(p1d, 'phi_potential'):
              self.im.fill_data(p1d.phi_potential, phi_potential, ind=i, name='phi_potential')
          # momentum_tor was renamed momentum_phi in DD 3.42.
          if self.im.dd_at_least(3, 42, 0):
              if hasattr(p1d, 'momentum_phi'):
                  self.im.fill_data(p1d.momentum_phi, momentum_phi, ind=i, name='momentum_phi')
          elif hasattr(p1d, 'momentum_tor'):
              self.im.fill_data(p1d.momentum_tor, momentum_phi, ind=i, name='momentum_tor')

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
              self.im.set_species_name(cp.profiles_1d[i].ion[ij], ion.label)
              arr = len(ion.aas)>i and len(ion.zns)>i and len(ion.zions)>i
              aa = ion.aas[i].astype(float) if arr else ion.a
              z_n = ion.zns[i].astype(float) if arr else ion.z_n
              zion = ion.zions[i].astype(float) if arr else ion.z_ion
              cp.profiles_1d[i].ion[ij].z_ion = zion
              cp.profiles_1d[i].ion[ij].element[0].a = aa
              self.im.set_nuclear_charge(cp.profiles_1d[i].ion[ij].element[0], z_n)

              #self.im.fill_data(cp.profiles_1d[i].ion[ij].z_ion_1d, ion.z_ion_1d, ind=i, name='z_ion_1d')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].temperature, ion.temp, ind=i, name='temp')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].density, ion.dens, ind=i, name='dens')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].density_fast, ion.dens_fast, ind=i, name='densf')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].density_thermal, ion.dens_th, ind=i, name='densth')
              # Thermal species pressure p=n*T, with T in eV and n in m^-3.
              ion_node = cp.profiles_1d[i].ion[ij]
              if hasattr(ion_node, 'pressure_thermal') and self.im.is_valid(ion.dens_th, ion.temp):
                  try:
                      pth = ion.dens_th * ion.temp * constants.e
                      self.im.fill_data(ion_node.pressure_thermal, pth, ind=i, name='ion_pressure_thermal')
                  except Exception:
                      pass
              # Radial NCLASS velocity is stored by TRANSP on the XB grid and has
              # the same number of points as the core-profile X grid.  Keep it for
              # now; a later cleanup can interpolate XB -> X explicitly.
              self.im.fill_data(cp.profiles_1d[i].ion[ij].velocity.radial,
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
              _vtor = ion.vel_tor[i] if self.im.is_valid(ion.vel_tor) and len(ion.vel_tor) > i else None
              _vpol = ion.vel_pol[i] if self.im.is_valid(ion.vel_pol) and len(ion.vel_pol) > i else None
              if _vtor is not None and not np.isscalar(_vtor) and len(_vtor) == _nr:
                  self.im.fill_data(cp.profiles_1d[i].ion[ij].velocity.toroidal,
                                    ion.vel_tor, ind=i, name='vel_tor')
              elif i == 0 and _vtor is not None and not np.isscalar(_vtor):
                  self.logger.info(f'\t\tskipping vel_tor: source grid has {len(_vtor)} points, core_profiles grid has {_nr}')

              if _vpol is not None and not np.isscalar(_vpol) and len(_vpol) == _nr:
                  self.im.fill_data(cp.profiles_1d[i].ion[ij].velocity.poloidal,
                                    ion.vel_pol, ind=i, name='vel_pol')
              elif i == 0 and _vpol is not None and not np.isscalar(_vpol):
                  self.logger.info(f'\t\tskipping vel_pol: source grid has {len(_vpol)} points, core_profiles grid has {_nr}')

              # OMEGA/OMEGF/OMEGB are TRANSP X-grid quantities and therefore map
              # directly to the IMAS radial coordinate.
              self.im.fill_data(cp.profiles_1d[i].ion[ij].rotation_frequency_tor,
                                ion.rot_freq_tor, ind=i, name='rftor', chk=len(rho_tor_norm[i]))
              #self.im.fill_data(cp.profiles_1d[i].ion[ij].velocity.parallel, ion.vel_par, ind=i, name='vel_par')

              #states
              nstates = len(ion.states)
              cp.profiles_1d[i].ion[ij].state.resize(nstates)
              for ist, istate in enumerate(ion.states):
                  self.im.fill_data(cp.profiles_1d[i].ion[ij].state[ist].density, istate.sdens, ind=i, name='sdens')
                  self.im.set_species_name(cp.profiles_1d[i].ion[ij].state[ist], istate.label)

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
              self.im.set_species_name(cp.profiles_1d[i].neutral[j], neut.label)
              cp.profiles_1d[i].neutral[j].element[0].a = neut.a
              self.im.set_nuclear_charge(cp.profiles_1d[i].neutral[j].element[0], neut.z_n)

              self.im.fill_data(cp.profiles_1d[i].neutral[j].temperature, neut.temp, ind=i, name='n0temp')
              self.im.fill_data(cp.profiles_1d[i].neutral[j].density, neut.dens, ind=i, name='n0dens')

          #
          self.im.fill_data_nc(cp.profiles_1d[i].t_i_average, 'TI', ind=i)
          self.im.fill_data(cp.profiles_1d[i].q, qprof, ind=i, name='qprof')
          self.im.fill_data(cp.profiles_1d[i].zeff, zeff, ind=i, name='zeff')
          self.im.fill_data(cp.profiles_1d[i].n_i_total_over_n_e, n_i_total_over_n_e, ind=i,name='ni_o_ne')

          self.im.fill_data(cp.profiles_1d[i].j_bootstrap, j_bootstrap, ind=i, name='j_bootstrap_parallel')
          self.im.fill_data(cp.profiles_1d[i].j_ohmic, j_ohmic, ind=i, name='j_ohmic_parallel')
          self.im.fill_data(cp.profiles_1d[i].j_total, j_total, ind=i, name='j_total_parallel')
          self.im.fill_data(cp.profiles_1d[i].j_non_inductive, j_nonind, ind=i, name='j_nonind_parallel')
          # Total toroidal current density was renamed j_tor -> j_phi
          # in late DD3.  DD3.39 uses j_tor; DD4.1.1 uses j_phi.
          if self.im.dd_at_least(3, 42, 0):
              self.im.fill_data(cp.profiles_1d[i].j_phi, cur_phi, ind=i, name='j_phi')
          else:
              self.im.fill_data(cp.profiles_1d[i].j_tor, cur_phi, ind=i, name='j_tor')

          # e_field.radial intentionally left empty: VRPOT is preserved as
          # phi_potential, but converting -grad(phi) requires an explicit
          # radial metric convention.

        self.im.putdb_imas(cp)
        #self.im.putdb_imas_profile(cp, occur)
        self.logger.info(f'  Done core_profiles') #: occur {occur}')

if __name__ == '__main__':
    print("Please use: python convert.py ")
