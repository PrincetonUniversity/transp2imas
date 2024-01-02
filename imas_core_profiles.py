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
            ion.vel_tor = self.im.get_ncvar(f'VTOR{short}_NC', 0.01)
            self.logger.info(f'\t{symb} vel_tor {ion.vel_tor.shape}')
            #ion.vel_par = [] #UPLLHATI/(BZXR/RAXIS)*1.0D-2
            #ion.z_ion_1d = []
            # RAD/SEC
            ion.rot_freq_tor = self.im.get_ncvar(f'OMEGA') #CMG_H CMG_D

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
            ion.rot_freq_tor = self.im.get_ncvar('OMGX')

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
            ion.rot_freq_tor = self.im.get_ncvar(f'{freqs}_{short}')

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
        cp = self.im.ids.get('core_profiles')
        cm2m = 0.01
        aconv = 1.0e-4
        vconv = 1.0e-6
        pi = constants.pi  #scipy
        sgn = 1 #-1 if self.im.iter_ else 1 TODO

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
        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', vconv)
        #Surface area of the toroidal flux surface
        #FLUX SURFACE AREA (TIME3,X) cm^2
        surface = self.im.get_ncvar('SURF', aconv)
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
        self.im.fill_data_nc(cp.global_quantities.ip, 'PCURC', mult=sgn, interp=True)

        # current_non_inductive. Toroidal component.
        # CUR(TIME3,X) AMPS/CM2. TOTAL PLASMA CURRENT
        cur_tot = self.im.get_ncvar('CUR', aconv)

        #CUROH(TIME3, X)
        j_ohmic = self.im.get_ncvar('CUROH', aconv)
        j_nonind = []
        if self.im.is_valid(cur_tot, j_ohmic, grid_area):
            j_nonind = self.im.operate_arrays(cur_tot, '-', j_ohmic)
            curt =  self.im.np_multiply_sum(j_nonind, grid_area, axis=1)
            self.im.fill_data(cp.global_quantities.current_non_inductive, curt, mult=sgn)

        #Bootstrap current (toroidal component). CURBS(TIME3, X) AMPS/CM2
        # check surface integration
        curbs = self.im.get_ncvar('CURBS', aconv)
        if len(curbs) > 0:
            curbsxa = sgn * self.im.np_multiply_sum(curbs, grid_area, axis=1)
            self.im.fill_data(cp.global_quantities.current_bootstrap, curbsxa)

        #vloop VSURC(TIME)
        self.im.fill_data_nc(cp.global_quantities.v_loop, 'VSURC', interp=True)

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
            currents = self.im.get_ncvar('PCURC')
            a = np.sqrt(parea/pi)
            beta_tor_norm = 100.*btor*a*bz/currents
        except:
            pass

        #beta_pol BETAT(TIME)
        self.im.fill_data_nc(cp.global_quantities.beta_pol, 'BETAT', interp=True)

        #Plasma diamagnetic energy content = 3/2*integral over the plasma volume of
        #  the total perpendicular pressure [J].
        #PPLAS(TIME3, X) PASCALS PLASMA PRESSURE
        pplas = self.im.get_ncvar('PPLAS') #TODO
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
        nion = self.im.get_ncvar('NI')
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
            rho_tor = np.sqrt(bflux_tor/(pi*np.tile(bz[None,:].T, bflux_tor.shape[1])))
        except:
            self.im.logger.info('error calculating rho_tor')

        #(TIME3, XB)
        psi = self.im.get_ncvar('PLFLX2PI')
        if not self.im.is_valid(psi):
            psi = self.im.get_ncvar('PLFLX', 2.*pi) #(TIME3, XB)

        #PSI0_TR(TIME)  Wb/rad
        psi_axis = self.im.get_ncvar('PSI0_TR')
        if ntimep != ntimes:
            psi_axis = self.im.interp1d_time(self.im.timep, self.im.timesc, psi_axis)

        #Normalised poloidal flux coordinate
        # = sqrt((psi(rho)-psi(magnetic_axis)) / (psi(LCFS)-psi(magnetic_axis)))
        #rho_pol_norm
        plflx = psi
        rho_pol_norm = []
        try:
            rho_pol_norm = np.sqrt(abs((psi - psi[0])/(psi[-1]-psi[0])))
        except:
            pass

        nrho = rho_tor_norm.shape[0]

        self.logger.info(f'  ntimep {ntimep} nrho {nrho}')
        self.logger.info(f'\tbz {len(bz)} rho_tor {len(rho_tor)} bflux_tor {len(bflux_tor)}')

        n_i_total_over_n_e = []
        if self.im.is_valid(nion, nel):
            n_i_total_over_n_e = nion/nel

        #ZEFFP(TIME3, X) PLASMA COMPOSITION ZEFF PROFILE
        zeff = self.im.get_ncvar('ZEFFP')
        
        #VRPOT(TIME3, XB); V; RADIAL ELECTRICAL POTENTIAL
        vrpot = self.im.get_ncvar('VRPOT')
        # GXI(TIME3, X); 1/CM; <GRAD(XI)> FLUX SURF VOL.AVG
        gxi = self.im.get_ncvar('GXI', scale=1/0.01)
        efield_rad = []
        try:
            n = vrpot.shape[1]
            vrpot = vrpot[1:n] - vrpot[0:n-1]
            vol = np.cumsum(gvol)
            area = np.cumsum(darea)
            rminor = surface*area/(2*self.im.pi*vol)
            efield_rad[1:n] = vrpot/(rminor*gxi)  # m to cm ?
            efield_rad[0] = zefld[1]
        except:
            self.logger.info(f'\tfailed to calculate efield_rad '
                              f'sizes:: vrpot: {len(vrpot)}')

        # Q profile at XB ; Q(TIME3, XB)
        qtmp = self.im.get_ncvar('Q') 
        qprof = []
        for i in range(ntimep):
            #interpolate onto centers of zones
            qprof.append(np.interp(coord_x[i], coord_xb[i], qtmp[i]))

        self.logger.info(f'  Filling profile1d time slices: {ntimep}')
        cp.profiles_1d.resize(ntimep)

        #fill profiles_1d[*].grid
        #self.im.fill_profiles1d_grid(cp) #TODO use this

        for i in range(ntimep):
          #grid
          self.im.fill_data(cp.profiles_1d[i].grid.rho_tor_norm, rho_tor_norm, ind=i, name='tho_tor_norm')
          self.im.fill_data(cp.profiles_1d[i].grid.rho_tor, rho_tor, ind=i, name='rho_tor')
          self.im.fill_data(cp.profiles_1d[i].grid.rho_pol_norm, rho_pol_norm, ind=i, name='rho_pol_norm')
          self.im.fill_data(cp.profiles_1d[i].grid.volume, gvol, ind=i, name='gvol')
          self.im.fill_data(cp.profiles_1d[i].grid.psi, psi, ind=i, name='psi')
          self.im.fill_data(cp.profiles_1d[i].grid.area, grid_area, ind=i, name='grid_area')
          self.im.fill_data(cp.profiles_1d[i].grid.surface, surface, ind=i, name='surface')

          if len(psi_axis) > i:
              cp.profiles_1d[i].grid.psi_magnetic_axis = psi_axis[i].astype(float)
          if len(psi) > i:
              cp.profiles_1d[i].grid.psi_boundary = psi[i,-1].astype(float)

          #electrons
          self.im.fill_data(cp.profiles_1d[i].electrons.temperature, tel, ind=i, name='tel')
          self.im.fill_data(cp.profiles_1d[i].electrons.density, nel, ind=i, name='nel')

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
              cp.profiles_1d[i].ion[ij].label = ion.label
              arr = len(ion.aas)>i and len(ion.zns)>i and len(ion.zions)>i
              aa = ion.aas[i].astype(float) if arr else ion.a
              z_n = ion.zns[i].astype(float) if arr else ion.z_n
              zion = ion.zions[i].astype(float) if arr else ion.z_ion
              cp.profiles_1d[i].ion[ij].z_ion = zion
              cp.profiles_1d[i].ion[ij].element[0].a = aa
              cp.profiles_1d[i].ion[ij].element[0].z_n = z_n

              #self.im.fill_data(cp.profiles_1d[i].ion[ij].z_ion_1d, ion.z_ion_1d, ind=i, name='z_ion_1d')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].temperature, ion.temp, ind=i, name='temp')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].density, ion.dens, ind=i, name='dens')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].density_fast, ion.dens_fast, ind=i, name='densf')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].density_thermal, ion.dens_th, ind=i, name='densth')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].velocity.radial, ion.vel_rad, ind=i, name='vel_rad')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].velocity.toroidal, ion.vel_tor, ind=i, name='vel_tor')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].velocity.poloidal, ion.vel_pol, ind=i, name='vel_pol')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].rotation_frequency_tor, ion.rot_freq_tor, ind=i, name='rftor')
              #self.im.fill_data(cp.profiles_1d[i].ion[ij].velocity.parallel, ion.vel_par, ind=i, name='vel_par')

              #states
              nstates = len(ion.states)
              cp.profiles_1d[i].ion[ij].state.resize(nstates)
              for ist, istate in enumerate(ion.states):
                  self.im.fill_data(cp.profiles_1d[i].ion[ij].state[ist].density, istate.sdens, ind=i, name='sdens')
                  cp.profiles_1d[i].ion[ij].state[ist].label = istate.label

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
              cp.profiles_1d[i].neutral[j].label = neut.label
              cp.profiles_1d[i].neutral[j].element[0].a = neut.a
              cp.profiles_1d[i].neutral[j].element[0].z_n = neut.z_n

              self.im.fill_data(cp.profiles_1d[i].neutral[j].temperature, neut.temp, ind=i, name='n0temp')
              self.im.fill_data(cp.profiles_1d[i].neutral[j].density, neut.dens, ind=i, name='n0dens')

          #
          self.im.fill_data_nc(cp.profiles_1d[i].t_i_average, 'TI', ind=i)
          self.im.fill_data(cp.profiles_1d[i].q, qprof, ind=i, name='qprof')
          self.im.fill_data(cp.profiles_1d[i].zeff, zeff, ind=i, name='zeff')
          self.im.fill_data(cp.profiles_1d[i].n_i_total_over_n_e, n_i_total_over_n_e, ind=i,name='ni_o_ne')

          self.im.fill_data(cp.profiles_1d[i].j_bootstrap, curbs, ind=i, name='curbs')
          self.im.fill_data(cp.profiles_1d[i].j_ohmic, j_ohmic, ind=i, name='johmic')
          self.im.fill_data(cp.profiles_1d[i].j_total, cur_tot, ind=i, name='cur_tot')
          self.im.fill_data(cp.profiles_1d[i].j_non_inductive, j_nonind, ind=i, name='Jnonind')
          #self.im.fill_data(cp.profiles_1d[i].j_tor, j_tor, ind=i)

          self.im.fill_data(cp.profiles_1d[i].e_field.radial, efield_rad, ind=i, name='efield_rad')

        self.im.putdb_imas(cp)
        #self.im.putdb_imas_profile(cp, occur)
        self.logger.info(f'  Done core_profiles') #: occur {occur}')

if __name__ == '__main__':
    print("Please use: python convert.py ")
