
from nc2imas import *

class ion_data(atom_base):
    def __init__(self, symb):
        super().__init__(symb)
        #arrays for times
        self.dens = []
        self.temp = []
        self.vel_rad = []
        self.vel_tor = []
        self.vel_pol = []
        self.vel_par = []
        self.z_ion_1d = []
        self.rotation_frequency_tor = []

    def is_valid(self):
        return (self.a>0. and self.atoms_n>0 and
               self.z_ion>0 and len(self.dens)>0 )

class neutrals_data(atom_base):
     def __init__(self, symb):
        super().__init__(symb)
        #arrays for times
        self.dens = []
        self.temp = []
        self.pres = []

class imas_core_profiles:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()

    #ions data
    def get_ion_data(self) -> list:
        ionDict = atom_base('').ions
        ions = []
        rmaj = np.average(self.im.get_ncvar('RMJSYM'))  # RMJSYM RMJM

        for short, label in ionDict.items():
            symb = 'HE' + short if short in ['3', '4'] else short
            dens = self.im.get_ncvar(f'N{symb}', 1.0e6)
            if len(dens) == 0:
                continue
            self.logger.info(f'ion: {symb} => {label}')
            ion = ion_data(symb)
            ion.dens = dens
            #FUSN HE4 ION DENSITY'%'N/CM**3'%FDENS_4
            #EFUSN_4  average fusion helium 4
            self.logger.info(f'{symb}, dens {len(ion.dens)}')
            ion.temp = self.im.get_ncvar(f'TI')
            self.logger.info(f'TI temp {len(ion.temp)}')
            ion.vel_rad = self.im.get_ncvar(f'VNDNC_{symb}', 0.01)
            self.logger.info(f'{symb} vel_rad {len(ion.vel_rad)}')
            ion.vel_pol = self.im.get_ncvar(f'VPOL{short}_NC', 0.01)
            self.logger.info(f'{symb} vel_pol {len(ion.vel_pol)}')
            ion.vel_tor = self.im.get_ncvar(f'VTOR{short}_NC', 0.01)
            self.logger.info(f'{symb} vel_tor {len(ion.vel_tor)}')
            ion.vel_par = [] #UPLLHATI/(BZXR/RAXIS)*1.0D-2
            ion.z_ion_1d = []
            ang_vel = self.im.get_ncvar(f'OMG_{short}') #CMG_H CMG_D

            if self.im.is_valid(ang_vel, rmaj):
                ion.rotation_frequency_tor = ang_vel / rmaj
            if ion.is_valid():
               ions.append(ion)

        #impurity = average in TRANSP, not only an element.
        symb = 'X'
        label = 'X'
        dens = self.im.get_ncvar('NIMP', 1.0e6) #XDENS.  NIMP or NIMP_NC
        if len(dens) == 0:
           return ions
        self.logger.info(f'impurity ion: {label}  => {symb}')
        ion = ion_data(symb)
        ion.dens = dens
        ion.temp = self.im.get_ncvar(f'TXNC')                   # TX  TXNC
        self.logger.info(f'{symb} temp {len(ion.temp)}')
        ion.vel_rad = self.im.get_ncvar('VELIM', 0.01)          # VELIM
        self.logger.info(f'{symb} vel_rad {len(ion.vel_rad)}')
        ion.vel_pol = self.im.get_ncvar('VPOLX_NC', 0.01)        # VPOLX_NC
        self.logger.info(f'{symb} vel_pol {len(ion.vel_pol)}')
        ion.vel_tor = self.im.get_ncvar(f'VTORX_NC', 0.01)       # VTORX  VTORX_NC
        self.logger.info(f'{symb} vel_tor {len(ion.vel_tor)}')
        ion.vel_par = []
        ion.z_ion_1d = []
        ang_vel = self.im.get_ncvar('OMGX') # OMGX CMGX
        if self.im.is_valid(ang_vel, rmaj):
            ion.rotation_frequency_tor = ang_vel/rmaj
        #AZIMP
        zn = -1.0 #not filled in IDS #TODO
        aa = np.average(self.im.get_ncvar('AIMP'))
        zion = np.average(self.im.get_ncvar('XZIMP'))   #XSIMP XZIMPD
        if self.im.is_valid(aa, zion, zn):
            ion.z_n = zn #zn.astype(float)
            ion.a = aa.astype(float)
            ion.z_ion = zion.astype(float)
        if ion.is_valid():
            ions.append(ion)
        return ions

    #volume/wall neutral density
    def get_neutral_data(self, wall:str=None):
        dst = "DN0W" if wall else "DN0V" #volume neutral
        neutDict = atom_base('').neutrals
        neuts = []
        for symb, label in neutDict.items():
            dens = self.im.get_ncvar(f'{dst}{symb}', 1.0e6)
            if len(dens)==0:
                continue
            self.logger.info(f'neutral: {symb} => {label}')
            neut = neutrals_data(label)
            neut.dens = dens
            self.logger.info(f'{symb}, dens {len(neut.dens)}')
            neut.temp = []
            neut.pres = []
            neuts.append(neut)
        return neuts


    def create_ids(self):
        self.im.logger.info('\ncore_profiles')
        cp = self.im.ids.core_profiles
        cm2m = 0.01
        aconv = 1.0e-4
        vconv = 1.0e-6
        pi = constants.pi  #scipy
        sgn = 1 #-1 if self.im.iter_ else 1 TODO

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        timep = self.im.get_ncvar('TIME3')
        ntimep = len(timep)
        coord_x = self.im.get_ncvar('X')
        coord_xb = self.im.get_ncvar('XB')
        self.im.checkErrExit(timesc, timep)

        cp = self.im.fill_properties(cp)
        self.im.fill_code_properties(cp)
        cp.ids_properties.homogeneous_time = 1
        cp.time.resize(ntimes)
        cp.time = timesc

        #vacuum toroidal field
        self.im.fill_vac_tor_field(cp)

        #global_quantities
        #------------------------------------------------------------------
        self.im.fill_data_nc(cp.global_quantities.ip, 'PCURC', mult=sgn)

        #current_non_inductive
        # grid area. DAREA (TIME3, X) CM**2 ZONE CROSS SECTIONAL AREA
        grid_area = self.im.get_ncvar('DAREA', aconv)

        # current_non_inductive. Toroidal component.
        # CUR(TIME3,X) AMPS/CM2. TOTAL PLASMA CURRENT
        cur_tot = self.im.get_ncvar('CUR', aconv)
        j_ohmic = self.im.get_ncvar('CUROH', aconv)
        if self.im.is_valid(cur_tot, j_ohmic, grid_area):
            curt = self.im.operate_arrays(cur_tot, '-', j_ohmic)
            j_nonind =  self.im.np_multiply_sum(curt, grid_area, axis=1)
            self.im.fill_data(cp.global_quantities.current_non_inductive, j_nonind, mult=sgn)

        #Bootstrap current (toroidal component). CURBS(TIME3, X) AMPS/CM2
        # check surface integration
        curbs = self.im.get_ncvar('CURBS', aconv)
        if len(curbs) > 0:
            curbsxa = sgn * self.im.np_multiply_sum(curbs, grid_area, axis=1)
            self.im.fill_data(cp.global_quantities.current_bootstrap, curbsxa)

        #vloop VSURC(TIME)
        self.im.fill_data_nc(cp.global_quantities.v_loop, 'VSURC')

        #li  LI_3(TIME)
        self.im.fill_data_nc(cp.global_quantities.li, 'LI_3')

        #beta_tor BTEQ(TIME) EQUILIBRIUM BETA(TOROIDAL)
        self.im.fill_data_nc(cp.global_quantities.beta_tor, 'BTEQ')

        #beta_tor_norm
        bz = self.im.get_ncvar('BZ') #(TIME) T. BZ outside plasma
        currents = self.im.get_ncvar('PCURC') #(TIME)
        btor = self.im.get_ncvar('BTEQ') #toroidal beta
        # PAREA(TIME) - plasma CS area[cm^2]
        parea = self.im.get_ncvar('PAREA', aconv)
        beta_tor_norm = []
        try:
            # 100 * beta_tor * a[m] * B0 [T] / ip [MA]   ???
            a = np.sqrt(parea/pi)
            beta_tor_norm = 100.*btor*a*bz/currents
        except:
            pass

        #beta_pol BETAT(TIME)
        self.im.fill_data_nc(cp.global_quantities.beta_pol, 'BETAT')

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', vconv)

        #Plasma diamagnetic energy content = 3/2*integral over the plasma volume of
        #  the total perpendicular pressure [J].
        #PPLAS(TIME3, X) PASCALS PLASMA PRESSURE
        pplas = self.im.get_ncvar('PPLAS') #TODO
        if len(pplas)> 0:
            edia = self.im.np_multiply_sum(pplas, gvol, axis=1, scale=1.5)
            self.im.fill_data(cp.global_quantities.energy_diamagnetic, edia)

        #z_eff_resistive

        #t_e_peaking
        te0 = self.im.get_ncvar('TE0') #magnetic axis
        tel = self.im.get_ncvar('TE')  #TE(TIME3,X) eV
        tev = []
        #volume PVOL(TIME) CM^3 : PLASMA VOLUME
        pvol = self.im.get_ncvar('PVOL', vconv)
        try:
            tev = np.sum(np.multiply(tel, gvol), axis=1) / pvol
            self.im.fill_data(cp.global_quantities.t_e_peaking, te0/tev)
        except:
            pass

        #t_i_average_peaking
        ti0 = self.im.get_ncvar('TI0') #magnetic axis
        tion = self.im.get_ncvar('TI')  #TE(TIME3,X)EV
        tiv = []
        if self.im.is_valid(ti0, tion, pvol):
            ti_ = self.im.np_multiply_sum(tion, gvol, axis=1)
            ti_ = self.im.operate_arrays(ti_, '/', pvol)
            ti0_ti = self.im.operate_arrays(ti0, '/', ti_)
            self.im.fill_data(cp.global_quantities.t_i_average_peaking, ti0_ti, name='t_i_average_peak')

        #resistive_psi_losses
        #ejima

        #t_e_volume_average
        if len(tev) >0:
          self.im.fill_data(cp.global_quantities.t_e_volume_average, tev)

        #n_e_volume_average
        #NE(TIME3, X) N/CM3
        nel = self.im.get_ncvar('NE', self.im.vol_inv)
        if len(nel) >0:
          ne_ = self.im.np_multiply_sum(nel, gvol, axis=1)
          ne_ = self.im.operate_arrays(ne_, '/', pvol)
          self.im.fill_data(cp.global_quantities.n_e_volume_average, ne_)

        #IONS
        #which ion ??
        #t_i_volume_average
        if len(tiv) >0:
          self.im.fill_data(cp.global_quantities.t_e_volume_average, tiv)

        #n_i_volume_average
        nion = self.im.get_ncvar('NI')
        if len(nion) >0:
          ni_ = self.im.np_multiply_sum(nion, gvol, axis=1)
          ni_ = self.im.operate_arrays( ni_, '/', pvol)
          self.im.fill_data(cp.global_quantities.n_e_volume_average, ni_)

        #ion_time_slice
        # ----------end of global -----------------------------------------

        #profiles_1d
        # X(TIME3,X) is sqrt(norm. toroidal flux)
        rho_tor_norm = coord_x #np.multiply(coord_x, coord_x)

        #TRFLX(TIME3, XB)
        bflux_tor = self.im.get_ncvar('TRFLX')
        rho_tor = []
        try:
            rho_tor = np.sqrt(bflux_tor/(pi*np.tile(bz[None,:].T, bflux_tor.shape[1])))
        except:
            self.im.logger.info('error calculating rho_tor')

        #(TIME3, XB)
        psi = self.im.get_ncvar('PLFLX2PI')
        if not self.im.is_valid:
            psi = self.im.get_ncvar('PLFLX', 2.*pi) #(TIME3, XB)

        #Surface area of the toroidal flux surface
        #FLUX SURFACE AREA (TIME3,X) cm^2
        surface = self.im.get_ncvar('SURF', 1.0e-4)

        psi_axis = self.im.get_ncvar('PSI0_TR') # Wb/rad
        #Normalised poloidal flux coordinate
        # = sqrt((psi(rho)-psi(magnetic_axis)) / (psi(LCFS)-psi(magnetic_axis)))
        #rho_pol_norm
        plflx = self.im.get_ncvar('PLFLX') #(TIME3, XB)
        rho_pol_norm = []
        try:
            rho_pol_norm = np.sqrt(abs((psi - psi[0])/(psi[-1]-psi[0])))
        except:
            pass

        nrho = rho_tor_norm.shape[1]

        self.logger.info(f'ntimep {ntimep} nrho {nrho}')
        self.logger.info(f'bz {len(bz)} rho_tor {len(rho_tor)} bflux_tor {len(bflux_tor)}')

        #ions
        ions = self.get_ion_data()

        n_i_total_over_n_e = []
        if self.im.is_valid(nion, nel):
          n_i_total_over_n_e = nion/nel

        neutrals = [] #self.get_neutral_data()

        #ZEFFP(TIME3, X) PLASMA COMPOSITION ZEFF PROFILE
        zeff = self.im.get_ncvar('ZEFFP')

        qprof = self.im.get_ncvar('Q')
        curbs = self.im.get_ncvar('CURBS', 1.e4)

        # SUMCUR Summed toroidal currents
        pcureq = self.im.get_ncvar('PCUREQ')
        jpfc = self.im.get_ncvar('CURSUM_PFC_MHD')
        jpassive = self.im.get_ncvar('CURSUM_PASS_MHD')
        j_tor = []
        if self.im.is_valid(pcureq, jpfc, jpassive):
          try :
            j_tor = pcureq + jpfc + jpassive
            #TODO convert to density
          except :
            self.logger.info(f'pcureq {len(pcureq)} jpfc {len(jpfc)} jpassive {len(jpassive)}')


        #CUROH(TIME3, X)AMPS/CM2
        j_ohmic = self.im.get_ncvar('CUROH', 1.0e4)

        # current_non_inductive. Toroidal component.
        # CUR(TIME3,X) AMPS/CM2. TOTAL PLASMA CURRENT
        cur_tot = self.im.get_ncvar('CUR', 1.0e4)
        j_nonind = cur_tot - j_ohmic

        self.logger.info(f'Filling profile1d time slices: {ntimep}')
        cp.profiles_1d.resize(ntimep)

        #fill profiles_1d[*].grid
        #self.im.fill_profiles1d_grid(cp) #TODO use this

        for i in range(ntimep):

          #grid
          self.im.fill_data(cp.profiles_1d[i].grid.rho_tor_norm, rho_tor_norm, ind=i)
          self.im.fill_data(cp.profiles_1d[i].grid.rho_tor, rho_tor, ind=i)
          self.im.fill_data(cp.profiles_1d[i].grid.rho_pol_norm, rho_pol_norm, ind=i)
          self.im.fill_data(cp.profiles_1d[i].grid.volume, gvol, ind=i)
          self.im.fill_data(cp.profiles_1d[i].grid.psi, psi, ind=i)
          self.im.fill_data(cp.profiles_1d[i].grid.area, grid_area, ind=i)
          self.im.fill_data(cp.profiles_1d[i].grid.surface, surface, ind=i)

          if len(psi_axis) > i:
              cp.profiles_1d[i].grid.psi_magnetic_axis = psi_axis[i].astype(float)
          if len(psi) > i:
              cp.profiles_1d[i].grid.psi_boundary = psi[i,-1].astype(float)

          #electrons
          self.im.fill_data(cp.profiles_1d[i].electrons.temperature, tel, ind=i)
          self.im.fill_data(cp.profiles_1d[i].electrons.density, nel, ind=i)

          #ions
          cp.profiles_1d[i].ion.resize(len(ions))
          if i==0:
              self.logger.info(f' number of ions {len(ions)}')
          for ij in range(len(ions)):
              ion = ions[ij]
              if i==0:
                  self.logger.info(f'* filling ion {ion.label}')
              cp.profiles_1d[i].ion[ij].element.resize(1)
              cp.profiles_1d[i].ion[ij].element[0].a = ion.a
              cp.profiles_1d[i].ion[ij].element[0].z_n = ion.z_n
              cp.profiles_1d[i].ion[ij].element[0].atoms_n = ion.atoms_n
              cp.profiles_1d[i].ion[ij].z_ion = ion.z_ion
              cp.profiles_1d[i].ion[ij].label = ion.label

              self.im.fill_data(cp.profiles_1d[i].ion[ij].z_ion_1d, ion.z_ion_1d, ind=i, name='z_ion_1d')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].temperature, ion.temp, ind=i, name='temp')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].density, ion.dens, ind=i, name='dens')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].velocity.radial, ion.vel_rad, ind=i, name='vel_rad')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].velocity.toroidal, ion.vel_tor, ind=i, name='vel_tor')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].velocity.poloidal, ion.vel_pol, ind=i, name='vel_pol')
              self.im.fill_data(cp.profiles_1d[i].ion[ij].velocity.parallel, ion.vel_par, ind=i, name='vel_par')

          #neutrals
          n_neutrals = len(neutrals)
          cp.profiles_1d[i].neutral.resize(n_neutrals)
          if i==0:
              self.logger.info(f'number of neutrals {n_neutrals}')
          for ij in range(n_neutrals):
              neut = neutrals[ij]
              if i==0:
                  self.logger.info(f'* filling neutral {neut.label}')
              cp.profiles_1d[i].neutral[ij].label = neut.label
              if neut.a > 0. and neut.z_n > 0. and neut.atoms_n > 0:
                  cp.profiles_1d[i].neutral[ij].element.resize(1)
                  cp.profiles_1d[i].neutral[ij].element[0].a = neut.a
                  cp.profiles_1d[i].neutral[ij].element[0].z_n = neut.z_n
                  cp.profiles_1d[i].neutral[ij].element[0].atoms_n = neut.atoms_n
                  self.im.fill_data(cp.profiles_1d[i].neutral[ij].density, neut.dens, ind=i, name='dens')

          #
          self.im.fill_data_nc(cp.profiles_1d[i].t_i_average, 'TI', ind=i)
          self.im.fill_data(cp.profiles_1d[i].q, qprof, ind=i)
          self.im.fill_data(cp.profiles_1d[i].zeff, zeff, ind=i)
          self.im.fill_data(cp.profiles_1d[i].n_i_total_over_n_e, n_i_total_over_n_e, ind=i)

          self.im.fill_data(cp.profiles_1d[i].j_bootstrap, curbs, ind=i)
          self.im.fill_data(cp.profiles_1d[i].j_ohmic, j_ohmic, ind=i)
          self.im.fill_data(cp.profiles_1d[i].j_total, cur_tot, ind=i)
          self.im.fill_data(cp.profiles_1d[i].j_non_inductive, j_nonind, ind=i)
          self.im.fill_data(cp.profiles_1d[i].j_tor, j_tor, ind=i)

        self.im.putdb_imas(cp)
        self.logger.info('Done core_profiles')

if __name__ == '__main__':
    print("Please use: python convert.py ")
