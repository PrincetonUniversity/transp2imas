from nc2imas import *

class imas_summary:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.src = None
        self.create_ids()

    #from TRANSP output CDF file
    def create_ids(self):
        self.im.logger.info('\nsummary')
        sid = self.im.ids.summary
        iter_ = self.im.iter_
        sgn = -1 if iter_ else 1
        src = 'Calculated'

        #X(TIME3, X); XB(TIME3, XB)
        coord_x = self.im.get_ncvar('X')
        coord_xb = self.im.get_ncvar('XB')
        naxis_x = coord_x.shape[1]
        naxis_xb = coord_xb.shape[1]

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)

        #profile time
        timep = self.im.get_ncvar('TIME3')
        ntimep = len(timep)
        self.logger.info(f'ntime_prof {ntimep} ntime_sc {ntimes} nx {naxis_x} nxb {naxis_xb}')
        self.logger.info(f'ntimes {ntimes} {timesc[:3]} ..{timesc[-3:]}')

        # both TIME and TIME3 from the NETCDF file
        sid.ids_properties.homogeneous_time = 0
        sid.time.resize(ntimes)
        #scalar time
        sid.time = timesc
        sid = self.im.fill_properties(sid)

        #global_quantities from CDF file
        #ip : Total plasma current (toroidal component).
        #PCURC (TIME)
        # TODO current and magnetic field consistent  with COCOS definition
        self.im.fill_value_nc(sid.global_quantities.ip, 'PCURC', mult=sgn)
        #b0  BZ(TIME) - BZ @R=RMAJOR OUTSIDE PLASMA
        bz = self.im.get_ncvar('BZ')
        self.im.fill_value(sid.global_quantities.b0, bz, mult=sgn)

        #r0 FLT0.
        bzxR = self.im.get_ncvar('BZXR')
        if self.im.is_valid(bzxR):
            rvalues = bzxR/bz
            # I do not understand the line below
            r0s = 0.5*(rvalues[0]+rvalues[-1])*1.0e-2
            self.im.fill_value(sid.global_quantities.r0, r0s, rsz=0)

        #li  LI_3(TIME)
        self.im.fill_value_nc(sid.global_quantities.li, 'LI_3')
        #vloop VSURC(TIME)
        self.im.fill_value_nc(sid.global_quantities.v_loop, 'VSURC')

        #beta_pol BETAT(TIME)
        self.im.fill_value_nc(sid.global_quantities.beta_pol, 'BETAT')

        #beta_tor BTEQ(TIME) EQUILIBRIUM BETA(TOROIDAL)
        self.im.fill_value_nc(sid.global_quantities.beta_tor, 'BTEQ')

        #beta_tor_nor . 100*beta_tor*a[m]*B0[T]/ip[MA]
        btor = sid.global_quantities.beta_tor.value
        b0 = sid.global_quantities.b0.value
        ip = sid.global_quantities.ip.value
        # PAREA(TIME) - plasma CS area[cm^2]
        area = self.im.get_ncvar('PAREA')*1.e-4
        try:
            a = np.sqrt(area/self.im.pi)
            data = 1.0e8*btor*a*b0/ip
        except:
            self.im.logger.info('Error: beta_tor_norm')
        else:
            self.im.fill_value(sid.global_quantities.beta_tor_norm, data)

        #beta_tor_thermal_norm

        #TODO Question: scalar TIME is filled in sid.time. Can profile time
        # be different in CDF file? If so interpolation needed ?

        #energy_total FLT_1D; from UTOTL(TIME3,X) JLES/CM3
        utot = self.im.get_ncvar('UTOTL')
        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL')

        if self.im.is_valid(utot, gvol):
            utotv = self.im.np_multiply_sum(utot, gvol, axis=1)
            etot = self.im.operate_arrays(utotv, '*', 1.5)
            self.im.fill_value(sid.global_quantities.energy_total, etot)

        #PMHD(TIME3, X). Pressure from equilibrium solver 
        pmhd = self.im.get_ncvar('PMHD')
        if self.im.is_valid(pmhd, gvol):
            emhd = self.im.np_multiply_sum(pmhd, gvol, axis=1, scale=1.5)
            self.im.fill_value(sid.global_quantities.energy_mhd, emhd)

        #volume PVOL(TIME) CM^3 : PLASMA VOLUME
        pvol = self.im.get_ncvar('PVOL')
        self.im.fill_value_nc(sid.global_quantities.volume, pvol, mult=1.0e-6)

        #energy_thermal UTHRM(TIME3, X) JLES/CM3 THERMAL ENERGY DENSITY
        utherm = self.im.get_ncvar('UTHRM')
        if self.im.is_valid(utherm):
            etherm = self.im.np_multiply_sum(utherm, gvol, axis=1, scale=1.5)
            self.im.fill_value(sid.global_quantities.energy_thermal, etherm)

        #BPOL(TIME3, XB) Tesla
        bpol = self.im.get_ncvar('BPOL')
        if self.im.is_valid(bpol):
            try:
                for ib, bp in enumerate(bpol):
                    bpol[ib] = np.interp(coord_x[ib], coord_xb[ib], bp)
                bpol2 = np.multiply(bpol, bpol)
                ebpol = 1.0e6 * 1.0/(2.0*constants.mu_0) * np.sum(np.multiply(bpol2, gvol), axis=1)
            except:
                self.im.logger.info('error calculating ebpol')
            self.im.fill_value(sid.global_quantities.energy_b_field_pol, ebpol)

        # UFASTPP(TIME3, X) JLES/CM3 FAST ION PERP ENERGY DENSITY
        efastpp = self.im.get_ncvar('UFASTPP')
        if self.im.is_valid(efastpp):
            efastpp = self.im.np_multiply_sum(efastpp, gvol, axis=1, scale=1.5)
            self.im.fill_value(sid.global_quantities.energy_fast_perpendicular, efastpp)

        # UFASTPA(IME3, X) JLES/CM3 FAST ION PLL ENERGY DENSITY
        efastpl = self.im.get_ncvar('UFASTPA')
        if self.im.is_valid(efastpl):
            efastpl = self.im.np_multiply_sum(efastpl, gvol, axis=1, scale=1.5)
            self.im.fill_value(sid.global_quantities.energy_fast_parallel, efastpl)


        #Plasma diamagnetic energy content = 3/2*integral over the plasma volume of
        #  the total perpendicular pressure [J].
        #PPLAS(TIME3, X) PASCALS PLASMA PRESSURE
        pplas = self.im.get_ncvar('PPLAS') #TODO
        if self.im.is_valid(pplas):
            edia = self.im.np_multiply_sum(pplas, gvol, axis=1, scale=1.5*1.0e-6) #cm3 to m3
            self.im.fill_value(sid.global_quantities.energy_diamagnetic, edia)


        #denergy_diamagnetic_dt #TODO


        #current_ohm
        # grid area. DAREA (TIME3, X) CM**2 ZONE CROSS SECTIONAL AREA
        grid_area = self.im.get_ncvar('DAREA')
        #CUROH(TIME3, X)AMPS/CM2
        j_ohmic = self.im.get_ncvar('CUROH')
        if self.im.is_valid(grid_area, j_ohmic):
            current_ohm = self.im.np_multiply_sum(j_ohmic, grid_area, axis=1)
            self.im.fill_value(sid.global_quantities.current_ohm, current_ohm, mult=sgn)

        # current_non_inductive. Toroidal component.
        # CUR(TIME3,X) AMPS/CM2. TOTAL PLASMA CURRENT
        cur_tot = self.im.get_ncvar('CUR')
        if self.im.is_valid(cur_tot):
            curm = self.im.operate_arrays(cur_tot, '-', j_ohmic)
            j_nonind =  self.im.np_multiply_sum(curm, grid_area, axis=1)
            self.im.fill_value(sid.global_quantities.current_non_inductive, j_nonind, mult=sgn)

        #tau_energy. TAUA1(TIME), SECONDS, ENERGY CONFINEMENT (TOTAL)
        self.im.fill_value_nc(sid.global_quantities.tau_energy, 'TAUA1')

        #h_98; H98Y2 (TIME); TauE98y,2 confinement Hfactor
        self.im.fill_value_nc(sid.global_quantities.h_98, 'H98Y2')

        #h_mode, LHMODE(TIME); H-Mode indicator
        self.im.fill_value_nc(sid.global_quantities.h_mode, 'LHMODE')

        #q95
        self.im.fill_value_nc(sid.global_quantities.q_95, 'Q95')

        #Bootstrap current (toroidal component). CURBS(TIME3, X) AMPS/CM2
        # check surface itnegration
        curbs = self.im.get_ncvar('CURBS')
        if self.im.is_valid(curbs):
            curbsxa = self.im.np_multiply_sum(curbs, grid_area, axis=1, scale=sgn)
            self.im.fill_value(sid.global_quantities.current_bootstrap, curbsxa)


        #fusion_gain


        # TAUES(TIME3, XB) PLASMA ENERGY CONFINEMENT
        # this section is incorrect, I am commenting it out for now
        #taues = self.im.get_ncvar('TAUES')
        #if taues is not None:
        #    for it,ta in enumerate(taues) :
        #        taues[it] = np.interp(coord_x[0], coord_xb[0], ta)
        #    self.im.fill_value(sid.global_quantities.tau_energy, taues)

        #TAUE98Y2(TIME) TauE98y,2 confinement scaling
        self.im.fill_value_nc(sid.global_quantities.tau_energy_98, 'TAUE98Y2')

        #ETA_USE(TIME3, X)OHM*CM. RESISTIVITY USED OR INFERRED
        #resist = self.im.get_ncvar('ETA_USE')
        #if resist is not None:
        #    resist = self.im.np_multiply_sum(resist, coord_x, axis=1) #TODO
        #    self.im.fill_value(sid.global_quantities.resistance, resist)

        #POH(TIME3, X) WATTS/CM3, OHMIC HEATING POWER
        poh = self.im.get_ncvar('POH')
        if self.im.is_valid(poh):
            poh = self.im.np_multiply_sum(poh, gvol, axis=1)
        self.im.fill_value(sid.global_quantities.power_ohm, poh)

        #PRAD(TIME3, X) WATTS/CM3 NET RADIATED POWER USED
        prad = self.im.get_ncvar('PRAD')
        if self.im.is_valid(prad):
            prad = self.im.np_multiply_sum(prad, gvol, axis=1)
        self.im.fill_value(sid.global_quantities.power_radiated, prad)

        #PRAD_LI(TIME3, X) WATTS/CM3 , LINE RADIATION
        prad_li = self.im.get_ncvar('PRAD_LI')
        if self.im.is_valid(prad_li):
            prad_li = self.im.np_multiply_sum(prad_li, gvol, axis=1)
        self.im.fill_value(sid.global_quantities.power_line, prad_li)

        #PRAD_BR(TIME3, X) WATTS/CM3  BREMSSTRAHLUNG RADIATION
        prad_br = self.im.get_ncvar('PRAD_BR')
        if self.im.is_valid(prad_br):
            prad_br = self.im.np_multiply_sum(prad_br, gvol, axis=1)
        self.im.fill_value(sid.global_quantities.power_bremsstrahlung, prad_br)

        #PRAD_CY(TIME3, X) WATTS/CM3 CYCLOTRON RADIATION
        prad_cy = self.im.get_ncvar('PRAD_CY')
        if self.im.is_valid(prad_cy):
            prad_cy = self.im.np_multiply_sum(prad_cy, gvol, axis=1)
        self.im.fill_value(sid.global_quantities.power_synchrotron, prad_cy)

        #power_loss through separatrix [W]

        #fusion_fluence


        #si_external_average


        # X and XB represent two different radial coordinates.
        # Need to make sure that the radial coordinate is consistent with the profiles
        # X is sqrt(norm. toroidal flux) defined in zone center
        # XB is sqrt(norm. toroidal flux) defined in zone boundary
        #rho_tor_norm. XB(TIME3, XB).

        xcoord = self.im.operate_arrays(coord_x, '*', coord_x) #multiply ?
        # only one coordinate at one time ?
        #self.im.fill_value(sid.local.magnetic_axis.rho_tor_norm, xcoord[:,-1]) #TODO


        #magnetic axis position
        self.im.fill_data_nc(sid.local.magnetic_axis.position.r, 'RAXIS', mult=1.e-2)
        self.im.fill_data_nc(sid.local.magnetic_axis.position.z, 'YAXIS', mult=1.e-2)

        #magnetic axis n_e, NE(TIME3,X) N/CM**3
        n_e = self.im.get_ncvar('NE')
        if self.im.is_valid(n_e):
            self.im.fill_value(sid.local.magnetic_axis.n_e, n_e[:,0], mult=1.0e6)

        #magnetic axis  t_e.  TE(TIME3,X); EV,  ELECTRON TEMPERATURE
        self.im.fill_value_nc(sid.local.magnetic_axis.t_e, 'TE0')

        #n_i_total
        n_i = self.im.get_ncvar('NI')
        if self.im.is_valid(n_i):
            self.im.fill_value(sid.local.magnetic_axis.n_i_total, n_i[:,0], mult=1.0e6)

        #t_i_average  over ion species TODO
        self.im.fill_value_nc(sid.local.magnetic_axis.t_i_average, 'TI0')
        #zeff
        # if only the value on-axis is needed, we can take ZEFFC
        #zeff = self.im.get_ncvar('ZEFFP')
        self.im.fill_value_nc(sid.local.magnetic_axis.zeff, 'ZEFFC')
        #self.im.fill_value(sid.local.magnetic_axis.zeff, zeff[:,0])
        #q
        self.im.fill_value_nc(sid.local.magnetic_axis.q, 'Q0')

        #magnetic_shear=rho_tor/q*dq/drho_tor; SHAT(TIME3, X)=(r/q)*(dq/dr)
        # not sure this is correct, I am commenting this out
        #mag = self.im.get_ncvar('SHAT')
        #if mag is not None:
        #    magsh = mag.copy()
        #    #average of consecutive columns, except for the last
        #    magsh[:,:-1] += mag[:,1:] #add next column
        #    magsh[:,-1] = 3*mag[:,-1] - mag[:,-2] # rdum
        #    magsh /= 2
        #    self.im.fill_value(sid.local.magnetic_axis.magnetic_shear, magsh[:,0])

        #Ion densities
        #DEUTERIUM  ND(TIME3, X); N/CM**3
        self.im.fill_ion_densities(sid.local.magnetic_axis.n_i, ind2=0)

        #VTORD_NC(TIME3, RMAJM) CM/SEC ; NC D+ toroidal velocity

        #OMEGA(TIME3, X)
        #sid.local.magnetic_axis.velocity_tor.hydrogen


        #runaways.particles. Number of runaway electrons [-]. FLT_1D 1..time.
        # SERUN(TIME3, X) N/CM3/SEC . RUNAWAY ELEC SOURCE RATE
        serun = self.im.get_ncvar('SERUN')
        if self.im.is_valid(serun):
            run_ptcls = self.im.np_multiply_sum(serun, gvol, axis=1)
            self.im.fill_value(sid.runaways.particles, run_ptcls)

        #runaways.current. Parallel current driven by the runaway electrons [A]

        #EC
        #Set of ECRH/ECCD launchers. 6 Max in MDS+.
        #PEECH+1..n gyrotrons. (TIME3,X) WATTS/CM3
        # please, fix this. The number of gyrotrons should not be hardcoded
        ec_ngyro = self.im.get_var_names_re('PECIN?', nd1=1, nd2=2, length=True)

        self.logger.info(rf'# EC gyrotrons {ec_ngyro}')
        peech, pecin = 0., 0.
        sid.heating_current_drive.ec.resize(ec_ngyro)
        for i in range(ec_ngyro):
            peechi = self.im.get_ncvar('PEECH'+str(i+1))
            if self.im.is_valid(peechi):
                peech += peechi
                power = self.im.np_multiply_sum(peechi, gvol, axis=1)
                self.im.fill_value(sid.heating_current_drive.ec[i].power, power)

            # launched or injected is PECINn WATTS
            pecini = self.im.get_ncvar('PECIN'+str(i+1))
            if self.im.is_valid(pecini):
                pecin += pecini
                self.im.fill_value(sid.heating_current_drive.ec[i].power_launched, pecini)

            #ECCURn(TIME3, X) AMPS/CM2  ;ECRH CURRENT (GYROTRON n)
            eccur = self.im.get_ncvar('ECCUR'+str(i+1))
            if self.im.is_valid(eccur):
                current = self.im.np_multiply_sum(eccur, grid_area, axis=1)
                self.im.fill_value(sid.heating_current_drive.ec[i].current, current)

            #Toroidal angle of ECRH at resonance [rad]
            #PHAIECHn(TIME) DEGREES; Toroidal aiming, GYROTRON n
            #TRANSP defines angles in the laboratory frame. IDS defines angles in the EC frame.
            angle = self.im.get_ncvar('PHAIECH'+str(i+1))
            if self.im.is_valid(angle):
                angle = angle*self.im.pi/180.  # (angle-180.0) ec_launchers
                self.im.fill_value(sid.heating_current_drive.ec[i].angle_tor, angle)

            #Poloidal angle of ECRH at resonance [rad]
            #THETECHn(TIME) DEGREES Poloidal aiming, GYROTRON n
            angle = self.im.get_ncvar('THETECH'+str(i+1))
            if self.im.is_valid(angle):
                angle = -(angle+90.0)*self.im.pi/180. #(angle-90.0) #ec_launchers
                self.im.fill_value(sid.heating_current_drive.ec[i].angle_pol, angle)

        # Total EC power, into vacuum vessel [W]. PECIN(TIME), ECRF INPUT POWER
        if self.im.is_valid(pecin):
            self.im.fill_value(sid.heating_current_drive.power_launched_ec, pecin)
        # Total EC power coupled to plasma. PEECH(TIME3,X) WATTS/CM3 ;ECRH ELECTRON HEATING
        if self.im.is_valid(peech):
            power_ec = self.im.np_multiply_sum(peech, gvol, axis=1)
            self.im.fill_value(sid.heating_current_drive.power_ec, power_ec)


        #IC
        # do not hardcode the number of IC antennae
        nic_ant = self.im.get_var_names_re('PICHA?', nd1=1, nd2=2, length=True)
        self.logger.info(rf'# IC Antenna {nic_ant}')
        picha, pichant = 0., 0.
        sid.heating_current_drive.ic.resize(nic_ant)
        for i in range(nic_ant):
            #PICHAn(TIME)
            pichai = self.im.get_ncvar('PICHA'+str(i+1))
            if self.im.is_valid(pichai):
                picha += pichai
                self.im.fill_value(sid.heating_current_drive.ic[i].power_launched, pichai)
            #RFPWRABn  WATTS. Power ABSORBED, Antenna #n  (PICF01N01 WATTS/CM3)
            #  I will doublecheck this
            self.im.fill_value_nc(sid.heating_current_drive.ic[i].power, 'RFPWRAB'+str(i+1))
            #FREQAn(TIME)
            self.im.fill_value_nc(sid.heating_current_drive.ic[i].frequency, 'FREQA'+str(i+1))

            pichanti = self.im.get_ncvar('PICHANT'+str(i+1))
            if self.im.is_valid(pichanti):
                pichant += pichanti

            #PRFFI  WATTS ICRF POWER to FAST IONS

        #  PICHA = 1 ;
        #  PICHANT = 0 ;

        # PICHA(TIME) WATTS. TOTAL ICRF ANTENNA POWER
        if np.any(picha > 0.):
            self.im.fill_value(sid.heating_current_drive.power_launched_ic, picha)
        # PICHANT (TIME3,X)  WATTS/CM3 . ICRF HEATING by ANTENNA
        if np.any(pichant > 0.):
            pichant = self.im.np_multiply_sum(pichant, gvol, axis=1)
            self.im.fill_value(sid.heating_current_drive.power_ic, pichant)

        #heating_current_drive.nbi




        #Boundary
        #RMJSYM CM,  MAJOR RADII (DATA MAPPING), radial grid for symmetrization;
        #RMAJDIFF(TIME) ; CM ;MAX DIFF BOUND AND FLUX MIDPLANE
        #RMAJBLIM(TIME) ;CM;  RMAJ EXTENT OF THE BOUNDARY   '
        #RBOUN(TIME3, XB)

        bzxR = self.im.get_ncvar('BZXR')
        try:
            rvalues = bzxR/bz
            rcent = 0.5*(rvalues[0]+rvalues[-1])*1.0e-2
            rcent = np.tile(np.atleast_1d(rcent), ntimep)
        except:
            self.im.logger.info('Error: geometric_axis_r')
        else:
            self.im.fill_value(sid.boundary.geometric_axis_r, rcent, name='rcent')

        #sid.boundary.geometric_axis_z

        self.im.fill_value_nc(sid.boundary.magnetic_axis_r, 'RAXIS', mult=1.e-2)
        self.im.fill_value_nc(sid.boundary.magnetic_axis_z, 'YAXIS', mult=1.e-2)

        #boundary.minor_radius
        # 1/XB last * radius  ?

        #ELONG(TIME3, XB)
        elong = self.im.get_ncvar('ELONG')
        self.im.fill_value(sid.boundary.elongation, elong[-1])

        #TRIANGU(TIME3, XB)
        triang = self.im.get_ncvar('TRIANGU')
        self.im.fill_value(sid.boundary.triangularity_upper, triang[-1])

        #TRIANGL(TIME3, XB)
        triang = self.im.get_ncvar('TRIANGL')
        self.im.fill_value(sid.boundary.triangularity_lower, triang[-1])

        #boundary.strike_point_inner_r
        #boundary.strike_point_inner_z
        #boundary.strike_point_outer_r
        #boundary.strike_point_outer_z
        #boundary.strike_point_configuration
        #boundary.gap_limiter_wall
        #boundary.distance_inner_outer_separatrices
        #boundary.x_point_main

        if not self.im.is_valid(timep, ind=len(timesc)-1):
            self.im.logger.info(f'timesc {timesc.shape} timep {timep.shape}')

        # volume_average
        def fill_volume_average(name:str):
            dat = self.im.get_ncvar(name)
            if self.im.is_valid(dat, gvol, pvol):
                dat = self.im.np_multiply_sum(dat, gvol, axis=1)
                dat = self.im.operate_arrays(dat, '/', pvol)
                if name=='TE':
                    self.im.fill_value(sid.volume_average.t_e, dat, name='t_e')
                elif name=='NE':
                    self.im.fill_value(sid.volume_average.n_e, dat, name='n_e')
                elif name=='NI':
                    self.im.fill_value(sid.volume_average.n_i_total, dat, name='n_i')
                elif name == 'DNEDT':
                    self.im.fill_value(sid.volume_average.dn_e_dt, dat, name='dnedt')
                elif name=='ZEFFP':
                    self.im.fill_value(sid.volume_average.zeff, dat, name='zeffp')

        #t_e.  TE(TIME3,X); EV,  ELECTRON TEMPERATURE
        fill_volume_average('TE')

        #n_e average, NE(TIME3,X) N/CM**3
        fill_volume_average('NE')

        #NI(TIME3, X) N/CM**3
        fill_volume_average('NI')

        #n_i
        nsp = self.im.fill_ion_densities(sid.volume_average.n_i, pv=pvol, gv=gvol)

        #t_i_average
        #TI(TIME3, X) EV ION TEMPERATURE
        ntot = self.im.get_ncvar('TI')
        if self.im.is_valid(ntot, nsp, pvol):
            ntot /= nsp
            ntot = self.im.np_multiply_sum(ntot, gvol, axis=1)
            ntot = self.im.operate_arrays(ntot, '/', 1.0e-6*pvol)
            self.im.fill_value(sid.volume_average.t_i_average, ntot)

        #DNEDT(TIME3, X) N/CM3/SEC D/DT(ELECTRON DENSITY)
        fill_volume_average('DNEDT')

        #zeff
        #ZEFFP(TIME3, X) "PLASMA COMPOSITION ZEFF PROFILE
        fill_volume_average('ZEFFP')

        #meff_hydrogenic [amu]. Effective mass of the hydrogenic species
        #NH(TIME3, X) ;
        nh = self.im.get_ncvar('NH')
        nd = self.im.get_ncvar('ND')
        nt = self.im.get_ncvar('NT')
        umu2mev = scipy_phyc['atomic mass constant energy equivalent in MeV']
        #umu2mev is (931.49410242, 'MeV', 2.8e-07), last is uncertainty
        umu2mev = umu2mev[0]
        mh = atom_base('H').a * umu2mev # Mev
        md = atom_base('D').a * umu2mev # Mev
        mt = atom_base('T').a * umu2mev #
        if self.im.is_valid(nh, nd, nt, mh, md, mt, gvol):
            nhv = self.im.np_multiply_sum(nh, gvol, axis=1)
            nH = self.im.operate_arrays(nhv, '/', 1.0e-6*pvol)

            ndv = self.im.np_multiply_sum(nd, gvol, axis=1)
            nD = self.im.operate_arrays(ndv, '/', 1.0e-6*pvol)

            ntv = self.im.np_multiply_sum(nt, gvol, axis=1)
            nT = self.im.operate_arrays(ntv, '/', 1.0e-6*pvol)

            meffh = self.im.operate_arrays(mh*nH + md*nD + mt*nT, '/', nH+nD+nT)
            self.im.fill_value(sid.volume_average.meff_hydrogenic, meffh)

            #isotope_fraction_hydrogen. nH/(nH+nD+nT))
            fh = self.im.operate_arrays(nh, '/', nH+nD+nT)
            self.im.fill_value(sid.volume_average.isotope_fraction_hydrogen, fh)


        #MBTOR
        #TOTAL BETA TOROIDAL
        #BTTOT

        self.im.putdb_imas(sid)


if __name__ == '__main__':
    import argparse
    #a = get_args(sys.argv)
    ncimas = nc2imas(db='D3D', shot=147634, backend=imas.imasdef.MDSPLUS_BACKEND,
          run=1, occur=0, user=os.getcwd(), ncfile='147634G14.CDF', new_db=True)
    imas_summary(ncimas)
    im.close_imas_database()

