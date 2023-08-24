
from nc2imas import *

class imas_core_sources:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.vol_inv = self.im.vol_inv
        self.area_inv = self.im.area_inv
        self.filled = {}
        self.create_ids()

    def fill_element_data(self, ids, sp:str):
        self.im.fill_element_data(ids, sp)

    def check_entry(self, name_id, sp):
        if name_id in self.filled:
            print(f'Error: Skipping Existing Entry {name_id} {sp}')
            return False
        self.filled[name_id] = sp
        return True

    def fill_src_properties(self, sci, iname:str, iind:int, idesc:str,
            tname:str, tind:int, tdesc:str):
        sci.identifier.name = iname
        sci.identifier.index = iind
        sci.identifier.description = idesc
        #species
        sci.species.type.name = tname
        sci.species.type.index = tind
        sci.species.type.description = tdesc


    #sci source_ids, ntimes:n_time_frames
    def fill_nbi_sources(self, sci, ntimes:int, sp:str):
        valid = self.im.check_var_name(f'PBTH_{sp}', length=True)
        if not valid:
            self.im.logger.info(f'core source: nbi species error: {sp}')
            return

        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]

        #NBI global_quantities
        #global_quantities.power  coupled to plasma
        #BPTI (TIME) [WATTS] BEAM POWER TO IONS. =SUMINT(PBI)

        #BPTE (TIME) [WATTS] BEAM POWER TO ELECTRONS, =SUMINT(PBE)
        bpte_glob = self.im.get_ncvar('BPTE')
        #BPTH(TIME) WATTS; FAST ION POWER THERMALIZED ; BPTH_D
        bpth_glob = self.im.get_ncvar('BPTH')
        #BPTH=intgeral(BPTH_D +...)  (BHEAT_D..)

        #BPOH(TIME); WATTS ; POWER: OH CIRCUIT TO FAST IONS


        #particles
        #BSTH (TIME) N/SEC FAST ION THERMALIZATIONS
        # = BSTH_D + ....
        #WNMCTOT_D
        total_ion_particles = self.im.get_ncvar('BSTH')

        gptcls, ptcls_prof, pbth_sp = 0., 0., 0.

        bpti_glob = self.im.get_ncvar('BPST_D')

        #total= BPTI + BPTE + BPTH (BHEAT)
        global_power_tot = self.im.sum_arrays(bpti_glob, bpte_glob, bpth_glob)

        #BSTH_D (TIME) N/SEC D BEAM THERMALIZATIONS
        gptcls += self.im.get_ncvar(f'BSTH_{sp}')

        #SBTH_D(TIME3, X) N/CM3/SEC; D BEAM THERMALIZATION SOURCE; SBTH
        ptcls_prof += self.im.get_ncvar(f'SBTH_{sp}')
        #same for e_ and ions ?????????

        # PBTH_D (TIME3, X) WATTS/CM3 "D BEAM THERMALIZATION POWER
        pbth_sp += self.im.get_ncvar('PBTH_'+sp, self.vol_inv)


        #profile
        #PBI(TIME3, X) W/CM3; BEAM HEATING OF IONS
        ion_energy = self.im.get_ncvar('PBI', self.vol_inv)

        #PBTH(TIME3, X) WATTS/CM3 "FAST ION THERMALIZATION POWER
        total_ion_energy = self.im.get_ncvar('PBTH', self.vol_inv)


        #PBE(TIME3, X) WATTS/CM3  BEAM HEATING OF ELECTRONS. PBE_D  PBE01_TOT ??
        el_energy = self.im.get_ncvar('PBE', self.vol_inv)

        #SINJ_D PINJ_D

        #torque TQTHNB BTQJXBS TQB01
        #TQBTH(TIME3, X)Nt-M/CM3 ; BEAM THERMALIZATION TORQUE ; TQTH{02..14} (TIME3, X)
        #or TQBCO : BEAM COLLISIONAL TORQUE   ?
        torque = self.im.get_ncvar('TQBTH', self.vol_inv)

        #CURB(TIME3, X) ; AMPS/CM2; BEAM DRIVEN CURRENT
        j_parallel = self.im.get_ncvar('CURB', self.area_inv)

        self.fill_src_properties(sci, iname='nbi', iind=2, idesc="Source from Neutral Beam Injection",
                    tname='neutral', tind=4, tdesc='Neutral species in a single/average state')

        self.fill_element_data(sci.species.neutral, sp)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        #PBTH_D (TIME3, X) WATTS/CM3 "D BEAM THERMALIZATION POWER
        pbth = self.im.get_ncvar(f'PBTH_{sp}', self.im.vol_inv)

        #TODO : fusion has separate source entry in core_sources ???
        #PFTH_sp (TIME3, X); FUSN sp THERMALIZATION POWER
        #pfth = self.im.get_ncvar(f'PFTH_{sp}', self.im.vol_inv)
        #ion_energy = self.im.operate_arrays(pbth, '+', pfth)
        ion_energy = pbth

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if self.im.is_valid(global_power_tot, ind=i):
                glob.power = global_power_tot[i].astype(float)
            if self.im.is_valid(bpti_glob, ind=i):
                glob.total_ion_power = bpti_glob[i].astype(float)
            if self.im.is_valid(total_ion_particles, ind=i):
                glob.total_ion_particles = total_ion_particles[i].astype(float)
            if self.im.is_valid(gptcls, ind=i):
                glob.electrons.particles = gptcls[i].astype(float)
            if self.im.is_valid(bpte_glob, ind=i):
                glob.electrons.power = bpte_glob[i].astype(float)

            self.im.fill_data(profi.electrons.energy, el_energy, ind=i, chk=nprof, name='prof_el_energy')
            self.im.fill_data(profi.electrons.particles, ptcls_prof, ind=i, chk=nprof, name='prof_el_ptcls')

            self.im.fill_data(profi.torque_tor_inside, torque, ind=i, chk=nprof, name='torque')
            self.im.fill_data(profi.total_ion_energy, ion_energy, ind=i, chk=nprof, name='tot_ion_en')

            self.im.fill_data(profi.j_parallel, j_parallel, ind=i, chk=nprof, name='j_parallel')

            profi.ion.resize(1)
            j = 0
            #name from ['T', '3', '4', 'P', 'P14']
            spn = 'HE'+sp if sp in ['3','4'] else sp
            self.im.fill_element_data(profi.ion[j], spn)

            self.im.fill_data(profi.ion[j].energy, ion_energy, i, name='ion_energy')

            #BSTH_D (TIME) N/SEC D BEAM THERMALIZATIONS
            #self.im.fill_data_nc(profi.ion[j].particles, f'BSTH_{sp}', self.im.vol_inv) #, i)
            #neutral_index


    #sci source_ids, ntimes:n_time_frames
    def fill_ic_sources(self, sci, ntimes:int):
        #ions and electrons are the resulting sources due to ICRF heating.
        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]
        #species
        self.fill_src_properties(sci, iname='ic', iind=5,
                    idesc="Sources from heating at the ion cyclotron range of frequencies",
                    tname='unspecified', tind=0, tdesc='N/A')

        #PICHTOT(TIME) ; WATTS; TOTAL ICRF ANTENNA POWER
        #ant_power = self.im.get_ncvar('PICHTOT')

        #RFPWRAB1(TIME) ;WATTS; Power ABSORBED, Antenna #1
        ants =  self.im.get_var_names_re('RFPWRAB?', nd1=1, nd2=2)
        power_coupled = self.im.get_ncvars_lst(ants)

        #PRFB_D(TIME) WATTS; RF POWER -> D BEAM IONS

        #PRFFI POWER TO FAST IONS
        #PICF01N01


        #QICHA(TIME3, X) WATTS/CM3; TOTAL ICH HEATING
        #QICHE(TIME3, X) WATTS/CM3; ICH DIRECT ELECTRON HEATING
        #QICHI(TIME3, X) WATTS/CM3; ICH DIRECT TH.ION HEATING
        #QRFBAL

        #PIICH(TIME3, X) WATTS/CM3 ICRF ION HEATING
        el_heat_prof = self.im.get_ncvar('PEICH', self.im.vol_inv) #QICHE
        #PEICH(TIME3, X) WATTS/CM3 ICRF E- HEATING
        ion_heat_prof = self.im.get_ncvar('PIICH', self.im.vol_inv) #QICHI
        #tot_heat_prof = self.im.get_ncvar('QICHA', self.im.vol_inv)
        tot_heat_prof = ion_heat_prof + el_heat_prof

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)
        global_power_tot = self.im.np_multiply_sum(tot_heat_prof, gvol, axis=1)

        ion_power_glob = self.im.np_multiply_sum(ion_heat_prof, gvol, axis=1)

        el_power_glob = self.im.np_multiply_sum(el_heat_prof, gvol, axis=1)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(power_coupled)>i:
                glob.power = power_coupled[i].astype(float)
            if len(ion_power_glob)>i:
                glob.total_ion_power = ion_power_glob[i].astype(float)
            if len(el_power_glob)>i:
                glob.electrons.power = el_power_glob[i].astype(float)

            self.im.fill_data(profi.electrons.energy, el_heat_prof, ind=i, chk=nprof, name='prof_el_energy')
            #self.im.fill_data(profi.electrons.particles, el_particles_prof, ind=i, chk=nprof, name='prof_el_ptcls')

            #inside ???
            #self.im.fill_data(profi.torque_tor_inside, torque, ind=i, chk=nprof, name='torque')
            #self.im.fill_data(profi.total_ion_energy, total_ion_energy, ind=i, chk=nprof, name='tot_ion_en')

            profi.ion.resize(1)
            self.im.fill_data(profi.ion[0].energy, tot_heat_prof, ind=i, chk=nprof, name='pbi')
            #self.im.fill_data(profi.ion[0].particles, ion_particles_prof, ind=i, chk=nprof, name='ion_ptcls')



    #sci source_ids, ntimes:n_time_frames
    def fill_ec_sources(self, sci, ntimes:int):
        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]

        self.fill_src_properties(sci, iname='ec', iind=3,
                    idesc="Sources from electron cyclotron heating and current drive",
                    tname='unspecified', tind=0, tdesc='N/A')

        #PECHT(TIME) ; WATTS; ECRF ELECTRON HEATING
        power_glob = self.im.get_ncvar('PECHT')

        #PEECH(TIME3, X) WATTS/CM3; ECRH ELECTRON HEATING
        el_heat_prof = self.im.get_ncvar('PEECH', self.im.vol_inv)

        gvol = self.im.get_ncvar('DVOL', self.im.vconv)
        el_power_glob = self.im.np_multiply_sum(el_heat_prof, gvol, axis=1)

        #ECCUR(TIME3, X) ; AMPS/CM2; ECRH CURRENT
        curr = self.im.get_ncvar('ECCUR', self.im.area_inv)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(power_glob)>i:
                glob.power = power_glob[i].astype(float)
                glob.electrons.power = el_power_glob[i].astype(float)

            self.im.fill_data(profi.electrons.energy, el_heat_prof, ind=i, chk=nprof, name='el_energy')
            self.im.fill_data(profi.j_parallel, curr, ind=i, chk=nprof, name='curr')


    #sci source_ids, ntimes:n_time_frames
    def fill_lh_sources(self, sci, ntimes:int):
        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]

        self.fill_src_properties(sci, iname='lh', iind=4,
                    idesc="Sources from lower hybrid heating and current drive",
                    tname='unspecified', tind=0, tdesc='N/A')

        #PELH(TIME3, X) WATTS/CM3;  LH ELECTRON HEATING
        el_heat_prof = self.im.get_ncvar('PELH', self.im.vol_inv)
        #PILH(TIME3, X) WATTS/CM3; LH ION HEATING
        ion_heat_prof = self.im.get_ncvar('PILH', self.im.vol_inv)

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)

        ion_power_glob = self.im.np_multiply_sum(ion_heat_prof, gvol, 1, 'ion_heat_prof', 'gvol')
        el_power_glob = self.im.np_multiply_sum(el_heat_prof, gvol, 1, 'el_heat_prof', 'gvol')
        global_power_tot = self.im.operate_arrays(ion_power_glob, '+', el_power_glob,
                                                  'ion_power_glob', 'el_power_glob')

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if self.im.is_valid(global_power_tot, ion_power_glob, el_power_glob, ind=i, name='gptot'):
                glob.power = global_power_tot[i].astype(float)
                glob.total_ion_power = ion_power_glob[i].astype(float)
                glob.electrons.power = el_power_glob[i].astype(float)

            self.im.fill_data(profi.electrons.energy, el_heat_prof, ind=i, chk=nprof, name='prof_el_energy')

            profi.ion.resize(1)
            self.im.fill_data(profi.ion[0].energy, ion_heat_prof, ind=i, chk=nprof, name='pbi')



    #sci source_ids, ntimes:n_time_frames
    def fill_fusion_sources(self, sci, ntimes:int):
        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]

        self.fill_src_properties(sci, iname='fusion', iind=6,
                    idesc="Sources from fusion reactions, e.g. alpha particle heating",
                    tname='unspecified', tind=0, tdesc='N/A')

        #self.fill_element_data(sci.species.ion, sp)

        #global
        #FPTI (TIME) W ; FUSION ION POWER TO IONS ; sum of BFTI_
        gpow_i = self.im.get_ncvar('FPTI')
        #FPTE (TIME) W;  FUSION ION POWER TO ELECTRONS ; sum of BFTE_
        gpow_e = self.im.get_ncvar('FPTE')

        #(P14 is proton(14.7MeV) fusion product fast ion)
        # ['T', '3', '4', 'P', 'P14']; 3=>HE3, 4=>HE4
        sp_ions = self.im.get_species_names('BFTH_?', cut='BFTH_')

        #Total power coupled to the plasma.
        #BFTH_T  (TIME) W ; FUSN T POWER THERMALIZED ; PFION ?
        gpower = self.im.get_ncvar_sum(f'BFTH_?')

        #torque FPHTH_3(TIME) ; NT-M  ; FUSN HE3: THERMALIZATION ;  FPHXB_ ??
        gtorq = self.im.get_ncvar_sum(f'FPHTH_?')
        #FSTH_3(TIME)  N/SEC FUSN HE3 THERMALIZATIONS ; SFION ?
        gptcls = self.im.get_ncvar_sum(f'FSTH_?')
        #PFE_3(TIME3, X) ; WATTS/CM^3; FUSN HE3 ELECTRON HEATING
        prof_powe =  self.im.get_ncvar_sum(f'PFE_?', scale=self.im.vol_inv)

        #PFI_3(TIME3, X) ; WATTS/CM; FUSN HE3 ION HEATING
        prof_powi = self.im.get_ncvar_sum(f'PFI_?', scale=self.im.vol_inv)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(gpower)>i:
                glob.power = gpower[i].astype(float)
            if len(gpow_i)>i:
                glob.total_ion_power = gpow_i[i].astype(float)
            if len(gptcls)>i:
                glob.total_ion_particles = gptcls[i].astype(float)

            #glob.electrons.particles = 0
            if len(gpow_e)>i:
                glob.electrons.power = gpow_i[i].astype(float)

            if len(gtorq)>i:
                glob.torque_tor = gtorq[i].astype(float)

            self.im.fill_data(profi.electrons.energy, prof_powe, ind=i, chk=nprof, name='prof_el_energy')

            #self.im.fill_data(profi.electrons.particles,

            self.im.fill_data(profi.total_ion_energy, prof_powi, ind=i, chk=nprof, name='tot_ion_en')

            profi.ion.resize(len(sp_ions))
            for j in range(len(sp_ions)):
                profion = profi.ion[j]
                name = sp_ions[j]
                self.fill_element_data(profion, name)
                #PFTH_P (TIME3, X) WATTS/CM3  FUSN (P) THERMALIZATION POWER
                self.im.fill_data_nc(profion.energy, f'PFTH_{name}', self.im.vol_inv, i)
                #SFTH_P(TIME3, X) N/CM3/SEC  FUSN (P) THERMALIZATION SOURCE
                self.im.fill_data_nc(profion.particles, f'SFTH_{name}', self.im.vol_inv, i)

                # FDENS_T(TIME3, X) ; N/CM**3; FUSN (T) ION DENSITY


    #sci source_ids, ntimes:n_time_frames
    def fill_ohmic_sources(self, sci, ntimes:int):
        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]

        self.fill_src_properties(sci, iname='ohmic', iind=7, idesc="Source from ohmic heating",
                    tname='unspecified', tind=0, tdesc='N/A')

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)
        #BZ(TIME) ; Tesla; BZ @R=RMAJOR OUTSIDE PLASMA
        bz = self.im.get_ncvar('BZ')

        #PLJBXTR ! <J.B>ext from solution of Ohms law equation
        #PLJBXTR(TIME3, X) ; AMP*TESLA/CM2; <J.B> DRIVEN, FROM RESISTIVITY

        #parallel current. average(J.B)/B0 ; B0=vacuum_toroidal_field/b0
        #PLJBOH(TIME3, X) ;AMP*TESLA/CM2 ; <J.B> OHMIC
        jdotb = self.im.get_ncvar('PLJBOH', self.im.area_inv)

        #zeros in BZ ?
        j_parallel = []
        try:
            j_parallel = np.divide(jdotb, np.tile(bz[None,:].T, nprof))
        except:
            pass
        #beams_dict = atom_base().beams
        #beams = self.im.get_var_names_re('BPOH_?', beams_dict.keys())

        #BPOH(TIME); WATTS ; POWER: OH CIRCUIT TO FAST IONS
        power_fi_glob = self.im.get_ncvar('BPOH')

        #POHB(TIME3, X) W/CM3; POWER: OH CIRCUIT TO FAST IONS
        power_fi = self.im.get_ncvar('POHB', self.im.vol_inv)

        #BFOH_3(TIME) ; WATTS ; POWER: OH CIRCUIT TO FUSN HE3
        fusion_power = self.im.get_ncvar_sum(f'BFOH_?')
        #BPOH_D(TIME) ; POWER: OH CIRCUIT TO D BEAM
        beam_power = self.im.get_ncvar_sum(f'BPOH_?')

        self.im.logger.info(f'{fusion_power.shape}, {beam_power.shape}, {power_fi_glob.shape}')
        power_ions = self.im.sum_arrays(fusion_power, beam_power, power_fi_glob)

        #power_ions_glob = self.im.np_multiply_sum(power_ions, gvol, 1, 'power_ions', 'gvol')

        #POH(TIME3, X) WATTS/CM3; OHMIC HEATING POWER
        power_prof = self.im.get_ncvar('POH', self.im.vol_inv)
        power_glob = self.im.np_multiply_sum(power_prof, gvol, 1, 'power_prof', 'gvol')

        #conductivity from resistivity ?

        #TQOHB(TIME3, X) ;Nt-M/CM3; TORQUE OH CIRCUIT --> BEAM

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if self.im.is_valid(power_glob, power_ions):
                glob.power = power_glob[i].astype(float)
                glob.total_ion_power = power_ions[i].astype(float)
                #glob.electrons.power = el_power_glob[i].astype(float)

            #self.im.fill_data(profi.electrons.energy, el_heat_prof, ind=i, chk=nprof, name='prof_el_energy')

            self.im.fill_data(profi.j_parallel, j_parallel, ind=i, chk=nprof, name='j_par')
            profi.ion.resize(1)
            self.im.fill_data(profi.ion[0].energy, power_fi, ind=i, chk=nprof, name='power_ions')


    #sci source_ids, ntimes:n_time_frames
    def fill_bremsstrahlung_sources(self, sci, ntimes:int):
        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]

        self.fill_src_properties(sci, iname='bremsstrahlung', iind=8,
            idesc="Source from  bremsstrahlung; radiation losses are negative",
                    tname='electron', tind=1, tdesc='N/A')

        #self.fill_element_data(sci.species.ion, sp)

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)

        #PRAD_BR(TIME3, X) ; W/CM3; BREMSSTRAHLUNG RADIATION
        #negative for radiation power loss
        rad_prof = -1.0*self.im.get_ncvar('PLJBOH', self.im.vol_inv)
        rad_glob = self.im.np_multiply_sum(rad_prof, gvol, axis=1)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(rad_glob)>i:
                glob.power = rad_glob[i].astype(float)
                glob.electrons.power = rad_glob[i].astype(float)

            self.im.fill_data(profi.electrons.energy, rad_prof, ind=i, chk=nprof, name='el_energy')



    #sci source_ids, ntimes:n_time_frames
    def fill_impurity_radiation_sources(self, sci, ntimes:int):
        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]

        self.fill_src_properties(sci, iname='impurity_radiation', iind=203,
            idesc="Line radiation and Bremsstrahlung source; radiation losses are negative",
                    tname='unspecified', tind=0, tdesc='N/A')

        #self.fill_element_data(sci.species.ion, sp)

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)

        # PRBS_TOK(TIME3, X) ; W/CM3; TOK Impurity Brem Radiation
        #negative for radiation power loss
        rad_prof = -1.0*self.im.get_ncvar('PRBS_TOK', self.im.vol_inv)
        rad_glob = self.im.np_multiply_sum(rad_prof, gvol, axis=1)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(rad_glob)>i:
                glob.power = rad_glob[i].astype(float)
                glob.total_ion_power = rad_glob[i].astype(float)

            profi.ion.resize(1)
            self.im.fill_data(profi.ion[0].energy, rad_prof, ind=i, chk=nprof, name='power_ions')



    #sci source_ids, ntimes:n_time_frames
    def fill_line_radiation_sources(self, sci, ntimes:int):
        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]

        self.fill_src_properties(sci, iname='line_radiation', iind=10,
            idesc="Source from  line_radiation; radiation losses are negative sources",
                    tname='unspecified', tind=0, tdesc='N/A')

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)

        #PRAD_LI(TIME3, X) ; WATTS/CM3; LINE RADIATION
        #power loss due to line radiation is negative source.
        rad_prof = -1.0 * self.im.get_ncvar('PRAD_LI', self.im.vol_inv)
        rad_glob = self.im.np_multiply_sum(rad_prof, gvol, axis=1)

        # PRLS_TOK(TIME3, X) ?
        # PRADS_TOK(TIME3, X) ?

        #self.fill_element_data(sci.species.ion, sp)
        #self.fill_element_data(sci.species.neutral, sp)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(rad_glob)>i:
                glob.power = rad_glob[i].astype(float)

            #profi.ion.resize(1)
            #self.im.fill_data(profi.total_ion_power_inside, rad_glob, ind=i, chk=nprof, name='power_ions')
            #self.im.fill_data(profi.ion[0].energy, rad_glob, ind=i, chk=nprof, name='power_ions')

    def fill_collisional_equipartition_sources(self, sci, ntimes:int):

        #NUSTI/E(TIME3, X) ; ION/e- COLLISIONALITY
        #dimensionless. No place to fill this?
        return


    #sci source_ids, ntimes:n_time_frames
    def fill_runaways_sources(self, sci, ntimes:int):
        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]

        sp = 'electron'
        ind = 1 if sp=='electron' else 2 if sp=='ion' else -1
        self.fill_src_properties(sci, iname='runaways', iind=501,
            idesc="Source from run-away processes; includes both electron and ion run-away",
                    tname=sp, tind=ind, tdesc=sp)

        #self.fill_element_data(sci.species.ion, sp)

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)

        #SERUN(TIME3, X) ; N/CM3/SEC; RUNAWAY ELEC SOURCE RATE
        ptcl_prof = -1.0*self.im.get_ncvar('SERUN', self.im.vol_inv)
        ptcl_glob = self.im.np_multiply_sum(ptcl_prof, gvol, 1, 'ptcl_prof', 'gvol')

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(ptcl_glob)>i:
                glob.electrons.particles = ptcl_glob[i].astype(float)

            self.im.fill_data(profi.electrons.particles, ptcl_glob, ind=i, chk=nprof, name='ep_run')


    #sci source_ids, ntimes:n_time_frames
    def fill_pellet_sources(self, sci, ntimes:int):
        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]

        self.fill_src_properties(sci, iname='pellet', iind=14, idesc="Sources from injection",
                    tname='ion', tind=2, tdesc='Ion species in a single/average state')
        #sp = ''
        #self.fill_element_data(sci.species.ion, sp)

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)

        #PLABD(TIME3, X) ; N/CM**3; PELLET ABLATION (DATA), MEASURED DELTA(NE) BEFORE/AFTER PELLET

        #Error: particles flux needed. PLABD is density difference
        self.im.logger.info(f'Error: Core source-pellet : density filled instead of particle flux')
        ptcl_prof = -1.0*self.im.get_ncvar('PLABD', self.im.vol_inv)
        ptcl_glob = self.im.np_multiply_sum(ptcl_prof, gvol, axis=1)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(ptcl_glob)>i:
                glob.electrons.particles = ptcl_glob[i].astype(float)

            self.im.fill_data(profi.electrons.particles, ptcl_prof, ind=i, chk=nprof, name='pell_ptcl')


    #sci source_ids, ntimes:n_time_frames
    def fill_gas_puff_sources(self, sci, ntimes:int):
        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]

        self.fill_src_properties(sci, iname='gas_puff', iind=108, idesc="Gas puff",
                    tname='neutral', tind=4, tdesc='Neutral species in a single/average state')
        sp = ''
        #self.fill_element_data(sci.species.ion, sp)

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)
        #DN0OUT: NEUT. DENSITY OUTSIDE THE PLASMA
        #GASIN (TIME); ATOMS/SEC ; H+D+T GAS SOURCE AT BDY. RECYCD+GASFLD (Recycl+GASFLOW RATE DATA)






        #PLABD(TIME3, X) ; N/CM**3; PELLET ABLATION (DATA), MEASURED DELTA(NE) BEFORE/AFTER PELLET

        #Error: particles flux needed. PLABD is density difference
        self.im.logger.info(f'Error: Core source-pellet : density filled instead of particle flux')
        ptcl_prof = -1.0*self.im.get_ncvar('PLABD', self.im.vol_inv)
        ptcl_glob = self.im.np_multiply_sum(ptcl_prof, gvol, axis=1)

        #self.fill_element_data(sci.species.neutral, sp)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(ptcl_glob)>i:
                glob.electrons.particles = ptcl_glob[i].astype(float)

            self.im.fill_data(profi.electrons.particles, ptcl_prof, ind=i, chk=nprof, name='pell_ptcl')




    #sci source_ids, ntimes:n_time_frames
    def fill_ic_nbi_sources(self, sci, ntimes:int):
        #QICHFAST(TIME3, X). WATTS/CM3. ICH Heating of Beam & Fusn Ions
        heat_prof = self.im.get_ncvar('QICHFAST', self.im.vol_inv)


    def fill_total_sources(self, sci, ntimes:int):
        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]
        coord_xb = self.im.get_ncvar('XB')

        self.fill_src_properties(sci, iname='total', iind=1, idesc="Total source; combines all sources",
                    tname='unspecified', tind=0, tdesc='unspecified')

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)


        #EETR_OBS(TIME3, X) WATTS/CM3. Div(elec energy flux) observed
        #profiles_1d.electrons.energy W.m^-3. Source term for the electron energy equation
        een_flux_dens = self.im.get_ncvar('EETR_OBS', scale=self.im.vol_inv)
        #global_quantities.electrons.power (W)
        elpow = self.im.np_multiply_sum(een_flux_dens, gvol, axis=1)

        # EPTR_OBS(TIME3, X) N/CM3/SEC.  Div(electron flux)observed
        #profiles_1d.electrons.particles [m^-3.s^-1]. Source term for electron density equation
        eflux_dens = self.im.get_ncvar('EPTR_OBS', scale=self.im.vol_inv)
        #global_quantities.electrons.particles. Electron particle source [s^-1]
        elflux = self.im.np_multiply_sum(eflux_dens, gvol, axis=1)

        #UTOTL(TIME3, X) JLES/CM3 TOTAL ENERGY DENSITY

        #IETR_OBS(TIME3, X) WATTS/CM3. Div(ion energy flux) observed
        #profiles_1d.total_ion_energy
        ion_pow = self.im.get_ncvar('IETR_OBS', scale=self.im.vol_inv)
        #global_quantities.total_ion_power (W). Total power coupled to all ion species
        ion_pow_tot = self.im.np_multiply_sum(ion_pow, gvol, axis=1) #TODO

        # IPTR_OBS(TIME3, X) N/CM3/SEC.  Div(total ion flux)observed
        #profiles_1d.ion.particles [s^-1.m^-3]. Source term for ion density equation
        #WHICH ION ??
        iflux_dens = self.im.get_ncvar('IPTR_OBS', scale=self.im.vol_inv)
        #global_quantities.total_ion_particles. [(ions).s^-1]. Total ion particle source summed over ion species
        total_ion_particles = self.im.np_multiply_sum(iflux_dens, gvol, axis=1)

        nprof = coord_x.shape[1]

        sci.profiles_1d.resize(ntimes)
        self.im.fill_profiles1d_grid(sci)
        for i in range(ntimes):
            #glob =
            profi = sci.profiles_1d[i]

            #electrons
            profi.electrons.particles.resize(nprof)
            #DNEDT(TIME3, X),  N/CM3/SEC,  D/DT(ELECTRON DENSITY)
            self.im.fill_data_nc(profi.electrons.particles, 'DNEDT', self.im.vol_inv, ind=i)  #???

            #Electron source inside flux surface, cumulative integral of source term
            ptcls = self.im.get_ncvar('DNEDT', self.im.vol_inv, ind=i)       #???
            try:
                einside = np.cumsum(np.multiply(ptcls, gvol[i]))
                self.im.fill_data(profi.electrons.particles_inside, einside)
            except:
                pass

            #energy. Source term for electron energy equation [W.m^-3]
            #EETR_OBS, EETR_MOD(TIME3, X); WATTS/CM3, Div(elec energy flux) (observed, model)
            self.im.fill_data(profi.electrons.energy, 'EETR_OBS', self.im.vol_inv, ind=i)   #???

            #power_inside. Power coupled to e- inside flux surface.
            #Cumulative volume integral of source term for electron energy equation
            energy = self.im.get_ncvar('EETR_OBS', self.im.vol_inv, ind=i)    #???
            try:
                pinside = np.cumsum(np.multiply(energy, gvol[i]))
                self.im.fill_data(profi.electrons.power_inside, pinside)
            except:
                pass


            #total_ion_energy
            #Source term for the total (summed over ion species) energy equation [W.m^-3]

            #total_ion_power_inside


            #momentum_tor


            # DNIDT(TIME3, X),  N/CM3/SEC,  D/DT(TOTAL ION DENSITY)
            profi.ion.resize(1)
            self.im.fill_data_nc(profi.ion[0].particles, 'DNIDT', self.im.vol_inv, ind=i)

            #IETR_OBS, IETR_MOD(TIME3, X) ; WATTS/CM3 Div(ion energy flux) (model)
            self.im.fill_data(profi.ion[0].energy, 'IETR_OBS', self.im.vol_inv, ind=i)

            #power_inside
            energy = self.im.get_ncvar('IETR_OBS', self.im.vol_inv, ind=i)
            try:
                pinside = np.cumsum(np.multiply(energy, gvol[i]))
                self.im.fill_data(profi.total_ion_power_inside, pinside)
            except:
                pass



    def get_sources_info(self, sources:list)->dict:
        '''
        Return: dict  {str : [name, species]}
        '''
        srcs = {}
        for name in sources:
            if name == 'nbi':
                species1 = self.im.get_var_names('BSTH_?', cut='BSTH_')
                species2 = self.im.get_var_names('FSTH_?', cut='FSTH_')
                species = species1 + species2
                #add HE to beam if '3' or '4'
                species = ['HE' + s if (s=='3' or s=='4') else s for s in species]
                name = 'nbi'
                #srcs dictionary of key, value=[name, sp]
                src_ = { name+str(i+1):[name, species[i]] for i in range(len(species))}
            else:
                src_ = {name:[name, '']}

            srcs = {**srcs, **src_}
        return srcs


    #edge_sources
        #DN0OUT  Hydrogenic neutral density in scrape off region

    def create_ids(self):
        self.im.logger.info('\ncore_sources')

        csc = self.im.ids.core_sources
        #exact imas source term names in core_sources
        src_lst = ['nbi', 'ec', 'lh', 'ic', 'fusion', 'ohmic', 'bremsstrahlung', 'ic_nbi',
                  'line_radiation', 'impurity_radiation', 'collisional_equipartition',
                  'runaways', 'pellet', 'gas_puff', 'total']
        #src_dict{name_id : [name, species]}
        src_dict = self.get_sources_info(src_lst)
        csc.source.resize(len(src_dict))

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        timep = self.im.get_ncvar('TIME3')
        ntimep = len(timep)
        self.im.checkErrExit(timesc, timep)
        nfill = 0

        #fill source specific.
        for i, (name_id, srcs) in enumerate(src_dict.items()):
          name, sp = srcs[0], srcs[1]
          if not self.check_entry(name_id, sp):
              continue

          sci = csc.source[i]
          fname = f'fill_{name}_sources'
          #all functions fill_<source_name>_sources() in this class called
          if hasattr(self, fname):
              self.im.logger.info(f'\t {i+1} {fname} {sp}')
              funct = getattr(self, fname)
              if name == 'nbi':
                  funct(sci, ntimes, sp)
              else:
                  funct(sci, ntimes)
              nfill += 1

        if nfill > 0:
            csc = self.im.fill_properties(csc)
            self.im.fill_code_properties(csc)
            csc.ids_properties.homogeneous_time = 1
            csc.time.resize(ntimep)
            csc.time = timep
            self.im.putdb_imas(csc)
        self.logger.info('Done')


if __name__ == '__main__':
    print("Please use: python convert.py ")
