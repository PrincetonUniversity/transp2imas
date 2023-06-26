
from nc2imas import *

class imas_core_sources:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()

    def fill_element_data(self, ids, sp:str):
        self.im.fill_element_data(ids, sp)

    #name from ['T', '3', '4', 'P', 'P14']
    def fill_fusion_element_info(self, ids, name:str):
        if not name:
            return

        sp = 'HE'+name if name in ['3','4'] else name
        #sp is TRANSP symbol defined in atoms_class.entries
        atom_info = atom_base(sp)
        # no label in ids
        ids.element.resize(1)
        ids.element[0].a = atom_info.a
        ids.element[0].z_n = atom_info.z_n
        ids.element[0].atoms_n = 1

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
    def fill_nbi_sources(self, sci, ntimes:int, sp:str, sp:str):
        self.im.logger.info(f'core source: nbi')
        coord_x = self.im.get_ncvar('X')
        nprof = len(coord_x)

        #NBI global_quantities
        #global_quantities.power  coupled to plasma
        #BPTI (TIME) [WATTS] BEAM POWER TO IONS. =SUMINT(PBI)

        #BPTE (TIME) [WATTS] BEAM POWER TO ELECTRONS, =SUMINT(PBE)
        bpte_glob = self.im.get_ncvar('BPTE')
        #BPTH(TIME) WATTS; FAST ION POWER THERMALIZED ; BPTH_D
        bpth_glob = self.im.get_ncvar('BPTH')
        #total= BPTI + BPTE + BPTH (BHEAT)
        global_power_tot = bpti_glob + bpte_glob + bpth_glob
        #BPTH=intgeral(BPTH_D +...)  (BHEAT_D..)

        #BPOH(TIME); WATTS ; POWER: OH CIRCUIT TO FAST IONS


        #particles
        #BSTH (TIME) N/SEC FAST ION THERMALIZATIONS
        # = BSTH_D + ....
        #WNMCTOT_D
        total_ion_particles = self.im.get_ncvar('BSTH')

        gptcls, ptcls_prof, pbth_sp = 0., 0., 0.

        bpti_glob += self.im.get_ncvar('BPST_D')
        #BSTH_D (TIME) N/SEC D BEAM THERMALIZATIONS
        gptcls += self.im.get_ncvar(f'BSTH_{sp}').astype(float)

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
        pfth = self.im.get_ncvar(f'PFTH_{sp}', self.im.vol_inv)
        ion_energy = pbth + pfth

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(global_power_tot)>i:
                glob.power[i] = global_power_tot[i].astype(float)
                glob.total_ion_power[i] = bpti_glob[i].astype(float)
                glob.total_ion_particles[i] = total_ion_particles[i].astype(float)
                glob.electrons.particles[i] = gptcls[i]
                glob.electrons.power[i] = bpte_glob[i].astype(float)

            self.im.fill_data(profi.electrons.energy, el_energy, ind=i, chk=nprof, name='prof_el_energy')
            self.im.fill_data(profi.electrons.particles, ptcls_prof, ind=i, chk=nprof, name='prof_el_ptcls')

            #inside ???
            self.im.fill_data(profi.torque_tor_inside, torque, ind=i, chk=nprof, name='torque')
            self.im.fill_data(profi.total_ion_energy, ion_energy, ind=i, chk=nprof, name='tot_ion_en')

            self.im.fill_data(profi.j_parallel, j_parallel, ind=i, chk=nprof, name='j_parallel')

            profi.ion.resize(len(src_names))


            self.fill_fusion_element_info(profi, name)
            profi.ion[j].label = beams_dict[name]

            self.im.fill_data(profi.ion[j].energy, ion_energy, i, name='ion_energy')

            #BSTH_D (TIME) N/SEC D BEAM THERMALIZATIONS
            self.im.fill_data_nc(profi.ion[j].particles, f'BSTH_{sp}', self.im.vol_inv, i)
            #neutral_index


    #sci source_ids, ntimes:n_time_frames
    def fill_ic_sources(self, sci, ntimes:int, sp:str):
        self.im.logger.info(f'core source: IC')
        coord_x = self.im.get_ncvar('X')
        nprof = len(coord_x)
        #species
        self.fill_src_properties(sci, iname='ic', iind=5,
                    idesc="Sources from heating at the ion cyclotron range of frequencies",
                    tname='unspecified', tind=0, tdesc='ICRF antenna')
        #sp is TRANSP symbol defined in atoms_class.entries
        if sp:
            self.fill_element_data(sci.species.ion, sp)

        #PICHTOT(TIME) ; WATTS; TOTAL ICRF ANTENNA POWER
        #ant_power = self.im.get_ncvar('PICHTOT').astype(float)

        power_coupled = []
        #RFPWRAB1(TIME) ;WATTS; Power ABSORBED, Antenna #1
        ants =  self.im.get_var_names_re('RFPWRAB?', nd1=1, nd2=2)
        for i in range(len(ants)):
            name = ants[i]
            power_coupled += self.im.get_ncvar(name).astype(float)

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
        global_power_tot = np.sum(np.multiply(tot_heat_prof, gvol), axis=1)

        ion_power_glob = np.sum(np.multiply(ion_heat_prof, gvol), axis=1)

        el_power_glob = np.sum(np.multiply(el_heat_prof, gvol), axis=1)


        self.fill_species_neutral(sci, sp)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(power_coupled)>i:
                glob.power[i] = power_coupled[i]
                glob.total_ion_power[i] = ion_power_glob[i].astype(float)
                glob.electrons.power[i] = el_power_glob[i].astype(float)

            self.im.fill_data(profi.electrons.energy, el_heat_prof, ind=i, chk=nprof, name='prof_el_energy')
            #self.im.fill_data(profi.electrons.particles, el_particles_prof, ind=i, chk=nprof, name='prof_el_ptcls')

            #inside ???
            self.im.fill_data(profi.torque_tor_inside, torque, ind=i, chk=nprof, name='torque')
            #self.im.fill_data(profi.total_ion_energy, total_ion_energy, ind=i, chk=nprof, name='tot_ion_en')

            profi.ion.resize(1)
            self.im.fill_data(profi.ion[0].energy, tot_heat_prof, ind=i, chk=nprof, name='pbi')
            #self.im.fill_data(profi.ion[0].particles, ion_particles_prof, ind=i, chk=nprof, name='ion_ptcls')



    #sci source_ids, ntimes:n_time_frames
    def fill_ec_sources(self, sci, ntimes:int, sp:str):
        self.im.logger.info(f'core source: ec')
        coord_x = self.im.get_ncvar('X')
        nprof = len(coord_x)

        self.fill_src_properties(sci, iname='ec', iind=3,
                    idesc="Sources from electron cyclotron heating and current drive",
                    tname='unspecified', tind=0, tdesc='EC gyrotron')

        #sp is TRANSP symbol defined in atoms_class.entries
        if sp:
            self.fill_element_data(sci.species.ion, sp)

        #PECHT(TIME) ; WATTS; ECRF ELECTRON HEATING
        power_glob = self.im.get_ncvar('PECHT').astype(float)

        #PEECH(TIME3, X) WATTS/CM3; ECRH ELECTRON HEATING
        el_heat_prof = self.im.get_ncvar('PEECH', self.im.vol_inv)

        gvol = self.im.get_ncvar('DVOL', self.im.vconv)
        el_power_glob = np.sum(np.multiply(el_heat_prof, gvol), axis=1)

        #ECCUR(TIME3, X) ; AMPS/CM2; ECRH CURRENT
        curr = self.im.get_ncvar('ECCUR', self.im.area_inv)

        #
        self.fill_element_data(sci.species.neutral, sp)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(power_glob)>i:
                glob.power[i] = power_glob[i]
                glob.electrons.power[i] = el_power_glob[i].astype(float)

            self.im.fill_data(profi.electrons.energy, el_heat_prof, ind=i, chk=nprof, name='el_energy')
            self.im.fill_data(profi.j_parallel, curr, ind=i, chk=nprof, name='curr')


    #sci source_ids, ntimes:n_time_frames
    def fill_lh_sources(self, sci, ntimes:int, sp:str):
        self.im.logger.info(f'core source: lh')
        coord_x = self.im.get_ncvar('X')
        nprof = len(coord_x)

        self.fill_src_properties(sci, iname='lh', iind=4,
                    idesc="Sources from lower hybrid heating and current drive",
                    tname='unspecified', tind=0, tdesc='LH antenna')
        #sp is TRANSP symbol defined in atoms_class.entries
        if sp:
            self.fill_element_data(sci.species.ion, sp)

        #PELH(TIME3, X) WATTS/CM3;  LH ELECTRON HEATING
        el_heat_prof = self.im.get_ncvar('PELH', self.im.vol_inv)
        #PILH(TIME3, X) WATTS/CM3; LH ION HEATING
        ion_heat_prof = self.im.get_ncvar('PILH', self.im.vol_inv)

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)

        ion_power_glob = np.sum(np.multiply(ion_heat_prof, gvol), axis=1)

        el_power_glob = np.sum(np.multiply(el_heat_prof, gvol), axis=1)

        global_power_tot = ion_power_glob + el_power_glob

        #
        self.fill_element_data(sci.species.neutral, sp)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(global_power_tot)>i:
                glob.power[i] = global_power_tot[i].astype(float)
                glob.total_ion_power[i] = ion_power_glob[i].astype(float)
                glob.electrons.power[i] = el_power_glob[i].astype(float)

            self.im.fill_data(profi.electrons.energy, el_heat_prof, ind=i, chk=nprof, name='prof_el_energy')

            profi.ion.resize(1)
            self.im.fill_data(profi.ion[0].energy, ion_heat_prof, ind=i, chk=nprof, name='pbi')



    #sci source_ids, ntimes:n_time_frames
    def fill_fusion_sources(self, sci, ntimes:int, sp:str):
        self.im.logger.info(f'core source: fusion')
        coord_x = self.im.get_ncvar('X')
        nprof = len(coord_x)

        self.fill_src_properties(sci, iname='fusion', iind=6,
                    idesc="Sources from fusion reactions, e.g. alpha particle heating",
                    tname='unspecified', tind=0, tdesc='fusion')
        #sp is TRANSP symbol defined in atoms_class.entries
        if sp:
            self.fill_element_data(sci.species.ion, sp)

        #global
        #FPTI (TIME) W ; FUSION ION POWER TO IONS ; sum of BFTI_
        gpow_i = self.im.get_ncvar('FPTI').astype(float)
        #FPTE (TIME) W;  FUSION ION POWER TO ELECTRONS ; sum of BFTE_
        gpow_e = self.im.get_ncvar('FPTE').astype(float)

        #(P14 is proton(14.7MeV) fusion product fast ion)
        # ['T', '3', '4', 'P', 'P14']
        fions_dict = atom_base().fusion_ions
        src_names = self.im.get_var_names_lst('BFTH_?', fions_dict.keys())

        gpower, gtorq, prof_pow, gptcls, prof_powe, prof_powi = [], [], [], [], [], []
        fusion_power = []

        for i in range(len(src_names)):
            #selected from fions[..]
            name = src_names[i]
            #Total power coupled to the plasma.
            #BFTH_T  (TIME) W ; FUSN T POWER THERMALIZED ; PFION ?
            gpower += self.im.get_ncvar(f'BFTH_{name}').astype(float)

            #torque FPHTH_3(TIME) ; NT-M  ; FUSN HE3: THERMALIZATION ;  FPHXB_ ??
            gtorq += self.im.get_ncvar(f'FPHTH_{name}').astype(float)
            #FSTH_3(TIME)  N/SEC FUSN HE3 THERMALIZATIONS ; SFION ?
            gptcls += self.im.get_ncvar(f'FSTH_{name}').astype(float)
            #PFE_3(TIME3, X) ; WATTS/CM^3; FUSN HE3 ELECTRON HEATING
            prof_powe +=  self.im.get_ncvar(f'PFE_{name}', self.im.vol_inv)

            #PFI_3(TIME3, X) ; WATTS/CM; FUSN HE3 ION HEATING
            prof_powi += self.im.get_ncvar(f'PFI_{name}', self.im.vol_inv)

        #
        self.fill_element_data(sci.species.neutral, sp)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(gpower)>i:
                glob.power[i] = gpower[i]
            if len(gpow_i)>i:
                glob.total_ion_power[i] = gpow_i[i]
            if len(gptcls)>i:
                glob.total_ion_particles[i] = gptcls[i]

            #glob.electrons.particles[i] = 0
            if len(gpow_e)>i:
                glob.electrons.power[i] = gpow_i[i]

            if len(gtorq)>i:
                glob.torque[i] = gtorq[i]

            if len(prof_powe) >i:
                self.im.fill_data(profi.electrons.energy, prof_powe, ind=i, chk=nprof, name='prof_el_energy')

            #self.im.fill_data(profi.electrons.particles,

            self.im.fill_data(profi.total_ion_energy, prof_powi, ind=i, chk=nprof, name='tot_ion_en')

            profi.ion.resize(len(src_names))

            for j in range(len(src_names)):
                # existing from fions_dict.keys['T', '3', '4', 'P', 'P14']
                name = src_names[j]
                self.fill_fusion_element_info(profi, name)
                profi.ion[j].label = fions_dict[name]

                #PFTH_3 (TIME3, X) WATTS/CM3  FUSN HE3 THERMALIZATION POWER
                self.im.fill_data_nc(profi.ion[j].energy, f'PFTH_{name}', self.im.vol_inv, i)
                #SFTH_3(TIME3, X) N/CM3/SEC  FUSN HE3 THERMALIZATION SOURCE
                self.im.fill_data_nc(profi.ion[j].particles, f'SFTH_{name}', self.im.vol_inv, i)

                # FDENS_T(TIME3, X) ; N/CM**3; FUSN (T) ION DENSITY



    #sci source_ids, ntimes:n_time_frames
    def fill_ohmic_sources(self, sci, ntimes:int, sp:str):
        self.im.logger.info(f'Core source: ohmic')
        coord_x = self.im.get_ncvar('X')
        nprof = len(coord_x)

        self.fill_src_properties(sci, iname='ohmic', iind=7, idesc="Source from ohmic heating",
                    tname='unspecified', tind=0, tdesc='Ohmic Coil')
        #sp is TRANSP symbol defined in atoms_class.entries
        if sp:
            self.fill_element_data(sci.species.ion, sp)

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
        j_parallel = np.divide(jdotb, np.tile(bz[None,:].T, nprof))

        beams_dict = atom_base().beams
        beams = self.im.get_var_names_re('BPOH_?', beams_dict.keys())

        #BPOH(TIME); WATTS ; POWER: OH CIRCUIT TO FAST IONS
        power_fi_glob = self.im.get_ncvar('BPOH')

        #POHB(TIME3, X) W/CM3; POWER: OH CIRCUIT TO FAST IONS
        power_fi = self.im.get_ncvar('POHB', self.im.vol_inv)
        fusion_power, beam_power = [], []
        for i in range(len(beams)):
            #BFOH_3(TIME) ; WATTS ; POWER: OH CIRCUIT TO FUSN HE3
            fusion_power += self.im.get_ncvar(f'BFOH_{name}').astype(float)
            #BPOH_D(TIME) ; POWER: OH CIRCUIT TO D BEAM
            beam_power += self.im.get_ncvar(f'BPOH_{name}').astype(float)

        power_ions = fusion_power + beam_power + power_fi

        power_ions_glob = np.sum(np.multiply(power_ions, gvol), axis=1)

        #POH(TIME3, X) WATTS/CM3; OHMIC HEATING POWER
        power_prof = self.im.get_ncvar('POH', self.im.vol_inv).astype(float)
        power_glob = np.sum(np.multiply(power_prof, gvol), axis=1)

        #conductivity from resistivity ?

        #TQOHB(TIME3, X) ;Nt-M/CM3; TORQUE OH CIRCUIT --> BEAM


        self.fill_element_data(sci.species.neutral, sp)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(power_glob)>i:
                glob.power[i] = power_glob[i]
                glob.total_ion_power[i] = power_ions_glob[i]
                #glob.electrons.power[i] = el_power_glob[i].astype(float)

            #self.im.fill_data(profi.electrons.energy, el_heat_prof, ind=i, chk=nprof, name='prof_el_energy')

            profi.ion.resize(1)
            self.im.fill_data(profi.ion[0].j_parallel, j_parallel, ind=i, chk=nprof, name='j_par')
            self.im.fill_data(profi.ion[0].energy, power_ions, ind=i, chk=nprof, name='power_ions')




    #sci source_ids, ntimes:n_time_frames
    def fill_bremsstrahlung(self, sci, ntimes:int, sp:str):
        self.im.logger.info(f'Core source: bremsstrahlung')
        coord_x = self.im.get_ncvar('X')
        nprof = len(coord_x)

        self.fill_src_properties(sci, iname='bremsstrahlung', iind=8,
            idesc="Source from  bremsstrahlung; radiation losses are negative",
                    tname='electron', tind=1, tdesc='Electron')
        #sp is TRANSP symbol defined in atoms_class.entries
        if sp:
            self.fill_element_data(sci.species.ion, sp)

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)

        #PRAD_BR(TIME3, X) ; W/CM3; BREMSSTRAHLUNG RADIATION
        #negative for radiation power loss
        rad_prof = -1.0*self.im.get_ncvar('PLJBOH', self.im.vol_inv)
        rad_glob = np.sum(np.multiply(rad_prof, gvol), axis=1).astype(float)

        self.fill_element_data(sci.species.neutral, sp)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(rad_glob)>i:
                glob.power[i] = rad_glob[i]
                glob.electrons.power[i] = rad_glob[i]

            self.im.fill_data(profi.electrons.energy, rad_prof, ind=i, chk=nprof, name='el_energy')


    #sci source_ids, ntimes:n_time_frames
    def fill_line_radiation(self, sci, ntimes:int, sp:str):
        self.im.logger.info(f'Core source: line_radiation')
        coord_x = self.im.get_ncvar('X')
        nprof = len(coord_x)

        self.fill_src_properties(sci, iname='line_radiation', iind=10,
            idesc="Source from  line_radiation; radiation losses are negative sources",
                    tname='unspecified', tind=0, tdesc='unspecified')
        #sp is TRANSP symbol defined in atoms_class.entries
        if sp:
            self.fill_element_data(sci.species.ion, sp)

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)

        #PRAD_LI(TIME3, X) ; WATTS/CM3; LINE RADIATION
        #power loss due to line radiation is negative source.
        rad_prof = -1.0 * self.im.get_ncvar('PRAD_LI', self.im.vol_inv)
        rad_glob = np.sum(np.multiply(rad_prof, gvol), axis=1)

        self.fill_element_data(sci.species.neutral, sp)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(rad_glob)>i:
                glob.power[i] = rad_glob[i]

            #profi.ion.resize(1)
            #self.im.fill_data(profi.total_ion_power_inside, rad_glob, ind=i, chk=nprof, name='power_ions')
            #self.im.fill_data(profi.ion[0].energy, rad_glob, ind=i, chk=nprof, name='power_ions')

    def fill_collisional_equipartition(self, sci, ntimes:int, sp:str):

        #NUSTI/E(TIME3, X) ; ION/e- COLLISIONALITY
        #dimensionless. No place to fill this?
        return


    #sci source_ids, ntimes:n_time_frames
    def fill_runaways(self, sci, ntimes:int, sp:str):
        self.im.logger.info(f'Core source: runaways')
        coord_x = self.im.get_ncvar('X')
        nprof = len(coord_x)

        if sp != 'electron':
            self.im.logger.info(f'Core source: runaways for {sp} not available')
            return

        ind = 1 if sp=='electron' else 2 if sp=='ion' else -1
        self.fill_src_properties(sci, iname='runaways', iind=501,
            idesc="Source from run-away processes; includes both electron and ion run-away",
                    tname=sp, tind=ind, tdesc=sp)
        #sp is TRANSP symbol defined in atoms_class.entries
        if sp:
            self.fill_element_data(sci.species.ion, sp)

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)

        #SERUN(TIME3, X) ; N/CM3/SEC; RUNAWAY ELEC SOURCE RATE
        ptcl_prof = -1.0*self.im.get_ncvar('SERUN', self.im.vol_inv)
        ptcl_glob = np.sum(np.multiply(run_ptcl_prof, gvol), axis=1).astype(float)

        self.fill_element_data(sci.species.neutral, sp)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(ptcl_glob)>i:
                glob.electrons.particles[i] = ptcl_glob[i]

            self.im.fill_data(profi.electrons.particles, ptcl_glob, ind=i, chk=nprof, name='ep_run')


    #sci source_ids, ntimes:n_time_frames
    def fill_pellet(self, sci, ntimes:int, sp:str):
        self.im.logger.info(f'Core source: pellet')
        coord_x = self.im.get_ncvar('X')
        nprof = len(coord_x)

        self.fill_src_properties(sci, iname='pellet', iind=14, idesc="Sources from injection",
                    tname='electron', tind=1, tdesc='electron')
        #sp is TRANSP symbol defined in atoms_class.entries
        if sp:
            self.fill_element_data(sci.species.ion, sp)

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)

        #PLABD(TIME3, X) ; N/CM**3; PELLET ABLATION (DATA), MEASURED DELTA(NE) BEFORE/AFTER PELLET

        #Error: particles flux needed. PLABD is density difference
        self.im.logger.info(f'Error: Core source-pellet : density filled instead of particle flux')
        ptcl_prof = -1.0*self.im.get_ncvar('PLABD', self.im.vol_inv)
        ptcl_glob = np.sum(np.multiply(ptcl_prof, gvol), axis=1).astype(float)

        self.fill_element_data(sci.species.neutral, sp)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(ptcl_glob)>i:
                glob.electrons.particles[i] = ptcl_glob[i]

            self.im.fill_data(profi.electrons.particles, ptcl_prof, ind=i, chk=nprof, name='pell_ptcl')

    #sci source_ids, ntimes:n_time_frames
    def fill_gas_puff(self, sci, ntimes:int, sp:str):
        self.im.logger.info(f'Core source: gas_puff')
        coord_x = self.im.get_ncvar('X')
        nprof = len(coord_x)

        self.fill_src_properties(sci, iname='gas_puff', iind=108, idesc="Gas puff",
                    tname='neutral', tind=4, tdesc=r'Neutral species in a single/average state')
        #sp is TRANSP symbol defined in atoms_class.entries
        if sp:
            self.fill_element_data(sci.species.ion, sp)

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)
        #DN0OUT: NEUT. DENSITY OUTSIDE THE PLASMA
        #GASIN (TIME); ATOMS/SEC ; H+D+T GAS SOURCE AT BDY. RECYCD+GASFLD (Recycl+GASFLOW RATE DATA)






        #PLABD(TIME3, X) ; N/CM**3; PELLET ABLATION (DATA), MEASURED DELTA(NE) BEFORE/AFTER PELLET

        #Error: particles flux needed. PLABD is density difference
        self.im.logger.info(f'Error: Core source-pellet : density filled instead of particle flux')
        ptcl_prof = -1.0*self.im.get_ncvar('PLABD', self.im.vol_inv)
        ptcl_glob = np.sum(np.multiply(ptcl_prof, gvol), axis=1).astype(float)

        self.fill_element_data(sci.species.neutral, sp)

        sci.profiles_1d.resize(ntimes)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(sci)
        sci.global_quantities.resize(ntimes)

        for i in range(ntimes):
            glob = sci.global_quantities[i]
            profi = sci.profiles_1d[i]
            if len(ptcl_glob)>i:
                glob.electrons.particles[i] = ptcl_glob[i]

            self.im.fill_data(profi.electrons.particles, ptcl_prof, ind=i, chk=nprof, name='pell_ptcl')




    #sci source_ids, ntimes:n_time_frames
    def fill_ic_nbi_sources(self, sci, ntimes:int, sp:str):
        self.im.logger.info(f'core source: IC+NBI')

        #QICHFAST(TIME3, X). WATTS/CM3. ICH Heating of Beam & Fusn Ions
        heat_prof = self.im.get_ncvar('QICHFAST', self.im.vol_inv)




    def get_sources_info(self, sources:list)->dict:
        '''
        Return: dict  {str : [name, species]}
        '''
        srcs = {}
        #NBI
        #H,D,T,3,4,N,A,K,X,Z ; Z=Hi_Z(heavy), 3=He4, 4=He4, N=neon, A=argon, K=krypton, X=xenon
        #'3','4' is replaced by 'HE3', 'HE4' in the next line
        sp_look = ['H', 'D', 'T', '3', '4', 'N', 'A', 'K', 'X', 'Z', 'P', 'P14']
        #selcts species from sp_look if exists in BSTH_ and FSTH_
        species = self.im.get_var_names_lst('BSTH_?', sp_look)
        species2 = self.im.get_var_names_lst('FSTH_?', sp_look)
        species = species + species2
        #add HE to beam if '3' or '4'
        species = ['HE' + s if (s=='3' or s=='4') else s for s in species]
        name = 'nbi'
        #srcs dictionary of key, value=[name, sp]
        srcs = { name+str(i+1):[name, species[i]] for i in range(len(species))}
        return srcs


    #edge_sources
        #DN0OUT  Hydrogenic neutral density in scrape off region

    def create_ids(self):
        self.im.logger.info('\ncore_sources')

        csc = self.im.ids.core_sources

        src_lst = ['total', 'nbi', 'ec', 'lh', 'ic', 'fusion', 'ohmic', 'bremsstrahlung', 'ic_nbi']
        #src_dict{name_id : [name, species]}
        src_dict = get_sources_info(src_lst)
        csc.source.resize(len(src_dict))

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        timep = self.im.get_ncvar('TIME3')
        ntimep = len(timep)
        self.im.checkErrExit(timesc, timep)

        #fill source specific
        for name_id, srcs in src_dict.items():
          sci = csc.source[i]
          name, sp= srcs[0], srcs[1]
          if name == 'nbi':
              self.fill_nbi_sources(sci, ntimep, sp)
          if name == 'ec':
              self.fill_ec_sources(sci, ntimep, sp)
          if name == 'lh':
              self.fill_lh_sources(sci, ntimep, sp)
          if name == 'ic':
              self.fill_ic_sources(sci, ntimep, sp)
          if name == 'fusion':
              self.fill_fusion_sources(sci, ntimep, sp)
          if name == 'ohmic':
              self.fill_ohmic_sources(sci, ntimep, sp)
          if name == 'bremsstrahlung':
              self.fill_bremsstrahlung(sci, ntimep, sp)
          if name == 'line_radiation':
              self.fill_line_radiation(sci, ntimep, sp)
          if name == 'ic_nbi':
              self.fill_ic_nbi_sources(sci, ntimep, sp)






        coord_x = self.im.get_ncvar('X')
        coord_xb = self.im.get_ncvar('XB')

        csc = self.im.fill_properties(csc)
        self.im.fill_code_properties(csc)
        csc.ids_properties.homogeneous_time = 1
        csc.time.resize(ntimep)
        csc.time = timep

        #vacuum toroidal field
        self.im.fill_vac_tor_field(csc)

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)


        #EETR_OBS(TIME3, X) WATTS/CM3. Div(elec energy flux) observed
        #profiles_1d.electrons.energy W.m^-3. Source term for the electron energy equation
        een_flux_dens = self.im.get_ncvar('EETR_OBS', mult=self.im.vol_inv)
        #global_quantities.electrons.power (W)
        elpow = np.sum(np.multiply(een_flux_dens, gvol), axis=1)

        # EPTR_OBS(TIME3, X) N/CM3/SEC.  Div(electron flux)observed
        #profiles_1d.electrons.particles [m^-3.s^-1]. Source term for electron density equation
        eflux_dens = self.im.get_ncvar('EPTR_OBS', mult=self.im.vol_inv)
        #global_quantities.electrons.particles. Electron particle source [s^-1]
        elflux = np.sum(np.multiply(eflux_dens, gvol), axis=1)

        #UTOTL(TIME3, X) JLES/CM3 TOTAL ENERGY DENSITY

        #IETR_OBS(TIME3, X) WATTS/CM3. Div(ion energy flux) observed
        #profiles_1d.total_ion_energy
        ion_pow = self.im.get_ncvar('IETR_OBS', mult=self.im.vol_inv)
        #global_quantities.total_ion_power (W). Total power coupled to all ion species
        ion_pow_tot = np.sum(np.multiply(ion_pow, gvol), axis=1) ?

        # IPTR_OBS(TIME3, X) N/CM3/SEC.  Div(total ion flux)observed
        #profiles_1d.ion.particles [s^-1.m^-3]. Source term for ion density equation
        #WHICH ION ??
        iflux_dens = self.im.get_ncvar('IPTR_OBS', mult=self.im.vol_inv)
        #global_quantities.total_ion_particles. [(ions).s^-1]. Total ion particle source summed over ion species
        total_ion_particles = np.sum(np.multiply(iflux_dens, gvol), axis=1)




        nprof =

        #Total sources
        csc.sources[0].identifier.name = "total"
        csc.sources[0].identifier.index = 1
        csc.sources[0].identifier.description = "combines all sources"
        csc.sources[0].profiles_1d.resize(ntimes)
        for i in range(ntimes):
            glob =
            profi = csc.sources[0].profiles_1d[i]

            #electrons
            profi.electrons.particles.resize(nprof)
            #DNEDT(TIME3, X),  N/CM3/SEC,  D/DT(ELECTRON DENSITY)
            self.im.fill_data_nc(profi.electrons.particles, 'DNEDT', self.im.vol_inv, ind=i)  #???

            #Electron source inside flux surface, cumulative integral of source term
            ptcls = self.im.get_ncvar('DNEDT', self.im.vol_inv, ind=i)       #???
            if self.im.is_valid(ptcls, gvol):
                einside = np.cumsum(np.multiply(ptcls, gvol[i]))
            self.im.fill_data(profi.electrons.particles_inside, einside)

            #energy. Source term for electron energy equation [W.m^-3]
            #EETR_OBS, EETR_MOD(TIME3, X); WATTS/CM3, Div(elec energy flux) (observed, model)
            self.im.fill_data(profi.electrons.energy, 'EETR_OBS', self.im.vol_inv, ind=i)   #???

            #power_inside. Power coupled to e- inside flux surface.
            #Cumulative volume integral of source term for electron energy equation
            energy = self.im.get_ncvar('EETR_OBS', self.im.vol_inv, ind=i)    #???
            if self.im.is_valid(energy, gvol):
                pinside = np.cumsum(np.multiply(energy, gvol[i]))
            self.im.fill_data(profi.electrons.power_inside, pinside)


            #total_ion_energy
            #Source term for the total (summed over ion species) energy equation [W.m^-3]

            #total_ion_power_inside


            #momentum_tor


            # DNIDT(TIME3, X),  N/CM3/SEC,  D/DT(TOTAL ION DENSITY)
            profi.ions.resize(1)
            self.im.fill_data_nc(profi.ions[0].particles, 'DNIDT', self.im.vol_inv, ind=i)

            #IETR_OBS, IETR_MOD(TIME3, X) ; WATTS/CM3 Div(ion energy flux) (model)
            self.im.fill_data(profi.ions[0].energy, 'IETR_OBS', self.im.vol_inv, ind=i)

            #power_inside
            energy = self.im.get_ncvar('IETR_OBS', self.im.vol_inv, ind=i)
            if self.im.is_valid(energy, gvol):
                pinside = np.cumsum(np.multiply(energy, gvol[i]))
            self.im.fill_data(profi.ions[0].power_inside, pinside)




        self.im.putdb_imas(csc)
        self.logger.info('Done')


if __name__ == '__main__':
    print("Please use: python convert.py ")