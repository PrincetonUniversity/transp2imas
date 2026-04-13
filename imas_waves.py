from nc2imas import *

class imas_waves:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\nwaves')
        ids = self.im.ids.get('waves')
        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        timesp = self.im.get_ncvar('TIME3')
        ntimep = len(timesc)
        self.im.checkErrExit(timesc, timesp)

        atomsDef = atom_base()
        beams = atomsDef.beams

        #ion species
        nspec = 0 #ng+nsbeam+nrhis+nmini

        gvol = self.im.get_ncvar('DVOL', self.im.vconv)
        ants = self.im.get_var_names_re('PICHA?', cut='PICHA', nd1=1, nd2=2)
        nicha = len(ants)
        self.logger.info(f'nich_antenna {nicha}')

        ## ????? not per antenna ??
        #QRFBTOT(TIME3, X) ; WATTS/CM3 ; ALL FAST ION HEATING by RF;
        fpdth_prof  = self.im.get_ncvar('QRFBTOT', self.im.vol_inv)
        #PICHFAST ?
        fpdth_glob = self.im.np_multiply_sum(fpdth_prof, gvol, axis=1)
        #PIICH(TIME3, X) ;  WATTS/CM3 ; ICRF ION HEATING
        ipdth_prof  = self.im.get_ncvar('PIICH', self.im.vol_inv)
        #PICHI  ?
        ipdth_glob = self.im.np_multiply_sum(ipdth_prof, gvol, axis=1)
        #PEICH(TIME3, X) ; WATTS/CM3 ; ICRF ELECTRON HEATING
        epdth_prof  = self.im.get_ncvar('PEICH', self.im.vol_inv)
        #PICHE?
        epdth_glob = self.im.np_multiply_sum(epdth_prof, gvol, axis=1)
        #total profile
        pow_dens = self.im.sum_arrays(fpdth_prof, ipdth_prof, epdth_prof)



        ids.coherent_wave.resize(nicha)
        for iw,astr in enumerate(ants):
            wav = ids.coherent_wave[iw]
            wav.identifier.type.name = 'IC'
            wav.identifier.type.index = 3
            wav.identifier.type.description = 'Wave field for IC heating and current drive'
            wav.wave_solver_type.name = 'TORIC'
            wav.wave_solver_type.index = 2
            wav.wave_solver_type.description = "full wave"

            freq = self.im.get_ncvar(f'FREQA{astr}')

            #RFPWRAB1(TIME) WATTS; Power ABSORBED, Antenna #1. PIC_F1 ??
            gpower = self.im.get_ncvar(f'RFPWRAB{astr}')

            #global
            wav.global_quantities.resize(ntimes)
            for i in range(ntimes):
                glob = wav.global_quantities[i]
                if self.im.is_valid(freq, ind=i):
                    glob.frequency = freq[i].astype(float)
                if self.im.is_valid(gpower, ind=i):
                    glob.power = gpower[i].astype(float)

                #XGRID_NPHI(TIME3, XGRID_NPHI); (numbers) ;Nphi Grid
                self.im.fill_data_nc(glob.n_tor, 'XGRID_NPHI', ind=i)

                if False: #data is for all antennas
                    #electrons
                    if self.im.is_valid(epdth_glob, ind=i):
                        glob.electrons.power_thermal = epdth_glob[i].astype(float)
                    if self.im.is_valid(fpdth_glob, ind=i):
                        glob.electrons.power_fast = fpdth_glob[i].astype(float)
                    #ions
                    glob.ion.resize(nspec)
                    for ii in range(nspec):
                        glob.ion.ion[ii].element.resize(1)
                        glob.ion.ion[ii].distribution_assumption = 0
                        #glob.ion.ion[ii].power_thermal =

            if False: #data is for al antennas
                wav.profiles_1d.resize(ntimep)
                #fill for all time slices
                self.im.fill_profiles1d_grid(wav)
                for i in range(ntimep):
                    iprof = wav.profiles_1d[i]
                    self.im.fill_data(iprof.power_density, pow_dens, ind=i, name='pow_dens')
                    #cur_pll =
                    #self.im.fill_data(iprof.current_parallel_density, cur_pll, ind=i, name='cur_pll')
                    self.im.fill_data(prof.electrons.power_density_thermal,
                        epdth_prof, ind=i, name='epdth_prof')
                    self.im.fill_profiles1d_grid(iprof)


        if len(ants)> 0:
            self.im.fill_vac_tor_field(ids)
            ids = self.im.fill_properties(ids)
            self.im.fill_code_properties(ids)
            ids.ids_properties.homogeneous_time = 1
            ids.time.resize(ntimes)
            ids.time = timesc
            self.im.putdb_imas(ids)
            self.logger.info('Done')

if __name__ == '__main__':
    print("Please use: python convert.py ")
