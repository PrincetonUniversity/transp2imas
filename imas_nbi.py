from nc2imas import *

class imas_nbi:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()

    def find_beam_times(self, stn, times, nlines):
        val = 0
        limit = 0
        if nlines>0:
           einje1 = self.im.get_ncvar(f'EINJ{stn}_E1')
           einje2 = self.im.get_ncvar(f'EINJ{stn}_E2')
           einje3 = self.im.get_ncvar(f'EINJ{stn}_E3')
           if self.im.is_valid(einje1, einje2, einje3):
              check = np.logical_and(einje1>limit, einje2>limit, einje3>limit)
              times = np.where(check, times, val)
        else:
            einjs = self.im.get_ncvar(f'EINJ')
            if self.im.is_valid(einjs):
               times = np.where(einj>0, times, val)
        #use original times if no data found?
        return times


    def create_ids(self):
        self.im.logger.info('\nnbi')
        nbi = self.im.ids.nbi
        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        self.im.checkErrExit(timesc)

        slist = self.im.get_var_names_re('PINJ?', nd1=1, nd2=2)
        nunits = len(slist)
        nlines = nunits

        #find species info
        spec_list = self.im.get_var_names('PINJ_?', cut= 'PINJ_')
        nspec_list = len(spec_list)

        if nunits==0 and nspec_list>0:
            nunits = nspec_list
        self.logger.info(f'NBeams nunits {nunits}; PINJ {slist}; species: {spec_list}')
        if nspec_list>1 and nlines>0:
            #beams numbers PINJ<n> are not matched for different species
            self.logger.info(f'Warning: {spec_list} not matched to {nlines} lines')

        nbi.unit.resize(nunits)
        for i in range(nunits):
            symb = spec_list[0]
            sname = symb
            if nspec_list>i:
                symb = spec_list[i]
                sname = spec_list[i]

            atomsDef = atom_base()
            sp_zn, sp_a, sp_label = atomsDef.find_beam_species_data(symb)

            nbu = nbi.unit[i]
            ib = i+1
            stn = f'0{ib}' if ib<=9 else f'{ib}'
            nbu.name = f"unit {ib}"

            if self.im.is_valid(sp_a, sp_zn):
                nbu.species.a = sp_a
                nbu.species.z_n = sp_zn
                nbu.species.label = sp_label

            beam_times = self.find_beam_times(stn, timesc, nlines)

            einj = []
            if nlines>0:
               #fill even if zeros
               #energy EINJ02_E1
               einj = self.im.get_ncvar(f'EINJ{stn}_E1')
            else:
               #EINJ(TIME)EV; MAX INITIAL BEAM ENERGY
               einj = self.im.get_ncvar('EINJ')

            nbu.energy.data.resize(len(einj))
            nbu.energy.data = einj

            #times
            #self.im.fill_data(nbu.energy.time, beam_times, name='energy_times')
            nbu.energy.time.resize(ntimes)
            nbu.energy.time = beam_times
            if i==0:
               print('NOTE: Scalar time filled')

            if nlines>0:
               #injected power PINJ02
               pinj = self.im.get_ncvar(f'PINJ{stn}')
            else:
               pinj = self.im.get_ncvar(f'PINJ_{stn}')

            nbu.power_launched.data.resize(len(pinj))
            nbu.power_launched.data = pinj

            nbu.power_launched.time.resize(ntimes)
            nbu.power_launched.time = beam_times


            if i==4:
               print(pinj)


        if nunits > 0:
            nbi = self.im.fill_properties(nbi)
            self.im.fill_code_properties(nbi)
            nbi.ids_properties.homogeneous_time = 1
            nbi.time.resize(ntimes)
            nbi.time = timesc
            self.im.putdb_imas(nbi)
            self.logger.info('Done NBI')

if __name__ == '__main__':
    print("Please use: python convert.py ")
