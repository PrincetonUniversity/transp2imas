from nc2imas import *

class imas_nbi:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\nnbi')
        nbi = self.im.ids.nbi
        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        self.im.checkErrExit(timesc)

        nunits = self.im.get_var_names_re('PINJ?', nd1=1, nd2=2, length=True)
        self.logger.info(f'NBeams nunits {nunits}')

        #find species info
        spec_list = self.im.get_var_names('PINJ_?', cut= 'PINJ_')
        symb = '' if len(spec_list)==0 else spec_list[0]
        atomsDef = atom_base()
        sp_zn, sp_a, sp_label = atomsDef.find_beam_species_data(symb)

        if len(spec_list)>1:
            #beams numbers PINJ<n> are not matched for different species
            raise Exception(f'Multiple species of NBI not handled: {spec_list}')
        if (symb and nunits<1) or (nunits>0 and not symb):
            raise Exception(f'NBI species/beams Error:: units:{nunits} species:{symb}')


        nbi.unit.resize(nunits)
        for i in range(nunits):
            nbu = nbi.unit[i]
            ib = i+1
            stn = f'0{ib}' if ib<=9 else f'{ib}'
            nbu.name = f"unit {ib}"

            if self.im.is_valid(sp_a, sp_zn):
                nbu.species.a = sp_a
                nbu.species.z_n = sp_zn
                nbu.species.label = sp_label

            #energy EINJ02_E1
            self.im.fill_data_nc(nbu.energy.data, f'EINJ{stn}_E1')
            self.im.fill_data(nbu.energy.time, timesc)

            #injected power PINJ02
            self.im.fill_data_nc(nbu.power_launched.data, f'PINJ{stn}')
            self.im.fill_data(nbu.power_launched.time, timesc)

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
