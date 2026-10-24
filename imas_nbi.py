#!/usr/bin/python
"""Populate the IMAS nbi IDS from TRANSP/NUBEAM injector configuration and time traces."""

from nc2imas import *
from writer_base import IdsWriter


class imas_nbi(IdsWriter):
    def find_beam_times(self, stn, times, nlines):
        val = 0
        limit = 0
        if nlines>0:
           einje1 = self.cdf.get_ncvar(f'EINJ{stn}_E1')
           einje2 = self.cdf.get_ncvar(f'EINJ{stn}_E2')
           einje3 = self.cdf.get_ncvar(f'EINJ{stn}_E3')
           if self.fill.is_valid(einje1, einje2, einje3):
              check = np.logical_and(einje1>limit, einje2>limit, einje3>limit)
              times = np.where(check, times, val)
        else:
            einj = self.cdf.get_ncvar(f'EINJ')
            if self.fill.is_valid(einj):
               times = np.where(einj>0, times, val)
        #use original times if no data found?
        return times


    def build(self):
        self.logger.info('\nnbi')
        nbi = self.dd.new_ids('nbi')
        #scalar time
        timesc = self.cdf.get_ncvar('TIME')
        ntimes = len(timesc)
        self.fill.checkErrExit(timesc)

        slist = self.cdf.get_var_names_re('PINJ?', nd1=1, nd2=2)
        nunits = len(slist)
        nlines = nunits

        #find species info
        spec_list = self.cdf.get_var_names('PINJ_?', cut= 'PINJ_')
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
            self.fill.set_nbi_unit_identity(
                nbu, f'unit_{ib}', f'TRANSP NBI unit {ib}'
            )

            if self.fill.is_valid(sp_a, sp_zn):
                nbu.species.a = float(sp_a)
                # z_n is INT_0D in DD4 (it was FLT_0D in older DDs).
                self.fill.set_nuclear_charge(nbu.species, sp_zn)
                # 'label' was renamed to 'name' in DD 3.42.
                # Strip the TRANSP bookkeeping suffix from names such as D_BEAM.
                sp_name = str(sp_label).replace('_BEAM', '')
                self.fill.set_species_name(nbu.species, sp_name)

            beam_times = self.find_beam_times(stn, timesc, nlines)

            einj = []
            if nlines>0:
               #fill even if zeros
               #energy EINJ02_E1
               einj = self.cdf.get_ncvar(f'EINJ{stn}_E1')
            else:
               #EINJ(TIME)EV; MAX INITIAL BEAM ENERGY
               einj = self.cdf.get_ncvar('EINJ')

            nbu.energy.data.resize(len(einj))
            nbu.energy.data = einj

            #times
            #self.fill.fill_data(nbu.energy.time, beam_times, name='energy_times')
            nbu.energy.time.resize(ntimes)
            nbu.energy.time = beam_times
            if i==0:
               self.logger.info('NOTE: Scalar time filled')

            if nlines>0:
               #injected power PINJ02
               pinj = self.cdf.get_ncvar(f'PINJ{stn}')
            else:
               pinj = self.cdf.get_ncvar(f'PINJ_{stn}')

            nbu.power_launched.data.resize(len(pinj))
            nbu.power_launched.data = pinj

            nbu.power_launched.time.resize(ntimes)
            nbu.power_launched.time = beam_times

        if nunits > 0:
            nbi = self.fill.fill_properties(nbi)
            self.fill.fill_code_properties(nbi)
            self.fill.fill_external_model_inputs(nbi, 'nbi')
            nbi.ids_properties.homogeneous_time = 1
            nbi.time.resize(ntimes)
            nbi.time = timesc
            self.logger.info('Done NBI')
            return [nbi]
        return []

if __name__ == '__main__':
    print("Please use: python convert.py ")
