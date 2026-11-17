#!/usr/bin/python
"""Populate the IMAS lh_antennas IDS from TRANSP lower-hybrid antenna data."""


from nc2imas import *
from writer_base import IdsWriter


class imas_lh_antennas(IdsWriter):
    def build(self):
        self.logger.info('\nlh_antennas')
        lha = self.dd.new_ids('lh_antennas')
        sgn = -1 if self.cdf.iter_ else 1

        #scalar time
        timesc = self.cdf.get_ncvar('TIME')
        self.fill.checkErrExit(timesc)

        lha = self.fill.fill_properties(lha)
        self.fill.fill_code_properties(lha)
        lha.ids_properties.homogeneous_time = 1
        lha.time = timesc

        nant = self.cdf.get_var_names_re('PLHANT?', nd1=1, nd2=2, length=True)
        self.logger.info(f'lh antenna {nant}')

        lha.antenna.resize(nant)
        for i in range(nant):
            ia = i+1
            sta = f'{ia}'

            #frequency
            freq = self.namelist.get_namelist_var('FREQLH', float, ind=ia-1)
            if self.fill.is_valid(freq):
                lha.antenna[i].frequency = freq

            # TRANSP does not know the difference between launched and forward power.
            # for consistency with what written in the CDF file documentation,
            # What users give to TRANSP is usually the forward power,
            # which TRANSP will interpret as generic input power
            # we need to write the 'input' power that is in the CDF file to both the launched and forward power.
            #power launched
            self.logger.info(f'Ant: {ia} PLHANT{ia}')
            self.fill.fill_data_nc(lha.antenna[i].power_launched.data, f'PLHANT{sta}')

            #power_forward. Forward power arriving to the back of the antenna
            self.logger.info(f'Ant: {ia} PLHANT{ia}')
            self.fill.fill_data_nc(lha.antenna[i].power_forward.data, f'PLHANT{sta}')

        if nant > 0:
            #power launched
            self.fill.fill_data_nc(lha.power_launched.data, 'PLH')

            self.logger.info('Done')
            return [lha]
        return []


if __name__ == '__main__':
    print("Please use: python convert.py ")
