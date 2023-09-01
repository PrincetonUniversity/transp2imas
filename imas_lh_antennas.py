
from nc2imas import *

class imas_lh_antennas:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\nlh_antennas')
        lha = self.im.ids.lh_antennas
        sgn = -1 if self.im.iter_ else 1

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        self.im.checkErrExit(timesc)

        lha = self.im.fill_properties(lha)
        self.im.fill_code_properties(lha)
        lha.ids_properties.homogeneous_time = 1
        lha.time.resize(ntimes)
        lha.time = timesc

        nant = self.im.get_var_names_re('PLHANT?', nd1=1, nd2=2, length=True)
        self.logger.info(f'lh antenna {nant}')

        lha.antenna.resize(nant)
        for i in range(nant):
            ia = i+1
            sta = f'{ia}'

            #frequency
            freq = self.im.get_namelist_var('FREQLH', float, ind=ia-1)
            if self.im.is_valid(freq):
                lha.antenna[i].frequency = freq

            # TRANSP does not know the difference between launched and forward power.
            # for consistency with what written in the CDF file documentation,
            # What users give to TRANSP is usually the forward power,
            # which TRANSP will interpret as generic input power
            # we need to write the 'input' power that is in the CDF file to both the launched and forward power.
            #power launched
            self.logger.info(f'Ant: {ia} PLHANT{ia}')
            self.im.fill_data(lha.antenna[i].power_launched.time, timesc)
            self.im.fill_data_nc(lha.antenna[i].power_launched.data, f'PLHANT{sta}')

            #power_forward. Forward power arriving to the back of the antenna
            self.logger.info(f'Ant: {ia} PLHANT{ia}')
            self.im.fill_data(lha.antenna[i].power_forward.time, timesc)
            self.im.fill_data_nc(lha.antenna[i].power_forward.data, f'PLHANT{sta}')

        if nant > 0:
            #power launched
            self.im.fill_data(lha.power_launched.time, timesc)
            self.im.fill_data_nc(lha.power_launched.data, 'PLH')

            self.im.putdb_imas(lha)
            self.logger.info('Done')


if __name__ == '__main__':
    print("Please use: python convert.py ")
