
from nc2imas import *

class imas_lh_antennas:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\nlh_antennas')
        lha = self.im.ids.lh_antennas
        cm2m = 0.01
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

            #power launced
            self.logger.info(f'Ant: {ia} PLHANT{ia}')
            self.im.fill_data(lha.antenna[i].power_launched.time, timesc)
            self.im.fill_data_nc(lha.antenna[i].power_launched.data, f'PLHANT{sta}')

            #power_forward. Forward power arriving to the back of the antenna


            #power_reflected = launched - absorbed
            plaunch = self.im.get_ncvar(f'PLHANT{sta}')
            self.logger.info(f'Ant: {ia} PLHABS{sta}')
            pabs = self.im.get_ncvar(f'PLHABS{sta}')
            if self.im.is_valid(plaunch, pabs):
                prefl = plaunch - pabs
                self.im.fill_data(lha.antenna[i].power_reflected.data, prefl)
                self.im.fill_data(lha.antenna[i].power_reflected.time, timesc)
            else:
                self.logger.info(f'Ant: {ia} Error: power_reflected. invalid: plaunch or pabs')

            #surface_current on the antenna straps and on passive components, for every time slice

            #module.  Set of antenna modules (each module is fed by a single transmission line)


        if nant > 0:
            #power launched
            self.im.fill_data(lha.power_launched.time, timesc)
            self.im.fill_data_nc(lha.power_launched.data, 'PLH')

            self.im.putdb_imas(lha)
            self.logger.info('Done')


if __name__ == '__main__':
    print("Please use: python convert.py ")
