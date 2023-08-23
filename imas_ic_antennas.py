
from nc2imas import *

class imas_ic_antennas:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\nic_antennas')
        ica = self.im.ids.ic_antennas
        sgn = -1 if self.im.iter_ else 1

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        self.im.checkErrExit(timesc)

        ica = self.im.fill_properties(ica)
        self.im.fill_code_properties(ica)
        ica.ids_properties.homogeneous_time = 1
        ica.time.resize(ntimes)
        ica.time = timesc

        nant = self.im.get_var_names_re('PICHA?', nd1=1, nd2=2, length=True)
        self.logger.info(f'ic_antenna {nant}')

        power = []

        ica.antenna.resize(nant)
        for i in range(nant):
            ant = ica.antenna[i]
            ia = i+1
            sta = f'{ia}'

            ant.module.resize(1)
            mod = ant.module[0]
            si = 0 # only 1 strap

            #frequency
            self.im.fill_data(mod.frequency.time, timesc)
            self.im.fill_data_nc(mod.frequency.data, f'FREQA{sta}')

            #power launched
            self.im.fill_data(mod.power_launched.time, timesc)
            self.im.fill_data_nc(mod.power_launched.data, f'PICHA{sta}')

            #power_forward. Forward power arriving to the back of the antenna
            # in TRANSP there is only launched power, but users fill the UFILES with the forward power.
            # when TRANSP runs the IC code, it assumes that injected power is the forward power. 
            self.im.fill_data(mod.power_forward.time, timesc)
            self.im.fill_data_nc(mod.power_forward.data, f'PICHA{sta}')

            #power_reflected. Don't use launched - absorbed

        if nant > 0:
            #power launched
            self.im.fill_data(ica.power_launched.time, timesc)
            self.im.fill_data_nc(ica.power_launched.data, 'PICHTOT')

            self.im.putdb_imas(ica)
            self.logger.info('Done')


if __name__ == '__main__':
    print("Please use: python convert.py ")
