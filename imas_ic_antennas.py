#!/usr/bin/python
"""Populate the IMAS ic_antennas IDS from TRANSP ion-cyclotron antenna data."""


from nc2imas import *

class imas_ic_antennas:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\nic_antennas')
        ica = self.im.new_ids('ic_antennas')

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
            # waves/coherent_wave/identifier/antenna_name must reference the
            # corresponding IC antenna name.
            if hasattr(ant, 'name'):
                ant.name = f'antenna{ia}'
            if hasattr(ant, 'description'):
                ant.description = f'TRANSP ion-cyclotron antenna {ia}'

            # frequency
            self.im.fill_data_nc(ant.frequency.data, f'FREQA{sta}')
            ant.frequency.time = np.asarray(timesc, dtype=np.float64)

            # power launched
            self.im.fill_data_nc(ant.power_launched.data, f'PICHA{sta}')
            ant.power_launched.time = np.asarray(timesc, dtype=np.float64)

            # Forward power arriving at the antenna.  TRANSP uses PICHA as the
            # injected/forward power supplied to the IC calculation.
            self.im.fill_data_nc(ant.power_forward.data, f'PICHA{sta}')
            ant.power_forward.time = np.asarray(timesc, dtype=np.float64)

        if nant > 0:
            self.im.fill_data_nc(ica.power_launched.data, 'PICHTOT')
            ica.power_launched.time = np.asarray(timesc, dtype=np.float64)
            self.im.putdb_imas(ica)
            self.logger.info('Done')
        else:
            self.logger.info('ic_antennas: no PICHAj output found; IDS not written')


if __name__ == '__main__':
    print("Please use: python convert.py ")
