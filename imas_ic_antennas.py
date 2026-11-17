#!/usr/bin/python
"""Populate the IMAS ic_antennas IDS from TRANSP ion-cyclotron antenna data."""


from nc2imas import *
from writer_base import IdsWriter


class imas_ic_antennas(IdsWriter):
    def build(self):
        self.logger.info('\nic_antennas')
        ica = self.dd.new_ids('ic_antennas')

        #scalar time
        timesc = self.cdf.get_ncvar('TIME')
        ntimes = len(timesc)
        self.fill.checkErrExit(timesc)

        ica = self.fill.fill_properties(ica)
        self.fill.fill_code_properties(ica)
        ica.ids_properties.homogeneous_time = 1
        ica.time.resize(ntimes, refcheck=False)
        ica.time = timesc

        nant = self.cdf.get_var_names_re('PICHA?', nd1=1, nd2=2, length=True)
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
            self.fill.fill_data_nc(ant.frequency.data, f'FREQA{sta}')
            ant.frequency.time = np.asarray(timesc, dtype=np.float64)

            # power launched
            self.fill.fill_data_nc(ant.power_launched.data, f'PICHA{sta}')
            ant.power_launched.time = np.asarray(timesc, dtype=np.float64)

            # Forward power arriving at the antenna.  TRANSP uses PICHA as the
            # injected/forward power supplied to the IC calculation.
            self.fill.fill_data_nc(ant.power_forward.data, f'PICHA{sta}')
            ant.power_forward.time = np.asarray(timesc, dtype=np.float64)

        if nant > 0:
            self.fill.fill_data_nc(ica.power_launched.data, 'PICHTOT')
            ica.power_launched.time = np.asarray(timesc, dtype=np.float64)
            self.logger.info('Done')
            return [ica]
        self.logger.info('ic_antennas: no PICHAj output found; IDS not written')
        return []


if __name__ == '__main__':
    print("Please use: python convert.py ")
