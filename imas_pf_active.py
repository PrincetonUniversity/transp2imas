
from nc2imas import *

class imas_pf_active:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\npf_active')
        pfc = self.im.ids.pf_active
        cm2m = 0.01
        vol_inv = 1.e6
        sgn = -1 if self.im.iter_ else 1

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        timep = self.im.get_ncvar('TIME3')
        ntimep = len(timep)
        self.im.checkErrExit(timesc, timep)

        pfc = self.im.fill_properties(pfc)
        self.im.fill_code_properties(pfc)
        pfc.ids_properties.homogeneous_time = 1
        pfc.time.resize(ntimes)
        pfc.time = timesc

        #coil names from VV_PF1AU(PF_ + PF1AU)
        coils = self.im.get_var_names('VV_?', cut='VV_')

        # not needed - we can parse the name from the names in the CDF file.
        # please, use this case for testing: /p/transpusers/fpoli/current/transp2imas_tests/ITER/38265R95/

        ncoils = len(coils)
        self.logger.info(f'ncoils {ncoils}')

        pfc.coil.resize(ncoils)

        #coils
        for i in range(ncoils):
            coil = coils[i]
            pfc.coil[i].name = coil
            self.logger.info(f'PF/OH coil: {coil}')

            #Current fed in the coil (for 1 turn)
            #KK_PF1AU(TIME) Current in PF1AU circuit; CC_PF.. from equil ; PF_PF1AU zeros
            self.im.fill_data_nc(pfc.coil[i].current.data, 'KK_'+coil, mult=sgn)
            self.im.fill_data(pfc.coil[i].current.time, timesc)

            #Voltage on the coil terminals
            #VV_PF1AU(TIME) Coil terminal voltage for pf1au
            self.im.fill_data_nc(pfc.coil[i].voltage.data, 'VV_'+coil)
            self.im.fill_data(pfc.coil[i].voltage.time, timesc)


        #supply
        pfc.supply.resize(ncoils)
        for i in range(ncoils):
            coil = coils[i]
            pfc.supply[i].name = coil
            self.logger.info(f'PF/OH supply: {coil}')

            #Current at the supply output
            #KT_PF2U(TIME) Target current in pf2u circuit. VS_PF2U is desired current
            self.im.fill_data_nc(pfc.supply[i].voltage.data, 'KT_'+coil)
            self.im.fill_data(pfc.supply[i].voltage.time, timesc)

            #Voltage at the supply output
            #VS_PF1AU(TIME) Coil source voltage for pf1au
            self.im.fill_data_nc(pfc.supply[i].voltage.data, 'VS_'+coil)
            self.im.fill_data(pfc.supply[i].voltage.time, timesc)


        if ncoils > 0:
            self.im.putdb_imas(pfc)
            self.logger.info('Done')


if __name__ == '__main__':
    print("Please use: python convert.py ")
