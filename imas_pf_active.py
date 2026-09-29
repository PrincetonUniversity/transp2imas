#!/usr/bin/python
"""Populate the IMAS pf_active IDS from TRANSP active poloidal-field circuit data."""


from nc2imas import *

class imas_pf_active:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\npf_active')
        pfc = self.im.new_ids('pf_active')
        cm2m = 0.01
        vol_inv = 1.e6

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        timep = self.im.get_ncvar('TIME3')
        ntimep = len(timep)
        self.im.checkErrExit(timesc, timep)

        pfc = self.im.fill_properties(pfc)
        self.im.fill_code_properties(pfc)
        pfc.ids_properties.homogeneous_time = 1
        # Assign through the numeric-array interface with an explicit float64
        # dtype.  DD4 FLT_1D nodes are double precision and direct structure
        # assignment can trigger noisy automatic-conversion diagnostics.
        # Direct assignment of a float64 ndarray is the cleanest DD4 path for
        # this top-level time coordinate.
        pfc.time = np.asarray(timesc, dtype=np.float64)

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
            # PF coil-current sign is an independent circuit/winding convention;
            # do not reuse the plasma-current orientation sign.
            self.im.fill_data_nc(pfc.coil[i].current.data, 'KK_'+coil)

            #Voltage on the coil terminals
            #VV_PF1AU(TIME) Coil terminal voltage for pf1au
            self.im.fill_data_nc(pfc.coil[i].voltage.data, 'VV_'+coil)


        #supply
        pfc.supply.resize(ncoils)
        for i in range(ncoils):
            coil = coils[i]
            pfc.supply[i].name = coil
            self.logger.info(f'PF/OH supply: {coil}')

            #Current at the supply output
            #KT_PF2U(TIME) Target current in pf2u circuit. VS_PF2U is desired current
            self.im.fill_data_nc(pfc.supply[i].voltage.data, 'KT_'+coil)

            #Voltage at the supply output
            #VS_PF1AU(TIME) Coil source voltage for pf1au
            self.im.fill_data_nc(pfc.supply[i].voltage.data, 'VS_'+coil)


        if ncoils==0:
            # use this to exclude pf_coils, since size of coils=1 if empty in trdat_ids !
            pfc.time = np.array([-999.0], dtype=np.float64)
        #if ncoils > 0: # error if empty; so filling with invalid times, to use it to skip this
        self.im.putdb_imas(pfc)
        self.logger.info('Done')


if __name__ == '__main__':
    print("Please use: python convert.py ")
