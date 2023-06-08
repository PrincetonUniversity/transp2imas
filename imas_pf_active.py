
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
        coils_pf = self.im.get_var_names_re('VV_PF?', cut='VV_', nd=1, nch1=1, nch2=2)

        #coil names from VV_OH
        coils_oh = self.im.get_var_names_re('VV_OH?', cut='VV_', nd=0, nch=0)

        coils = coils_pf + coils_oh
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
    a = get_args(sys.argv)

    sh, rn, mach = 142202, 'U01', 'NSTX'

    ncf = str(sh)+rn+'.CDF'
    nml = str(sh)+rn+'TR.DAT'
    nci = nc2imas(db=mach, shot=sh, backend='MDS', run=6, occur=0, user=os.getenv('USER'),
          new_db=True, ncfile=ncf, nml_file=nml, log='debug')
    imas_pf_active(nci)