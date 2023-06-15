
from nc2imas import *

class imas_pf_passive:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\npf_passive')
        pfc = self.im.ids.pf_passive
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

        #coil names from PS_1 .. PS_999
        coils = self.im.get_var_names_re('PS_?', cut='PS_', nd1=1, nd2=3)

        ncoils = len(coils)
        self.logger.info(f'ncoils {ncoils}')

        pfc.loop.resize(ncoils)

        #coils
        for i in range(ncoils):
            nst = coils[i]
            name = nst  #TODO get from description "..passive structure 1 <bvfb>"

            self.logger.info(f'passive coil loop: {nst}')

            pfc.loop[i].name = name

            #PS_1(TIME); A; Current in passive structure 1 <bvfb>
            self.im.fill_data_nc(pfc.loop[i].current, 'PS_'+nst, mult=sgn)

            #self.im.fill_data_nc(pfc.loop[i].resistance, '??_'+nst)

            self.im.fill_data(pfc.loop[i].time, timesc)

            nel = 0
            pfc.loop[i].element.resize(nel)





        if ncoils > 0:
            self.im.putdb_imas(pfc)
            self.logger.info('Done')


if __name__ == '__main__':
    #a = get_args(sys.argv)

    sh, rn, mach = 142202, 'U01', 'NSTX'

    ncf = str(sh)+rn+'.CDF'
    nml = str(sh)+rn+'TR.DAT'
    nci = nc2imas(db=mach, shot=sh, backend='MDS', run=6, occur=0, user=os.getenv('USER'),
          new_db=True, ncfile=ncf, nml_file=nml, log='debug')
    imas_pf_passive(nci)