#!/usr/bin/python
"""Populate the IMAS pf_passive IDS when passive-structure current and geometry data are available."""


from nc2imas import *

class imas_pf_passive:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()


    def create_ids(self):
        self.im.logger.info('\npf_passive')
        pfc = self.im.new_ids('pf_passive')

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

        # Passive-current channels saved by some TRANSP runs.  TRANSP2IMAS
        # deliberately does not reconstruct hardware geometry from run namelist
        # settings: authoritative geometry/resistance/connectivity belongs to the
        # machine-description/IDS database.
        coils = self.im.get_var_names_re('PS_?', cut='PS_', nd1=1, nd2=3)

        ncoils = len(coils)
        self.logger.info(f'ncoils {ncoils}')
        if ncoils == 0:
            self.logger.info('pf_passive: no PS_n passive-current channels; IDS not written')
            return

        self.logger.warning(
            'pf_passive: writing TRANSP passive-loop names/currents only; '
            'geometry, resistance, inductance, and circuit connectivity must come '
            'from the authoritative machine-description/IDS database')
        pfc.loop.resize(ncoils)

        #coils
        for i in range(ncoils):
            nst = coils[i]

            # PS_1:long_name  'Current in passive structure n <bvfb>'
            lname = self.im.get_ncattr(f'PS_{nst}:long_name')
            lname = str(lname).strip() if lname is not None else ''
            name = self.im.extract_str_between(lname, '<', '>', str(nst)) if lname else str(nst)
            if i==0 or i==ncoils-1:
                self.logger.info(f'passive coil loop: {nst}')

            pfc.loop[i].name = name

            #PS_1(TIME); A; Current in passive structure 1 <bvfb>
            # Passive-loop sign is an independent circuit orientation; do not
            # force it to follow the plasma-current direction.
            self.im.fill_data_nc(pfc.loop[i].current, f'PS_{nst}')

            #self.im.fill_data_nc(pfc.loop[i].resistance, '??_'+nst)

            nel = 0
            pfc.loop[i].element.resize(nel)





        if ncoils > 0:
            self.im.putdb_imas(pfc)
            self.logger.info('Done')


if __name__ == '__main__':
    print("Please use: python convert.py")
