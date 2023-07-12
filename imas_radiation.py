
from nc2imas import *

class imas_radiation:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()

    def fill_rad_properties(self, rid, iname:str, iind:int, idesc:str):
        rid.identifier.name = iname
        rid.identifier.index = iind
        rid.identifier.description = idesc


    def fill_rad_globals(self, rid, timesc, pwr, pwi=None, pwn=None,
                         pwe=None, vessel:bool=False):
        ntimes = len(timesc)
        if pwr is None or (ntimes != len(pwr)):
            self.im.logger.info(f'line_radiation : global ntime {ntimes} != size of pwr')
            return

        rid.global_quantities.resize(ntimes)
        for i in range(ntimes):
            glob = rid.global_quantities[i]
            glob.time = timesc[i].astype(float)
            if vessel: #not limited to within LCFS
                glob.inside_vessel.power = pwr[i].astype(float)
                if pwi is not None:
                    glob.inside_vessel.power_ion_total = pwi[i].astype(float)
                if pwn is not None:
                    glob.inside_vessel.power_neutral_total = pwn[i].astype(float)
                if pwe is not None:
                    glob.inside_vessel.power_electrons = pwe[i].astype(float)
            else:
                #Total power emitted by all species [W]
                glob.inside_lcfs.power = pwr[i].astype(float)
                if pwi is not None:
                    glob.inside_lcfs.power_ion_total = pwi[i].astype(float)
                if pwn is not None:
                    glob.inside_lcfs.power_neutral_total = pwn[i].astype(float)
                if pwe is not None:
                    glob.inside_lcfs.power_electrons = pwe[i].astype(float)


    #rid:radiation process
    def fill_bremsstrahlung(self, rid, timep, timesc):
        self.im.logger.info(f'radiation: bremsstrahlung')
        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]

        #PRAD_BR(TIME3, X) ; W/CM3; BREMSSTRAHLUNG RADIATION
        #negative for radiation power loss
        rad_prof = -1.0*self.im.get_ncvar('PRAD_BR', self.im.vol_inv)
        if not self.im.is_valid(rad_prof):
            self.im.logger.info(f'rad_prof empty')
            return

        self.fill_rad_properties(rid, 'bremsstrahlung', 8, 'Emission from bremsstrahlung')

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)
        radxv = self.im.operate_arrays(rad_prof, '*', gvol)
        rad_glob = self.im.np_multiply_sum(radxv, gvol, axis=1)

        self.fill_rad_globals(rid, timesc, pwr=rad_glob, pwi=rad_glob)

        ntime = len(timep)
        rid.profiles_1d.resize(ntime)
        #fill profiles_1d[*].grid
        self.im.fill_profiles1d_grid(rid)
        for i in range(ntime):
            profi = rid.profiles_1d[i]
            profi.time = timep[i].astype(float)

            if i==0:
                self.im.logger.info(f'PRAD_BR time {i}')
            self.im.fill_data(profi.electrons.power_inside, rad_prof, name='pwr_br')


    #ids:radiation process, ntimes:n_time_frames
    def fill_line_radiation(self, rid, timep, timesc):
        self.im.logger.info(f'radiation: line_radiation')

        #DVOL(TIME3, X) CM**3; ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)
        #PRAD_LI(TIME3, X) W/CM3: LINE RADIATION
        prad = self.im.get_ncvar('PRAD_LI', self.im.vol_inv)

        self.fill_rad_properties(rid, 'line_radiation', 10, 'Emission from line radiation')

        #PRLS_BE(TIME3, X) WATTS/CM3 ; Be Impurity Line Radiation
        neut_el = self.im.get_var_names_re('PRLS_?', cut='PRLS_', nch1=1, nch2=3)
        pwr_neut_names = ['PRLS_'+str(s) for s in  neut_el]
        self.im.logger.info(f'pwr_neut_names: {pwr_neut_names}')
        #total; separate is done in profiles
        imp_npwr = self.im.get_ncvars_lst(pwr_neut_names, self.im.vol_inv)

        #PRLX_BE_4(TIME3, X) WATTS/CM3 ; Be+4 Impurity Line Radiation
        pwx_ion_names = self.im.get_var_names_re('PRLX_?', nd1=1, nd2=2,
                                                nch1=1, nch2=3, dig_last=True)
        ions_dict = {s: s.split('_')[1:] for s in pwx_ion_names if len(s.split('_'))>1}
        self.im.logger.info(f'pwx_ion_names: {pwx_ion_names}')
        #total; separate is done in profiles
        imp_ipwr = self.im.get_ncvars_lst(pwx_ion_names, self.im.vol_inv)
        imp_ipwr_prof = self.im.operate_arrays(imp_ipwr, '*', gvol)
        imp_ipwr = self.im.np_multiply_sum(imp_ipwr, gvol, axis=1)

        #Total = PRAD_LI + PRLS_
        #prad = prad + imp_npwr  #Does PRAD_LI include PRLS_ ?
        prad_noname = self.im.operate_arrays(prad, '-', imp_npwr, 'prad', 'imp_npwr')

        pwn = self.im.np_multiply_sum(prad, gvol, axis=1)
        pwr = self.im.operate_arrays(pwn, '+', imp_ipwr)
        #global quantities
        self.fill_rad_globals(rid, timesc, pwr=pwr, pwn=pwn, pwi=imp_ipwr)

        ntimep = len(timep)
        rid.profiles_1d.resize(ntimep)

        for i in range(ntimep):
            #rid.profiles_1d.grid is filled in this
            self.im.fill_profiles1d_grid(rid)
            profi = rid.profiles_1d[i]

            profi.time = timep[i].astype(float)

            #PRAD_LI is total neutral  ;  unnamed neutral = PRAD_LI - PRLS_
            #add neutral impurities separately
            profi.neutral.resize(1+len(neut_el))
            #total neutrals, no name. Not filling species info
            self.im.fill_data(profi.neutral[0].power_inside, prad_noname, ind=i, name='prad_no')

            for ni, (name, sp) in enumerate(zip(pwr_neut_names, neut_el)):
                ni = ni+1 #
                if i==0:
                    self.im.logger.info(f' power for {name}')
                pwri = self.im.get_ncvar(name, self.im.vol_inv, ind=i)
                if self.im.is_valid(pwri, name=name):
                    self.im.fill_element_data(profi.neutral[ni], sp)
                    pwri = self.im.np_multiply_sum(pwri, gvol, axis=1)
                    self.im.fill_data(profi.neutral[ni].power_inside, pwri, name='pwr_n')

            nii = len(ions_dict)
            profi.ion.resize(nii)
            for ii, (name, lst) in enumerate(ions_dict.items()):
                sp, zn = lst
                if i==0:
                    self.im.logger.info(f' power for {name}')
                #sp is TRANSP symbol defined in atoms_class.entries
                self.im.fill_element_data(profi.ion[ii], sp, zn)
                pwri = self.im.get_ncvar(name, self.im.vol_inv, ind=i)
                self.im.fill_data(profi.ion[ii].power_inside, pwri, name='pwr_x')



    #rid:radiation process, ntimes:n_time_frames
    def fill_recombination(self, rid, timep, timesc):
        self.im.logger.info(f'radiation: recombination')

        #P0RECO(TIME3, X) WATTS/CM3 RECOMB NEUTRAL SCE POWER : PFLX0RECO ?
        #H isotopes only ?
        hpwr = self.im.get_ncvar('P0RECO', self.im.vol_inv)
        if not self.im.is_valid(hpwr):
            self.im.logger.info(f'Recombination: hpwr empty')
            return

        self.fill_rad_properties(rid, 'recombination', 11, 'Emission from recombination')

        #DVOL(TIME3, X) CM**3; ZONE VOLUME
        gvol = self.im.get_ncvar('DVOL', self.im.vconv)
        pwr = self.im.np_multiply_sum(hpwr, gvol, axis=1)

        #filling pwn:neutral ?
        self.fill_rad_globals(rid, timesc, pwr=pwr, pwn=pwr)

        ntimep = len(timep)
        rid.profiles_1d.resize(ntimep)

        for i in range(ntimep):
            #rid.profiles_1d.grid filled
            self.im.fill_profiles1d_grid(rid)
            idp = rid.profiles_1d[i]
            idp.time = timep[i].astype(float)

            #emissivity_neutral_total

            #Total power from ion species
            self.im.fill_data(idp.power_inside_neutral_total, hpwr, ind=i, name='hpwr_tot')
            sp = 'H'
            #idp.ion.resize(1) #H- only
            #self.im.fill_element_data(idp.ion[0], sp, '+') #label= H+
            idp.neutral.resize(1) #H- only
            self.im.fill_element_data(idp.neutral[0], sp) #label= H+

            self.im.fill_data(idp.neutral[0].power_inside, hpwr, ind=i, name='hpwr')


    def create_ids(self):
        self.im.logger.info('\n radiation')
        ids = self.im.ids.radiation
        sgn = -1 if self.im.iter_ else 1

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        timep = self.im.get_ncvar('TIME3')
        ntimep = len(timep)
        self.im.checkErrExit(timesc, timep)

        ids = self.im.fill_properties(ids)
        self.im.fill_code_properties(ids)
        ids.ids_properties.homogeneous_time = 1
        ids.time.resize(ntimes)
        ids.time = timesc

        self.im.fill_vac_tor_field(ids)

        #IMPURITY radiation ?????????
        #process
        ids.process.resize(3)  #TODO

        self.fill_line_radiation(ids.process[0], timep, timesc)
        self.fill_recombination(ids.process[1], timep, timesc)
        self.fill_bremsstrahlung(ids.process[2], timep, timesc)

        #if nnn > 0:
        self.im.putdb_imas(ids)
        self.logger.info('Done')


if __name__ == '__main__':
    print("Please use: python convert.py")