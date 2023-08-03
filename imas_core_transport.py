
from nc2imas import *

class imas_core_transport:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\ncore transport')
        ct = self.im.ids.core_transport
        cm2m = 0.01
        sgn = -1 if self.im.iter_ else 1

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        timep = self.im.get_ncvar('TIME3')
        ntimep = len(timep)
        self.im.checkErrExit(timesc, timep)

        ct = self.im.fill_properties(ct)
        self.im.fill_code_properties(ct)
        ct.ids_properties.homogeneous_time = 1
        ct.time.resize(ntimes)
        ct.time = timesc

        self.im.fill_vac_tor_field(ct)

        coord_x = self.im.get_ncvar('X')
        nprof = coord_x.shape[1]


        # check if this is a predictive run
        try:
            _ = self.im.get_ncvar('COOTRABSP')
            lpredict = 1
            ct.model.resize(3)
            ctp = ct.model[2]
            ctp.identifier.name = "combined"
            ctp.identifier.index = 1
            ctp.identifier.description = "PT_SOLVER for transport equations in" \
                  " tokamaks with anomalous and neoclassical transport: Combined output"
        except:
            lpredict = 0
            ct.model.resize(2)
            
        cti = ct.model[0]
        cti.identifier.name = "combined"
        cti.identifier.index = 1
        cti.identifier.description = "Interpetive results for diffusivity and convective " \
                  "velocity profiles"

        ctn = ct.model[1]
        ctn.identifier.name = "combined"
        ctn.identifier.index = 1
        ctn.identifier.description = "Neoclassical computations of diffusivity and convective " \
                  "velocity profiles"

        #electrons energy.flux
        #eetran_model(j) = ss%pe_trans(izon)/dvol(j,2)
        #C%EETRAN_MODEL(#1)%'Div(elec energy flux) (model)'%'WATTS/CM3'%EETR_MOD#+EETR
        # EETR_MOD(TIME3, X); WATTS/CM3; Div(elec energy flux) (model)
        # EETR_MOD * DVOL /  PAREA
        eetr_mod = self.im.get_ncvar('EETR_MOD', self.im.vol_inv)
        dvol = self.im.get_ncvar('DVOL', self.im.vconv)
        parea = self.im.get_ncvar('PAREA', self.im.aconv)
        eeflux = []
        if self.im.is_valid(eetr_mod, dvol, parea):
            eeflux = eetr_mod*dvol/parea

        #IETR_MOD(TIME3, X); WATTS/CM3; Div(ion energy flux) (model)
        ietr_mod = self.im.get_ncvar('IETR_MOD', self.im.vol_inv)
        ieflux = []
        if self.im.is_valid(ietr_mod, dvol, parea):
            ieflux = ietr_mod*dvol/parea

        #EPTR_MOD(TIME3, X); N/CM3/SEC; Div(electron flux) (model)
        eptr_mod = self.im.get_ncvar('EPTR_MOD', self.im.vol_inv)
        epflux = []
        if self.im.is_valid(eptr_mod, dvol, parea):
            epflux = eptr_mod*dvol/parea

        #IPTR_MOD(TIME3, X) ; N/CM3/SEC; Div(total ion flux) (model)
        iptr_mod = self.im.get_ncvar('IPTR_MOD', self.im.vol_inv)
        ipflux = []
        if self.im.is_valid(iptr_mod, dvol, parea):
            ipflux = iptr_mod*dvol/parea

        #ions = self.im.get_var_names_re('DFI_?', ....  )

        ct.profiles_1d.resize(ntimep)

        for i in range(ntimep):
            ctp = ct.profiles_1d[i]
            # --electron energy--
            #CONDE(TIME3, XB) ; CM**2/SEC; ELECTRON HEAT DIFFUSIVITY
            self.im.fill_data_nc(cti.electrons.energy.d, var='CONDE', mult=self.im.aconv, ind=i)
            self.im.fill_data_nc(ctp.electrons.energy.d, var='CONDEPR', mult=self.im.aconv, ind=i)
            #VELE(TIME3, XB); cm/sec; ELECTRON RADIAL VELOCITY
            self.im.fill_data_nc(cti.electrons.energy.v, var='VELE', mult=self.im.cm2m, ind=i)
            # --electron particles--
            #DIFFE(TIME3, XB) ; CM**2/SEC; ELEC PTCL DIFFUSIVITY 
            self.im.fill_data_nc(cti.electrons.particles.d, var='DIFFE', mult=self.im.aconv, ind=i)
            self.im.fill_data_nc(ctp.electrons.particles.d, var='DIFFNE', mult=self.im.aconv, ind=i)

            self.im.fill_data(cti.electrons.energy.flux, eeflux, ind=i, var='eeflux')


            #
            #self.im.fill_data_nc(ctp.electrons.particles.v,

            #
            self.im.fill_data(ctp.electrons.particles.epflux, epflux, ind=i, var='epflux')

            # --ions energy--
            #CONDI(TIME3, XB) ; CM**2/SEC ; ION HEAT DIFFUSIVITY
            self.im.fill_data_nc(ctp.total_ion_energy.d, 'CONDI', mult=1.0e4, ind=i, var='ie_d')

            #self.im.fill_data(ctp.total_ion_energy.v,

            #ions flux
            self.im.fill_data(ctp.total_ion_energy.flux, ipflux, ind=i, var='ipflux')


            # --toroidal momentum--
            #CHPHM(TIME3, XB) ; CM**2/SEC; MOMENTUM CHI(PHI) MODEL
            self.im.fill_data_nc(ctp.momentum_tor.d, 'CHPHM', mult=self.im.aconv, ind=i, var='chphm')

            #VMO_THMOD(TIME3, XB) ; CM/SEC ; Momentum v_rad, transport model
            self.im.fill_data_nc(ctp.momentum_tor.v, 'VMO_THMOD', mult=self.im.cm2m, ind=i, var='mtor_v')

            #MCONV(TIME3, X) ; Nt-M/CM3 ; CONVECTIVE TRANSPORT
            #self.im.fill_data_nc(ctp.momentum_tor.flux,

            ctp.ion.resize(len(ions))
            #DFI_D(TIME3, XB) ;CM**2/SEC; D+ ION DIFFUSIVITY (NMODEL=4)
            for ii, iname in enumerate(ions):
                self.im.fill_data_nc(ctp.ion[ii].energy.d, iname, mult=self.im.area_inv, ind=i, var=iname)

                #CONDWNCLI  CONDWNC4


if __name__ == '__main__':
    print("Please use: python convert.py")
