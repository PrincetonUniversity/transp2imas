from nc2imas import *

class imas_nbi:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\nnbi')
        nbi = self.im.ids.nbi
        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        self.im.checkErrExit(timesc)

        nunits = self.im.get_var_names_re('PINJ?', nd1=1, nd2=2, length=True)
        nbeams = self.im.get_namelist_var('NBEAM')
        na = self.im.get_namelist_nvars('ABEAMA?')
        self.logger.info(f'NBeams nunits {nunits} ABEAMA? {na} nbeam {nbeams}')
        #if nunits != na:
        #    self.logger.error(f'NBeams mismatch {nunits} != {na}')
        #    raise AssertionError(f'NBeams mismatch {nunits} != {na}')

        atomsDef = atom_base()
        beams = atomsDef.beams

        #D3D: input power, geometry, energy fractions

        nbi.unit.resize(nunits)
        for i in range(nunits):
            nbu = nbi.unit[i]
            ib = i+1
            cm2m = 0.01
            d2rad = constants.pi/180.0
            stn = f'0{ib}' if ib<=9 else f'{ib}'

            nbu.name = f"unit {ib}"

            #find NBeam species info from namelist
            sp_a = self.im.get_namelist_var(f'ABEAMA({ib})', float)
            sp_zn = self.im.get_namelist_var(f'XZBEAMA({ib})', float)

            self.logger.info(f' species: a {sp_a} zn {sp_zn} # {ib}' )
            sp_label = atomsDef.find_species_label(sp_a, sp_zn, beams)
            if self.im.is_valid(sp_a, sp_zn, sp_label):
                nbu.species.a = sp_a
                nbu.species.z_n = sp_zn
                nbu.species.label = sp_label

            #energy EINJ02_E1
            self.im.fill_data_nc(nbu.energy.data, f'EINJ{stn}_E1')
            self.im.fill_data(nbu.energy.time, timesc)

            #injected power PINJ02
            self.im.fill_data_nc(nbu.power_launched.data, f'PINJ{stn}')
            self.im.fill_data(nbu.power_launched.time, timesc)


            # please, ask Marina how this should be handled.
            # TRANSP outputs the energy at full, half and 1/3 of energy
            # It does not output the fraction of power
            # the fraction can be indicated in teh TRANSP namelist, but some namelists do not have it setup.
            # For example, the DIII-D namelists generated with OMFIT have FFULLA=0, FHALFA=0
            # so the fraction must be taken as the default value in TRANSP (1.0 for full fraction)
            # The namelists from JET have these setup
            # For this reason, we cannot take the fraction from the namelist either.
            # one way of doing this is to use the values in the namelist if they are >0
            # and assume that all power is at full energy if the value of FFULLA=0
            # please, note that the IDS wants a fraction value, not the actual power, 
            # so I will need to change the trdat_ids.f90 subroutine to be consistent with what TRANSP expects 
            # (whatever this is)
            #beam_power_fraction
            nbu.beam_power_fraction.data.resize((3, ntimes))
            self.im.fill_data(nbu.beam_power_fraction.time, timesc)
            self.im.fill_data_nc(nbu.beam_power_fraction.data[0], f'PINJ{stn}', resize=False)

            #self.im.fill_data_nc(nbu.beam_power_fraction.data[1],
            #self.im.fill_data_nc(nbu.beam_power_fraction.data[2],

        if nunits > 0:
            nbi = self.im.fill_properties(nbi)
            self.im.fill_code_properties(nbi)
            nbi.ids_properties.homogeneous_time = 1
            nbi.time.resize(ntimes)
            nbi.time = timesc
            self.im.putdb_imas(nbi)
            self.logger.info('Done NBI')

if __name__ == '__main__':
    print("Please use: python convert.py ")
