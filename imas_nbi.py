from nc2imas import *
from atoms_class import atom_base

class imas_nbi:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.im.print_debug()
        self.create_ids()

    def valid(self, *vars_) -> bool:
        return self.im.checkValid(vars_)

    def create_ids(self):
        nbi = self.im.ids.nbi
        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        self.im.checkErrExit(timesc)

        nbi = self.im.fill_properties(nbi)
        self.im.fill_code_properties(nbi)
        nbi.ids_properties.homogeneous_time = 1
        nbi.time.resize(ntimes)
        nbi.time = timesc

        nunits = self.im.get_nvar_names('PINJ?')
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
            if self.valid(sp_a, sp_zn, sp_label):
                nbu.species.a = sp_a
                nbu.species.z_n = sp_zn
                nbu.species.label = sp_label

            #energy EINJ02_E1
            self.im.fill_data_nc(nbu.energy.data, f'EINJ{stn}_E1')
            self.im.fill_data(nbu.energy.time, times)

            #injected power PINJ02
            self.im.fill_data_nc(nbu.power_launched.data, f'PINJ{stn}')
            self.im.fill_data(nbu.power_launched.time, times)

            #beam_power_fraction
            nbu.beam_power_fraction.data.resize((3, ntimes))
            self.im.fill_data(nbu.beam_power_fraction.time, times)
            self.im.fill_data_nc(nbu.beam_power_fraction.data[0], f'PINJ{stn}', resize=False)

            #energy/duration.  nml(TBOFFA - TBONA) ?
            #self.im.fill_data_nc(nbu.beam_power_fraction.data[1],
            #self.im.fill_data_nc(nbu.beam_power_fraction.data[2],

            #beam_current_fraction
            nbu.beam_current_fraction.data.resize((3, ntimes))
            self.im.fill_data(nbu.beam_current_fraction.time, times)
            #for j in range(3):
                #self.im.fill_data_nc(nbu.beam_current_fraction.data[j],

            #APERTURE
            nbu.aperture.resize(1)
            #beam aperture shape
            nbu.aperture[0].geometry_type = self.im.get_namelist_var(f'NBSHAPA({ib})', int)

            #see: sharepoint.iter.org ... nbi.html
            # (x1,x2,x3)=(e1,e2,e3); e1(X,Y,Z) https://.. utilities/detector_aperture_coordinates.svg

            #CENTER of 1st aperture ?
            #XRAPOFFA : horizontal offset of beam ctr rel. to ctr of aperture
            #dist of ctr of beam 1st aperture from midplane. Add  XRAPOFFA, XZAPOFFA ?
            nbu.aperture[0].centre.r = self.im.get_namelist_var(f'XYBAPA({ib})', float, cm2m)
            #dist of ctr of beam SCE above/below midplane.
            nbu.aperture[0].centre.z = self.im.get_namelist_var(f'XYBSCA({ib})', float, cm2m)
            #toroidal angle of beam source
            nbu.aperture[0].centre.phi = self.im.get_namelist_var(f'XBZETA({ib})', float, d2rad)

            nbu.aperture[0].x1_unit_vector.x
            nbu.aperture[0].x1_unit_vector.y
            nbu.aperture[0].x1_unit_vector.z
            nbu.aperture[0].x2_unit_vector.x
            nbu.aperture[0].x2_unit_vector.y
            nbu.aperture[0].x2_unit_vector.z
            nbu.aperture[0].x3_unit_vector.x
            nbu.aperture[0].x3_unit_vector.y
            nbu.aperture[0].x3_unit_vector.z
            nbu.aperture[0].x1_width = self.im.get_namelist_var(f'RAPEDGA({ib})', float, cm2m)
            nbu.aperture[0].x2_width = self.im.get_namelist_var(f'XZPEDGA({ib})', float, cm2m)

        self.im.putdb_imas(nbi)
        self.logger.info('Done NBI')

if __name__ == '__main__':
    a = get_args(sys.argv)

    #sh, rn, mach = 142202, 'U01', 'NSTX' #  142202U01.CDF
    sh, rn, mach = 147634, 'G14', 'D3D' # 147634G14.CDF
    sh, rn, mach = 175286, 'T57', 'D3D' # 175286T57.CDF

    ncf = str(sh)+rn+'.CDF'
    nml = str(sh)+rn+'TR.DAT'
    nci = nc2imas(db=mach, shot=sh, backend='MDS', run=4, occur=0, user=os.getenv('USER'),
          new_db=True, ncfile=ncf, nml_file=nml, log='debug')
    imas_nbi(nci)