
from nc2imas import *

class imas_lh_antennas:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.im.print_debug()
        self.create_ids()

    def create_ids(self):
        lha = self.im.ids.lh_antennas
        cm2m = 0.01
        sgn = -1 if self.im.iter_ else 1

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        timep = self.im.get_ncvar('TIME3')
        ntimep = len(timep)
        self.im.checkErrExit(timesc, timep)

        lha = self.im.fill_properties(lha)
        self.im.fill_code_properties(lha)
        lha.ids_properties.homogeneous_time = 1
        lha.time.resize(ntimes)
        lha.time = timesc

        #ref point to define the poloidal angle, e.g. geometrical centre of vacuum vessel
        raxis = self.im.get_ncvar('RAXIS', cm2m)
        zaxis = self.im.get_ncvar('YAXIS', cm2m)
        if self.im.checkValid(raxis, zaxis):
            lha.reference_point.r = raxis[0].astype(float)
            lha.reference_point.z = zaxis[0].astype(float)

        nant = self.im.get_nvar_names('PICHA?')
        power = []

        lha.antenna.resize(nant)
        for i in range(nant):
            ia = i+1
            sta = f'{ia}'

            lha.antenna[i].name = f'LH Antenna{sta}'
            lha.antenna[i].identifier = f'LH Antenna{sta}'

            #frequency
            freq = self.im.get_namelist_var('FREQLH', ia)
            if self.im.checkValid(freq):
                lha.antenna[i].frequency = freq

            #position
            pos = self.imget_namelist_var('RHOPSI0LH', ia)

            #power launced
            self.im.fill_data_nc(lha.antenna[i].power_launched.time, timesc)
            self.im.fill_data_nc(lha.antenna[i].power_launched.data, f'PLHANT{sta}')

            #power_forward. Forward power arriving to the back of the antenna


            #power_reflected = launched - absorbed
            plaunch = self.im.get_ncvar(f'PLHANT{sta}')
            pabs = self.im.get_ncvar(f'PLHABS{sta}')
            if self.im.checkValid(plaunch, pabs):
                prefl = plaunch - pabs
                self.im.fill_data(lha.antenna[i].power_reflected.data, prefl)
                self.im.fill_data(lha.antenna[i].power_reflected.time, times)
            else:
                self.logger.info(f'Error: power_reflected. invalid: plaunch or pabs')

            #surface_current on the antenna straps and on passive components, for every time slice

            #module.  Set of antenna modules (each module is fed by a single transmission line)



        #power launched
        self.im.fill_data(lha.power_launched.time, timesc)
        self.im.fill_data_nc(lha.power_launched.data, 'PLH')

        self.im.putdb_imas(lha)
        self.logger.info('Done')


if __name__ == '__main__':
    a = get_args(sys.argv)

    sh, rn, mach = 33618, 'Y19', 'ITER' # 33618Y19.CDF

    ncf = str(sh)+rn+'.CDF'
    nml = str(sh)+rn+'TR.DAT'
    nci = nc2imas(db=mach, shot=sh, backend='MDS', run=6, occur=0, user=os.getenv('USER'),
          new_db=True, ncfile=ncf, nml_file=nml, log='debug')
    imas_lh_antennas(nci)