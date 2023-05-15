
from nc2imas import *

class imas_ic_antennas:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.im.print_debug()
        self.create_ids()

    def create_ids(self):
        ica = self.im.ids.ic_antennas
        cm2m = 0.01
        sgn = -1 if self.im.iter_ else 1

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        timep = self.im.get_ncvar('TIME3')
        ntimep = len(timep)
        self.im.checkErrExit(timesc, timep)

        ica = self.im.fill_properties(ica)
        self.im.fill_code_properties(ica)
        ica.ids_properties.homogeneous_time = 1
        ica.time.resize(ntimes)
        ica.time = timesc

        #ref point to define the poloidal angle, e.g. geometrical centre of vacuum vessel
        raxis = self.im.get_ncvar('RAXIS', cm2m)
        zaxis = self.im.get_ncvar('YAXIS', cm2m)
        if self.im.checkValid(raxis, zaxis):
            ica.reference_point.r = raxis[0].astype(float)
            ica.reference_point.z = zaxis[0].astype(float)

        nant = self.im.get_nvar_names('PICHA?')
        power = []

        ica.antenna.resize(nant)
        for i in range(nant):
            ia = i+1
            sta = f'{ia}'

            ica.antenna[i].name = f'IC Antenna{sta}'
            ica.antenna[i].identifier = f'IC Antenna{sta}'

            #frequency
            self.im.fill_data(ica.antenna[i].frequency.time, timesc)
            self.im.fill_data_nc(ica.antenna[i].frequency.data, f'FREQA{sta}')

            #power launced
            self.im.fill_data_nc(ica.antenna[i].power_launched.time, timesc)
            self.im.fill_data_nc(ica.antenna[i].power_launched.data, f'PICHA{sta}')

            #power_forward. Forward power arriving to the back of the antenna


            #power_reflected = launched - absorbed
            plaunch = self.im.get_ncvar(f'PICHA{sta}')
            pabs = self.im.get_ncvar(f'RFPWRAB{sta}')
            if self.im.checkValid(plaunch, pabs):
                prefl = plaunch - pabs
                self.im.fill_data(ica.antenna[i].power_reflected.data, prefl)
                self.im.fill_data(ica.antenna[i].power_reflected.time, timesc)


            #surface_current on the antenna straps and on passive components, for every time slice

            #module.  Set of antenna modules (each module is fed by a single transmission line)



        #power launched
        self.im.fill_data(ica.power_launched.time, timesc)
        self.im.fill_data_nc(ica.power_launched.data, 'PICHTOT')

        self.im.putdb_imas(ica)
        self.logger.info('Done')


if __name__ == '__main__':
    a = get_args(sys.argv)

    sh, rn, mach = 89846, '', 'CMOD' # 89846.CDF

    ncf = str(sh)+rn+'.CDF'
    nml = str(sh)+rn+'TR.DAT'
    nci = nc2imas(db=mach, shot=sh, backend='MDS', run=5, occur=0, user=os.getenv('USER'),
          new_db=True, ncfile=ncf, nml_file=nml, log='debug')
    imas_ic_antennas(nci)