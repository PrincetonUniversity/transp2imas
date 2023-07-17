
from nc2imas import *

class imas_ic_antennas:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger

        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\nic_antennas')
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
        if self.im.is_valid(raxis, zaxis):
            ica.reference_point.r = raxis[0].astype(float)
            ica.reference_point.z = zaxis[0].astype(float)

        nant = self.im.get_var_names_re('PICHA?', nd1=1, nd2=2, length=True)
        self.logger.info(f'ic_antenna {nant}')

        power = []

        ica.antenna.resize(nant)
        for i in range(nant):
            ant = ica.antenna[i]
            ia = i+1
            sta = f'{ia}'

            ant.module.resize(1)
            mod = ant.module[0]
            si = 0 # only 1 strap

            # don't fill geometry from namelist
            if False:
              mod.strap.resize(1)
              #dist to Faraday shield
              dFar = self.im.get_namelist_var(rf'RFARTR_A({ia})', float, cm2m)
              # regex need escaped \(
              straps = self.im.get_namelist_keys_re(rf'RGEOANT_A\(?,{ia}\)', cut='RGEO', nd1=1, nd2=3)

              rs, ys, phis = [], [], []
              for sname in straps:
                  rs.append( self.im.get_namelist_var(rf'RGEO{sname}', float, cm2m) )
                  ys.append( self.im.get_namelist_var(rf'YGEO{sname}', float, cm2m) )
                  phis.append(0.)

              mod.strap[si].distance_to_conductor = dFar
              mod.strap[si].outline.r = np.array(rs, dtype=float)
              mod.strap[si].outline.z = np.array(rs, dtype=float)
              mod.strap[si].outline.phi = np.array(phis, dtype=float)

            ant.name = f'IC Antenna{sta}'
            ant.identifier = f'IC Antenna{sta}'

            #frequency
            self.im.fill_data(mod.frequency.time, timesc)
            self.im.fill_data_nc(mod.frequency.data, f'FREQA{sta}')

            #power launced
            self.im.fill_data(mod.power_launched.time, timesc)
            self.im.fill_data_nc(mod.power_launched.data, f'PICHA{sta}')

            #power_forward. Forward power arriving to the back of the antenna


            #power_reflected = launched - absorbed
            plaunch = self.im.get_ncvar(f'PICHA{sta}')
            pabs = self.im.get_ncvar(f'RFPWRAB{sta}')
            if self.im.is_valid(plaunch, pabs):
                prefl = plaunch - pabs
                self.im.fill_data(ant.power_reflected.data, prefl)
                self.im.fill_data(ant.power_reflected.time, timesc)


            #surface_current on the antenna straps and on passive components, for every time slice

            #module.  Set of antenna modules (each module is fed by a single transmission line)


        if nant > 0:
            #power launched
            self.im.fill_data(ica.power_launched.time, timesc)
            self.im.fill_data_nc(ica.power_launched.data, 'PICHTOT')

            self.im.putdb_imas(ica)
            self.logger.info('Done')


if __name__ == '__main__':
    print("Please use: python convert.py ")