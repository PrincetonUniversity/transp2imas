

from nc2imas import *

class imas_ec_launchers:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.im.print_debug()
        self.create_ids()

    def get_ec_launch_positions(self):
        #launching_position r. major radius of source. From namelist
        src_r = self.im.get_namelist_vars('XECECH?')
        # height above midplane of source
        src_z = self.im.get_namelist_vars('ZECECH?')
        #Phi location of source (degrees)
        src_phi = self.im.get_namelist_vars('PHECECH?')
        if len(src_r)>0 and len(src_z)>0 and len(src_phi)>0:
            print('type src_r', len(src_r))
            src_r = 1.0e-2 * src_r
            src_z *= 1.0e-2
            src_phi *= scipy.constants.pi/180.
        return src_r, src_z, src_phi

    def create_ids(self):
        ec = self.im.ids.ec_launchers
        #get from namelist
        src_r, src_z, src_phi = self.get_ec_launch_positions()
        freq = self.im.get_namelist_vars('FREQECH?')

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)

        self.im.fill_properties(ec)
        self.im.fill_code_properties(ec)
        ec.ids_properties.homogeneous_time = 1
        ec.time.resize(ntimes)
        ec.time = timesc

        nbeams = self.im.get_nvar_names('PECIN?')
        ec.beam.resize(nbeams)
        #ec.latency =

        self.logger.info(f'#EC beams {nbeams}')
        for i in range(nbeams):
            ib = i+1
            #frequency
            if len(freq)>i:
                self.im.fill_data(ec.beam[i].frequency.data, freq[i])
                #self.im.fill_data(ec.beam[i].frequency.time, freq_time[i])

            #PECIN+1..n(TIME)
            self.im.fill_data_nc(ec.beam[i].power_launched.data, 'PECIN'+str(ib))
            #O-mode ?
            ec.beam[i].mode = 1

            #launching_position of source
            if (len(src_r)>i and len(src_z)>i and len(src_phi)>i):
                #r. major radius of source. XECECH in namelist
                self.im.fill_data(ec.beam[i].launching_position.r, src_r[i])
                # height above midplane of source. ZECECH
                self.im.fill_data(ec.beam[i].launching_position.z, src_z[i])
                #Phi location of source (degrees). PHECECH
                self.im.fill_data(ec.beam[i].launching_position.phi, src_phi[i])

            #ec.beam[i].launching_position.r_limit_min = r_limit_min
            #ec.beam[i].launching_position.r_limit_max = r_limit_max

            #steering_angle_tor [rad]. Azimuthal launch angle
            #PHAIECHn(TIME) DEGREES; Toroidal aiming, GYROTRON n
            #TRANSP defines angles in the laboratory frame. IDS defines angles in the EC frame.
            angle = self.im.get_ncvar('PHAIECH'+str(i+1))
            if angle is not None:
                angle = (angle-180.0)*constants.pi/180.
                self.im.fill_data(ec.beam[i].steering_angle_tor, angle)

            #steering_angle_pol [rad]
            #THETECHn(TIME) DEGREES Poloidal aiming, GYROTRON n
            angle = self.im.get_ncvar('THETECH'+str(i+1))
            if angle is not None:
                angle = (angle-90.0)*constants.pi/180.
                self.im.fill_data(ec.beam[i].steering_angle_pol, angle)


            #spot
            #phase
        self.im.putdb_imas(ec)


if __name__ == '__main__':
    a = get_args(sys.argv)
    nci = nc2imas(db='D3D', shot=147634, backend='MDS', run=1, occur=0, user=os.getenv('USER'),
          new_db=True, ncfile='147634G14.CDF', nml_file='147634G14TR.DAT', log='debug')
    imas_ec_launchers(nci)

