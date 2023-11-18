

from nc2imas import *

class imas_ec_launchers:
    def __init__(self, ncimas:nc2imas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()

    def create_ids(self):
        self.im.logger.info('\nec_launchers')

        ec = self.im.ids.ec_launchers
        cm2m = 0.01
        pi = self.im.pi

        #scalar time
        timesc = self.im.get_ncvar('TIME')
        ntimes = len(timesc)
        self.im.checkErrExit(timesc)

        ec = self.im.fill_properties(ec)
        ec = self.im.fill_code_properties(ec)
        ec.ids_properties.homogeneous_time = 1
        ec.time.resize(ntimes)
        ec.time = timesc

        nbeams = self.im.get_var_names_re('PECIN?', nd1=1, nd2=2, length=True)
        ec.beam.resize(nbeams)

        src_r = np.zeros(nbeams)
        src_z = np.zeros(nbeams)
        src_phi = np.zeros(nbeams)


        # TO BE REMOVED ---------
        if False:
           #launching_position r. major radius of source. From namelist
           src_r = self.im.get_namelist_vars('XECECH', float, scale=cm2m)
           # height above midplane of source
           src_z = self.im.get_namelist_vars('ZECECH', float, scale=cm2m)
           #Phi location of source (degrees). If single value, resize it to array
           src_phi = self.im.get_namelist_vars('PHECECH', float,
                     scale=constants.pi/180.0, size=nbeams, resize=True)
           # remove code above this line (X, Z, phi)
           if not self.im.is_valid(src_r, src_z, src_phi, ind=nbeams-1):
               src_r = src_z = src_phi = np.zeros(nbeams)
           # ---------TO BE REMOVED

        self.logger.info(f'#EC beams {nbeams}')
        for i in range(nbeams):
            ib = i+1
            #frequency. search either FREQECH or FREQECH(i)
            freq = self.im.get_namelist_vars(f'FREQECH', float, ind=i)
            if freq is None:
               freq = 0.0

            freq = np.tile(np.atleast_1d(freq), ntimes)
            self.im.fill_data(ec.beam[i].frequency.data, freq)
#           self.im.fill_data(ec.beam[i].frequency.time, timesc)

            #PECIN+1..n(TIME)
            self.im.fill_data_nc(ec.beam[i].power_launched.data, 'PECIN'+str(ib))
#            self.im.fill_data(ec.beam[i].power_launched.time, timesc)

            if self.im.is_valid(src_r, src_z, src_phi, ind=i):
                #r. major radius of source. XECECH in namelist
                self.im.fill_data(ec.beam[i].launching_position.r, src_r[i])
                # height above midplane of source. ZECECH
                self.im.fill_data(ec.beam[i].launching_position.z, src_z[i])
                #Phi location of source (degrees). PHECECH
                self.im.fill_data(ec.beam[i].launching_position.phi, src_phi[i])

            #O-mode ?
#            ec.beam[i].mode = 1

            # fmp needs to check the code below for a number of cases with different B, I geometry
            #steering_angle_tor [rad]. Azimuthal launch angle
            #TRANSP defines angles in the laboratory frame +R to -R. 90 deg along +Phi
            #IDS defines angles in the EC frame, increasing towards the positive phi axis.
            # the formula for conversion can be found in the paper by Prater, Eq. 
            # theta=pi-acos(cos(beta)*sin(alpha));phi=pi+asin(sin(beta)/sin(theta))

            eps = 1.0e-10
            #PHAIECHn(TIME) DEGREES; Toroidal aiming, GYROTRON n
            phi = self.im.get_ncvar('PHAIECH'+str(i+1), scale=pi/180.0)
            theta = self.im.get_ncvar('THETECH'+str(i+1), scale=pi/180.0)
            if phi is None:
               phi = 0.0
            if theta is None:
               theta = 0.0
            # check for None in case the above is removed
            if phi is not None and theta is not None:
               beta = np.arcsin(np.sin(theta)*np.sin(phi-pi))
               cosbeta = np.cos(beta)
               alpha = np.where(abs(cosbeta)>eps, np.arcsin(np.cos(pi-theta)/cosbeta), 0)
               self.im.fill_data(ec.beam[i].steering_angle_tor, beta)
               self.im.fill_data(ec.beam[i].steering_angle_pol, alpha)

        self.im.putdb_imas(ec)
        self.logger.info('Done')


if __name__ == '__main__':
    print("Please use: python convert.py ")


