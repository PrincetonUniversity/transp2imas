

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

        freq = self.im.get_namelist_vars('FREQECH?', float)

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
        #ec.latency =

        self.logger.info(f'#EC beams {nbeams}')
        for i in range(nbeams):
            ib = i+1
            #frequency
            if self.im.is_valid(freq, ind=i):
                self.im.fill_data(ec.beam[i].frequency.data, freq[i])
                self.im.fill_data(ec.beam[i].frequency.time, timesc)

            #PECIN+1..n(TIME)
            self.im.fill_data_nc(ec.beam[i].power_launched.data, 'PECIN'+str(ib))
            self.im.fill_data(ec.beam[i].power_launched.time, timesc)

            #O-mode ?
            ec.beam[i].mode = 1

            #ec.beam[i].launching_position.

            # fmp needs to check the code below for a number of cases with different B, I geometry
            #steering_angle_tor [rad]. Azimuthal launch angle
            #PHAIECHn(TIME) DEGREES; Toroidal aiming, GYROTRON n
            #TRANSP defines angles in the laboratory frame +R to -R. 90 deg along +Phi
            #IDS defines angles in the EC frame ?? increasing towards the positive phi axis.
            angle = self.im.get_ncvar('PHAIECH'+str(i+1))
            if angle is not None:
                angle = (angle-180.0)*constants.pi/180.
                self.im.fill_data(ec.beam[i].steering_angle_tor, angle)
            #If angle-180.0 then inconsistency between tr_imas_mod.f90:468, 470
            #deg2rad*phaiech_use(iant), deg2rad*thetech_use(iant)

            #steering_angle_pol [rad]
            #THETECHn(TIME) DEGREES Poloidal aiming, GYROTRON n
            #TRANSP: calculated with respect to the positive Z axis. THETECH=90.0 points horizontally,
            # THETECH=180.0 points downwards.
            #IMAS: in the R,Z plane (from the -R axis towards the -Z axis). Not in detector_aperture_coordinates
            angle = self.im.get_ncvar('THETECH'+str(i+1))
            if angle is not None:
                angle = (angle-90.0)*constants.pi/180.
                #angle = -(angle+90.0)*constants.pi/180.
                self.im.fill_data(ec.beam[i].steering_angle_pol, angle)


            #spot
            #phase
        self.im.putdb_imas(ec)


if __name__ == '__main__':
    print("Please use: python convert.py ")


