

from nc2imas import *

class imas_ec_launchers:
    def __init__(self, ncimas:nc2imas, ec_nml_launch=False):
        self.im = ncimas
        self.logger = ncimas.logger
        self.ec_nml_launch = ec_nml_launch
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

        def get_namelist_data(name, tp, ind, scale=1, ntimes=None, fill=False):
            dat = self.im.get_namelist_vars(name=name, tp=tp, scale=scale, ind=ind)
            # resize option is suitable when a scalar stand for array, not for filling times
            if dat is None:
               dat = 0.0
            if fill:
               arr = False
               try:
                  n = len(dat)
                  if n>0:
                     arr = True
               except:
                  pass

               if arr:
                  print(f'Error in resizing namelist {name}; already array {dat}')
                  raise Exception('failed to get {name} from namelist')
               if not arr:
                  dat = np.tile(np.atleast_1d(dat), ntimes)
            return dat


        self.logger.info(f'#EC beams {nbeams}')
        for i in range(nbeams):
            ib = i+1
            #frequency. search either FREQECH or FREQECH(i)
            freq = get_namelist_data('FREQECH', float, i, ntimes=ntimes, fill=True)
            self.im.fill_data(ec.beam[i].frequency.data, freq)
#           self.im.fill_data(ec.beam[i].frequency.time, timesc)

            #PECIN+1..n(TIME)
            self.im.fill_data_nc(ec.beam[i].power_launched.data, 'PECIN'+str(ib))
#            self.im.fill_data(ec.beam[i].power_launched.time, timesc)

            # Only for testing. Disabled by default
            if self.ec_nml_launch:
               if i==0:
                  print('\t..using launch data from namelist')
               #r. major radius of source. XECECH in namelist
               launchx = get_namelist_data('XECECH', float, i, cm2m, ntimes, fill=True)
               self.im.fill_data(ec.beam[i].launching_position.r, launchx, name='launchX')
               # height above midplane of source. ZECECH
               launchz = get_namelist_data('ZECECH', float, i, cm2m, ntimes, fill=True)
               self.im.fill_data(ec.beam[i].launching_position.z, launchz, name='launchZ')
               #Phi location of source (degrees). PHECECH
               launch_phi = get_namelist_data('PHECECH', float, ind=-1, scale=constants.pi/180.0,
                                                ntimes=ntimes, fill=True)
               self.im.fill_data(ec.beam[i].launching_position.phi, launch_phi, name='launchPHI')

               angle = 0.0
               angles = np.full(ntimes, angle)
               self.im.fill_data(ec.beam[i].spot.angle, angles)
               si, sj = 0, 1
               if math.cos(angle)**2>0.5 :
                  si, sj = 1, 0
               # initial principal horizontal beam width
               spot_hor = get_namelist_data('ECB_WIDTH_HOR', float, i, cm2m, ntimes, fill=True)
               spot_vert = get_namelist_data('ECB_WIDTH_VERT', float, i, cm2m, ntimes, fill=True)
               curv_hor = get_namelist_data('ECB_CURV_HOR', float, i, cm2m, ntimes, fill=True)
               curv_vert = get_namelist_data('ECB_CURV_VERT', float, i, cm2m, ntimes, fill=True)
               if self.im.is_valid(spot_hor, spot_vert, ind=ntimes-1):
                  ec.beam[i].spot.size.resize((2,ntimes));
                  ec.beam[i].spot.size[si,:] = spot_hor
                  ec.beam[i].spot.size[sj,:] = spot_vert

               if self.im.is_valid(curv_hor, curv_vert, ind=ntimes-1) and \
                  np.all(curv_hor) and np.all(curv_vert) :
                  curv_hor = -1.0/curv_hor
                  curv_vert = -1.0/curv_vert
                  ec.beam[i].phase.curvature.resize((2,ntimes));
                  ec.beam[i].phase.curvature[si,:] = curv_hor
                  ec.beam[i].phase.curvature[sj,:] = curv_vert
            # --- above for testing only ---

            # using mode from namelist
            rfmode = get_namelist_data('RFMODECH', float, i)
            if rfmode is not None:
               # extra ordinary mode -1; ordinary mode 1
               ec.beam[i].mode = -1 if (rfmode<0.01) else 1

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


