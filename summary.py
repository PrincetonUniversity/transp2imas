import os
import sys
import numpy as np
import scipy
import re

import imas
from imas_class import Convertor

from scipy.interpolate import interp1d

from pathlib import Path
from netCDF4 import Dataset

class CdfReader:
    def __init__(self, ncfile:str=None):
        self.ds = None
        self.vars_ = []
        self.set_ncfile(ncfile)
        self.open_ncfile()

    def set_ncfile(self, ncfile):
        print('setting cdf file')
        self.ncfile = ncfile

    def open_ncfile(self):
        self.ds = Dataset(self.ncfile,'r')
        self.vars_.clear()
        self.expand_vars(self.ds)

    def print_sample(self):
        for k,v in self.ds.dimensions.items():
            print(k,'\n',v)
            break
        for k,v in self.ds.variables.items():
            print(k,'\n',v)
            break
        print('\n',self.ds.variables['TIME'])

        print('\n', self.ds.variables['TIME'].dimensions)

        print('DIMENSIONS')
        for dim in self.ds.dimensions:
            print(dim)

    def expand_vars(self, gp):
        for key, _ in gp.variables.items():
            self.vars_.append(key)
        for _, sub in gp.groups.items():
            expand_vars(self.vars_, sub)

    def get_nvar_names(self, name, deb:bool=False):
        pattern = r'\b'+name+r'\d{1,2}'
        matches = [self.vars_[i] for i,s in enumerate(self.vars_) if re.search(pattern, s)]
        if deb:
            print(matches)
        return len(matches)

    def get_ncvar(self, var:str):
        try:
            return np.array(self.ds.variables[var][:])
        except:
            return None

    def get_ncattr(name):
        try:
            return self.ncfile.getncattr(name)
        except:
            print(f'{name} not found in {self.ncfile}')
            return None

    def get_ncvarx(self, var:str, xvar:str=None, n:int=-1):
        yv = self.ds.variables[var][:]
        if(numpy.ma.isMaskedArray(yv)):
            yv = np.array(yv)
        if xvar is None:
            return yv

        ydims = self.ds.variables[var].dimensions
        xv = self.ds.variables[xvar][:]
        if(numpy.ma.isMaskedArray(xv)):
            xv = np.array(xv)
        if len(ydims) == 1 or n < 0:
            return yv, xv

        if n > -1:
            if len(ydims) == 2:
                if xvar in ydims[0]:
                    ydimn = 0
                    zvar = ydims[1]
                else:
                    ydimn = 1
                    zvar = ydims[0]

                zz = self.ds.variables[zvar][:]
                xdims = self.ds.variables[xvar].dimensions
                if len(xdims) == 2:
                    if ydimn == 0:
                        xv1 = xv[:][n]
                        yv1 = yv[:][n]
                    else:
                        xv1 = xv[n][:]
                        yv1 = yv[n][:]
                else:
                    xv1 = xv[:]
                    yv1 = yv[n][:]

                return xv1, yv1

            else:
                return xv, yv
        return 0, 0


class Summary:
    def __init__(self, conv:Convertor, ncfile:str=None):
        self.conv = conv
        self.ncfile = ncfile
        self.nc = None
        self.src = None
        self._deb = conv._deb
        self.do_convert()

    def set_resize_len(self, rsz):
        self.nresize = rsz

    def fill_value(self, ids, dat, rsz:int=-1, mult:float=1.0, name:str=None):
        name_ = name if name else 'data'
        if dat is None:
            if self._deb:
                print(name, ' is None')
            return

        nresize = len(dat) if rsz==-1 else rsz

        if nresize > 0:
            ids.value.resize(nresize, refcheck=False)
        ids.value = mult*dat
        #if self.src:
        #    ids.source = self.src

    def fill_value_nc(self, ids, var:str, rsz:int=-1, mult:float=1.0):
        dat = self.nc.get_ncvar(var)
        if dat is None:
            if self._deb:
                print(f'{var} not found')
        else:
            self.fill_value(ids, dat, rsz=rsz, mult=mult, name=var)

    def fill_data_nc(self, ids, var:str, mult:float=1.0):
        dat = self.nc.get_ncvar(var)
        if dat is None:
            if self._deb:
                print(f'{var} not found')
        else:
            ids.resize(len(dat), refcheck=False)
            ids = mult*dat


    #from TRANSP output CDF file
    def convert_summary(self):
        sid = self.conv.outdb.summary
        iter_ = conv.iter_
        sgn = -1 if iter_ else 1
        self.nc = CdfReader(self.ncfile)
        src = 'Calculated'

        coord_x = self.nc.get_ncvar('X')
        coord_xb = self.nc.get_ncvar('XB')
        naxis_x = coord_x.shape[1]
        naxis_xb = coord_xb.shape[1]

        #scalar time
        timesc = self.nc.get_ncvar('TIME')
        ntimes = len(timesc)
        self.set_resize_len(ntimes)

        #profile time
        timep = self.nc.get_ncvar('TIME3')
        ntimep = len(timep)
        if self._deb:
            print(f'ntime_prof {ntimep} ntime_sc {ntimes} nx {naxis_x} nxb {naxis_xb}')
            print('ntimes', ntimes, timesc[:3], '..', timesc[-3:])

        # both TIME and TIME3 from the NETCDF file
        sid.ids_properties.homogeneous_time = 0
        sid.time.resize(ntimes)
        #scalar time
        sid.time = timesc
        sid = self.conv.fill_properties(sid)

        #global_quantities from CDF file
        #ip : Total plasma current (toroidal component).
        #PCURC (TIME)
        # TODO current and magnetic field consistent  with COCOS definition
        self.fill_value_nc(sid.global_quantities.ip, 'PCURC', mult=sgn)
        #b0  BZ(TIME) - BZ @R=RMAJOR OUTSIDE PLASMA
        bz = self.nc.get_ncvar('BZ')
        self.fill_value(sid.global_quantities.b0, bz, mult=sgn)

        #r0 FLT0.
        bzxR = self.nc.get_ncvar('BZXR')
        if bzxR is not None:
            rvalues = bzxR/bz
            # I do not understand the line below
            r0s = 0.5*(rvalues[0]+rvalues[-1])*1.0e-2
            self.fill_value(sid.global_quantities.r0, r0s, rsz=0)

        #li  LI_3(TIME)
        self.fill_value_nc(sid.global_quantities.li, 'LI_3')
        #vloop VSURC(TIME)
        self.fill_value_nc(sid.global_quantities.v_loop, 'VSURC')

        #beta_pol BETAT(TIME)
        self.fill_value_nc(sid.global_quantities.beta_pol, 'BETAT')

        #beta_tor BTEQ(TIME) EQUILIBRIUM BETA(TOROIDAL)
        self.fill_value_nc(sid.global_quantities.beta_tor, 'BTEQ')

        #beta_tor_nor . 100*beta_tor*a[m]*B0[T]/ip[MA]
        btor = sid.global_quantities.beta_tor.value
        b0 = sid.global_quantities.b0.value
        ip = sid.global_quantities.ip.value
        # PAREA(TIME) - plasma CS area[cm^2]
        area = self.nc.get_ncvar('PAREA')*1.e-4
        if area is not None:
            a = np.sqrt(area/scipy.pi)
            data = 1.0e8*btor*a*b0/ip
            self.fill_value(sid.global_quantities.beta_tor_norm, data)

        #beta_tor_thermal_norm

        #energy_total FLT_1D; from UTOTL(TIME3,X) JLES/CM3
        utot = self.nc.get_ncvar('UTOTL')
        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.nc.get_ncvar('DVOL')

        if utot is not None and gvol is not None:
            etot = 3.0/2.0 * np.sum(np.multiply(utot, gvol), axis=1)
            self.fill_value(sid.global_quantities.energy_total, etot)

        #PMHD(TIME3, X). Pressure from equilibrium solver 
        pmhd = self.nc.get_ncvar('PMHD')
        if pmhd is not None and gvol is not None:
            emhd = 1.5 * np.sum(np.multiply(pmhd, gvol), axis=1)
            self.fill_value(sid.global_quantities.energy_mhd, emhd)

        #volume PVOL(TIME) CM^3 : PLASMA VOLUME
        self.fill_value_nc(sid.global_quantities.volume, 'PVOL', mult=1.0e-6)

        #energy_thermal UTHRM(TIME3, X) JLES/CM3 THERMAL ENERGY DENSITY
        utherm = self.nc.get_ncvar('UTHRM')
        if utherm is not None:
            etherm = 1.5 * np.sum(np.multiply(utherm, gvol), axis=1)
            self.fill_value(sid.global_quantities.energy_thermal, etherm)

        #energy_b_field_pol
        #BPOL(TIME3, XB) Tesla
        bpol = self.nc.get_ncvar('BPOL')
        if bpol is not None:
            #TODO interpolate to X
            for ib, bp in enumerate(bpol):
                bpol[ib] = np.interp(coord_x[0], coord_xb[0], bp)
            bpol2 = np.multiply(bpol, bpol)
            ebpol = 1.0e6 * 1.0/(2.0*scipy.constants.mu_0) * np.sum(np.multiply(bpol2, gvol), axis=1)
            self.fill_value(sid.global_quantities.energy_b_field_pol, ebpol)

        # UFASTPP(TIME3, X) JLES/CM3 FAST ION PERP ENERGY DENSITY
        efastpp = self.nc.get_ncvar('UFASTPP')
        if bpol is not None:
            efastpp = 1.5 * np.sum(np.multiply(efastpp, gvol), axis=1)
            self.fill_value(sid.global_quantities.energy_fast_perpendicular, efastpp)

        # UFASTPA(IME3, X) JLES/CM3 FAST ION PLL ENERGY DENSITY
        efastpl = self.nc.get_ncvar('UFASTPA')
        if bpol is not None:
            efastpl = 1.5 * np.sum(np.multiply(efastpl, gvol), axis=1)
            self.fill_value(sid.global_quantities.energy_fast_parallel, efastpl)


        #Plasma diamagnetic energy content = 3/2*integral over the plasma volume of
        #  the total perpendicular pressure [J].
        #PPLAS(TIME3, X) PASCALS PLASMA PRESSURE
        pplas = self.nc.get_ncvar('PPLAS') #TODO
        if pplas is not None:
            edia = 3.0/2.0 * np.sum(np.multiply(pplas, gvol), axis=1)*1.0e-6 #cm3 to m3
            self.fill_value(sid.global_quantities.energy_diamagnetic, edia)


        #denergy_diamagnetic_dt #TODO


        #current_ohm
        # grid area. DAREA (TIME3, X) CM**2 ZONE CROSS SECTIONAL AREA
        grid_area = self.nc.get_ncvar('DAREA')
        #CUROH(TIME3, X)AMPS/CM2
        j_ohmic = self.nc.get_ncvar('CUROH')
        if grid_area is not None and j_ohmic is not None:
            current_ohm = np.sum(np.multiply(j_ohmic, grid_area), axis=1)
            self.fill_value(sid.global_quantities.current_ohm, current_ohm, mult=sgn)

        # current_non_inductive. Toroidal component.
        # CUR(TIME3,X) AMPS/CM2. TOTAL PLASMA CURRENT
        cur_tot = self.nc.get_ncvar('CUR')
        if cur_tot is not None:
            j_nonind =  np.sum(np.multiply(cur_tot - j_ohmic, grid_area), axis=1)
            self.fill_value(sid.global_quantities.current_non_inductive, j_nonind, mult=sgn)

        #tau_energy. TAUA1(TIME), SECONDS, ENERGY CONFINEMENT (TOTAL)
        self.fill_value_nc(sid.global_quantities.tau_energy, 'TAUA1')

        #h_98; H98Y2 (TIME); TauE98y,2 confinement Hfactor
        self.fill_value_nc(sid.global_quantities.h_98, 'H98Y2')

        #h_mode, LHMODE(TIME); H-Mode indicator
        self.fill_value_nc(sid.global_quantities.h_mode, 'LHMODE')

        #q95
        self.fill_value_nc(sid.global_quantities.q_95, 'Q95')

        #Bootstrap current (toroidal component). CURBS(TIME3, X) AMPS/CM2
        # check surface itnegration
        curbs = self.nc.get_ncvar('CURBS')
        if curbs is not None:
            curbsxa = sgn * np.sum(np.multiply(curbs, grid_area), axis=1)
            self.fill_value(sid.global_quantities.current_bootstrap, curbsxa)


        #fusion_gain


        # TAUES(TIME3, XB) PLASMA ENERGY CONFINEMENT
        # this section is incorrect, I am commenting it out for now
        #taues = self.nc.get_ncvar('TAUES')
        #if taues is not None:
        #    for it,ta in enumerate(taues) :
        #        taues[it] = np.interp(coord_x[0], coord_xb[0], ta)
        #    self.fill_value(sid.global_quantities.tau_energy, taues)

        #TAUE98Y2(TIME) TauE98y,2 confinement scaling
        self.fill_value_nc(sid.global_quantities.tau_energy_98, 'TAUE98Y2')

        #ETA_USE(TIME3, X)OHM*CM. RESISTIVITY USED OR INFERRED
        #resist = self.nc.get_ncvar('ETA_USE')
        #if resist is not None:
        #    resist = np.sum(np.multiply(resist, coord_x), axis=1) #TODO
        #    self.fill_value(sid.global_quantities.resistance, resist)

        #POH(TIME3, X) WATTS/CM3, OHMIC HEATING POWER
        poh = self.nc.get_ncvar('POH')
        if poh is not None:
            poh = np.sum(np.multiply(poh, gvol), axis=1)
        self.fill_value(sid.global_quantities.power_ohm, poh)

        #PRAD(TIME3, X) WATTS/CM3 NET RADIATED POWER USED
        prad = self.nc.get_ncvar('PRAD')
        if prad is not None:
            prad = np.sum(np.multiply(prad, gvol), axis=1)
        self.fill_value(sid.global_quantities.power_radiated, prad)

        #PRAD_LI(TIME3, X) WATTS/CM3 , LINE RADIATION
        prad_li = self.nc.get_ncvar('PRAD_LI')
        if prad_li is not None:
            prad_li = np.sum(np.multiply(prad_li, gvol), axis=1)
        self.fill_value(sid.global_quantities.power_line, prad_li)

        #PRAD_BR(TIME3, X) WATTS/CM3  BREMSSTRAHLUNG RADIATION
        prad_br = self.nc.get_ncvar('PRAD_BR')
        if prad_br is not None:
            prad_br = np.sum(np.multiply(prad_br, gvol), axis=1)
        self.fill_value(sid.global_quantities.power_bremsstrahlung, prad_br)

        #PRAD_CY(TIME3, X) WATTS/CM3 CYCLOTRON RADIATION
        prad_cy = self.nc.get_ncvar('PRAD_CY')
        if prad_cy is not None:
            prad_cy = np.sum(np.multiply(prad_cy, gvol), axis=1)
        self.fill_value(sid.global_quantities.power_synchrotron, prad_cy)

        #power_loss through separatrix [W]

        #fusion_fluence


        #si_external_average


        # X and XB represent two different radial coordinates.
        # Need to make sure that the radial coordinate is consistent with the profiles
        # X is sqrt(norm. toroidal flux) defined in zone center
        # XB is sqrt(norm. toroidal flux) defined in zone boundary
        #rho_tor_norm. XB(TIME3, XB).

        xcoord = np.multiply(coord_x, coord_x) #?
        # only one coordinate at one time ?
        #self.fill_value(sid.local.magnetic_axis.rho_tor_norm, xcoord[:,-1]) #TODO


        #magnetic axis position
        self.fill_data_nc(sid.local.magnetic_axis.position.r, 'RAXIS', mult=1.e-2)
        self.fill_data_nc(sid.local.magnetic_axis.position.z, 'YAXIS', mult=1.e-2)

        #magnetic axis n_e, NE(TIME3,X) EV.
        n_e = self.nc.get_ncvar('NE')
        if n_e is not None:
            self.fill_value(sid.local.magnetic_axis.n_e, n_e[:,0], mult=1.0e6)

        #magnetic axis  t_e.  TE(TIME3,X); EV,  ELECTRON TEMPERATURE
        self.fill_value_nc(sid.local.magnetic_axis.t_e, 'TE0')

        #n_i_total
        n_i = self.nc.get_ncvar('NI')
        if n_i is not None:
            self.fill_value(sid.local.magnetic_axis.n_i_total, n_i[:,0], mult=1.0e6)

        #t_i_average  over ion species TODO
        self.fill_value_nc(sid.local.magnetic_axis.t_i_average, 'TI0')
        #zeff
        # if only the value on-axis is needed, we can take ZEFFC
        #zeff = self.nc.get_ncvar('ZEFFP')
        self.fill_value_nc(sid.local.magnetic_axis.zeff, 'ZEFFC')
        #self.fill_value(sid.local.magnetic_axis.zeff, zeff[:,0])
        #q
        self.fill_value_nc(sid.local.magnetic_axis.q, 'Q0')

        #magnetic_shear=rho_tor/q*dq/drho_tor; SHAT(TIME3, X)=(r/q)*(dq/dr)
        # not sure this is correct, I am commenting this out
        #mag = self.nc.get_ncvar('SHAT')
        #if mag is not None:
        #    magsh = mag.copy()
        #    #average of consecutive columns, except for the last
        #    magsh[:,:-1] += mag[:,1:] #add next column
        #    magsh[:,-1] = 3*mag[:,-1] - mag[:,-2] # rdum
        #    magsh /= 2
        #    self.fill_value(sid.local.magnetic_axis.magnetic_shear, magsh[:,0])

        #Ion densities
        #DEUTERIUM  ND(TIME3, X); N/CM**3
        dens_h = self.nc.get_ncvar('ND')
        if dens_h is not None:
            dens_h = 1.0e6 * dens_h
            self.fill_value(sid.local.magnetic_axis.n_i.deuterium, dens_h[:,0])
        #VTORD_NC(TIME3, RMAJM) CM/SEC ; NC D+ toroidal velocity

        #OMEGA(TIME3, X)
        #sid.local.magnetic_axis.velocity_tor.hydrogen


        #runaways.particles. Number of runaway electrons [-]. FLT_1D 1..time.
        # SERUN(TIME3, X) N/CM3/SEC . RUNAWAY ELEC SOURCE RATE
        serun = self.nc.get_ncvar('SERUN')
        if serun is not None:
            run_ptcls = np.sum(np.multiply(serun, gvol), axis=1)
            self.fill_value(sid.runaways.particles, run_ptcls)

        #runaways.current. Parallel current driven by the runaway electrons [A]


        #EC
        #Set of ECRH/ECCD launchers. 6 Max in MDS+.
        #PECIN1 .. PEECHn for n gyrotrons. (TIME3,X) WATTS/CM3
        # please, fix this. The number of gyrotrons should not be hardcoded
        ec_ngyro = self.nc.get_nvar_names('PECIN')
        if self._deb:
            print('# EC gyrotrons', ec_ngyro)
        peech, pecin = 0., 0.
        sid.heating_current_drive.ec.resize(ec_ngyro)
        for i in range(ec_ngyro):
            peechi = self.nc.get_ncvar('PEECH'+str(i+1))
            if peechi is not None:
                peech += peechi
                power = np.sum(np.multiply(peechi, gvol), axis=1)
                self.fill_value(sid.heating_current_drive.ec[i].power, power)

            # launched or injected is PECINn
            pecini = self.nc.get_ncvar('PECIN'+str(i+1))
            if pecini is not None:
                pecin += pecini
                self.fill_value(sid.heating_current_drive.ec[i].power_launched, pecini)

            #ECCURn(TIME3, X) AMPS/CM2  ;ECRH CURRENT (GYROTRON n)
            eccur = self.nc.get_ncvar('ECCUR'+str(i+1))
            if eccur is not None:
                current = np.sum(np.multiply(eccur, grid_area), axis=1)
                self.fill_value(sid.heating_current_drive.ec[i].current, current)

            #Toroidal angle of ECRH at resonance [rad]
            #PHAIECHn(TIME) DEGREES; Toroidal aiming, GYROTRON n
            #TRANSP defines angles in the laboratory frame. IDS defines angles in the EC frame.
            angle = self.nc.get_ncvar('PHAIECH'+str(i+1))
            if angle is not None:
                angle = (angle-180.0)*scipy.pi/180.
                self.fill_value(sid.heating_current_drive.ec[i].angle_tor, angle)

            #Poloidal angle of ECRH at resonance [rad]
            #THETECHn(TIME) DEGREES Poloidal aiming, GYROTRON n
            angle = self.nc.get_ncvar('THETECH'+str(i+1))
            if angle is not None:
                angle = (angle-90.0)*scipy.pi/180.
                self.fill_value(sid.heating_current_drive.ec[i].angle_pol, angle)

        # Total EC power, into vacuum vessel [W]. PECIN(TIME), ECRF INPUT POWER
        if np.any(pecin > 0.):
            self.fill_value(sid.heating_current_drive.power_launched_ec, pecin)
        # Total EC power coupled to plasma. PEECH(TIME3,X) WATTS/CM3 ;ECRH ELECTRON HEATING
        if np.any(peech > 0.):
            power_ec =  np.sum(np.multiply(peech, gvol), axis=1)
            self.fill_value(sid.heating_current_drive.power_ec, power_ec)


        #IC
        # do not hardcode the number of IC antennae
        nic_ant = self.nc.get_nvar_names('PICHA')
        if self._deb:
            print('# IC Antenna', nic_ant)
        picha, pichant = 0., 0.
        sid.heating_current_drive.ic.resize(nic_ant)
        for i in range(nic_ant):
            #PICHAn(TIME)
            pichai = self.nc.get_ncvar('PICHA'+str(i+1))
            if pichai is not None:
                picha += pichai
                self.fill_value(sid.heating_current_drive.ic[i].power_launched, pichai)
            #RFPWRABn  WATTS. Power ABSORBED, Antenna #n  (PICF01N01 WATTS/CM3)
            #  I will doublecheck this
            self.fill_value_nc(sid.heating_current_drive.ic[i].power, 'RFPWRAB'+str(i+1))
            #FREQAn(TIME)
            self.fill_value_nc(sid.heating_current_drive.ic[i].frequency, 'FREQA'+str(i+1))

            pichanti = self.nc.get_ncvar('PICHANT'+str(i+1))
            if pichanti is not None:
                pichant += pichanti

            #PRFFI  WATTS ICRF POWER to FAST IONS

        #  PICHA = 1 ;
        #  PICHANT = 0 ;

        # PICHA(TIME) WATTS. TOTAL ICRF ANTENNA POWER
        if np.any(picha > 0.):
            self.fill_value(sid.heating_current_drive.power_launched_ic, picha)
        # PICHANT (TIME3,X)  WATTS/CM3 . ICRF HEATING by ANTENNA
        if np.any(pichant > 0.):
            pichant = np.sum(np.multiply(pichant, gvol), axis=1)
            self.fill_value(sid.heating_current_drive.power_ic, pichant)




        #MBTOR
        #TOTAL BETA TOROIDAL
        #BTTOT

        return sid

    def do_convert(self):
        print('\nconvert summary')
        ids = self.convert_summary()
        #save IDS
        self.conv.putdb_imas(ids)
        print('done summary')

def getArgs(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-s', '--shot', type=int, help='integer shot number')
    parser.add_argument('-i', '--input-cdf', type=str, help='CDF TRANSP output file')
    parser.add_argument('-b', '--backend', type=str, help='backend: MDS or HDF')
    parser.add_argument('-d', '--db', type=str, help='tokamak or database')
    parser.add_argument('-u', '--user', type=str, help='username or path of database')
    parser.add_argument('-r', '--run', type=int, help='run number')
    parser.add_argument('-o', '--occur', type=int, help='occurrence')
    parser.add_argument('-n', '--new_db', help='new database', action="store_true")
    parser.add_argument('-f', '--default-val', help='use default', action="store_true")

    args = parser.parse_args()
    if args.default_val:
        args = defaults(args)
    if not (args.db or args.shot or args.run or args.occur or args.input_cdf):
        parser.print_help()
        exit(1)

    return args

def defaults(args):
    #fill with default values if not passed in
    if not args.default_val:
        return args
    if args.occur is None:
        args.occur = 0
    if args.verbose is None:
        args.verbose = False
    if args.new_db is None:
        args.new_db = True
    if args.run is None:
        args.run = 0 #2302 #W02 #this is occur
    if args.shot is None:
        args.shot = 147634 #200001 #
    if args.db is None:
        args.db = 'D3D' #NSTX-U
    if args.user is None:
        args.user = os.getcwd() # os.getenv('USER') # or '~/public/imasdb'
    if args.input_cdf is None:
        args.input_cdf = '147634T16.CDF'
    if args.backend == 'MDS':
        args.backend = imas.imasdef.MDSPLUS_BACKEND
    else:
        args.backend = imas.imasdef.HDF5_BACKEND
    return args


if __name__ == '__main__':
    import argparse
    a = getArgs(sys.argv)

    print(a.db, a.shot, a.run, a.occur, a.new_db, a.user, a.backend, a.verbose,a.input_cdf)

    conv = Convertor(a.db, a.shot, a.run, a.occur, a.new_db, a.user, a.backend, _deb=a.verbose)
    Summary(conv, ncfile=a.input_cdf)
    conv.close_imas_database()

