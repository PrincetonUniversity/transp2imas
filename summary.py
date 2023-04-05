import os
import sys
import numpy as np
import scipy

import imas
from imas_class import Convertor

from scipy.interpolate import interp1d

from pathlib import Path
from netCDF4 import Dataset

class CdfReader:
    def __init__(self, ncfile:str=None):
        self.ds = None
        self.set_ncfile(ncfile)
        self.open_ncfile()

    def set_ncfile(self, ncfile):
        print('setting cdf file')
        self.ncfile = ncfile

    def open_ncfile(self):
        self.ds = Dataset(self.ncfile,'r')


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
        for dim in self.ds.dimensionss():
            print(dim)

    def get_ncvar(self, var:str):
        try:
            return np.array(self.ds.variables[var][:])
        except:
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
    def __init__(self, conv:Convertor, occur:int=None, plot:bool=False,
      ncfile:str=None, iter_:bool=False, _deb:bool=False):
        self.conv = conv
        self.ncfile = ncfile
        self.nc = None
        self.src = None
        self.nresize = 0
        self.iter_ = iter_
        self._deb = _deb
        self._plot = plot
        self.occur = occur if occur is not None else conv.occur
        self.do_convert()

    # why do we need mds_tree and mds_node?
    def fill_ids_interp_mds(self, times, ids, mds_tree, mds_node):
        try:
            dat, vtimes = self.conv.get_mds_data(mds_tree, mds_node, dim=True)
        except:
            print(f'{mds_node} not found in {self.mds_shot}.{mds_tree}')
            return
        dat = np.interp(times, vtimes, dat)
        self.fill_value(ids, dat, len(times))

    def set_resize_len(self, rsz):
        self.nresize = rsz

    def fill_value(self, ids, dat, rsz:int=-1, mult:float=1.0, name:str=None):
        name_ = name if name else 'data'
        if dat is None:
            if self._deb:
                print(name, ' is None')
            return

        #nresize = self.nresize if rsz==-1 else rsz
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
            # wrong spelling in line below
            print('poistion', var, ids)

    #from TRANSP output CDF file
    def convert_summary(self):
        occur = 0
        sid = self.conv.outdb.summary
        sgn = -1 if iter_ else 1
        self.nc = CdfReader(self.ncfile)
        src = 'Calculated'

        naxis_x = self.nc.get_ncvar('X').shape[1]
        naxis_xb = self.nc.get_ncvar('XB').shape[1]

        #scalar time
        times = self.nc.get_ncvar('TIME')
        ntimes = len(times)
        self.set_resize_len(ntimes)

        #profile time
        timep = self.nc.get_ncvar('TIME3')
        ntimep = len(timep)
        if self._deb:
            print(f'ntime_prof {ntimep} ntime_sc {ntimes} nx {naxis_x} nxb {naxis_xb}')
            print('ntimes', ntimes, times)

        # this is incorrect if you are taking both TIME and TIME3 from the NETCDF file
        sid.ids_properties.homogeneous_time = 1
        sid.time.resize(ntimes)
        #scalar time
        sid.time = times
        sid = self.conv.fill_properties(sid)
        #global_quantities from CDF file
        #ip : Total plasma current (toroidal component).
        #PCURC (TIME)

        # need to check that the current and magnetic field are written consistently with the COCOS definition 
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
        area = self.nc.get_ncvar('PAREA')*1.e-4 # PAREA(TIME) - plasma CS area[cm^2]
        if area is not None:
            a = np.sqrt(area/scipy.pi)
            print(type )
            data = 1.0e8*btor*a*b0/ip
            self.fill_value(sid.global_quantities.beta_tor_norm, data)

        #energy_total FLT_1D; from UTOTL(TIME3,X) JLES/CM3
        utot = 1.0e6*self.nc.get_ncvar('UTOTL')
        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = 1.0e-6*self.nc.get_ncvar('DVOL')

        if utot is not None and gvol is not None:
            etot = np.sum(np.multiply(utot, gvol), axis=1)
            self.fill_value(sid.global_quantities.energy_total, etot)

        #volume PVOL(TIME) CM^3 : PLASMA VOLUME
        pvol = 1.0e-6 * self.nc.get_ncvar('PVOL')
        if pvol is not None:
            self.fill_value(sid.global_quantities.volume, pvol)

        #energy_thermal UTHRM(TIME3, X) JLES/CM3 THERMAL ENERGY DENSITY
        utherm = 1.0e6*self.nc.get_ncvar('UTHRM')
        if utherm is not None:
            etherm = np.sum(np.multiply(utherm, gvol), axis=1)
            self.fill_value(sid.global_quantities.energy_thermal, etherm)

        #Plasma diamagnetic energy content = 3/2*integral over the plasma volume of
        #  the total perpendicular pressure [J].
        #PPLAS(TIME3, X) PASCALS PLASMA PRESSURE
        pplas = self.nc.get_ncvar('PPLAS') #TODO
        if pplas is not None:
            edia = 3.0/2.0 * np.sum(np.multiply(pplas, gvol), axis=1)*1.0e-6 #cm3 to m3
            self.fill_value(sid.global_quantities.energy_diamagnetic, edia)

        #current_ohm
        # grid area. DAREA (TIME3, X) CM**2 ZONE CROSS SECTIONAL AREA
        grid_area = 1.0e-4 * self.nc.get_ncvar('DAREA')
        #CUROH(TIME3, X)AMPS/CM2
        j_ohmic = 1.0e4 * self.nc.get_ncvar('CUROH')
        if grid_area is not None and j_ohmic is not None:
            current_ohm = np.sum(np.multiply(j_ohmic, grid_area), axis=1)
            self.fill_value(sid.global_quantities.current_ohm, current_ohm, mult=sgn)

        #PLCURTOT(TIME3, XB) ;AMPS; TOTAL POLOIDAL CUR TO WALL
        #POLOIDAL CUR : PLCURPRP, PLCURPLL, PLCURTOT

        # current_non_inductive. Toroidal component.
        # CUR(TIME3,X) AMPS/CM2. TOTAL PLASMA CURRENT
        cur_tot = self.nc.get_ncvar('CUR')
        if cur_tot is not None:
            j_nonind = 1.0e4 * np.sum(np.multiply(cur_tot - j_ohmic, grid_area), axis=1)
            print(cur_tot.shape, grid_area.shape, j_nonind.shape)
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
            curbsxa = 1.0e4 * sgn * np.sum(np.multiply(curbs, grid_area), axis=1)
            self.fill_value(sid.global_quantities.current_bootstrap, curbsxa)


        #rho_tor_norm. XB(TIME3, XB). RHOVN
        # X and XB represent two different radial coordinates.
        # Need to make sure that the radial coordinate is consistent with the profiles
        # X is sqrt(norm. toroidal flux) defined in zone center
        # XB is sqrt(norm. toroidal flux) defined in zone boundary
        xb = self.nc.get_ncvar('XB')
        #self.fill_value(sid.local.magnetic_axis.rho_tor_norm, xb[:,-1]) #TODO
        # the line below should be removed
        #rho_tor= sqrt(b_flux_tor/(pi*b0)) ~ sqrt(pi*r^2*b0/(pi*b0)) ~ r [m].

        #self.fill_value_nc(sid.local.magnetic_axis.rho_tor, rho_tor)


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
            n_i_over_ne = 1.0e6 * np.divide(n_i, n_e)
            self.fill_value(sid.local.magnetic_axis.n_i_total, n_i_over_ne[:,0])

        #t_i_average
        self.fill_value_nc(sid.local.magnetic_axis.t_i_average, 'TI0')
        #zeff
        # if only the value on-axis is needed, we can take ZEFFC
        zeff = self.nc.get_ncvar('ZEFFP')
        self.fill_value(sid.local.magnetic_axis.zeff, zeff[:,0])
        #q
        self.fill_value_nc(sid.local.magnetic_axis.q, 'Q0')

        #magnetic_shear=rho_tor/q*dq/drho_tor; SHAT(TIME3, X)=(r/q)*(dq/dr)
        # not sure this is correct
        mag = self.nc.get_ncvar('SHAT')
        if mag is not None:
            magsh = mag.copy()
            #average of consecutive columns, except for the last
            magsh[:,:-1] += mag[:,1:] #add next column
            magsh[:,-1] = 3*mag[:,-1] - mag[:,-2] # rdum
            magsh /= 2
            self.fill_value(sid.local.magnetic_axis.magnetic_shear, magsh[:,0])

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
            run_ptcls = 1.0e6 * np.sum(np.multiply(serun, gvol), axis=1)
            self.fill_value(sid.runaways.particles, run_ptcls)

        #runaways.current. Parallel current driven by the runaway electrons [A]


        #EC
        #Set of ECRH/ECCD launchers. 6 Max. in MDS+.
        #PECIN1 .. PEECH6,  6 gyrotrons. (TIME3, X) WATTS/CM3
        # do not harcode the number of gyrotrons. The IDS needs to be allocated for the number of gyrotron used.
        ec_ngyro = 6
        peech, pecin = 0., 0.
        sid.heating_current_drive.ec.resize(ec_ngyro)
        for i in range(ec_ngyro):
            # the power launched is PECINx
            peechi = self.nc.get_ncvar('PEECH'+str(i))
            if peechi is not None:
                peech += peechi
                powerl = 1.0e6 * np.sum(np.multiply(peechi, gvol), axis=1)
                self.fill_value(sid.heating_current_drive.ec[i].power_launched, powerl)

            # this is wrong. PECIN is the input power, not the absorbed power.
            # We need to integrate PEECH over volume to get the total power absorbed.
            pecini = self.nc.get_ncvar('PECIN'+str(i))
            if pecini is not None:
                pecin += pecini
                self.fill_value_nc(sid.heating_current_drive.ec[i].power, pecini)

            #ECCURn(TIME3, X) AMPS/CM2  ;ECRH CURRENT (GYROTRON n)
            eccur = self.nc.get_ncvar('ECCUR'+str(i))
            if eccur is not None:
                current = 1.0e4 * np.sum(np.multiply(eccur, grid_area), axis=1)
                self.fill_value(sid.heating_current_drive.ec[i].current, current)

            #Toroidal angle of ECRH at resonance [rad]
            #PHAIECHn(TIME) DEGREES; Toroidal aiming, GYROTRON n
            # this is wrong. TRANSP defines angles in the laboratory frame. IDS defines angles in the EC frame.
            self.fill_value_nc(sid.heating_current_drive.ec[i].angle_tor, 'PHAIECH'+str(i), scipy.pi/180.0)

            #Poloidal angle of ECRH at resonance [rad]
            #THETECHn(TIME) DEGREES Poloidal aiming, GYROTRON n
            # this is wrong. TRANSP defines angles in the laboratory frame. IDS defines angles in the EC frame.
            self.fill_value_nc(sid.heating_current_drive.ec[i].angle_tor, 'THETECH'+str(i), scipy.pi/180.0)

        # Total EC power, into vacuum vessel [W]. PECIN(TIME), ECRF INPUT POWER
        if np.any(pecin > 0.):
            self.fill_value(sid.heating_current_drive.power_launched_ec, pecin)
        # Total EC power coupled to plasma. PEECH(TIME3,X) WATTS/CM3 ;ECRH ELECTRON HEATING
        # this is wrong. The power coupled to the plasma is PEECH, PECIN is the total injected power.
        # PECIN1, 2, ... is the power injected from individual gyrotrons
        if np.any(peech > 0.):
            power_ec =  1.0e6 * np.sum(np.multiply(peech, gvol), axis=1)
            self.fill_value(sid.heating_current_drive.power_ec, power_ec)


        #IC
        # do not hardcode the number of IC antennae
        nic_ant = 2
        picha, pichant = 0., 0.
        sid.heating_current_drive.ic.resize(nic_ant)
        for i in range(nic_ant):
            #PICHAn(TIME)
            pichai = self.nc.get_ncvar('PICHA'+str(i))
            if pichai is not None:
                picha += pichai
                self.fill_value(sid.heating_current_drive.ic[i].power_launched, pichai)
            #RFPWRABn  WATTS. Power ABSORBED, Antenna #n  (PICF01N01 WATTS/CM3)
            #  I will doublecheck this 
            self.fill_value_nc(sid.heating_current_drive.ic[i].power, 'RFPWRAB'+str(i))
            #FREQAn(TIME)
            self.fill_value_nc(sid.heating_current_drive.ic[i].frequency, 'FREQA'+str(i))

            pichanti = self.nc.get_ncvar('PICHANT'+str(i))
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
            pichant = 1.0e6 * np.sum(np.multiply(pichant, gvol), axis=1)
            self.fill_value(sid.heating_current_drive.power_ic, pichant)




        #MBTOR
        #TOTAL BETA TOROIDAL
        #BTTOT

        return sid

    def do_convert(self):
        print('\nconvert summary')
        ids = self.convert_summary()
        #save IDS
        self.conv.putdb_imas(ids, self.occur)
        print('done summary')

# not needed
#call this only to create separate instance of conversion
def run_convertor(db, shot, run, occur:int=0, mds_shot:int=None, new_db:bool=False,
      _deb:bool=False, user:str=None, backend:str=None, plot:bool=False,
       ncfile:str=None, iter_:bool=False):
    #create convertor
    conv = Convertor(db, shot, run, occur, new_db, user, backend)
    Summary(conv, ncfile=ncfile, iter_=iter_, _deb=_deb)
    # close datafile
    conv.close_imas_database()

if __name__ == '__main__':
    occur = 0
    if len(sys.argv) > 1:
        occur = int(sys.argv[1])
        print('summary occurrence', occur)
    _deb = True
    new_db = True
    plot = True
    mds_shot =  None #141711 #205430 # also 142000, 141711
    iter_ = False #ITER
    #mds_run = 2502 #Y02
    run = 0 #2302 #W02 #this is occur
    shot = 163303 #200001 #
    db = 'D3D' #NSTX-U
    user = os.getenv('USER') # or 'public'
    dbpath = '/public/imasdb/'
    ncfile = '163303W02.CDF'
    run_convertor(db, shot, run, occur, mds_shot, new_db, _deb=_deb, plot=plot, \
      ncfile=ncfile, iter_=iter_)

