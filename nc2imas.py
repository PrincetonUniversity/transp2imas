#!/usr/bin/python3
#
import os
import sys
import numpy as np
import inspect

from scipy import constants
from scipy.interpolate import interp1d
from datetime import datetime
import re
from pathlib import Path
import logging
import argparse
from netcdf_class import netcdf_connection
from imas_class import imas_connection
from atoms_class import *


def getVarName(var):
    localVars = inspect.currentframe().f_back.f_locals.items()
    return str([k for k, v in localVars if v is var][0])


class nc2imas(imas_connection, netcdf_connection):
    def __init__(self, run, backend, occur, user, new_db, ncfile,
         nml_file:str=None, log:str='error'):
        netcdf_connection.__init__(self, ncfile)
        imas_connection.__init__(self, self.device, self.shot, run, backend, occur, user, new_db)

        self.nml_name = nml_file
        self.logger = None
        #namelist contents dict
        self.nmlist = None
        self.set_logger(log)
        self.logger.info(f'creating nc2imas using {ncfile} for {self.device} {self.shot}')
        self.load_namelist()
        self.cm2m = 0.01
        self.area_inv = 1.0e4
        self.vconv = 1.e-6
        self.aconv = 1.0e-4
        self.vol_inv = 1.e6
        self.sgn = -1 if self.iter_ else 1
        self.pi = constants.pi  #scipy

    #scalar 0 is invalid. Don't check if it can be valid
    #Don't call this using varargs forwarded by another function
    def is_valid(self, *variables, ind:int=0, name:str='var'):
        ii = ind+1
        for var in variables:
          if var is None:
              self.logger.info(f'{name} is None')
              return False
          elif isinstance(var, np.ndarray):
              sh = var.shape
              if sh[0]<ii or (len(sh)>1 and sh[1]==0):
                  self.logger.info(f'{name} shpe {sh}')
                  return False
          elif isinstance(var, (list, tuple, set)):
              if len(var)<ii:
                  self.logger.info(f'{name} len: {len(var)}')
                  return False
          else:
              if var == 0:
                  self.logger.info(f'{name} = 0')
                  return False
        return True

    def operate(self, var1, symb:str, var2, name1:str='', name2:str=''):
        if not self.is_valid(var1, var2):
            self.logger.info(f' operate error: {name1} {symb} {name2}')
            return

        if symb == '-':
            return var1 - var2
        elif symb == '+':
            return var1 - var2
        elif symb == '/':
            return var1/var2
        elif symb == '*':
            return var1 * var2
        else:
            self.logger.info(f' operate error {symb}')
            return None


    def checkErrExit(self, *variables):
        for var in variables:
            name = getVarName(var)
            stat = self.is_valid(var)
            if stat:
                self.logger.exception(f'{name} is empty')


    #fill ids.value
    def fill_value(self, ids, dat, rsz:int=-1, mult:float=None, name:str=None):
        name_ = name if name else 'data'
        if not self.is_valid(dat):
            self.logger.info(f'{name} is None')
            return

        nresize = len(dat) if (not np.isscalar(dat) and rsz==-1) else rsz

        if nresize > 0:
            ids.value.resize(nresize, refcheck=False)
        if mult is not None:
            dat = mult*dat
        if np.isscalar(dat):
            ids.value = dat
        else:
            ids.value[:] = dat
        #if self.src:
        #    ids.source = self.src

    def fill_value_nc(self, ids, var:str, rsz:int=-1, mult:float=None):
        dat = self.get_ncvar(var)
        if dat is None:
            self.logger.info(f'{var} not found')
        else:
            self.fill_value(ids, dat, rsz=rsz, mult=mult)

    def fill_data_nc(self, ids, var:str, mult:float=None, ind:int=-1, resize:bool=True):
        if not isinstance(var, str):
            raise ValueError(f'variable {var} is not string')
        dat = self.get_ncvar(var)
        if not self.is_valid(dat):
            self.logger.info(f'{var} not-found')
        else:
            dat = dat if ind<0 else dat[ind]
            if resize:
                ids.resize(len(dat), refcheck=False)
            if mult is not None:
                dat = mult*dat
            ids[:] = dat

    #don't use this to fill scalar
    def fill_data(self, ids, dat, mult:float=None, ind:int=-1, name:str=None,
          resize:bool=True, chk:int=-1) -> bool:
        nm = name if name is not None else 'data'
        try:
          if self.is_valid(dat, name=nm):
            dat = dat[ind] if ind>= 0 else dat
            if chk>0 and not np.isscalar(dat) and len(dat)!=chk:
                self.logger.info(f'Error: size mismatch {len(dat)} != {chk}')
                return

            if (not np.isscalar(dat)) and resize:
                ids.resize(len(dat), refcheck=False)
            if mult is not None:
                dat = mult*dat
            ids[:] = dat
          else:
            self.logger.info(f'{nm} is invalid')
        except:
            self.logger.info(f'{nm} is invalid')

    def load_namelist(self):
        if self.nml_name is None:
            print('**Error: namelist file not set')
            return
        self.nmlist = {}
        with open(self.nml_name) as nmf:
            self.logger.info(f'loading {self.nml_name}')
            for line in nmf:
                lns = line.strip()
                if lns and (lns[0] not in ['!', '&', '/']):
                    lns = lns.split('!',1)[0]
                    kv = lns.split('=',1)
                    if len(kv) == 2:
                        val = kv[1].strip()
                        #cast to int or float
                        try:
                            val = int(val)
                        except ValueError:
                            try:
                                vv = val.upper().replace('D0','')
                                val = float(vv)
                            except ValueError:
                                pass
                        self.nmlist[kv[0].strip()] = val

    def reload_namelist(self, fname:str):
        self.nml_name = fname
        try:
            self.nmlist.clear()
        except:
            pass
        self.load_namelist()

    #name's last character is discarded. Eg. XECECH?
    def get_namelist_nvars(self, name:str):
        if name[-1] != '?':
            raise Exception(f'Error: {name} needs "?" as last char')
        name = name[:-1]
        pattern = fr'\b{name}(\d{1,2})\b'
        try:
            matches = [self.nmlist[k] for k in self.nmlist.keys() if re.search(pattern, k)]
            len(matches)
        except:
            self.logger.info(f'**Error: {name}: namelist not found')
        return 0


    def make_namelist_re_pattern(self, name:str, nd1:int=0, nd2:int=0, nch1:int=0,
          nch2:int=0, nd:int=0, nch:int=0) -> str:
        #only upto ?,? handled (two occurences)
        nn = name.count('?')
        st = 1 if nn>1 else -1
        pat = self.make_re_pattern(name, nd1, nd2, nch1, nch2, nd, nch, step=st)
        ii = name.find('?')

        if nn > 1:
            pat = self.make_re_pattern(fr'{pat}{name[ii+1:]}', nd1, nd2, nch1, nch2, nd, nch)
        return pat

    def get_namelist_keys_re(self, name:str, cut:str='', nd1:int=0, nd2:int=0, nch1:int=0,
          nch2:int=0, nd:int=0, nch:int=0) -> list:
        # name: regex need escaped \(
        pat = self.make_namelist_re_pattern(name, nd1, nd2, nch1, nch2, nd, nch)

        if pat is None:
            return []
        matches = [k for k in self.nmlist.keys() if re.search(pat, k)]
        if cut:
            matches = [s.split(cut,1)[1] for s in matches if len(s.split(cut,1))>1]
        return matches


    def get_namelist_vars(self, name:str, tp=None, scale=0):
        #return list. name's last character is discarded. Eg. XECECH?
        n = self.get_namelist_nvars(name)
        name0 = name[-1]
        try:
            val = [self.nmlist[f'{name0}({i})'] for i in range(1,n+1)]
        except:
            self.logger.info(f'** Failed to get {name} from name list')
            if self.nmlist is None:
                self.logger.info(f'** {self.nmlist} is not set')
        else:
            if self.is_valid(val):
                if tp in [int, float]:
                    val = [tp(x) for x in val]
                if scale > 0:
                    val = scale*val
                return val
        return []

    #returns scalar
    def get_namelist_var_(self, name:str, tp=None, scale=0, ind:int=-1):
        try:
            val = self.nmlist[name]
            if ind>=0 and isinstance(val, str):
                val = re.split(r'[\s,;]+', val)[ind]
        except:
            self.logger.info(f'** {name} NOT found in name list')
            if self.nmlist is None:
                self.logger.info(f'** {self.nmlist} is not set')
            return None

        if tp in [int, float]:
            val = tp(val)
            if scale > 0:
                val = scale*val
        return val


    def get_namelist_var(self, name:str, tp=None, scale=0, ind:int=-1):
        #name: don't escape (); only regex needs \(
        res = self.get_namelist_var_(name, tp, scale, ind)
        if res is None and ind>=0:
            name = f'{name}({ind+1})'
            res = self.get_namelist_var_(name, tp, scale, ind)
        return res


    def extract_str_between(self, name:str, pat1:str, pat2:str, default:str):
        str_ = None
        try:
            li = name.find(pat1)+1
            ri = name.find(pat2)
            str_ = name[li:ri] if ri>li+1 else default
        except:
            self.logger.debug(f'extract error: {name} doesnot have {pat1} and {pat2}')
        return str_

    def fluxBdryToCenter(self, data, dat0=None):
        dat = []
        if not self.valid(data):
            return dat

        dat = dat.copy()
        #average of consecutive columns, except first
        dat[:,1:] += dat[:,:-1] #add previous column
        if self.valid(dat0):
            dat[:,0] += dat0
        else:
            dat[:,0] *= 2.
        return dat/2.

    def set_logger(self, log:str = None):
        self.logger = logging.getLogger(__name__)
        levels = {'debug': logging.DEBUG, 'info': logging.INFO,
                  'error': logging.ERROR}
        log = log if log is not None else 'error'
        level = levels.get(log.lower())
        self.logger.setLevel(level=level)
        self.logger.addHandler(logging.FileHandler('log.out', mode='w'))

    #for all times
    def fill_vac_tor_field(self, ids):
        raxis = self.get_ncvar('RAXIS', self.cm2m)
        #(TIME) TESLA*CM. vac field
        bzxR = self.get_ncvar('BZXR', self.cm2m)
        if self.is_valid(raxis, bzxR):
            b_field_tor = bzxR / raxis
            rave = (raxis[0] + raxis[-1])/2.0
            self.logger.info('Vac_tor_field')
            ids.vacuum_toroidal_field.r0 = rave.astype(float)
            self.fill_data(ids.vacuum_toroidal_field.b0, b_field_tor)

    #fill ids.profiles_1d[*].grid. Resize if not laready done
    def fill_profiles1d_grid(self, ids):
        #scalar time
        timesc = self.get_ncvar('TIME')
        ntimes = len(timesc)
        coord_x = self.get_ncvar('X')
        #profiles_1d
        # X(TIME3,X) is sqrt(norm. toroidal flux)
        rho_tor_norm = np.multiply(coord_x, coord_x)
        #beta_tor_norm
        bz = self.get_ncvar('BZ') #(TIME) T. BZ outside plasma
        #TRFLX(TIME3, XB)
        bflux_tor = self.get_ncvar('TRFLX')
        rho_tor = []
        if self.is_valid(bz, bflux_tor):
            rho_tor = np.sqrt(bflux_tor/(self.pi*np.tile(bz[None,:].T, bflux_tor.shape[1])))

        psi = self.get_ncvar('PLFLX', 2.*self.pi) #(TIME3, XB)
        rho_pol_norm = np.sqrt(abs((psi - psi[0])/(psi[-1]-psi[0]))) #TODO

        #DVOL(TIME3, X) CM**3  ZONE VOLUME
        gvol = self.get_ncvar('DVOL', self.vconv)
        # grid area. DAREA (TIME3, X) CM**2 ZONE CROSS SECTIONAL AREA
        grid_area = self.get_ncvar('DAREA', self.aconv)
        #Surface area of the toroidal flux surface
        #FLUX SURFACE AREA (TIME3,X) cm^2
        surface = self.get_ncvar('SURF', 1.0e-4)

        psi_axis = self.get_ncvar('PSI0_TR') # Wb/rad

        if len(ids.profiles_1d)==0:
            ids.profiles_1d.resize(ntimes)

        for i in range(ntimes):
            #grid
            gd = ids.profiles_1d[i].grid
            self.fill_data(gd.rho_tor_norm, rho_tor_norm, ind=i)
            self.fill_data(gd.rho_tor, rho_tor, ind=i)
            self.fill_data(gd.rho_pol_norm, rho_pol_norm, ind=i)
            self.fill_data(gd.psi, psi, ind=i)
            self.fill_data(gd.volume, gvol, ind=i)
            self.fill_data(gd.area, grid_area, ind=i)
            self.fill_data(gd.surface, surface, ind=i)
            if len(psi_axis) > i:
                gd.psi_magnetic_axis = psi_axis[i].astype(float)
            if len(psi) > i:
                gd.psi_boundary = psi[i,-1].astype(float)

    def fill_element_data(self, ids, sp:str, ilab:str='', zn:float=None):
        #sp is TRANSP symbol
        if not sp:
            return
        #sp is TRANSP symbol defined in atoms_class.entries
        atom_info = atom_base(sp, ilab)
        ids.label = sp
        ids.element.resize(1)
        ids.element[0].a = atom_info.a
        if zn:
            ids.element[0].z_n = zn
        else:
            ids.element[0].z_n = atom_info.z_n
        ids.element[0].atoms_n = 1

    #returns n species. fills volume average if plasma and grid volumes are passed in
    def fill_ion_densities(self, ids, ind1:int=-1, ind2:int=-1, pv=None, gv=None):
        #N/CM**3
        mult = 1.0e6
        #ions
        #calculate total ion densities
        nsp = 0
        idss = {'NH':ids.hydrogen, 'ND':ids.deuterium, 'NT':ids.tritium, 'NHe3':ids.helium_3,
               'NHe4':ids.helium_4, 'NBe':ids.beryllium, 'NLith':ids.lithium, 'NC':ids.carbon,
               'NNi':ids.nitrogen, 'NNe':ids.neon, 'NAr':ids.argon, 'NXe':ids.xenon,
               'NO2':ids.oxygen, 'NW':ids.tungsten, 'NFe':ids.iron, 'NKr':ids.krypton}

        for name, ids in idss.items():
            dat = self.get_ncvar(name)
            if not self.is_valid(dat):
                self.logger.info(f'{name} not found')
                continue
            if ind1 >= 0:
                dat = dat[ind1,:]
            elif ind2 >= 0:
                dat = dat[:, ind2]
            elif pv is not None:
                 #volume average if plasma volume is passed in, and NO index passed
                 dat = np.divide(np.sum(np.multiply(dat, gv), axis=1), 1.0e-6*pv)

            ids.value.resize(len(dat), refcheck=False)
            ids.value[:] = mult*dat
            nsp += 1
        return nsp


    def fill_properties(self, ids):
        runid, device, year, cdf_date, build_date, tglf, glf = self.get_nc_global_atts()
        comment =  f'TRANSP results for run {runid} {device} year: {year} CDF date: {cdf_date}'
        ids.ids_properties.comment = comment + f' build: {build_date} tglf: {tglf} glf23: {glf}'
        ids.ids_properties.creation_date = f'{str(datetime.now())} CDF: {cdf_date}'
        #ids.ids_properties.version_put = str(os.getenv('IMAS_VERSION'))
        ids.ids_properties.data_dictionary = str(os.getenv('IMAS_VERSION'))
        ids.ids_properties.access_layer = str(os.getenv('IMAS_VERSION'))
        #ids_properties/provenance.node.resize()
        return ids

    def fill_code_properties(self, ids):
        ids.code.name = 'transp2imas'
        ids.code.repository = 'https://github.com/PrincetonUniversity/transp2imas'
        return ids


def get_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-cdf', type=str, help='CDF TRANSP output file, with or without .CDF')
    parser.add_argument('-b', '--backend', type=str, help='backend: MDS or HDF')
    parser.add_argument('-r', '--run', type=int, help='run number of o/p IDS')
    parser.add_argument('-new', '--new_db', help='create new database', action="store_true")
    parser.add_argument('-o', '--occur', type=int, help='occurrence')
    parser.add_argument('-u', '--user', type=str, help='username or path of database')

    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument("-l", "--log", help="Provide logging level=debug,info,warn,error,critical")

    args = parser.parse_args()

    if not args.input_cdf:
        print("** Error: Failed to get command line arguments**\n")
        parser.print_help()
        print("Example: \npython convert.py -i 147634G14 [ -r 1 -b MDS -new -o 0 -l debug]\n")
        sys.exit(1)

    runid = args.input_cdf
    runid = runid[:-4] if runid[-4:].upper()=='.CDF' else runid

    ii = 0
    for i in range(len(runid)):
      if runid[i].isdigit():
        ii = i+1
      else:
        break
    shot, rid = int(runid[:ii]), runid[ii:]
    ncfile = str(shot)+rid+'.CDF'
    nml = str(shot)+rid+'TR.DAT'

    run = 1
    backend='MDS'
    newdb=True
    occur=0
    log='debug'  #for logger

    if args.run:
        run = int(args.run)
    if args.backend:
        backend = args.backend
    if args.new_db:
        new_db = args.new_db
    if args.occur:
        occur = int(args.occur)
    if args.log:
        log = args.log

    user = os.path.expanduser('~') +'/public/imasdb' if backend=='MDS' else os.getenv('USER')
    #new_db = True if occur == 0 else False

    return backend, run, occur, user, newdb, ncfile, nml, log


if __name__ == '__main__':
    print("Please use: python convert.py")

