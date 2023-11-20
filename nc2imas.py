#!/usr/bin/python3
#
import os
import sys
import math
import numpy as np
import inspect

from scipy import constants
from scipy.interpolate import interp1d
from datetime import datetime
import re
import logging
from netcdf_class import netcdf_connection
from imas_class import imas_connection
from atoms_class import *
from utilities import *


class LineNumber:
    def __str__(self):
        return str(inspect.currentframe().f_back.f_lineno)

__LINE__ = LineNumber()


def getVarName(var):
    localVars = inspect.currentframe().f_back.f_locals.items()
    return str([k for k, v in localVars if v is var][0])


class nc2imas(imas_connection, netcdf_connection):
    def __init__(self, run, backend, occur, user, machine, new_db, ncfile,
         nml_file:str=None, log:str='error', debug:int=0):
        netcdf_connection.__init__(self, ncfile, machine)
        imas_connection.__init__(self, self.device, self.shot, run, backend, occur, user, new_db)
        self.nml_name = nml_file
        self.logger = None
        #namelist contents dict
        self.nmlist = None
        self.set_logger(log)
        self.logger.info(f'creating nc2imas using {ncfile} for {self.device} {self.shot} {self.run}')
        self.load_namelist()
        self.cm2m = 0.01
        self.area_inv = 1.0e4
        self.vconv = 1.e-6
        self.aconv = 1.0e-4
        self.vol_inv = 1.e6
        self.timesc = []
        self.timep = []
        self.ntimes = 0
        self.ntimep = 0
        self.sgn = -1 if self.iter_ else 1
        self.pi = constants.pi  #scipy

        self.debug = False if debug<1 else True
        self.debug_max = False if debug<2 else True
        self.get_nc_times()


    def get_nc_times(self):
        self.timesc = self.get_ncvar('TIME')
        self.ntimes = len(self.timesc)
        self.timep = self.get_ncvar('TIME3')
        self.ntimep = len(self.timep)
        if self.ntimes>1 and self.ntimep>1:
            coords_x = self.get_ncvar('X')
            sh = coords_x.shape[1] if coords_x.ndim>1 else coords_x.shape
            self.logger.info(f' nprof {sh}: ntimes {self.ntimes}'
               f' {self.timesc[0]} .. {self.timesc[-1]}')
            self.logger.info(f' ntimep {self.ntimep} {self.timep[0]} .. {self.timep[-1]}')
            self.logger.info(f' Dtimep {self.timep[1]-self.timep[0]};'
               f' Dtimesc {self.timesc[1]-self.timesc[0]}')
        else:
            raise Exception(f'Invalid times prof:{self.ntimep} scalar {self.ntimes}')


    def check_1dim(self, dat):
         if isinstance(dat, np.ndarray):
            return dat.ndim==1
         elif isinstance(dat, (list, tuple, set)):
            return np.array(dat).ndim==1


    #Don't call this using varargs forwarded by another function
    def is_valid(self, *variables, ind:int=0, name:str=''):
        ii = ind+1
        for var in variables:
          if var is None:
              if ind==1:
                  self.logger.info(f'{name} is None')
              return False
          elif isinstance(var, np.ndarray):
              #ndim=0; eg:np.array(5)
              if var.ndim==0:
                  return False
              sh = var.shape
              if sh[0]<ii or (len(sh)>1 and sh[1]==0):
                  if ind==1:
                      self.logger.info(f'{name} shape {sh}')
                  return False
          elif isinstance(var, (list, tuple, set)):
              if len(var)<ii or not any(var):
                  if ind==1:
                      self.logger.info(f'{name} len: {len(var)}')
                  return False
        return True

    def sum_arrays(self, *arrays):
        res = []
        init = False
        try:
            for arr in arrays:
                if len(arr)==0:
                    continue
                if init:
                    res = res+arr
                else:
                    res = arr
                    init = True
        except:
            self.logger.info(f' sum arrays error')
        return res


    def operate_arrays(self, arr1, op:str, arr2=None, name1:str='', name2:str=''):
        if not self.is_valid(arr1, arr2):
            self.logger.info(f' operate_arrays error: {name1} {op} {name2}')
            return np.array([])

        err = False
        try:
            if op == '-':
                return arr1 - arr2
            elif op == '+':
                return arr1 - arr2
            elif op == '/':
                return np.divide(arr1, arr2, out=np.zeros_like(arr1), where=arr2!=0)
            elif op == '*':
                return arr1 * arr2
            elif op == 'sqrt':
                if isinstance(arr1, np.ndarray):
                    return np.sqrt(arr1)
                elif isinstance(arr1, list):
                    return [math.sqrt(x) for x in arr1]
                elif isinstance(arr1, set):
                    return {math.sqrt(x) for x in arr1}
                elif isinstance(arr1, tuple):
                    return tuple(math.sqrt(x) for x in arr1)
            else:
                err = True
        except:
            err = True

        if err:
            self.logger.info(f' operate_arrays error:: {name1} {op} {name2}')
            return np.array([])

    def np_multiply_sum(self, arr1, arr2, axis:int, name1:str='', name2:str='', scale=None):
        if not self.is_valid(arr1, arr2):
            mes = ' np_multiply_sum error'
            s1, s2 = len(arr1), len(arr2)
            if isinstance(arr1, np.ndarray) and isinstance(arr2, np.ndarray):
               s1, s2 = arr1.shape, arr2.shape
            self.logger.info(f'{mes}: {name1}:{s1} {name2}: {s2}')
            return np.array([])
        try:
            if scale is not None:
                return scale * np.sum(np.multiply(arr1, arr2), axis=axis)
            else:
                return np.sum(np.multiply(arr1, arr2), axis=axis)
        except:
            self.logger.info(f' np_multiply_sum error:: {name1} {name2}')
            return np.array([])


    def checkErrExit(self, *variables, name:str=''):
        for var in variables:
            stat = self.is_valid(var)
            if stat is False:
               self.logger.exception(f'var {name} is invalid')
               raise Exception(f'**ERROR**: {name} is invalid: len={len(var)}: {var}')


    def interp1d_time(self, time_new, time, dat, ind1:int=None, ind2:int=None):
        r'''
        interpolate dat from time onto time_new
        time_new: 1d times; output dat has this size.
        dat: 2d (time, 1d)
        time: 1d times (x axis); Same size of 1st dim or values or f(x) of dat
        '''
        #if dat is 2D, 1st dim should have len of len(time)
        if self.is_valid(time_new) and self.is_valid(dat, time, ind=len(time)-1):
           if isinstance(dat, np.ndarray) and dat.ndim==2:
               if ind1 is not None and ind2 is None:
                  dat = dat[ind1, :]
               if ind1 is None and ind2 is not None:
                  dat = dat[:, ind2]
           if len(time_new) == len(time):
               return dat
           else:
               return np.interp(time_new, time, dat)
        else:
           return []


    #fill ids.value
    def fill_value(self, ids, dat, rsz:int=-1, mult:float=None, name:str=None):
        name_ = name if name else 'data'
        if not self.is_valid(dat):
            #self.logger.info(f'{name} is None')
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
            pass
            #self.logger.info(f'{var} not found')
        else:
            self.fill_value(ids, dat, rsz=rsz, mult=mult)

    # if var is 2D, ind must be valid. If ind<0, use this to fill ONLY 1D data
    def fill_data_nc(self, ids, var:str, mult:float=None, ind:int=-1, resize:bool=True,
         interp:bool=False):
        if not isinstance(var, str):
            raise ValueError(f'variable {var} is not string')
        dat = self.get_ncvar(var, ind=ind)
        if not self.is_valid(dat):
            if ind==0:
                self.logger.info(f'{var} not-found')
        else:
            if mult is not None:
                dat = mult*dat
            if interp is True and self.ntimes != self.ntimep:
                dat = np.interp(self.timep, self.timesc, dat)
            if resize:
                ids.resize(len(dat), refcheck=False)
            ids[:] = dat


    #If ind<0, use this to fill ONLY 1D data. if var is 2D, ind must be valid.
    #Don't use this to fill scalar
    def fill_data(self, ids, dat, mult:float=None, ind:int=-1, ind2:int=-1,
          interp:bool=False, name:str=None, resize:bool=True, chk:int=-1) -> bool:
        nm = name if name is not None else 'data'

        try:
          if self.is_valid(dat, name=nm):
            if (ind>= 0 and ind2>=0):
               dat = dat[ind:ind2]
            elif ind>= 0:
               dat = dat[ind]

            if chk>0 and not np.isscalar(dat) and len(dat)!=chk:
               if ind==0:
                  self.logger.info(f'Error: size mismatch {len(dat)} != {chk}')
               return

            if np.isscalar(dat):
               raise Exception('Cannot fill scalar value')
            if mult is not None:
                dat = mult*dat
            if interp is True and self.ntimes != self.ntimep:
                dat = np.interp(self.timep, self.timesc, dat)
            if (not np.isscalar(dat)) and resize:
                ids.resize(len(dat), refcheck=False)
            ids[:] = dat
          elif ind==0:
            self.logger.info(f'{nm} is invalid')
        except:
            if ind==0:
                self.logger.info(f'{nm} is invalid')

        #if not self.check_1dim(dat):
            #print('Check dim')
            #raise Exception('cannot fill data that is not 1D')

    def load_namelist(self):
        if self.nml_name is None:
            print('**Error: namelist file not set')
            return
        self.nmlist = {}
        if not check_file_exists(self.nml_name):
            print(f'**WARNING: namelist file {self.nml_name} not found')
            return

        with open(self.nml_name, encoding = 'ISO-8859-1') as nmf:
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


    #TODO move to utilities
    def parse_number(self, ss:str, tp):
        if not isinstance(ss, str):
            return tp(ss)

        #replace exponential 10** if any
        ss0 = ss
        res= ss.find('10**')
        if res >=0:
            ss = ss.replace('10**', 'E')
        opls = ['*', '/', '+', '-']
        nn = []
        ops = []
        p = 0

        #splitting into separate number strings
        for i,a in enumerate(ss):
            if a in opls and i>0 and ss[i-1].upper() !='E':
                nn.append(ss[p:i])
                ops.append(a)
                p = i+1

        if p<len(ss):
            nn.append(ss[p:])
        if len(ops)!=len(nn)-1:
            raise Exception(f'Error in parsing string to number {ss0}')
        n = 0
        try:
            n = tp(nn[0]) if len(nn)>0 else None
        except:
            print(f'parse error: {ss0} extracted: {nn}')
        for i,s in enumerate(ops):
            d = tp(nn[i+1])
            if s=='*':
                n = n * d
            elif s=='/':
                n = n/d
            elif s=='+':
                n = n+d
            elif s=='-':
                n = n-d
        return n


    #name's last character is discarded. Eg. XECECH?. () added separately
    def get_namelist_nvars(self, name:str, ndig:int=2):
        if name[-1] != '?':
            raise Exception(f'Error: {name} needs "?" as last char')
        name = name[:-1]
        pattern = fr'\b{name}(\((d{0,{ndig}}))?\b'
        try:
            matches = [self.nmlist[k] for k in self.nmlist.keys()
                        if name==k or re.search(pattern, k)]
            return len(matches)
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
        # name: regex need to be escaped \(
        pat = self.make_namelist_re_pattern(name, nd1, nd2, nch1, nch2, nd, nch)

        if pat is None:
            return []
        matches = [k for k in self.nmlist.keys() if re.search(pat, k)]
        if cut:
            matches = [s.split(cut,1)[1] for s in matches if len(s.split(cut,1))>1]
        return matches


    #returns scalar or list
    def get_namelist_var_(self, name:str, tp=None, scale=1., ind:int=-1):
        try:
            val = self.nmlist[name]
            if isinstance(val, str):
               val = re.split(r'[\s,;]+', val)
            if ind>=0:
               val = val[ind]
        except:
            self.logger.info(f'** {name} NOT found in name list')
            if self.nmlist is None:
                self.logger.info(f'** {self.nmlist} is not set')
            return None

        narr = 0
        try:
            narr = len(val)
        except:
            pass

        if tp in [int, float]:
            if narr>0:
               val =  [self.parse_number(val[i], tp) for i in range(narr)]
            else:
               val = self.parse_number(val, tp)
            if scale != 1.0:
                val = scale*val
        return val

    # Some scalar can stand for arrays, e.g: FREQECH, which causes error when
    #  an array is expected. So, call for array should handle single variable also.
    # Returns scalar for single, list for array on RHS
    def get_namelist_vars(self, name:str, tp=None, scale=1.0, ind:int=-1,
          size:int=1, resize:bool=False):
        #name: don't escape (); only regex needs \(
        # search for single or without () for array on rhs
        res = self.get_namelist_var_(name, tp, scale)
        n = 0
        # search for index
        if ind>=0 and res is None:
           namei = f'{name}({ind+1})'
           res = self.get_namelist_var_(namei, tp, scale)
           n = 1
        # search for array
        elif res is None:
           n = self.get_namelist_nvars(f'{name}?')
           res = []
           for i in range(n):
              namei = f'{name}({i+1})'
              val = self.get_namelist_var_(namei, tp, scale)
              if val is not None:
                 res.append(val)

        check = True
        try:
            if resize and n!= size:
              n = size
              res = np.full(n, res)
            if n>0 and np.size(res)!=n:
               check = False
        except:
            pass
        if not check:
            raise Exception(f'size of {name} != {n}')
        return res

    def get_namelist_var(self, name:str, tp=None, scale=1.0, ind:int=-1, size:int=1):
        return get_namelist_vars(name, tp, scale, ind, size)


    def extract_str_between(self, name:str, pat1:str, pat2:str, default:str, tp=None):
        res = None
        try:
            li = name.find(pat1) + len(pat1) + 1
            #can have multiple occurences. eg: ','
            ris = [s.start() for s in re.finditer(pat2, name)]
            ri = [i for i in ris if i>li][0]
            res = name[li:ri] if ri>li+1 else default
            if tp is not None:
                res = tp(res)
        except:
            self.logger.debug(f'extract error: {name} doesnot have {pat1} and {pat2}')
        return res

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
        logname = 'log.out'
        if os.path.exists(logname):
            os.rename(logname, f'{logname}_')
        self.logger.addHandler(logging.FileHandler(logname, mode='w')) #'a'

    #for all times
    def fill_vac_tor_field(self, ids):
        #RAXIS(TIME)
        raxis = self.get_ncvar('RAXIS', self.cm2m)
        #(TIME) TESLA*CM. vac field
        bzxR = self.get_ncvar('BZXR', self.cm2m)
        if self.is_valid(raxis, bzxR):
            if self.ntimep != self.ntimes:
               raxis = np.interp(self.timep, self.timesc, raxis)
               bzxR = np.interp(self.timep, self.timesc, bzxR)
            b_field_tor = bzxR / raxis
            rave = (raxis[0] + raxis[-1])/2.0
            self.logger.info('Vac_tor_field')
            ids.vacuum_toroidal_field.r0 = rave.astype(float)
            self.fill_data(ids.vacuum_toroidal_field.b0, b_field_tor)

    #fill ids.profiles_1d[*].grid. Resize if not laready done
    def fill_profiles1d_grid(self, ids):
        #scalar time
        timesp = self.get_ncvar('TIME3')
        ntimep = len(timesp)
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
            bzf = self.pi*np.tile(bz[None,:].T, bflux_tor.shape[1])
            rho_tor = self.operate_arrays(bflux_tor, '/', bzf)
            rho_tor = self.operate_arrays(rho_tor, 'sqrt')

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
            ids.profiles_1d.resize(ntimep)

        for i in range(ntimep):
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
        ids.label = atom_info.label
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
                 datv = self.np_multiply_sum(dat, gv, axis=1)
                 dat = self.operate_arrays(datv, '/', 1.0e-6*pv)

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


if __name__ == '__main__':
    print("Please use: python convert.py")

