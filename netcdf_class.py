
import os
import sys
from pathlib import Path
import numpy as np
import scipy
import re
from atoms_class import atom_base

from scipy.interpolate import interp1d

from pathlib import Path
from netCDF4 import Dataset

class netcdf_connection:
    def __init__(self, ncfile:str):
        self.nc = None
        self.vars_ = []
        self.ncfile = ncfile
        self.device = None
        self.shot = 0
        self.open_ncfile()
        self.set_machine_runid()

    def __del__(self):
        self.close_ncfile()

    def set_ncfile(self, ncfile):
        print(f're-setting cdf file to {ncfile}')
        self.ncfile = ncfile

    def open_ncfile(self):
        print(f'{"-"*30} \n netcdf file {self.ncfile}')
        path_ = Path(self.ncfile).absolute()
        if path_.is_symlink():
            path_ = os.readlink(self.ncfile)
        self.nc = Dataset(path_, 'r')
        self.expand_vars(self.nc)

    def close_ncfile(self):
        try:
            self.nc.close()
        except:
            print(f'Error in closing ncfile {self.ncfile}')

    def expand_vars(self, gp):
        for key, _ in gp.variables.items():
            self.vars_.append(key)
        for _, sub in gp.groups.items():
            expand_vars(self.vars_, sub)

    def reload_ncfile(self, fname:str):
        self.ncfile = fname
        try:
            self.nc.close()
            self.vars_.clear()
        except:
            pass
        self.open_ncfile()

    def make_re_pattern(self, name:str, nd1:int=0, nd2:int=0, nch1:int=0,
          nch2:int=0, nd:int=0, nch:int=0, dig_last:bool=False, step:int=-1) -> str:
        r'''
          nd1=digit_begin nd2=digit_end. same for char:nch1,nch2.
          nd, nch for fixed number; override nd1,nd2,nch1,nch2
          Used for specified ranges of # of digits and chars
          step=-1 :  add ^ + pattern + $;
          step!=-1 string including 1st '?' after replacing ? by pattern.
          eg: get_var_names_re('VV_PF?', cut='VV_', nd=1, nch1=1, nch2=2).

          #step call
          p = make_re_pattern('NAME\(?,?\)', nd1=1, nd2=3, step=1) => NAME\(\d{1,3}[A-z]{0,0}
          #second of the step call
          make_re_pattern(fr'{p},?\)', nd1=1, nd2=3)  #p is pat from 1st call= \d{1,3}[A-z]{0,0}
                        => NAME\(\d{1,3}[A-z]{0,0},\d{1,3}[A-z]{0,0}\)
          NOTE: escape ().  Error: NAME(?,?)
          For non-regex don't do escape
        '''
        ind = name.find('?') #only 1st position of '?''
        if ind < 0:
            raise Exception(f'Error in pattern build of {name} : missing "?"')

        nd1 = nd if nd>0 else nd1
        nd2 = nd if nd>0 else nd2
        nch1 = nch if nch>0 else nch1
        nch2 = nch if nch>0 else nch2
        pat = f'{name[:ind]}'
        if dig_last:
            pat =  pat + r'[A-z]{' + f'{str(nch1)},{str(nch2)}'+ r'}' + \
                 r'\d{' + f'{str(nd1)},{str(nd2)}' + r'}'
        else:
            pat =  pat + r'\d{' + f'{str(nd1)},{str(nd2)}' + r'}' + \
                r'[A-z]{' + f'{str(nch1)},{str(nch2)}' + r'}'

        if step == -1:
            pat = pat + rf'{name[ind+1:]}' + '$'

        return pat

    def make_re_pattern_quest(self, name:str) -> str:
        r'''
         Requires at least one '?' in name
         All '?' are replaced by \w+
         BSTH_?4_P? converted to BSTH_\w+4_P\w+
        '''
        inds = [ind for ind, ch in enumerate(name.lower()) if '?' == ch]
        if len(inds) < 1:
            raise Exception(f'Error in pattern build of {name} : missing "?"')
        pat = rf'{name[:inds[0]]}\w+'
        for i in range(1,len(inds)):
            pat += rf'{name[inds[i-1]+1:inds[i]]}\w+'
        return pat + name[inds[-1]+1:] + '$'

    def get_var_names_re_(self, pat, cut:str='', length:bool=False):
        matches = [self.vars_[i] for i,s in enumerate(self.vars_) if re.match(pat, s)]
        if cut:
            matches = [str(s.split(cut,1)[1]) for s in matches if len(s.split(cut,1))>1]
        return matches if not length else len(matches)

    def get_var_names_re(self, name:str, cut:str='', nd1:int=0, nd2:int=0, nch1:int=0,
          nch2:int=0, nd:int=0, nch:int=0, dig_last:bool=False, length:bool=False):
        #use number of digits and characters to make pattern
        pat = self.make_re_pattern(name, nd1, nd2, nch1, nch2, nd, nd, dig_last)
        return self.get_var_names_re_(pat, cut, length)

    def get_var_names(self, name:str, cut:str='', length=False):
        #'?' in names are replaced by \w+'
        pat = self.make_re_pattern_quest(name)
        return self.get_var_names_re_(pat, cut, length)

    def check_var_name(self, name:str, length=False):
        pat = rf'{name}$'
        n = self.get_var_names_re_(pat, length=length)
        return n>0

    #names: list[str] or np.array
    def split_var_name(self, names, cut:str='') -> dict[str:list]:
        try:
            res = {s.split(cut)[0]: s.split(cut)[1:] for s in names if len(s.split(cut))>1}
        except:
            print(f'{name} split error {name.split(chars,1)}')


    #TODO use chars after '?'
    def get_var_names_lst(self, name:str, lst:list[str], crop=False, cut:str='') -> list[str]:
        r'''
          select existing variables matching pattern + passed in lst
          eg: get_var_names_lst('BSTH_?', sp_look); sp_look=['H', 'D', 'T',..]
          Only 1st position of '?' used to split. Others are treated as characters:
          BSTH_?_P? is used as BSTH_... while _P? discarded
        '''
        ind = name.find('?') #only 1st position of '?' replaced
        if ind < 0:
            raise Exception(f'Error in search of array variables {name}')
        name = name[:ind]
        look = [name+ str(s) for s in lst]
        matches = [self.vars_[i] for i,s in enumerate(self.vars_) if s in look]
        st = cut if cut else name
        if cut or crop:
            matches = [str(s.split(st,1)[1]) for s in matches if len(s.split(st,1))>1]
        return matches

    def get_species_names(self, name:str, cut:str) -> list[str]:
        #characters after 1st position of '?' returned
        ind = name.find('?') #only 1st position of '?' replaced
        if ind < 0:
            raise Exception(f'Error in search of array variables {name}')
        name = name[:ind]
        sps = self.get_var_names_re_(name, cut)
        sps = ['HE' + s if (s=='3' or s=='4') else s for s in sps]
        return sps


    def get_ncvar(self, var:str, scale:float=1., ind:int=-1) -> np.array:
        try:
            res = None
            if ind>=0:
               res = np.array(self.nc.variables[var][:][ind])
            else:
               res = np.array(self.nc.variables[var][:])

            if scale < 1. or scale > 1.:
                return scale*res
            else:
                return res
        except:
            return np.array([])

    def get_ncvars_lst(self, lst:list[str], scale:float=1., ind:int=-1) -> np.array:
        tot = []
        n = 0
        for name in lst:
            res = np.array(self.get_ncvar(name, scale, ind))
            if n>0 and len(res) == len(tot):
                tot += res
            if n==0 and len(res)>0:
                tot = res
                n = 1
        return tot if n>0 else np.array([])

    def get_ncvar_sum(self, name:str, scale:float=1.):
        r'''
         Requires at least one '?' in name
         All '?' are replaced by \w+
         BSTH_?4_P? converted to BSTH_\w+4_P\w+
        '''
        names = self.get_var_names(name)
        return self.get_ncvars_lst(names, scale=scale)


    def get_nc_global_atts(self):
        runid, cdf_date, build_date, device, year, tglf, glf = None, \
                                  None, None, None, None, None, None
        try:
            runid = self.nc.getncattr('Runid').strip()
            cdf_date = self.nc.getncattr('CDF_date').strip()
            build_date = self.nc.getncattr('BUILD_date').strip()
            device = self.nc.getncattr('device').strip()
            year = self.nc.getncattr('year').strip()
        except:
            pass
        try:
            tglf = self.nc.getncattr('TGLF_hash').strip()
            glf = self.nc.getncattr('GLF23_hash').strip()
        except:
            pass
        return runid, device, year, cdf_date, build_date, tglf, glf


    def print_ncattributes(self, key):
        nc_attrs = self.nc.ncattrs()
        if key in nc_attrs:
            print(self.nc.getncattr(key))
        else:
            try:
                print("\ttype:", repr(self.nc.variables[key].dtype) )
                for ncattr in self.nc.variables[key].ncattrs():
                    print(f'\t{ncattr} :  {repr(self.nc.variables[key].getncattr(ncattr))}' )
            except KeyError:
                print(f"\tWARNING: variable attribute not found: {key}")


    def set_machine_runid(self):
        self.device = self.get_ncattr('device').strip()
        self.shot = int(self.get_ncattr('shot'))

        try:
            runf = Path(self.ncfile).stem #.CDF
            runid = self.get_ncattr('Runid').strip()
            if runf != runid:
                print(f'Runid in nc file header {runid}')
        except:
            print(f'Could not verify Runid with filename {self.ncfile}')


    def get_ncattr(self, key):
        name = key
        keys = key.split(':',1)
        if len(keys) > 1:
            key, name = keys[0], keys[1]

        #some global attributes are ints, others are strings
        if key in self.nc.ncattrs():
            return self.nc.getncattr(key)

        try:
            if name in self.nc.variables[key].ncattrs():
                return repr(self.nc.variables[key].getncattr(name))
        except KeyError:
            print(f"\tWARNING: attribute not found: {key}")

        return None

    def get_ncvarx(self, var:str, xvar:str=None, n:int=-1):
        yv = self.nc.variables[var][:]
        if(numpy.ma.isMaskedArray(yv)):
            yv = np.array(yv)
        if xvar is None:
            return yv

        ydims = self.nc.variables[var].dimensions
        xv = self.nc.variables[xvar][:]
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

                zz = self.nc.variables[zvar][:]
                xdims = self.nc.variables[xvar].dimensions
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


    def print_sample(self):
        for k,v in self.nc.dimensions.items():
            print(k,'\n',v)
            break
        for k,v in self.nc.variables.items():
            print(k,'\n',v)
            break
        print('\n',self.nc.variables['TIME'])
        print('\n', self.nc.variables['TIME'].dimensions)
        print('DIMENSIONS')
        for dim in self.nc.dimensions:
            print(dim)