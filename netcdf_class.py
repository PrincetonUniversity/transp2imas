
import os
import sys
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
        self.open_ncfile()

    def set_ncfile(self, ncfile):
        print('re-setting cdf file')
        self.ncfile = ncfile

    def open_ncfile(self):
        print(f'netcdf file {self.ncfile}')
        path_ = os.readlink(self.ncfile)
        self.nc = Dataset(path_, 'r')
        self.expand_vars(self.nc)

    def close_ncfile(self):
        try:
            self.nc.close()
        except:
            print('Error in closing ncfile')

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

    #default pattern has 1 or 2 digits \d{1,2}
    #name's last character is discarded. Eg. PECIN?
    def get_nvar_names(self, name:str) -> int:
        if name[-1] != '?':
            raise Exception(f'Error in search of array variables {name}')
        name = name[:-1]
        pat = r'\b' + name + r'\d{1,2}\b'
        matches = [self.vars_[i] for i,s in enumerate(self.vars_) if re.search(pat, s)]
        return len(matches)


    def get_var_names_re(self, name:str, cut:str='', nd1:int=0, nd2:int=0, nch1:int=0,
          nch2:int=0, nd:int=0, nch:int=0) -> list:
        '''
          nd1=digit_begin nd2=digit_end. same for char:nch1,nch2.
          nd, nch for fixed number; override nd1,nd2,nch1,nch2
          Used for specified ranges of # of digits and chars
          eg: get_var_names_re('VV_PF?', cut='VV_', nd=1, nch1=1, nch2=2)
        '''
        if name[-1] != '?':
            raise Exception(f'Error in search of array variables {name}')
        name = name[:-1]

        nd1 = nd if nd>0 else nd1
        nd2 = nd if nd>0 else nd2
        nch1 = nch if nch>0 else nch1
        nch2 = nch if nch>0 else nch2
        pat = f'^{name}' + r'\d{' + f'{str(nd1)},{str(nd2)}' + r'}' + \
              r'[A-z]{' + f'{str(nch1)},{str(nch2)}' + r'}$'

        matches = [self.vars_[i] for i,s in enumerate(self.vars_) if re.search(pat, s)]
        if cut:
            matches = [s.split(cut,1)[1] for s in matches if len(s.split(cut,1))>1]
        return matches



    def get_var_names_lst(self, name:str, lst:list, crop=True, cut:str='') -> list:
        '''
          select existing variables matching pattern + passed in lst
          eg: get_var_names_lst('BSTH_?', sp_look); sp_look=['H', 'D', 'T',..]
        '''
        if name[-1] != '?':
            raise Exception(f'Error in search of array variables {name}')
        name = name[:-1]
        look = [name+ str(s) for s in lst]
        matches = [self.vars_[i] for i,s in enumerate(self.vars_) if s in look]
        st = cut if cut else name
        if cut or crop:
            matches = [str(s.split(st,1)[1]) for s in matches if len(s.split(st,1))>1]
        return matches


    def get_ncvar(self, var:str, scale:float=0., ind:int=-1) -> np.array:
        try:
            res = np.array(self.nc.variables[var][:])
            if ind>=0:
                res = res[ind]
            if scale > 0:
                return scale*res
            else:
                return res
        except:
            return np.array([])

    def get_ncattr(self, name):
        try:
            return self.nc.getncattr(name)
        except:
            print(f'{name} not found in {self.ncfile}')
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