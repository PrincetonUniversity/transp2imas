
import os
import sys
import numpy as np
import scipy
import re

from scipy.interpolate import interp1d

from pathlib import Path
from netCDF4 import Dataset

class netcdf_connection:
    def __init__(self, ncfile:str=None):
        self.nc = None
        self.vars_ = []
        self.ncfile = ncfile
        self.open_ncfile()

    def set_ncfile(self, ncfile):
        print('re-setting cdf file')
        self.ncfile = ncfile

    def open_ncfile(self):
        self.nc = Dataset(self.ncfile,'r')
        print(f'netcdf file {self.ncfile}')
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


    #name's last character is discarded. Eg. PECIN?
    def get_nvar_names(self, name:str) -> int:
        if name[-1] != '?':
            raise Exception(f'Error in search of array variables {name}')
        name = name[:-1]
        pattern = r'\b' + name + r'\d{1,2}\b'
        matches = [self.vars_[i] for i,s in enumerate(self.vars_) if re.search(pattern, s)]
        return len(matches)

    def get_ncvar(self, var:str):
        try:
            return np.array(self.nc.variables[var][:])
        except:
            return None

    def get_ncattr(name):
        try:
            return self.ncfile.getncattr(name)
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