#!/usr/bin/python3
#
import os
import sys
import numpy as np
from scipy import constants
from scipy.interpolate import interp1d
from datetime import datetime
import re
from pathlib import Path
import logging
import argparse
from netcdf_class import netcdf_connection
from imas_class import imas_connection
#import imas
      
class nc2imas(imas_connection, netcdf_connection):
    def __init__(self, db, shot, run, backend, occur, user, new_db, ncfile,
         nml_file:str=None, log:str='error'):
        imas_connection.__init__(self, db, shot, run, backend, occur, user, new_db)
        netcdf_connection.__init__(self, ncfile)

        self.nml_name = nml_file
        self.logger = None
        #namelist contents dict
        self.nmlist = None
        self.set_logger(log)
        self.logger.info('creating nc2imas')
        self.load_namelist()

    def checkErrExit(*variables):
        for var in variables:
            if isinstance(var, np.ndarray):
                sh = var.shape
                if sh[0]==0 or (len(sh)>1 and sh[1]==0):
                    self.logger.Exception('empty np.array')
            elif isinstance(var, (list, tuple, set)):
                if len(var) == 0:
                    self.logger.Exception('empty array')
            else:
                if var == 0:
                    self.logger.Exception(f'empty {var}')

    def fill_value(self, ids, dat, rsz:int=-1, mult:float=None, name:str=None):
        name_ = name if name else 'data'
        if dat is None:
            self.logger.info(f'{name} is None')
            return

        nresize = len(dat) if rsz==-1 else rsz

        if nresize > 0:
            ids.value.resize(nresize, refcheck=False)
        if mult is not None:
            dat = mult*dat
        ids.value[:] = dat
        #if self.src:
        #    ids.source = self.src

    def fill_value_nc(self, ids, var:str, rsz:int=-1, mult:float=None):
        dat = self.get_ncvar(var)
        if dat is None:
            self.logger.info(f'{var} not found')
        else:
            self.fill_value(ids, dat, rsz=rsz, mult=mult, name=var)

    def fill_data_nc(self, ids, var:str, mult:float=None):
        dat = self.get_ncvar(var)
        if dat is None:
            self.logger.info(f'{var} not found')
        else:
            ids.resize(len(dat), refcheck=False)
            if mult is not None:
                dat = mult*dat
            ids[:] = dat

    def fill_data(self, ids, dat, mult:float=None):
        if not(dat is None or np.isscalar(dat)):
            ids.resize(len(dat), refcheck=False)
            if mult is not None:
                dat = mult*dat
            ids[:] = dat

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
            if dat is None:
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
        ids.ids_properties.comment = 'TRANSP output to IMAS IDS'
        ids.ids_properties.creation_date = str(datetime.now())
        #ids.ids_properties.version_put = str(os.getenv('IMAS_VERSION'))
        ids.ids_properties.data_dictionary = str(os.getenv('IMAS_VERSION'))
        ids.ids_properties.access_layer = str(os.getenv('IMAS_VERSION'))
        #ids_properties/provenance.node.resize()
        return ids

    def fill_code_properties(self, ids):
        ids.code.name = 'transp2imas'
        ids.code.repository = 'https://github.com/PrincetonUniversity/transp2imas'
        return ids

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
                                val = float(val)
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
        except:
            print(f'**Error: {name}: namelist not found')
            return 0
        return len(matches)


    #return list. name's last character is discarded. Eg. XECECH?
    def get_namelist_vars(self, name:str, tp=None, scale=0):
        n = self.get_namelist_nvars(name)
        name0 = name[-1]
        try:
            val = [self.nmlist[f'{name0}({i})'] for i in range(1,n+1)]
        except:
            print(f'** Failed to get {name} from name list')
            if self.nmlist is None:
                print(f'** {self.nmlist} is not set')
            return []
        else:
            if tp in [int, float]:
                val = [tp(x) for x in val]
            if scale > 0:
                val = scale*val
            return val

    #returns single value
    def get_namelist_var(self, name:str, tp=None):
        try:
            val = self.nmlist[name]
        except:
            print(f'** Failed to get {name} from name list')
            if self.nmlist is None:
                print(f'** {self.nmlist} is not set')
        else:
            if tp in [int, float]:
                val = tp(val)
            return val

    def set_logger(self, log:str = None):
        self.logger = logging.getLogger(__name__)
        levels = {'debug': logging.DEBUG, 'info': logging.INFO,
                  'error': logging.ERROR}
        log = log if log is not None else 'error'
        level = levels.get(log.lower())
        self.logger.setLevel(level=level)
        self.logger.addHandler(logging.FileHandler('log.out', mode='w'))


def get_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument("-l", "--log", help="Provide logging level")
    parser.add_argument('-s', '--shot', type=int, help='integer shot number')
    parser.add_argument('-i', '--input-cdf', type=str, help='CDF TRANSP output file')
    parser.add_argument('-b', '--backend', type=str, help='backend: MDS or HDF')
    parser.add_argument('-d', '--db', type=str, help='tokamak or database')
    parser.add_argument('-u', '--user', type=str, help='username or path of database')
    parser.add_argument('-r', '--run', type=int, help='run number')
    parser.add_argument('-o', '--occur', type=int, help='occurrence')
    parser.add_argument('-n', '--new_db', help='new database', action="store_true")
    parser.add_argument('-df', '--default-val', help='use default', action="store_true")

    args = parser.parse_args()
    if not (args.db or args.shot or args.run or args.occur or args.input_cdf):
        pass
        #parser.print_help()
        #exit(1)
    return args


if __name__ == '__main__':
    a = get_args(sys.argv)
    nc = nc2imas(db='D3D', shot=147634, backend='MDS', run=1, occur=0, user=os.getenv('USER'),
           new_db=True, ncfile='147634G14.CDF', nml_file='147634G14TR.DAT', log='debug')

