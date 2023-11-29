import os
import sys
from datetime import datetime
# from types import NoneType
import numpy as np

import imas


class imas_connection:
    def __init__(self, db:str, shot:int, run:int=0, backend:str=None,
                 occur:int=0, user:str=None, new_db:bool=False):
        self.imas_dbname = db
        self.shot = shot
        self.run = run
        self.occur = occur
        self.create_new = new_db
        self.user = user
        self.backend = backend
        self.ids = None
        self.outdir = ''
        self.iter_ = True if db[:4].upper()=='ITER' else False
        self.init_imas_database()
        self.get_imas_database()

    def print_debug(self, msg:str='', force=False):
        print(msg)
        print(' user', self.user)
        print(' create_new:', self.create_new, '; occur:', self.occur)
        print(' machine/imas_db:', self.imas_dbname, ' backend:', self.backend)
        print('shot:', self.shot, '; converted run number:', self.run)

    #for MDS+ backend the directory needs to exist in ~/public/imasdb
    def set_outdir(self, public_db:bool=True):
        diry = self.user
        #if home dir, use default path
        if self.user == os.getenv('USER'):
            diry = os.path.expanduser('~') + '/public/imasdb'

        dbpath = diry + '/' + self.imas_dbname
        version = str(os.getenv('IMAS_VERSION')[0])
        dbpath = dbpath + '/' + version
        if self.backend == imas.imasdef.HDF5_BACKEND:
            dbpath = dbpath + '/' + str(self.shot) + '/' + str(self.run)
        elif self.backend == imas.imasdef.MDSPLUS_BACKEND:
            dbpath = dbpath + '/0/' #occur?
        self.outdir = dbpath


    def init_imas_database(self):
        if self.ids is not None:
            print('ids exists')
            return
        mkdir = False
        if self.backend=='HDF5' or self.backend=='HDF':
            self.backend = imas.imasdef.HDF5_BACKEND
        elif self.backend=='MDS' or self.backend=='MDSPLUS':
            self.backend = imas.imasdef.MDSPLUS_BACKEND
            mkdir = True

        if not self.user:
            self.user = os.getenv('USER')
        if self.create_new is None:
            self.create_new = True

        self.set_outdir()

        #TODO thread lock
        if self.ids is None:
            if self.shot > 214747 and self.backend==imas.imasdef.MDSPLUS_BACKEND:
                print('adjusting shot number for MDS backend')
                self.shot = self.shot%214747
            if mkdir and not os.path.exists(self.outdir):
                os.makedirs(self.outdir) #only for MDS+

            try:
                self.ids = imas.ids(self.shot, self.run, 0, 0)
            except:
                raise Exception('Error initializing IMAS database')

    def get_imas_database(self):
        if not self.ids:
            raise Exception('ids is not valid, call init_database first')

        if self.create_new:
            try:
                status, idx = self.ids.create_env_backend(self.user,self.imas_dbname,
                                            os.environ['IMAS_VERSION'], self.backend)
            except Exception:
                status = 1
                import traceback
                traceback.print_exc(file=sys.stdout)
            if status != 0:
                self.print_debug('Error creating new databse using input:')
                raise sys.exit(1)
        else:
            status, idx = self.ids.open_env_backend(self.user,self.imas_dbname,
                                       os.environ['IMAS_VERSION'], self.backend)
            self.print_debug(f'Attempting to open existing databse for run {self.run}')
            if status != 0:
                mes = f'\n ***Error opening database for run {self.run}.'
                mes = mes + ' Make sure it exists, else use -new'
                raise Exception(mes)


    def putdb_imas(self, ids):
        if ids:
            ids.put(self.occur)
        else:
            raise TypeError('ids not valid')


    def close_imas_database(self):
        if self.ids:
            self.ids.close()
            del self.ids
        else:
            raise TypeError('ids not valid')
        print('ids stored in', self.outdir)


if __name__ == '__main__':
    print("Please use: python convert.py")
