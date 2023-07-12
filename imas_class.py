import os
import sys
from datetime import datetime
from types import NoneType
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
        self.iter_ = True if db[:4].upper()=='ITER' else False
        self.init_imas_database()
        self.get_imas_database()

    def print_debug(self, msg:str='', force=False):
        print(msg)
        print(' create_new:', self.create_new, '; occur:', self.occur)
        print(' machine/imas_db:', self.imas_dbname, ' backend:', self.backend)
        print('shot:', self.shot, '; converted run number:', self.run)

    #for MDS+ backend the directory needs to exist in ~/public/imasdb
    def create_dir(self, public_db:bool=True):
        if self.backend != imas.imasdef.MDSPLUS_BACKEND:
            return
        hom = os.path.expanduser('~')
        dbpath = f'{hom}/public/imasdb/{self.imas_dbname}'
        version = str(os.getenv('IMAS_VERSION')[0])
        outdir = f"{dbpath}/{version}/0"
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    def init_imas_database(self):
        if self.ids is not None:
            print('ids exists')
            return
        if self.backend=='HDF5':
            self.backend = imas.imasdef.HDF5_BACKEND
        elif self.backend=='MDS':
            self.backend = imas.imasdef.MDSPLUS_BACKEND
        if not self.user:
            self.user = os.getenv('USER')
        if self.create_new is None:
            self.create_new = True

        #TODO thread lock
        if self.ids is None:
            self.create_dir() #only for MDS+
            try:
                self.ids = imas.ids(self.shot, self.run, 0, 0)
            except:
                raise Exception('Error initializing IMAS database')

    def get_imas_database(self):
        if not self.ids:
            raise Exception('ids is not valid, call init_database first')
        if self.create_new:
            status, idx = self.ids.create_env_backend(self.user,self.imas_dbname,
                                            os.environ['IMAS_VERSION'], self.backend)
            self.print_debug('create new databse')
            if status != 0:
                raise Exception('Error creating database')
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

if __name__ == '__main__':
    print("Please use: python convert.py")
