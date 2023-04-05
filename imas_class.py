import os
import sys
from datetime import datetime
from types import NoneType
import numpy as np

import imas



class Convertor:
    def __init__(self, db, shot, run, occur:int=0, new_db:bool=False,
      _deb:bool=False, user:str=None, backend:str=None):
        self.imas_dbname = db
        self.shot = shot
        self.run = run
        self.occur = occur
        self._deb = _deb
        self.create_new = new_db
        self.user = user
        self.backend = backend
        self.outdb = None
        self.init_imas_database()
        self.get_imas_database()

    def fill_properties(self, ids):
        ids.ids_properties.comment = 'TRANSP to IMAS'
        ids.ids_properties.creation_date = str(datetime.now())
        ids.code.name = 'transp2imas'
        ids.code.repository = 'https://github.com/PrincetonUniversity/transp2imas'
        return ids

    def print_debug(self, msg:str=''):
        if self._deb:
            print(msg)
            print(' backend:', self.backend, 'run:', self.run)
            print(' create_new:', self.create_new, 'user:', self.user)
            print(' imas_db:', self.imas_dbname, 'shot:', self.shot)

    #for MDS+ the directory needs to exist
    def create_dir(self, public_db:bool=True):
        if self.backend != imas.imasdef.MDSPLUS_BACKEND:
            return
        version = str(os.getenv('IMAS_VERSION')[0])
        hom = os.path.expanduser('~')
        if self.imas_dbname[0] == '~':
            n = 2 if self.imas_dbname[1] == '/' else 1
            dbpath = hom + self.imas_dbname[n:]
        elif self.imas_dbname[0] != '/':
            dbpath = hom + '/public/imasdb/' + self.imas_dbname
        else:
            dbpath = self.imas_dbname

        outdir = f"{dbpath}/{version}/{self.run}"
        if not os.path.exists(outdir):
            if self._deb:
                print('Creating', outdir)
            os.makedirs(outdir)

    def init_imas_database(self):
        if self.outdb is not None:
            print('outdb exists')
            return
        if not self.backend:
            self.backend = imas.imasdef.MDSPLUS_BACKEND
        if not self.user:
            self.user = os.getenv('USER')
        if self.create_new is None:
            self.create_new = True

        version = str(os.getenv('IMAS_VERSION')[0])
        #TODO thread lock
        if self.outdb is None:
            self.print_debug('* variables to be used:')
            self.create_dir()
            try:
                self.outdb = imas.ids(self.shot, self.run, 0, 0)
            except:
                raise Exception('Error initializing IMAS database')

    def get_imas_database(self, create:bool=False):
        if create:
            self.create_new = True
        if not self.outdb:
            raise Exception('outdb is not valid, call init_database first')

        if self.create_new:
            status, idx = self.outdb.create_env_backend(self.user,self.imas_dbname, \
                                                        os.environ['IMAS_VERSION'], \
                                                        self.backend) 
        else:
            status, idx = self.outdb.open_env_backend(self.user,self.imas_dbname, \
                                                      os.environ['IMAS_VERSION'], \
                                                      self.backend) 
            
        if status != 0:
                raise Exception('Error creating/openning the database')


    def putdb_imas(self, ids, occur:int=0):
        if ids:
            ids.put(occur)
        else:
            raise TypeError('outdb not valid')


    def close_imas_database(self):
        if self.outdb:
            self.outdb.close()
            del self.outdb
        else:
            raise TypeError('outdb not valid')


if __name__ == '__main__':
    occur = 0
    _deb = True
    new_db = True
    plot = True
    iter_ = False #ITER
    run = 0 #2302 #W02 #this is occur
    shot = 163303 #200001 #141711 #205430 # also 142000, 141711
    db = 'D3D' #NSTX-U
    user = os.getenv('USER') # or 'public'
    dbpath = '/public/imasdb/'
    ncfile = '163303W02.CDF'

    #create convertor
    conv = Convertor(db, shot, run, occur, new_db, _deb, user)
    #call IDS creator

    # close datafile
    conv.close_imas_database()
