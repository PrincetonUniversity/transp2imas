#!/usr/bin/python
"""IMAS database-entry helpers shared by the individual IDS writer modules."""

import os
import sys
from datetime import datetime
# from types import NoneType
import numpy as np

import imas


class imas_connection:
    def __init__(self, db:str, shot:int, run:int=0, backend:str=None,
                 occur:int=0, user:str=None, new_db:bool=False, dd_version:str=None):
        self.imas_dbname = db
        self.shot = shot
        self.run = run
        self.occur = occur
        self.create_new = new_db
        self.user = user
        self.backend = backend
        # Select the target Data Dictionary explicitly when requested.  This is
        # independent of the installed IMAS-Python version: modern IMAS-Python
        # can create/write IDS objects for several DD versions in one environment.
        requested_dd = None if dd_version in (None, '', 'auto') else str(dd_version)
        env_dd = os.getenv('IMAS_VERSION') or os.getenv('IMAS_DD_VERSION')
        if requested_dd:
            self.dd_version = requested_dd
        elif env_dd:
            self.dd_version = env_dd
        else:
            # Use the actual IMAS-Python default instead of assuming DD4.
            self.dd_version = str(imas.IDSFactory().version)
        self.dd_version_tuple = self._parse_dd_version(self.dd_version)
        self.dd_major_version = str(self.dd_version_tuple[0])
        self.ids = None
        # Updated only after IMAS DBEntry.put() accepts an IDS.
        self._written_ids = []
        self._active_writer = None
        self.outdir = ''
        self.iter_ = True if db[:4].upper()=='ITER' else False
        self.init_imas_database()
        self.get_imas_database()


    @staticmethod
    def _parse_dd_version(version):
        """Convert a DD version string such as 3.39.0 or 4.1.1 to a tuple."""
        vals = []
        for token in str(version).split('.')[:3]:
            digits = ''.join(ch for ch in token if ch.isdigit())
            vals.append(int(digits) if digits else 0)
        while len(vals) < 3:
            vals.append(0)
        return tuple(vals)

    def dd_at_least(self, major, minor=0, patch=0):
        return self.dd_version_tuple >= (int(major), int(minor), int(patch))

    def dd_before(self, major, minor=0, patch=0):
        return self.dd_version_tuple < (int(major), int(minor), int(patch))

    @property
    def is_dd4(self):
        return self.dd_version_tuple[0] >= 4

    @property
    def is_dd3(self):
        return self.dd_version_tuple[0] == 3

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
        dbpath = dbpath + '/' + self.dd_major_version
        if self.backend == imas.ids_defs.HDF5_BACKEND:
            dbpath = dbpath + '/' + str(self.shot) + '/' + str(self.run)
        elif self.backend == imas.ids_defs.MDSPLUS_BACKEND:
            dbpath = dbpath + '/0/' #occur?
        self.outdir = dbpath


    def init_imas_database(self):
        if self.ids is not None:
            print('ids exists')
            return
        mkdir = False
        if self.backend=='HDF5' or self.backend=='HDF':
            self.backend = imas.ids_defs.HDF5_BACKEND
        elif self.backend=='MDS' or self.backend=='MDSPLUS':
            self.backend = imas.ids_defs.MDSPLUS_BACKEND
            mkdir = True

        if not self.user:
            self.user = os.getenv('USER')
        if self.create_new is None:
            self.create_new = True

        self.set_outdir()

        #TODO thread lock
        if self.ids is None:
            if self.shot > 214747 and self.backend==imas.ids_defs.MDSPLUS_BACKEND:
                print('adjusting shot number for MDS backend')
                self.shot = self.shot%214747
            if mkdir and not os.path.exists(self.outdir):
                os.makedirs(self.outdir) #only for MDS+

            try:
                #self.ids = imas.ids(self.shot, self.run, 0, 0)
                # IMAS-Python API. The legacy pulse/run constructor remains supported
                # and is useful here because transp2imas already exposes database, pulse,
                # run and user independently. dd_version selects the Data Dictionary; when
                # unset IMAS-Python uses its configured default.
                kwargs = {}
                if self.dd_version:
                    kwargs["dd_version"] = self.dd_version
                self.ids = imas.DBEntry(
                    self.backend, self.imas_dbname, self.shot, self.run,
                    self.user, **kwargs
                )
                # Keep the converter's schema switches synchronized with the
                # actual factory used by the DBEntry.
                self.dd_version = str(self.ids.factory.version)
                self.dd_version_tuple = self._parse_dd_version(self.dd_version)
                self.dd_major_version = str(self.dd_version_tuple[0])
            except Exception as exc:
                raise RuntimeError(
                    f'Error initializing IMAS database for DD {self.dd_version}'
                ) from exc

    def get_imas_database(self):
        if not self.ids:
            raise Exception('ids is not valid, call init_database first')

        if self.create_new:
            try:
                # IMAS-Python create()/open() raise on failure; they do not return
                # the old Access-Layer (status, index) tuple.
                self.ids.create()
            except Exception as exc:
                self.print_debug('Error creating new database using input:')
                raise RuntimeError('Failed to create IMAS database entry') from exc
        else:
            self.print_debug(f'Attempting to open existing database for run {self.run}')
            try:
                self.ids.open()
            except Exception as exc:
                mes = f'\n ***Error opening database for run {self.run}.'
                mes += ' Make sure it exists, else use -new'
                raise RuntimeError(mes) from exc


    def new_ids(self, ids_name: str):
        """Create an empty IDS using the Data Dictionary attached to this DBEntry."""
        if not self.ids:
            raise TypeError('database entry is not valid')
        return self.ids.factory.new(ids_name)

    def _record_accepted_put(self, ids, occurrence):
        """Keep a readback manifest; never record a failed database write."""
        name = self._active_writer
        if name is None:
            metadata = getattr(ids, 'metadata', None)
            name = getattr(metadata, 'name', None) or getattr(ids, '_name', None)
        self._written_ids.append({'ids': str(name or 'unknown'),
                                  'occurrence': int(occurrence)})

    def putdb_imas(self, ids):
        if ids is None:
            raise TypeError('ids not valid')
        self.ids.put(ids, occurrence=self.occur)
        self._record_accepted_put(ids, self.occur)

    def putdb_imas_(self, ids):
        if ids is None:
            raise TypeError('ids not valid')
        for i in range(3):
            self.ids.put(ids, occurrence=i)
            self._record_accepted_put(ids, i)

    def putdb_imas_profile(self, ids, occur):
        """Write an explicit occurrence without changing self.occur."""
        if ids is None or not isinstance(occur, int):
            raise TypeError('ids/occur not valid')
        self.ids.put(ids, occurrence=occur)
        self._record_accepted_put(ids, occur)

    def close_imas_database(self):
        if self.ids:
            self.ids.close()
            del self.ids
        else:
            raise TypeError('ids not valid')
        print('ids stored in', self.outdir)


if __name__ == '__main__':
    print("Please use: python convert.py")
