#!/usr/bin/python
"""Example reader for inspecting IDS data through the legacy-style workflow."""

import imas
import numpy
import sys
import os


dbname='AUGD'
shot = 39220
run = 1301


def test_read_ids():
    conn = imas.ids(shot, run, 0, 0)
    conn.open_env(os.environ['USER'],dbname,os.environ['IMAS_VERSION'][:1])
    ids = conn.core_profiles


    ids.get(1)
    nprof1d = len(ids.profiles_1d)
    print('cp_id_len', nprof1d)
    print('time:', ids.time)

    for i in range(nprof1d):
      nions = len(ids.profiles_1d[i].ion)
      for j in range(nions):
         ion = ids.profiles_1d[i].ion[j]
         nel = len(ion.element)
         nstates = len(ion.state)
         print('time', i, 'ion', j, ion.label, ion.z_ion, nel, nstates, nstates)
         for iel in range(nel):
            elem = ion.element[iel]
         for ist in range(nstates):
            print(ion.state[ist].label)


    conn.close()


if __name__ == '__main__':
   test_read_ids()