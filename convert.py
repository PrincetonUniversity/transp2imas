

import os
import sys
import numpy as np

from nc2imas import nc2imas, get_args
from summary import imas_summary
from imas_equilibrium import imas_equilibrium
from imas_wall import imas_wall
from imas_core_profiles import imas_core_profiles
from imas_ec_launchers import imas_ec_launchers
from imas_lh_antennas import imas_lh_antennas
from imas_ic_antennas import imas_ic_antennas
from imas_nbi import imas_nbi
from imas_pf_active import imas_pf_active


#example use: python convert.py D3D 147634G14

machine, shot, rid = 'D3D', 147634, 'G14' #147634G14

#use machine, runid from command arguments
arg = get_args(sys.argv)

if len(sys.argv)>2:
    machine = sys.argv[1]
    runid = sys.argv[2]
    ii = 0
    for i in range(len(runid)):
      if runid[i].isdigit():
        ii = i+1
      else:
        break
    shot, rid = int(runid[:ii]), runid[ii:]

run = 111
backend='MDS'
newdb=True

ncfile = str(shot)+rid+'.CDF'
nml = str(shot)+rid+'TR.DAT'
occur=0
user=os.getenv('USER')
log='debug'  #for logger

all_ = True

#new_db = True if occur == 0 else False

user = os.path.expanduser('~') +'/public/imasdb'

#common
nci = nc2imas(db=machine, shot=shot, backend=backend, run=run, occur=occur, user=user,
          new_db=newdb, ncfile=ncfile, nml_file=nml, log=log)
nci.logger.info(f'\n{ncfile}, {nml}')

if all_:
    print('summary')
    imas_summary(nci)

    print('equilibrium')
    imas_equilibrium(nci)

    print('core_profiles')
    imas_core_profiles(nci)

    print('ec_launchers')
    imas_ec_launchers(nci)

    print('ic_antenna')
    imas_ic_antennas(nci)

    print('lh_antenna')
    imas_lh_antennas(nci)

    print('nbi')
    imas_nbi(nci)

    print('pf_active')
    imas_pf_active(nci)

    print('wall')
    imas_wall(nci)


# close datafile
nci.close_imas_database()


