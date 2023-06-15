

import os
import sys
import numpy as np

print('imas loading...')
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


#example use: python convert.py -m D3D -i 147634G14 -r 1 -b MDS -new -o 0 -l debug

#example. These will be replaced below
mach = 'D3D'
nput_cdf = '147634G14.CDF'
run = 1
backend = 'MDS'
newdb = True
occur = 0
user = os.getenv('USER')
log = 'debug'  #for logger


all_ = True

#get arguments from command line
mach, shot, backend, run, occur, user, newdb, ncfile, nml, log = get_args(sys.argv)

nci = nc2imas(db=mach, shot=shot, backend=backend, run=run, occur=occur, user=user,
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


