

import os
import sys
import numpy as np

print('imas loading...')
from nc2imas import nc2imas, get_args
from summary import *
from imas_equilibrium import *
from imas_wall import *
from imas_core_profiles import *
from imas_ec_launchers import *
from imas_lh_antennas import *
from imas_ic_antennas import *
from imas_nbi import *
from imas_pf_active import *
from imas_pf_passive import *
from imas_edge_sources import *
#from imas_core_sources import *
#from imas_core_transport import *

#example use: python convert.py -i 147634G14 -r 1 -b MDS -new -o 0 -l debug

#example. These will be replaced below
nput_cdf = '147634G14.CDF'
run = 1
backend = 'MDS'
newdb = True
occur = 0
user = os.getenv('USER')
log = 'debug'  #for logger


all_ = True

#get arguments from command line
backend, run, occur, user, newdb, ncfile, nml, log = get_args(sys.argv)

nci = nc2imas(backend=backend, run=run, occur=occur, user=user,
          new_db=newdb, ncfile=ncfile, nml_file=nml, log=log)
nci.logger.info(f'\n{ncfile}, {nml}')


imas_pf_passive(nci) # 142202U01 NSTX
all_ = False



if all_:
    print('summary')
    imas_summary(nci)

    print('equilibrium') #142202U01.CDF
    imas_equilibrium(nci)

    print('core_profiles')
    imas_core_profiles(nci)

    print('ec_launchers')
    imas_ec_launchers(nci)

    print('ic_antenna') #89846.CDF
    imas_ic_antennas(nci)

    print('lh_antenna') #33618Y19.CDF
    imas_lh_antennas(nci)

    print('nbi') # 175286T57.CDF
    imas_nbi(nci)

    print('pf_active') # 142202U01 NSTX
    imas_pf_active(nci)

    print('pf_passive')
    imas_pf_passive(nci)

    print('wall')
    imas_wall(nci)

    print('core transport') #89846.CDF
    #imas_core_transport(nci)

    print('edge sources') #112989Z09
    #imas_edge_sources(nci)

    print('imas_core_sources') # 112989Z09(NSTX), 175286T57
    #imas_core_sources(nci)

# close datafile
nci.close_imas_database()


