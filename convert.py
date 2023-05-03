
import os
import sys
import numpy as np
from subprocess import run

from summary import Summary
from imas_class import Convertor

import os
import sys
import numpy as np
from subprocess import run

from utils import *
from equilibrium import Equilibrium
from limiters import Limiters
from core_profiles import CoreProfiles
from thomson_scattering import Thomson
from coil_parameters import PFActive
from charge_exchange import CHERS
from nbi import NBI
from convert import Convertor

mds_shot = 141711
nstx = True
nstxu = False
deb = False
incr = 0
plot = True

thoms = False
wall = False
cprof = False
equil = False
pf_active = False
chers = False
nbi = False

all_ = True
if len(sys.argv) > 1:
    all_ = False
    var = sys.argv[1]
    if var == 'all':
        all_ = True
    elif var == 'pf':
        pf_active = True
    elif var == 'ts':
        thoms = True
    elif var == 'lim':
        wall = True
    elif var == 'cp':
        cprof = True
    elif var == 'eq':
        equil = True
    elif var == 'cx':
        chers = True
    elif var == 'nb':
        nbi = True
    else:
        print('option not defined')
        sys.exit()


if len(sys.argv) > 2:
    var = sys.argv[2]
    if var.upper() == 'Q':
        deb = False
    elif var == 'NSTX':
        nstx = True
        nstxu = False
    elif var == 'NSTXU' or var == 'NSTX-U':
        nstxu = True
        nstx = False

if len(sys.argv) > 3:
    var = sys.argv[3]
    if var == 'i' or var == 'incr':
        incr = 1

mds_run = 2502 #Y02
run = 1
occur = 0
new_db = True if occur == 0 else False
if nstx:
    db = 'NSTX'
    shot = 100001
elif nstxu:
    db = 'NSTX-U' #NSTX
    shot = 200001  #100001

user = os.path.expanduser('~') +'/public/imasdb'
backend=imas.imasdef.HDF5_BACKEND
#create convertor
conv = Convertor(db, shot, run, occur, mds_shot, new_db, _deb=deb,
   user=user, backend=backend, plot=plot)

#occur passed in Limiters(..) etc. override the occur in convertor.
if wall or all_:
    print('\nLimiters...')
    Limiters(conv, occur+0, mds_shot=mds_shot)

if cprof or all_:
    print('\nCore Profiles....')
    CoreProfiles(conv, occur+incr, mds_shot=mds_shot)

if equil or all_:
    print('\nEequilibrium...')
    Equilibrium(conv, occur+incr, mds_shot=mds_shot)

if thoms or all_:
    print('\nThomson scattering (MPTS)...')
    Thomson(conv, occur+incr, mds_shot=mds_shot)

if pf_active or all_:
    #valid shot 205430
    mds_shot = 205430
    print('\nCoil PFActive...', mds_shot)
    #override mds_shot set in Convertor
    PFActive(conv, occur+incr, mds_shot=mds_shot)

if chers or all_:
    mds_shot = 141711
    print('\n CHERS...', mds_shot)
    CHERS(conv, occur+incr, mds_shot=mds_shot)

if nbi or all_:
    mds_shot = 205430
    print('\n NBI...', mds_shot)
    NBI(conv, occur+incr, mds_shot=mds_shot)

# close datafile
conv.close_imas_database()








run = 1
occur = 0
new_db = True if occur == 0 else False
shot = 100001

db = 'D3D'
ncfile = '163303W02.CDF' #in pwd

user = os.path.expanduser('~') +'/public/imasdb'
backend=imas.imasdef.HDF5_BACKEND

#create convertor. Use this for all IDS
conv = Convertor(db, shot, run, occur, new_db, user=user, backend=backend,
       _deb=True, plot=False)

#make summary IDS
Summary(conv, ncfile)



# close datafile
conv.close_imas_database()