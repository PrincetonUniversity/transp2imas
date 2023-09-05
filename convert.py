

import os
import sys
import numpy as np

from nc2imas import nc2imas, get_args
from utilities import *
from imas_summary import *
from imas_equilibrium import *
from imas_wall import *
from imas_core_profiles import *
from imas_ec_launchers import *
from imas_lh_antennas import *
from imas_ic_antennas import *
from imas_nbi import *
from imas_tf import imas_tf
from imas_pf_active import *
from imas_pf_passive import *
from imas_edge_sources import *
from imas_sawteeth import *
from imas_radiation import *
from imas_core_sources import *
from imas_waves import *
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
def do_convert(argvs):
    #get arguments from command line
    #if -i all -d <dir> ;instead of CDF file name, then search for all CDF files in that dir,
    #and run in a loop
    backend, run, occur, user, machine, newdb, ncfile, nml, log = get_args(argvs) #(sys.argv)

    nci = nc2imas(backend=backend, run=run, occur=occur, user=user,
         machine=machine, new_db=newdb, ncfile=ncfile, nml_file=nml, log=log)
    nci.logger.info(f'{"-"*30} \n{ncfile}, {nml}')

    #test using 97781K31 different profile size

    if all_:
        print('wall')
        imas_wall(nci)

        print('lh_antenna') 
        imas_lh_antennas(nci)

        print('ic_antenna') 
        imas_ic_antennas(nci)

        # to be fixed: since homogeneous time is 1, all 1D waveform must have the same time base
        # to be fixed in pf_active: currents are COCOS-dependent
        print('pf_active') 
        imas_pf_active(nci)

#  O H I O => Only Handle It Once ...

# PRIORITY

        # COCOS consistency
        #  TRANSP uses COCOS=5 
        # for the time being I suggest to do the following for the Milestone of FY23: 
        # - since TRANSP is COCOS=5 => write the CDF output as is
        # - add a comment where relevant that the IDS has been written assuming COCOS=5 (used by TRANSP)

        # write q profile => need to convert from zone boundary to zone center
        # fill/fix vacuum_toroidal_field => this is COCOS dependent
        # check global quantities: ip and b0 are COCOS dependent
        # check radial grids => psi is COCOS dependent
        # check units Weber vs Weber/rad
        # check values of psi at center and at boundary
        # focus on electrons, ignore ions for now
        # need to check the size of the 1D and 2D time vector.
        # if the size is different, homogeneous_time should be 0, not 1  
        # and the arrays need to be filled consistently
        print('core_profiles')
        imas_core_profiles(nci)

        # need to write TF IDS => only r0 and b_field_tor_vacuum_r => COCOS dependent    
        print('tf')
        imas_tf(nci)

         # to be fixed: consistency with homogenous_time
#        print('nbi') # 175286T57.CDF
#        imas_nbi(nci)

         # to be fixed: steering angles/definition/consistency
         # consistency with homogenous_time
#        print('ec_launchers')
#        imas_ec_launchers(nci)


# IGNORE EVERYTHING LISTED BELOW THIS LINE. NEED TO COMPLETE EVERYTHING ABOVE THIS LINE FIRST


# everything below this line needs to be verified

        # to be fixed: everything
#        print('pf_passive')
#        imas_pf_passive(nci)


        # to be fixed: everything
#        print('radiation') #33618Y19.CDF  38530P94
#        imas_radiation(nci)

#        print('summary')
#        imas_summary(nci)

         # to be fixed: remove moment calculations
         # all global quantities
         # RBZ (and Ro, Bo)
         # boundary
         # handling of size of arrays
         # transpose/not transpose
         # COCOS transform
         # grids: chec consistency
         # rho_tor_norm is incorrect
         # psi is COCOS dependent (both 1D and 2D)
#        print('equilibrium') #142202U01.CDF
#        imas_equilibrium(nci)

        # write q profile => need to convert from zone boundary to zone center
        # fill/fix vacuum_toroidal_field => this is COCOS dependent
        # check global quantities: ip and b0 are COCOS dependent
        # check radial grids => psi is COCOS dependent
        # check units Weber vs Weber/rad
        # check values of psi at center and at boundary
#        print('core_profiles')
#        imas_core_profiles(nci)

        # need to write TF IDS => only r0 and b_field_tor_vacuum_r => COCOS dependent    

        # to be fixed: everything
#        print('imas_core_sources') # 112989Z09(NSTX), 175286T57
#        imas_core_sources(nci)

        # to be fixed: everything
#        print('imas_waves')
#        imas_waves(nci)


        #print('core transport') #89846.CDF
        #imas_core_transport(nci)

        #print('edge sources') #112989Z09
        #imas_edge_sources(nci)

        #print('sawteeth')
        #imas_sawteeth()

    # close datafile
    nci.close_imas_database()


#if -i all -d <dir> ; then run in a loop; otherwise run provided input file only
if __name__ == '__main__':
    argvs = sys.argv
    indf = [i for i,x in enumerate(argvs) if x=='all' and i>0 and argvs[i-1]=="-i"]
    if len(indf)==0:
        do_convert(argvs)
    else:
        i = indf[0]
        indd = [i+1 for i,x in enumerate(argvs) if x=='-d' or x=='--input-dir']
        if len(indd)==0:
            print('Couldnot get valid directory')
        else:
            dir_ = argvs[indd[0]]
            fls = find_all_files(".CDF", dir_)
            for f in fls:
                argvs[i] = f
                do_convert(argvs)
