#!/usr/bin/python
"""Command-line entry point for converting TRANSP CDF output into IMAS IDSs."""



import sys

from nc2imas import nc2imas, get_args, TRANSP2IMAS_VERSION
from utilities import find_all_files
from imas_summary import imas_summary
from imas_equilibrium import imas_equilibrium
from imas_wall import imas_wall
from imas_core_profiles import imas_core_profiles
from imas_ec_launchers import imas_ec_launchers
from imas_lh_antennas import imas_lh_antennas
from imas_ic_antennas import imas_ic_antennas
from imas_nbi import imas_nbi
from imas_tf import imas_tf
from imas_pf_active import imas_pf_active
from imas_pf_passive import imas_pf_passive
from imas_edge_sources import imas_edge_sources
from imas_radiation import imas_radiation
from imas_core_sources import imas_core_sources
from imas_waves import imas_waves
from imas_workflow import imas_workflow
from imas_dataset_fair import imas_dataset_fair
from imas_distribution_sources import imas_distribution_sources
from imas_core_transport import imas_core_transport

TRANSP2IMAS_DD4_PATCHLEVEL = TRANSP2IMAS_VERSION
all_ = True
def do_convert(argvs):
    print(f'transp2imas patchlevel {TRANSP2IMAS_DD4_PATCHLEVEL}')
    #get arguments from command line
    #if -i all -d <dir> ;instead of CDF file name, then search for all CDF files in that dir,
    #and run in a loop
    parsed_args = get_args(argvs)
    if len(parsed_args) == 16:
        backend, run, occur, user, machine, newdb, ncfile, nml, log, \
          moments, nzones, deb, ec_nml_launch, ip_direction, bt_direction, dd_version = parsed_args
    elif len(parsed_args) == 15:
        backend, run, occur, user, machine, newdb, ncfile, nml, log, \
          moments, nzones, deb, ec_nml_launch, ip_direction, bt_direction = parsed_args
        dd_version = 'auto'
    elif len(parsed_args) == 13:
        backend, run, occur, user, machine, newdb, ncfile, nml, log, \
          moments, nzones, deb, ec_nml_launch = parsed_args
        ip_direction = 'auto'
        bt_direction = 'auto'
        dd_version = 'auto'
    else:
        raise ValueError(
            f'Unexpected get_args() return length {len(parsed_args)}; expected 13, 15 or 16'
        )

    nci = nc2imas(backend=backend, run=run, occur=occur, user=user,
         machine=machine, new_db=newdb, ncfile=ncfile, nml_file=nml, log=log, debug=deb,
         ip_direction=ip_direction, bt_direction=bt_direction, dd_version=dd_version)
    nci.logger.info(f'{"-"*30} \n{ncfile}, {nml}')
    print(f'IMAS target DD: {nci.dd_version} (schema family DD{nci.dd_major_version})')


    if all_:
        print('wall')
        imas_wall(nci)

        print('tf')
        imas_tf(nci)

        print('lh_antenna') 
        imas_lh_antennas(nci)

        print('ic_antenna') 
        imas_ic_antennas(nci)

        print('ec_launchers')
        # get launch position from namelist
        imas_ec_launchers(nci, ec_nml_launch=ec_nml_launch)

        # NBI mapping has been audited; writer-specific time handling is kept local.
        print('nbi') 
        imas_nbi(nci)

        # NUBEAM birth/deposition source distribution. This is distinct from
        # the full fast-ion distributions IDS and is available even when FBM
        # phase-space output was not saved.
        print('distribution_sources')
        imas_distribution_sources(nci)

        # PF circuit-current signs are independent of plasma Ip/Bphi orientation.
        print('pf_active') 
        imas_pf_active(nci)

        # COCOS/sign consistency is handled centrally by nc2imas.
        # DD4 is written in COCOS=17.  Physical Ip/Bphi directions are obtained
        # from TRANSP NLJCCW/NLBCCW (or explicit command-line overrides), while
        # psi/q/phi and signed current profiles are transformed accordingly.

        # Core profiles, radial grids, units, and COCOS-dependent signs have
        # been audited against the regression CDFs.
        print('core_profiles')
        imas_core_profiles(nci)

        # Passive structures remain conditional on TRANSP PS_n channels and
        # need a dedicated passive-current regression case for full validation.
        print('pf_passive')
        imas_pf_passive(nci)


        # Dedicated spatial radiation output: bremsstrahlung, synchrotron,
        # and line/impurity radiation when the corresponding TRANSP profiles
        # are present.  P0RECO is neutral kinetic source power, not radiation.
        print('radiation')
        imas_radiation(nci)

# summary has been re-audited against TRANSP energy/current units and DD4.
# It is now part of the standard conversion.

        # Equilibrium geometry, grids, moments, units, and COCOS transforms
        # have been re-audited against the regression cases.
        print('equilibrium') #142202U01.CDF
        imas_equilibrium(nci, conv_moments=moments, nzones_in=nzones)

        print('summary')
        imas_summary(nci)

        # Audited source mappings are documented in README.md.
        print('core_sources')
        imas_core_sources(nci)

        # RF wave propagation/deposition.  The writer is conservative: it
        # creates waves only when supported absorbed/deposited EC or IC output
        # is present, and leaves launcher-only runs to the antenna IDSs.
        print('waves')
        imas_waves(nci)

        print('core_transport')
        imas_core_transport(nci)

        print('edge_sources') #112989Z09
        imas_edge_sources(nci)

        #print('sawteeth')
        #imas_sawteeth()

        # Dataset-level metadata and workflow provenance. These are written
        # last so they describe the complete converted entry. Whole-IDS source
        # provenance is also attached by fill_properties() to every IDS.
        print('workflow')
        imas_workflow(nci)

        print('dataset_fair')
        imas_dataset_fair(nci)

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
