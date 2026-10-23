#!/usr/bin/python
"""General utility functions shared across the TRANSP-to-IMAS converter."""

import os
import sys
from pathlib import Path
import math
import numpy as np
import inspect
import argparse

import mapping_status


def check_file_exists(name:str)->bool:
    path_ = Path(name).absolute()
    if path_.is_symlink():
        path_ = os.readlink(name)
    return os.path.isfile(path_)

def vacuum_field_reference(raxis_m, bzxr_tm, bt_sign, times, source_times):
    """Reference radius ``r0`` and vacuum toroidal field ``b0(times)`` for IMAS.

    IMAS requires ``vacuum_toroidal_field.r0 * b0(t)`` to equal the TF IDS vacuum field
    times radius (TRANSP ``BZXR`` with the toroidal-field sign applied).  ``r0`` is therefore
    one constant radius, the mean of the first and last magnetic-axis radius of the source
    time grid, and ``b0 = bt_sign * BZXR / r0``.  Dividing by the time-dependent axis radius
    instead breaks the product whenever the axis moves (up to 17% in the NSTX examples).

    ``raxis_m`` and ``bzxr_tm`` are on ``source_times``; ``b0`` is interpolated to ``times``.
    Returns ``(r0, b0)``.
    """
    raxis = np.asarray(raxis_m, dtype=np.float64)
    field_r = bt_sign * np.asarray(bzxr_tm, dtype=np.float64)
    r0 = 0.5 * (raxis[0] + raxis[-1])
    if not np.array_equal(times, source_times):
        field_r = np.interp(times, source_times, field_r)
    return float(r0), field_r / r0


def append_successful_entries(new_ids, ids, array_name, builders):
    """Fill an IDS array of structures with only the entries whose builder succeeds.

    IMAS-Python's ``resize(n)`` defaults to ``keep=False`` and recreates every element, so
    allocating a larger array, filling it and then shrinking with ``resize(used)`` erases
    all the filled entries (even a same-size resize does).  Filling slots in place also lets
    a builder that fails half way leave partial data in a slot the next builder reuses.
    Each builder therefore fills an empty element of a scratch IDS; only elements whose
    builder returned a true value are appended to ``ids.<array_name>``.

    ``new_ids`` is a callable returning an empty IDS of the same type (for example
    ``lambda: im.new_ids('radiation')``), ``builders`` an iterable of callables taking an
    empty element.  Returns the number of entries appended.
    """
    target = getattr(ids, array_name)
    used = 0
    for build in builders:
        scratch = getattr(new_ids(), array_name)
        scratch.resize(1)
        if build(scratch[0]):
            target.append(scratch[0])
            used += 1
    return used


def find_all_files(pat, path_):
    result = []
    for root, dir_, files in os.walk(path_):
        for fshort in files:
            if fshort.endswith(pat):
                result.append(os.path.join(root,fshort))
    return result


def find_file(fname, path_):
   result = []
   for root, dir_, files in os.walk(path_):
      if fname in files:
         result.append(os.path.join(root, fname))
   return result


#runid is also CDF filename with or without extension
def find_runid_shot_code(runid):
    cdf_name = Path(runid).stem
    #rids = [a for a in cdf_name if (a.isdigit() or a.isalpha())]
    #numbers = ''.join(rids)
    i0 = 0
    ii = 0
    found = False
    for i in range(len(cdf_name)):
      if cdf_name[i].isdigit():
        ii = i+1
        if not found:
           found = True
           i0 = i
      elif found:
        break

    try:
      shot = int(cdf_name[i0:ii])
      rid = cdf_name[ii:]
      ris = [i for i,x in enumerate(rid) if x=='_']
      ri = ris[0] if len(ris)>0 else len(rid)
      rid = rid[:ri]
      part = [str(ord(x)-96) if x.isalpha() else x for x in rid.lower()]
    except (ValueError, TypeError, IndexError) as exc:
      mes = f'Error getting shot, run-id from filename stem {cdf_name}. Please use <runid>'
      raise ValueError(mes) from exc

    run = 0
    try:
      run = int(''.join(part))
    except (ValueError, TypeError):
      print(' Error getting run-number:', part, 'shot:', shot,
            'code:', rid, 'runid:', runid)
    return run, shot, rid


def print_progress(tot, at):
    if not sys.stdout.isatty():
      return
    bar = 50
    progr = int(at*bar/tot)
    comp = str(int(at*100/tot)) + '%'
    if at==tot:
      sys.stdout.write("\033[K")
    else:
      print('[', chr(9618)*progr,' ',comp, '.'*(bar-progr), '] ',
         'time:' + str(at)+'/'+str(tot), sep='', end='\r', flush=True)


import cProfile
import functools
import pstats

def profile_it(funct):
    @functools.wraps(funct)
    def func(*args, **kwargs):
        prof = cProfile.Profile()
        prof.enable()
        try:
            res = funct(*args, **kwargs)
        finally:
            prof.disable()
            with open('prof_out.txt', 'w') as pfile:
                pst = pstats.Stats(prof, stream=pfile)
                pst.sort_stats('cumulative').print_stats(50)
        return res
    return func


def get_args(argv=None):
    desc = f' (To run all CDF files in <dir_path> and its subdirs: -i all -d <dir_path>)'
    parser = argparse.ArgumentParser(description='Convert TRANSP output CDF to IMAS IDS' + desc)
    parser.add_argument('-i', '--input-cdf', type=str,
        help='CDF TRANSP output file, with or without .CDF')
    parser.add_argument('-d', '--input-dir', type=str,
        help='input directory (file can be in its subdir), if file not in pwd')
    parser.add_argument('-b', '--backend', type=str, default='MDS',
        choices=['MDS','MDSPLUS', 'HDF', 'HDF5'], help='IMAS backend')
    parser.add_argument('-r', '--run', type=int,
        help='run number of o/p IDS (overrides default converted number: eg: Z10 -> 2610)')
    parser.add_argument('-new', '--new-db', help='create new database', action='store_true')
    parser.add_argument('-o', '--occur', type=int, help='occurrence')
    parser.add_argument('-u', '--user', type=str, help='username or path of database')
    parser.add_argument('-m', '--machine', type=str, help='machine. eg: ITER')
    parser.add_argument('--ip-direction', choices=['auto','ccw','cw'], default='auto',
        help='plasma current direction viewed from above; auto uses NLJCCW from TR.DAT')
    parser.add_argument('--bt-direction', choices=['auto','ccw','cw'], default='auto',
        help='toroidal field direction viewed from above; auto uses NLBCCW from TR.DAT')
    parser.add_argument('--dd-version', type=str, default='auto',
        help='target IMAS Data Dictionary version (e.g. 3.39.0 or 4.1.1); auto uses IMAS environment/default')
    parser.add_argument('--pf-passive-geometry', type=str, default=None,
        help='authoritative passive conductor geometry JSON; maps PS_n channels to physical loops')
    parser.add_argument('--pf-active-geometry', type=str, default=None,
        help='authoritative active PF coil/supply/circuit JSON; validates winding signs and channels')
    parser.add_argument('--ids', default=None,
        help='comma-separated IDSs to convert; default all (e.g. workflow,equilibrium,core_profiles)')
    parser.add_argument('--lenient-mapping', action='store_true',
        help='record a CDF variable that exists but cannot be mapped (bad shape, unreadable) as a skipped '
             'field with its reason in the --report-json journal instead of aborting; default is strict')
    parser.add_argument('--report-json', default=None,
        help='write per-IDS conversion results and successful IMAS puts as JSON (including on failure)')
    parser.add_argument('-mom', '--moments', type=str, choices=['yes', 'no'], default='no',
         help='convert moments to equilibrium, default is no')
    parser.add_argument('-zones', '--nzones', type=int, default=0,
         help='num poloidal grid points per flux surface = max(65,min(251,(2*nzones+1))). Default is NR')
    parser.add_argument('-ecnml', '--ec-nml-launch', type=str, choices=['yes', 'no'], default='no',
         help='get EC launch position from namelist (for testing)')
    # -h is added by default for help
    parser.add_argument('--debug', help='extra debug prints', action='store_true')
    parser.add_argument('--debug-max', help='All possible debug output including plots', action='store_true')
    parser.add_argument('-v', '--verbose', action='count', default=0, help='increase output verbosity')
    parser.add_argument("-l", "--log", choices=['debug','info','warn','error','critical'],
        help="Provide logging level")

    args = parser.parse_args(argv[1:] if argv is not None else None)
    mapping_status.set_strict(not args.lenient_mapping)

    if not args.input_cdf:
        print("** Error: Failed to get command line arguments**\n")
        parser.print_help(sys.stderr)
        print("Example: \npython convert.py -i 147634G14 [ -r 1 -b MDS -new -l debug]\n")
        sys.exit(1)

    #get runid, DAT filename from input cdf filename
    runno, shot, rid = find_runid_shot_code(args.input_cdf)

    #ncfile and nml file
    nml = str(shot)+rid+'TR.DAT'
    extn = Path(args.input_cdf).suffix
    ncfile = args.input_cdf if 'CDF' in extn.upper() else args.input_cdf + '.CDF'
    fdir = Path(ncfile).parent
    cwd = Path().absolute()
    #if cdf file is not in cwd, then use DAT from dir of cdf file
    if fdir != cwd and fdir:
        nml = os.path.join(fdir, nml)

    #use input_dir if present
    if args.input_dir:
        fname = Path(ncfile).name
        loc = Path(args.input_dir)
        res = find_file(fname, loc)
        if len(res)==0:
            raise Exception(f'{fname} not found in {loc}')
        ncfile = res[0]
        loc = Path(ncfile).parent
        nml = Path(nml).name
        res = find_file(nml, loc)
        nml = ''
        if len(res)>0:
            nml = res[0]

    if not check_file_exists(ncfile):
        raise Exception(f'NC file {ncfile} not found.')
    if not check_file_exists(nml):
        print(f'*** WARNING: namelist {nml} not found.')

    run = 0
    if args.run:
        run = int(args.run)
    else:
        run = runno

    backend = 'HDF5'
    newdb = False
    occur = 0
    log = 'debug'  #for logger
    machine = None
    moments = False
    nzones = 0
    ec_nml = False
    debug = 0
    if args.backend:
        backend = args.backend
    if args.occur:
        occur = int(args.occur)
    if args.log:
        log = args.log
    if args.new_db:
        newdb = True
    if args.machine:
        machine = args.machine
    if args.moments and args.moments=='yes':
        moments = True
    if args.nzones:
        nzones = int(args.nzones)
    if args.ec_nml_launch and args.ec_nml_launch=='yes':
        ec_nml = True
    if args.debug:
        debug = 1
    if args.debug_max:
        debug = 2

    user = args.user if args.user else None
    #new_db = True if occur == 0 else False

    return backend, run, occur, user, machine, newdb, ncfile, nml, log, moments, nzones, debug, ec_nml, args.ip_direction, args.bt_direction, args.dd_version, args.pf_passive_geometry, args.pf_active_geometry, args.ids, args.report_json

