import os
import sys
from pathlib import Path
import math
import numpy as np
import inspect
import argparse


def check_file_exists(name:str)->bool:
    path_ = Path(name).absolute()
    if path_.is_symlink():
        path_ = os.readlink(name)
    return os.path.isfile(path_)

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
    except:
      mes = f'Error getting shot, run-id from filename stem {cdf_name}. Please use <runid>'
      raise Exception(mes)

    run = 0
    try:
      run = int(''.join(part))
    except:
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
    parser.add_argument('-new', '--new-db', help='create new database', action="store_true")
    parser.add_argument('-o', '--occur', type=int, help='occurrence')
    parser.add_argument('-u', '--user', type=str, help='username or path of database')
    parser.add_argument('-m', '--machine', type=str, help='machine. eg: ITER')
    # -h is added by default for help
    parser.add_argument('-v', '--verbose', action='count', default=0, help='increase output verbosity')
    parser.add_argument("-l", "--log", choices=['debug','info','warn','error','critical'],
        help="Provide logging level")

    args = parser.parse_args()

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

    if args.backend:
        backend = args.backend
    if args.new_db:
        new_db = args.new_db
    if args.occur:
        occur = int(args.occur)
    if args.log:
        log = args.log
    if args.new_db:
        newdb = True
    if args.machine:
        machine = args.machine

    user = os.path.expanduser('~') +'/public/imasdb' if \
            (backend=='MDS' or backend=='MDSPLUS') else os.getenv('USER')
    #new_db = True if occur == 0 else False

    return backend, run, occur, user, machine, newdb, ncfile, nml, log

