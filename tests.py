import os
import sys
import numpy as np
from matplotlib import pyplot as plt

import imas
import MDSplus
from plotter import *



class Tests:
    def __init__(self, pre:str='', save_plots:bool=True):
        self.save = save_plots
        self.pre = pre
        print('Running tests')

    def init_database(self, db, shot, run, backend=None, user=None):
        if backend is None:
            backend = imas.imasdef.MDSPLUS_BACKEND
        if user is None:
            user = os.getenv('USER')
        print('backend', backend)
        self.input = imas.DBEntry(backend,db,shot,run,user)

    def test_equilibrium(self):
        print('* testing equilibrium')
        eq = self.input.get('equilibrium')
        pre = self.pre+'_eq' if self.pre else 'eq'
        #ni for number of indices to plot for all times
        plot_equilibrium_ids(eq, ni=6, pre=pre, save=self.save)

    def test_core_profiles(self):
        print('* testing test_core_profiles')
        cp = self.input.get('core_profiles')
        pre = self.pre+'_cp' if self.pre else 'cp'
        plot_core_profiles(cp, pre=pre, save=self.save)


    def test_limiters(self):
        print('* testing limiters')
        wall = self.input.get('wall')
        pre = self.pre+'_wall' if self.pre else 'wall'
        plot_limiter_ids(wall, pre=pre, save=self.save)

    def test_thomson_scattering(self):
        print('*thomson_scattering')
        ts = self.input.get('thomson_scattering')
        pre = self.pre+'_thomson' if self.pre else 'thomson'
        plot_thomson_scattering_ids(ts, pre=pre, save=self.save)

    def test_PFActive(self):
        pf = self.input.get('pf_active')
        pre = self.pre+'_pf' if self.pre else 'pf'
        plot_PFActive(pf, pre=pre, save=self.save)

    def test_CHERS(self):
        cx = self.input.get('charge_exchange')
        pre = self.pre+'_chers' if self.pre else 'chers'
        plot_CHERS(cx, pre=pre, save=self.save)


    def test_NBI(self):
        nb = self.input.get('nbi')
        pre = self.pre+'_nbi' if self.pre else 'nbi'
        plot_NBI(nb, pre=pre, save=self.save)

    def do_test_eq(self, shot, run, db):
        backend = imas.imasdef.MDSPLUS_BACKEND
        self.init_database(db,shot,run,backend)
        self.input.open()
        self.test_equilibrium()
        self.input.close()

    def do_test(self, shot, run, db):
        #shot,run,db = 141711,2502,'NSTX'
        # 10000,99,'test'
        backend = None #imas.imasdef.HDF5_BACKEND
        self.init_database(db,shot,run,backend)
        self.input.open()
        self.test_NBI()
        self.test_CHERS()
        self.test_PFActive()
        self.test_thomson_scattering()
        self.test_core_profiles()
        self.test_limiters()
        self.test_equilibrium()

        self.input.close()
    #end of tests


if __name__ == '__main__':
    db = ''
    shot = ''
    run = ''
    if len(sys.argv) > 1 and sys.argv[1]!= '-':
        db = sys.argv[1]
    if len(sys.argv) > 2 and sys.argv[2]!= '-':
        shot = int(sys.argv[2])
    if len(sys.argv) > 3 and sys.argv[3]!= '-':
        run = int(sys.argv[3])

    if not (db or shot or run):
        #shot,run,db = 163303,8702,'D3Dr' #r'
        shot,run,db = 100001, 1, 'NSTX'
    print(f'testing stored database {db} {shot} {run}')

    pre = db+str(run) # +'_h5'
    save_plots = True
    tt = Tests(pre, save_plots)
    tt.do_test_eq(shot, run, db)
    #tt.do_test(shot,run,db)
