#!/usr/bin/python
"""Exercise plot_ids argument parsing, figure creation and file output without IMAS."""

import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import numpy as np

SCRIPT = Path(__file__).with_name('plot_ids.py')
spec = importlib.util.spec_from_file_location('plot_ids_under_test', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

assert module.parse_args(['38530P94', 'ITER']).field == 'psirz'
assert module.parse_args(['38530P94', 'ITER', '--run', '1']).run == 1
assert module.parse_args(['38530P94', 'ITER', 'HDF5', 'q', '--time-index', '2']).time_index == 2
help_result = subprocess.run([sys.executable, str(SCRIPT), '--help'],
                             capture_output=True, text=True)
assert help_result.returncode == 0 and '--save' in help_result.stdout

r = np.linspace(2.0, 3.0, 6)
z = np.linspace(-0.6, 0.6, 8)
base = np.add.outer(r, z)
profile = SimpleNamespace(grid=SimpleNamespace(dim1=r, dim2=z), psi=base,
                          b_field_r=base*0.01, b_field_tor=base*0+2.0,
                          b_field_z=base*0.02)
p1 = SimpleNamespace(rho_tor_norm=np.linspace(0, 1, 6), q=np.arange(6)+1,
                     psi=np.linspace(0, 1, 6), pressure=np.ones(6), j_tor=np.ones(6))
ts = SimpleNamespace(profiles_2d=[profile], profiles_1d=p1,
                     boundary=SimpleNamespace(outline=SimpleNamespace(r=r, z=np.zeros(6))))
eq = SimpleNamespace(time_slice=[ts], time=np.array([0.3]))
with tempfile.TemporaryDirectory() as temp:
    for field in module.FIELDS:
        fig = module.make_figure(eq, field, 0)
        path = Path(temp) / f'{field}.png'
        fig.savefig(path)
        assert path.stat().st_size > 1000, field
        from matplotlib import pyplot as plt
        plt.close(fig)
print('PLOT IDS: PASS — all six plot types and CLI help (mock equilibrium)')
