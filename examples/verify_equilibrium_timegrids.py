#!/usr/bin/python
"""Verify equilibrium TIME/TIME3 alignment against every bundled CDF.

Exercises the actual writer's LCFS-moment methods and nc2imas time-coordinate
methods with small source-backed adapters; no IMAS or netCDF4 install needed.
Also checks the writer's equilibrium.time/time_slice and boundary-shaping
contracts.  The test does not replace a real IMAS DD4 write/readback.
"""

import argparse
import ast
import logging
from pathlib import Path
import re

import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def extract_method(filename, classname, methodname):
    """Run the source method, without importing modules that require IMAS."""
    path = ROOT / filename
    mod = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    cls = next(n for n in mod.body if isinstance(n, ast.ClassDef)
               and n.name == classname)
    func = next(n for n in cls.body if isinstance(n, ast.FunctionDef)
                and n.name == methodname)
    isolated = ast.Module(body=[func], type_ignores=[])
    namespace = {'np': np}
    exec(compile(ast.fix_missing_locations(isolated), str(path), 'exec'), namespace)
    return namespace[methodname]


class CaseReader:
    """Small SciPy-backed stand-in for an nc2imas CDF connection."""

    def __init__(self, nc):
        self.nc = nc
        self.timesc = self.get_ncvar('TIME')
        self.timep = self.get_ncvar('TIME3')
        self.ntimes = len(self.timesc)
        self.ntimep = len(self.timep)
        self.pi = np.pi
        self.theta_raw_from_target = -1.0
        self.debug_max = False
        self.logger = logging.getLogger('verify-equilibrium-timegrids')

    def get_ncvar(self, name, scale=1.0, ind=-1, interp=False):
        if name not in self.nc.variables:
            return np.array([])
        var = self.nc.variables[name]
        if ind >= 0:
            values = np.asarray(var[ind], dtype=np.float64).copy()
        else:
            values = np.asarray(var[:], dtype=np.float64).copy()
        if interp:
            values = np.interp(self.timep, self.timesc, values)
        return values * scale

    def get_var_names_re(self, expr, nd1=0, nd2=0, length=False, **kwargs):
        require(expr in ('RMCB?', 'RMC?'), f'unexpected moment pattern {expr}')
        regex = r'^RMCB[0-9]+$' if expr == 'RMCB?' else r'^RMC[0-9]+$'
        names = [name for name in self.nc.variables if re.fullmatch(regex, name)]
        return len(names) if length else names

    def is_valid(self, *values, ind=0, name=''):
        for value in values:
            if value is None:
                return False
            arr = np.asarray(value)
            if arr.ndim == 0 or len(arr) <= ind or arr.size == 0:
                return False
        return True


# Bind the real implementation under test, not an independently reimplemented
# reader or interpolator.  Only netCDF/IMAS allocation is mocked.
CaseReader.interp1d_time = extract_method('nc2imas.py', 'nc2imas', 'interp1d_time')
CaseReader.get_ncvar_on_time = extract_method('nc2imas.py', 'nc2imas', 'get_ncvar_on_time')


class EquilibriumAdapter:
    def __init__(self, reader):
        self.im = reader
        self.logger = reader.logger
        self.rmcb = []
        self.rmsb = []
        self.ymcb = []
        self.ymsb = []
        self.nmoments = 0
        self.dthetab = 2 * np.pi / 128
        self.dtheta = None
        self.debug_max = False


EquilibriumAdapter.get_moments_bdry = extract_method(
    'imas_equilibrium.py', 'imas_equilibrium', 'get_moments_bdry')
EquilibriumAdapter.convert_moments_bdry = extract_method(
    'imas_equilibrium.py', 'imas_equilibrium', 'convert_moments_bdry')


def writer_contract():
    path = ROOT / 'imas_equilibrium.py'
    tree = ast.parse(path.read_text(encoding='utf-8'))
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef)
               and n.name == 'imas_equilibrium')
    method = next(n for n in cls.body if isinstance(n, ast.FunctionDef)
                  and n.name == 'create_ids')
    lines = ast.unparse(method)
    require('eqm.time_slice.resize(ntimep)' in lines,
            'equilibrium.time_slice must use TIME3 length')
    require('eqm.time[:] = timep' in lines,
            'equilibrium.time must use TIME3 values')
    for name in ('elongb', 'triangb', 'triang_lob', 'triang_upb'):
        require(f'{name} = ' in lines and '[:, -1]' in lines,
                f'{name}: missing TIME3 boundary shaping')
    require('self.im.is_valid(psi, psi_axis, ind=ntimep - 1)' in lines,
            'psi_axis length validation must use ntimep, not scalar ntimes')
    require('self.convert_moments_bdry(nthetab, it)' in lines,
            'boundary moments must be converted for every TIME3 slice')
    require('if ntimep != ntimes' not in lines,
            'equilibrium may not gate interpolation on array length')
    require('get_ncvar_on_time' in ast.unparse(cls),
            'equilibrium writer must use dimension-aware time coordinates')


def synthetic():
    class Var:
        def __init__(self, dimension):
            self.dimensions = (dimension,)

    class Reader(CaseReader):
        def __init__(self):
            self.timesc = np.array([0., 1., 2.])
            self.timep = np.array([0., 0.5, 2.])
            self.nc = type('NC', (), {'variables': {'A': Var('TIME'),
                                                     'B': Var('TIME3')}})()
        def get_ncvar(self, name, scale=1.0, **kwargs):
            return scale * np.array([0., 10., 20.])

    reader = Reader()
    require(np.allclose(reader.get_ncvar_on_time('A'), [0., 5., 20.]),
            'equal-length, different-valued TIME/TIME3 need interpolation')
    require(np.allclose(reader.get_ncvar_on_time('B', target='TIME'), [0., 13.3333333333333, 20.]),
            'TIME3 -> TIME interpolation incorrect')
    require(np.allclose(reader.interp1d_time(reader.timep, reader.timesc,
                          np.array([0., 10., 20.])), [0., 5., 20.]),
            'interp1d_time must not skip equal-size, different-time coordinates')
    require(np.array_equal(reader.get_ncvar_on_time('B'), [0., 10., 20.]),
            'already TIME3-valued data must not be interpolated')
    try:
        reader.get_ncvar_on_time('A', 'OTHER')
    except ValueError:
        pass
    else:
        raise AssertionError('unrecognized time coordinate was accepted')
    require(reader.get_ncvar_on_time('Z').size == 0,
            'missing optional source history should return an empty array')


def verify_case(nc, label):
    reader = CaseReader(nc)
    nt = reader.ntimes
    nt3 = reader.ntimep
    require(nt > 1 and nt3 > 1, f'{label}: missing TIME or TIME3')
    require(np.all(np.diff(reader.timesc) > 0), f'{label}: TIME is not increasing')
    require(np.all(np.diff(reader.timep) > 0), f'{label}: TIME3 is not increasing')
    # Native coordinates matter; shapes alone are insufficient.
    for name in ('X', 'XB', 'Q', 'PLFLX', 'RMC00',
                 'ELONG', 'TRIANG', 'TRIANGU', 'TRIANGL'):
        if name in nc.variables:
            var = nc.variables[name]
            require(var.dimensions[0] == 'TIME3' and var.shape[0] == nt3,
                    f'{label}: {name} unexpectedly not on TIME3')
    for name in ('PSIRZ', 'B_FIELD', 'RGRID', 'ZGRID'):
        if name in nc.variables:
            var = nc.variables[name]
            require(var.dimensions[0] == 'TIME3' and var.shape[0] == nt3,
                    f'{label}: {name} unexpectedly not on TIME3')
    for name in ('PCURC', 'RAXIS', 'YAXIS', 'RMCB0'):
        if name in nc.variables:
            var = nc.variables[name]
            require(var.dimensions == ('TIME',) and var.shape == (nt,),
                    f'{label}: {name} unexpectedly not on TIME')
            mapped = reader.get_ncvar_on_time(name)
            require(mapped.shape == (nt3,),
                    f'{label}: {name} resampled to wrong number of slices')
            require(np.allclose(mapped, np.interp(reader.timep, reader.timesc,
                                                   reader.get_ncvar(name))),
                    f'{label}: {name} has wrong interpolated values')
    eqm = EquilibriumAdapter(reader)
    eqm.get_moments_bdry()
    if 'RMCB0' in nc.variables:
        require(np.shape(eqm.rmcb)[0] == nt3 and np.shape(eqm.ymcb)[0] == nt3,
                f'{label}: boundary coefficients not defined on TIME3')
        require(np.allclose(eqm.rmcb[:, 0],
                  reader.get_ncvar_on_time('RMCB0', scale=0.01)),
                f'{label}: boundary RMCB0 source/resampling incorrect')
        for index in sorted({0, nt3 // 2, nt3 - 1}):
            r, z = eqm.convert_moments_bdry(128, index)
            require(len(r) == len(z) == 128 and np.all(np.isfinite(r))
                    and np.all(np.isfinite(z)),
                    f'{label}: missing or nonfinite LCFS at profile slice {index}')
            require(r[0] == r[-1] and z[0] == z[-1],
                    f'{label}: LCFS outline is not closed at {index}')
        # The last boundary coefficient must correspond to the last profile
        # time even where that index is beyond the scalar TIME array.
        require(np.isclose(eqm.rmcb[-1, 0],
                  np.interp(reader.timep[-1], reader.timesc,
                            reader.get_ncvar('RMCB0') * 0.01)),
                f'{label}: final TIME3 boundary not defined')
    return {'TIME': nt, 'TIME3': nt3, 'moments': eqm.nmoments}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, default=HERE / 'data')
    parser.add_argument('--quick', action='store_true',
                        help='synthetic plus writer-contract tests only')
    args = parser.parse_args(argv)
    writer_contract()
    synthetic()
    if args.quick:
        print('equilibrium TIME/TIME3 regression PASS (source + synthetic)')
        return 0
    files = sorted(args.data_root.rglob('*.CDF'))
    require(files, f'No CDF examples found under {args.data_root}')
    for file in files:
        with netcdf_file(str(file), 'r', mmap=False) as nc:
            result = verify_case(nc, file.name)
        print(f'{file.stem}: TIME={result["TIME"]} TIME3={result["TIME3"]} '
              f'LCFS moments={result["moments"]}')
    print(f'equilibrium TIME/TIME3 regression PASS ({len(files)}/{len(files)} CDFs)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
