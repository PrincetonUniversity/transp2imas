#!/usr/bin/python
"""Verify core_profiles axis-temperature peaking on TRANSP TIME/TIME3 grids.

All CDF cases are checked when run directly; ``--quick`` runs a synthetic
mismatched-time fixture and checks that the writer interpolates TE0.
This is a source-side verification, not an IMAS IDS write/readback test.
"""

import argparse
import ast
from pathlib import Path

import numpy as np
from scipy.io import netcdf_file

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def verify_writer_contract():
    """Prevent reintroduction of the missing TE0 TIME -> TIME3 conversion."""
    tree = ast.parse((ROOT / 'imas_core_profiles.py').read_text(encoding='utf-8'))
    converter_class = next((node for node in tree.body
                            if isinstance(node, ast.ClassDef)
                            and node.name == 'imas_core_profiles'), None)
    require(converter_class is not None, 'imas_core_profiles class not found')
    method = next((node for node in converter_class.body
                   if isinstance(node, ast.FunctionDef) and node.name == 'create_ids'), None)
    require(method is not None, 'create_ids method not found')
    found = False
    for node in ast.walk(method):
        if not isinstance(node, ast.Assign) or not any(
                isinstance(target, ast.Name) and target.id == 'te0'
                for target in node.targets):
            continue
        call = node.value
        if (not isinstance(call, ast.Call)
                or not isinstance(call.func, ast.Attribute)
                or call.func.attr != 'get_ncvar'
                or not call.args
                or not isinstance(call.args[0], ast.Constant)
                or call.args[0].value != 'TE0'):
            continue
        found = any(kw.arg == 'interp' and isinstance(kw.value, ast.Constant)
                    and kw.value.value is True for kw in call.keywords)
    require(found, 'core_profiles TE0 must be read with interp=True onto TIME3')


def verify_arrays(time, timep, te0, ti0, pvol, te, ti, zone_vol):
    """Check both temperature peaking histories using TRANSP's native units.

    TE0/TI0/PVOL are stored on TIME, while TE/TI/DVOL are on TIME3.
    DVOL and PVOL are both in cm^3 here, so their units cancel.
    """
    t = np.asarray(time, dtype=np.float64)
    tp = np.asarray(timep, dtype=np.float64)
    volumes = np.asarray(zone_vol, dtype=np.float64)
    require(t.ndim == tp.ndim == 1 and t.size >= 2 and tp.size >= 2,
            'TIME/TIME3 must be one-dimensional with at least two points')
    require(np.all(np.isfinite(t)) and np.all(np.diff(t) > 0), 'TIME invalid')
    require(np.all(np.isfinite(tp)) and np.all(np.diff(tp) > 0), 'TIME3 invalid')
    require(volumes.ndim == 2 and volumes.shape[0] == len(tp),
            'DVOL is not on TIME3')
    total = np.asarray(pvol, dtype=np.float64)
    require(total.shape == t.shape, 'PVOL is not on TIME')
    interp_volume = np.interp(tp, t, total)
    require(np.all(interp_volume > 0), 'PVOL contains nonpositive values')
    result = {'TIME': int(len(t)), 'TIME3': int(len(tp))}
    for symbol, original_axis, profile in (('TE', te0, te), ('TI', ti0, ti)):
        a = np.asarray(original_axis, dtype=np.float64)
        p = np.asarray(profile, dtype=np.float64)
        require(a.shape == t.shape, f'{symbol}0 is not on TIME')
        require(p.shape == volumes.shape, f'{symbol} does not match DVOL/TIME3')
        a_on_profile_grid = np.interp(tp, t, a)
        average = np.sum(p * volumes, axis=1) / interp_volume
        require(np.all(np.isfinite(average)) and np.all(average > 0),
                f'{symbol} volume average is invalid')
        peaking = a_on_profile_grid / average
        require(peaking.shape == tp.shape and np.all(np.isfinite(peaking)),
                f'{symbol} peaking cannot be written on TIME3')
        result[f'{symbol}_peaking_min'] = float(peaking.min())
        result[f'{symbol}_peaking_max'] = float(peaking.max())
    return result


def verify_case(nc, label='CDF'):
    required = ('TIME', 'TIME3', 'TE0', 'TI0', 'PVOL', 'TE', 'TI', 'DVOL')
    missing = [key for key in required if key not in nc.variables]
    require(not missing, f'{label}: core_profiles test missing variables {missing}')
    vars_ = [np.asarray(nc.variables[key][:]) for key in required]
    try:
        return verify_arrays(*vars_)
    except (AssertionError, ValueError, TypeError) as exc:
        raise AssertionError(f'{label}: {exc}') from exc


def verify_synthetic():
    t = np.array([0.0, 1.0, 2.0])
    tp = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
    te0 = np.array([200.0, 300.0, 400.0])
    ti0 = np.array([100.0, 150.0, 200.0])
    volume = np.array([20.0, 20.0, 20.0])
    vzone = np.full((5, 2), 10.0)
    te = np.full((5, 2), 100.0)
    ti = np.full((5, 2), 50.0)
    result = verify_arrays(t, tp, te0, ti0, volume, te, ti, vzone)
    require(result['TIME3'] == 5, 'synthetic interpolation returned wrong length')
    require(np.isclose(result['TE_peaking_min'], 2.0)
            and np.isclose(result['TE_peaking_max'], 4.0),
            'synthetic electron peaking is incorrect')
    require(np.isclose(result['TI_peaking_min'], 2.0)
            and np.isclose(result['TI_peaking_max'], 4.0),
            'synthetic ion peaking is incorrect')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, default=HERE / 'data')
    parser.add_argument('--quick', action='store_true',
                        help='test only source contract and synthetic mismatch')
    args = parser.parse_args(argv)
    verify_writer_contract()
    verify_synthetic()
    if args.quick:
        print('core_profiles TIME/TIME3 regression: PASS (synthetic + writer contract)')
        return 0
    cdfs = sorted(args.data_root.rglob('*.CDF'))
    require(cdfs, f'No CDF examples found under {args.data_root}')
    for cdf in cdfs:
        with netcdf_file(str(cdf), 'r', mmap=False) as nc:
            result = verify_case(nc, cdf.name)
        print(f'{cdf.stem}: TIME={result["TIME"]} TIME3={result["TIME3"]}; '
              f'Te peak={result["TE_peaking_min"]:.3f}..{result["TE_peaking_max"]:.3f}; '
              f'Ti peak={result["TI_peaking_min"]:.3f}..{result["TI_peaking_max"]:.3f}')
    print(f'core_profiles TIME/TIME3 regression: PASS ({len(cdfs)}/{len(cdfs)} CDFs)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
