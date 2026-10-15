#!/usr/bin/python
"""Verify neoclassical particle/heat theory transport independently of neutrals."""
import argparse
import json
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from theory_transport_mapping import (discover_neoclassical_channels,
    validate_profile, divergence_to_surface_flux, build_theory_plans)
from transport_test_support import write_core_transport


def verify_case(nc, label):
    shape = np.asarray(nc.variables['X'][:]).shape
    volume = np.asarray(nc.variables['DVOL'][:], dtype=np.float64)
    surf = np.asarray(nc.variables['SURF'][:], dtype=np.float64)
    channels = discover_neoclassical_channels(nc.variables)
    result = {'case': str(label), 'channels': {}, 'family': 'neoclassical'}
    for kind, items in channels.items():
        result['channels'][kind] = []
        for species, name in items.items():
            div = validate_profile(nc, name, kind, shape)
            flux = divergence_to_surface_flux(div, volume, surf, kind)
            assert np.isfinite(flux).all(), (label, name)
            # Independent finite-volume closure, excluding zero area surfaces.
            enclosed = np.cumsum(div * volume, axis=1)
            valid = np.abs(surf) > 0
            err = np.max(np.abs(flux[valid] * surf[valid] / 1e4 - enclosed[valid])) if np.any(valid) else 0
            relative = err / max(1., float(np.max(np.abs(enclosed))))
            assert relative < 2e-6, (label, name, relative)
            result['channels'][kind].append(name)
    assert any(result['channels'].values()), f'{label}: expected neoclassical channels'
    return result


class TestNeoclassical(unittest.TestCase):
    def test_aggregate_not_fabricated(self):
        d = discover_neoclassical_channels(['GFLNC_I','GFLNC_X','GFLNC_D','QFLNCC_E'])
        self.assertEqual(d['particles']['D'], 'GFLNC_D')
        self.assertIn('I', d['particles'])
        self.assertIn('X', d['particles'])

    def test_divergence_not_direct_flux(self):
        div = np.array([[1.,2.]])
        vol = np.array([[10.,20.]])
        surf = np.array([[5.,10.]])
        np.testing.assert_allclose(divergence_to_surface_flux(div,vol,surf,'heat'),
                                   np.array([[20000.,50000.]]))
        with self.assertRaises(ValueError):
            divergence_to_surface_flux(div,vol[:,:1],surf,'heat')

    def test_real_writer_neoclassical(self):
        path = ROOT/'examples/data/d3d/147634T13/147634T13.CDF'
        if not path.is_file():
            self.skipTest('T13 fixture missing')
        with netcdf_file(path, mmap=False) as nc:
            ids = write_core_transport(nc)
            models = {m.identifier.name: m for m in ids.model}
            self.assertIn('neoclassical',models)
            self.assertIn('neoclassical_classical_heat',models)
            species = [n.name for n in models['neoclassical'].profiles_1d[0].ion]
            self.assertNotIn('I',species)
            self.assertNotIn('X',species)
            self.assertNotIn('TOK',species)
            im_flux=models['neoclassical'].profiles_1d[0].ion[0].particles.flux.values
            source=validate_profile(nc,'GFLNC_D','particles')
            expected=divergence_to_surface_flux(source, np.asarray(nc.variables['DVOL'][:]),
                                                np.asarray(nc.variables['SURF'][:]),'particles')
            np.testing.assert_allclose(im_flux,expected[0],rtol=1e-12)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--quick', action='store_true')
    p.add_argument('--data-root', type=Path, default=ROOT/'examples/data')
    p.add_argument('--report-dir', type=Path, default=ROOT/'build/neoclassical-transport')
    args=p.parse_args(argv)
    if not unittest.TextTestRunner(verbosity=1).run(
            unittest.defaultTestLoader.loadTestsFromTestCase(TestNeoclassical)).wasSuccessful():
        return 1
    if args.quick:
        return 0
    files=sorted(args.data_root.rglob('*.CDF'))
    if not files:
        p.error('No CDF files')
    results=[]
    errors=[]
    for path in files:
        try:
            with netcdf_file(path, mmap=False) as nc:
                result=verify_case(nc,path.relative_to(args.data_root))
            results.append(result)
            print('PASS',result['case'],'NC channels:',sum(map(len,result['channels'].values())))
        except Exception as exc:
            errors.append(f'{path}: {type(exc).__name__}: {exc}')
            print('FAIL',errors[-1])
    args.report_dir.mkdir(parents=True,exist_ok=True)
    (args.report_dir/'summary.json').write_text(json.dumps({'checked':len(files),'cases':results,'errors':errors},indent=2)+'\n')
    return int(bool(errors))

if __name__=='__main__':
    raise SystemExit(main())
