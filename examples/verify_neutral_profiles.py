#!/usr/bin/python
"""Verify TRANSP wall, volume and recombination neutral state profiles only."""
import argparse
import json
from pathlib import Path
import sys
import types
import unittest
import xml.etree.ElementTree as ET

import numpy as np
from scipy.io import netcdf_file

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from neutral_profiles_mapping import (discover_neutral_channels, validate_profile,
                                      attach_rotation_xml)
from transport_test_support import FakeConverter, load_writer


def verify_case(nc, label):
    shape = np.asarray(nc.variables['X'][:]).shape
    specs=discover_neutral_channels(nc.variables)
    result={'case':str(label),'neutral_populations':[], 'family':'neutral_profiles'}
    for spec in specs:
        dens=validate_profile(nc,spec['density'],'density',shape)
        assert np.all(dens >= 0), (label,spec['density'])
        for kind in ('temperature','rotation'):
            if spec[kind]:
                values=validate_profile(nc,spec[kind],kind,shape)
                if kind=='temperature':
                    assert np.all(values>=0), (label,spec[kind])
        result['neutral_populations'].append(spec)
    assert specs, f'{label}: no neutral state source'
    return result


class TestNeutrals(unittest.TestCase):
    def test_population_separation(self):
        names=['DN0VD','DN0WD','N0V0_D','T0VD','OM0VD']
        result=discover_neutral_channels(names)
        self.assertEqual({s['origin'] for s in result},
                         {'wall','volume','recombination'})

    def test_xml_imas_string_wrapper(self):
        class IDSString0D:
            value='<root><inputs /></root>'
        n=types.SimpleNamespace(label='D_VOLUME',rotation_source='OM0VD',
            rotation_sign=1.0,omega=np.ones((2,3)))
        root=ET.fromstring(attach_rotation_xml(IDSString0D(),[n],[0.,1.]))
        self.assertIsNotNone(root.find('inputs'))
        self.assertEqual(root.find('neutral_toroidal_rotation/profile').get('source_variable'), 'OM0VD')
        self.assertEqual(attach_rotation_xml(IDSString0D(), [], []), IDSString0D.value)

    def test_writer_neutral_si_units(self):
        path=ROOT/'examples/data/d3d/147634T13/147634T13.CDF'
        if not path.exists():
            self.skipTest('T13 fixture missing')
        with netcdf_file(path,mmap=False) as nc:
            writer_cls=load_writer('core_profiles')
            writer=writer_cls.__new__(writer_cls)
            writer.im=FakeConverter(nc)
            writer.logger=writer.im.logger
            neutrals=writer.get_neutral_data()
            vol=next(n for n in neutrals if n.origin=='volume')
            np.testing.assert_allclose(vol.dens,validate_profile(nc,'DN0VD','density')*1e6)
            np.testing.assert_allclose(vol.temp,validate_profile(nc,'T0VD','temperature'))
            np.testing.assert_allclose(vol.omega,-validate_profile(nc,'OM0VD','rotation'))


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--quick',action='store_true')
    p.add_argument('--data-root',type=Path,default=ROOT/'examples/data')
    p.add_argument('--report-dir',type=Path,default=ROOT/'build/neutral-profiles')
    args=p.parse_args(argv)
    if not unittest.TextTestRunner(verbosity=1).run(
            unittest.defaultTestLoader.loadTestsFromTestCase(TestNeutrals)).wasSuccessful():
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
            with netcdf_file(path,mmap=False) as nc:
                item=verify_case(nc,path.relative_to(args.data_root))
            results.append(item)
            print('PASS',item['case'],'neutral populations:',len(item['neutral_populations']))
        except Exception as exc:
            errors.append(f'{path}: {type(exc).__name__}: {exc}')
            print('FAIL',errors[-1])
    args.report_dir.mkdir(parents=True,exist_ok=True)
    (args.report_dir/'summary.json').write_text(json.dumps({'checked':len(files),'cases':results,'errors':errors},indent=2)+'\n')
    return int(bool(errors))

if __name__=='__main__':
    raise SystemExit(main())
