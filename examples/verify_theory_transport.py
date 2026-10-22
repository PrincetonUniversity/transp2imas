#!/usr/bin/python
"""Verify theory transport outputs (MMM07/GLF23/TGLF) on EVERY CDF.

MMM (Weiland, DRIBM/DRBM, ETG, MTM) produces predicted anomalous
transport coefficients, NOT measured experimental diagnostics.  Check typed
core_transport fields, SI units, native TIME3/XB vs TIME3/X grids, electron
and ion component closure, and preservation of remaining reduced-model output.
"""

import argparse
import json
from pathlib import Path
import sys
import types
import unittest
from unittest import mock
import xml.etree.ElementTree as ET

import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from theory_transport_mapping import (
    MODEL_SPECS, GLF23_MODEL_SPEC, THEORY_FAMILIES, build_theory_plans,
    build_mmm_plan, build_glf23_plan, check_mmm_closure, mmm_model_xml,
    theory_model_xml, mmm_summary,
)


class FakeVar:
    def __init__(self, data, dimensions, units):
        self.data = np.asarray(data, dtype=np.float32)
        self.dimensions = dimensions
        self.units = units

    def __getitem__(self, item):
        return self.data[item]


class FakeCDF:
    def __init__(self, corrupt=False):
        t, n = 2, 3
        x = np.tile([0.15, 0.45, 0.8], (t, 1))
        xb = np.tile([0.2, 0.5, 1.0], (t, 1))
        self.variables = {
            'TIME3': FakeVar([0.1, 0.2], ('TIME3',), 'S'),
            'X': FakeVar(x, ('TIME3', 'X'), '-'),
            'XB': FakeVar(xb, ('TIME3', 'XB'), '-'),
        }
        for name, data in {
            'XKEMMM07': [[9., 12., 15.], [9., 12., 15.]],
            'XKEW19': [[1., 2., 3.], [1., 2., 3.]],
            'XKEDRBM': [[2., 2., 2.], [2., 2., 2.]],
            'XKEETG': [[4., 5., 6.], [4., 5., 6.]],
            'XKEMTM': [[2., 3., 4.], [2., 3., 4.]],
            'XKIMMM07': [[9., 9., 9.], [9., 9., 9.]],
            'XKIW19': [[8., 8., 8.], [8., 8., 8.]],
            'XKIDRBM': [[1., 1., 1.], [1., 1., 1.]],
            'XKDMMM07': [[1e4, 2e4, 3e4], [1e4, 2e4, 3e4]],
            'XPTMMM07': [[100, 150, 200], [100, 150, 200]],
            'VPHIMMM07': [[10, 20, -5], [10, 20, -5]],
            'XPPMMM07': [[1, 2, 3], [1, 2, 3]],
            'XKZMMM07': [[4, 5, 6], [4, 5, 6]],
        }.items():
            unit = 'CM/SEC' if name == 'VPHIMMM07' else 'CM**2/SEC'
            self.variables[name] = FakeVar(data, ('TIME3', 'XB'), unit)
        for name, unit in [('GAMMAMTM','1/SEC'),('OMEGAMTM','RAD/SEC'),
                           ('KYRSMTM','-'),('XDBMTM','-'),
                           ('GAMMAETG','1/SEC'),('OMEGAETG','RAD/SEC'),
                           ('KYRSETG','-'),('KYRSDRBM','-')]:
            self.variables[name] = FakeVar(np.ones((t, n)), ('TIME3', 'X'), unit)
        if corrupt:
            self.variables['GAMMAMTM'].dimensions = ('TIME', 'X')


class TestMMMContracts(unittest.TestCase):
    def test_source_decomposition_and_units(self):
        plan = build_mmm_plan(FakeCDF())
        self.assertTrue(plan['active'])
        for item in plan['closure'].values():
            self.assertEqual(item['status'], 'PASS')
            self.assertEqual(item['max_abs_difference_cm2_s'], 0)
        self.assertEqual([spec['name'] for spec in MODEL_SPECS],
                         ['mmm07_total', 'mmm07_weiland', 'mmm07_dribm',
                          'mmm07_etg', 'mmm07_mtm'])
        self.assertEqual(len(set(spec['index'] for spec in MODEL_SPECS)), 5)
        self.assertTrue(all(spec['index'] < 0 for spec in MODEL_SPECS))

    def test_input_absent_or_all_zero_does_not_fabricate_mmm(self):
        cdf = FakeCDF()
        for name in list(cdf.variables):
            if name not in ('TIME3', 'X', 'XB'):
                cdf.variables[name].data[:] = 0
        self.assertFalse(build_mmm_plan(cdf)['active'])
        for name in list(cdf.variables):
            if name not in ('TIME3', 'X', 'XB'):
                del cdf.variables[name]
        self.assertFalse(build_mmm_plan(cdf)['active'])

    def test_source_time_grid_and_units_rejected_on_mismatch(self):
        with self.assertRaisesRegex(ValueError, 'expected'):
            build_mmm_plan(FakeCDF(corrupt=True))
        cdf = FakeCDF()
        cdf.variables['XKEMTM'].units = 'EV'
        with self.assertRaisesRegex(ValueError, 'unexpected units'):
            build_mmm_plan(cdf)
        cdf = FakeCDF()
        cdf.variables['XKEMMM07'].data[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, 'non-finite'):
            build_mmm_plan(cdf)

    def test_failed_closure_is_reported_not_hidden(self):
        plan = build_mmm_plan(FakeCDF())
        plan['fields']['XKEMTM'][0, 0] += 1
        closure = check_mmm_closure(plan['fields'])
        self.assertEqual(closure['electron_heat']['status'], 'FAIL')
        self.assertEqual(closure['ion_heat']['status'], 'PASS')

    def test_reduced_model_xml_references_native_imas_coordinates(self):
        plan = build_mmm_plan(FakeCDF())
        mtm = ET.fromstring(mmm_model_xml(plan, MODEL_SPECS[-1]))
        self.assertEqual(mtm.get('classification'), 'reduced_anomalous_transport_model_output')
        self.assertEqual(mtm.get('parent_model'), 'mmm07_total')
        sp = mtm.find('supplementary_model_outputs')
        self.assertIsNotNone(sp)
        arrays = {el.get('source'): el for el in sp.findall('quantity')}
        self.assertEqual(arrays['GAMMAMTM'].get('source_dimensions'), 'TIME3,X')
        self.assertEqual(arrays['GAMMAMTM'].get('radial_coordinate'),
                         'grid_flux/rho_tor_norm')
        self.assertEqual(arrays['XDBMTM'].get('units'), '1')
        self.assertEqual(np.fromstring(arrays['OMEGAMTM'].text, sep=' ').shape, (6,))
        self.assertFalse(sp.findall('coordinate'), 'XML must not duplicate coordinate arrays')
        refs = {ref.get('source'): ref for ref in sp.findall('coordinate_reference')}
        self.assertEqual(set(refs), {'TIME3', 'X'})
        self.assertEqual(refs['TIME3'].get('ids_path'), 'core_transport/time')
        self.assertIn('grid_flux/rho_tor_norm', refs['X'].get('ids_path'))
        self.assertIsNone(refs['X'].text)
        total = ET.fromstring(mmm_model_xml(plan, MODEL_SPECS[0]))
        self.assertEqual(total.find('typed_coefficients/coefficient').get('source'), 'XKDMMM07')
        total_refs = {p.get('source'): p for p in
                      total.findall('supplementary_model_outputs/coordinate_reference')}
        self.assertEqual(set(total_refs), {'TIME3', 'XB'})
        self.assertIn('grid_d/rho_tor_norm', total_refs['XB'].get('ids_path'))
        self.assertEqual(total.find('source_component_closure/check').get('status'), 'PASS')

    def test_no_implicit_mixing_with_standalone_models(self):
        source = (ROOT / 'theory_transport_mapping.py').read_text()
        for name in ('XKECDBM', 'XKEPALEO', 'XKINC'):
            self.assertNotIn("'" + name + "': (", source)


def time_dependent_cdfs(root):
    """TRANSP result CDFs under ``root``; skips files without the radial/time grids.

    Companion files such as ``<run>PH.CDF`` (no TIME3/X/XB variables) are not
    time-dependent TRANSP output and cannot be passed to the transport writers.
    """
    found = []
    for path in sorted(Path(root).rglob('*.CDF')):
        with netcdf_file(path, mmap=False) as nc:
            if {'TIME3', 'X', 'XB'} <= set(nc.variables):
                found.append(path)
    return found


class TestMMMWriter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with mock.patch.dict(sys.modules, {
            'imas': types.ModuleType('imas'),
            'netCDF4': types.ModuleType('netCDF4'),
        }):
            sys.modules['netCDF4'].Dataset = object
            from imas_core_transport import imas_core_transport
            from transport_test_support import FakeConverter
        cls.Writer = imas_core_transport
        cls.FakeConverter = FakeConverter

    def _write(self, nc):
        im = self.FakeConverter(nc)
        writer = self.Writer.__new__(self.Writer)
        writer.im = im
        writer.logger = im.logger
        writer.time = im.get_ncvar('TIME3')
        writer.x = im.get_ncvar('X')
        writer.xb = im.get_ncvar('XB')
        writer.dvol = im.get_ncvar('DVOL', im.vconv)
        writer.surf = im.get_ncvar('SURF', im.aconv)
        writer.nt = len(writer.time)
        writer.nx = writer.xb.shape[1]
        writer.create_ids()
        ids = im.written
        self.last_core_transport = ids
        return {ids.model[i].identifier.name: ids.model[i]
                for i in range(len(ids.model))}

    def test_real_t16_coefficients_are_independent_models(self):
        case = ROOT / 'examples/data/d3d/147634T16/147634T16.CDF'
        if not case.is_file():
            self.skipTest('T16 MMM source fixture missing')
        with netcdf_file(case, mmap=False) as nc:
            models = self._write(nc)
            for spec in MODEL_SPECS:
                self.assertIn(spec['name'], models)
                self.assertEqual(models[spec['name']].identifier.index, spec['index'])
                self.assertEqual(len(models[spec['name']].profiles_1d), 9)
            source = lambda name: np.asarray(nc.variables[name][:], dtype=np.float64)
            np.testing.assert_allclose(self.last_core_transport.time, source('TIME3'))
            np.testing.assert_allclose(
                models['mmm07_total'].profiles_1d[3].electrons.energy.d.values,
                source('XKEMMM07')[3] * 1e-4)
            np.testing.assert_allclose(
                models['mmm07_etg'].profiles_1d[3].electrons.energy.d.values,
                source('XKEETG')[3] * 1e-4)
            np.testing.assert_allclose(
                models['mmm07_mtm'].profiles_1d[3].electrons.energy.d.values,
                source('XKEMTM')[3] * 1e-4)
            np.testing.assert_allclose(
                models['mmm07_total'].profiles_1d[3].momentum_phi.d.values,
                source('XPTMMM07')[3] * 1e-4)
            np.testing.assert_allclose(
                models['mmm07_total'].profiles_1d[3].momentum_phi.v.values,
                source('VPHIMMM07')[3] * 1e-2)
            self.assertNotIn('XKECDBM', str(models['mmm07_total'].comment))
            for mname in ('mmm07_mtm', 'mmm07_etg', 'mmm07_dribm'):
                component = models[mname]
                for it, profile in enumerate(component.profiles_1d._items):
                    np.testing.assert_allclose(profile.grid_d.rho_tor_norm.values,
                                               source('XB')[it])
                    np.testing.assert_allclose(profile.grid_flux.rho_tor_norm.values,
                                               source('X')[it])
                    self.assertAlmostEqual(profile.time, float(source('TIME3')[it]))
                    self.assertNotIn('volume', profile.grid_flux.__dict__)
            # Non-reduced components and total retain the original XB flux grid.
            np.testing.assert_allclose(
                models['mmm07_total'].profiles_1d[3].grid_flux.rho_tor_norm.values,
                source('XB')[3])
            mtm = ET.fromstring(models['mmm07_mtm'].code.parameters)
            self.assertIn('GAMMAMTM',
                          [x.get('source') for x in mtm.findall('.//quantity')])

    def test_mmm_absent_from_other_real_cdfs(self):
        data_root = ROOT / 'examples/data'
        files = time_dependent_cdfs(data_root)
        if not files:
            self.skipTest('no source fixtures')
        for path in files:
            if path.stem == '147634T16':
                continue
            with self.subTest(case=path.name), netcdf_file(path,mmap=False) as nc:
                models = self._write(nc)
                self.assertTrue(all(not name.startswith('mmm07') for name in models))



class TestGenericTheory(unittest.TestCase):
    def test_registry_and_tglf_not_fabricated(self):
        self.assertEqual(THEORY_FAMILIES, ('neoclassical','mmm07','glf23','tglf'))
        plans=build_theory_plans(FakeCDF())
        self.assertTrue(plans['mmm07']['active'])
        self.assertFalse(plans['glf23']['active'])
        self.assertFalse(plans['tglf']['active'])
        # A standalone scalar instability-location channel does not establish
        # TGLF diffusivity/model execution.
        cdf=FakeCDF()
        cdf.variables['TGLF_DR']=FakeVar([0.2, 0.3], ('TIME',), '-')
        self.assertFalse(build_theory_plans(cdf)['tglf']['active'])

    def test_glf23_real_iter_and_writer_si(self):
        path=ROOT/'examples/data/iter/38530P94.CDF'
        if not path.exists():
            self.skipTest('ITER GLF23 fixture missing')
        with netcdf_file(path,mmap=False) as nc:
            plans=build_theory_plans(nc)
            self.assertTrue(plans['glf23']['active'])
            self.assertFalse(plans['mmm07']['active'])
            self.assertEqual(set(GLF23_MODEL_SPEC['coefficients']),
                             {'XKEGLF23','XKIGLF23','ETPHIGLF'})
            models=TestMMMWriter()._write(nc)
            self.assertIn('glf23',models)
            self.assertNotIn('mmm07_total',models)
            profile=models['glf23'].profiles_1d[0]
            np.testing.assert_allclose(profile.electrons.energy.d.values,
                np.asarray(nc.variables['XKEGLF23'][:])[0]*1e-4)
            np.testing.assert_allclose(profile.total_ion_energy.d.values,
                np.asarray(nc.variables['XKIGLF23'][:])[0]*1e-4)
            np.testing.assert_allclose(profile.momentum_phi.d.values,
                np.asarray(nc.variables['ETPHIGLF'][:])[0]*1e-4)
            xml=ET.fromstring(models['glf23'].code.parameters)
            self.assertEqual(xml.get('family'),'glf23')
            self.assertEqual(xml.find('supplementary_quantity').get('source'),'DIFFIGLF')
            self.assertEqual(xml.find('supplementary_quantity').get('radial_coordinate'),None)

    def test_mmm_and_glf23_can_coexist_without_double_counting(self):
        case=ROOT/'examples/data/d3d/147634T16/147634T16.CDF'
        if not case.is_file():
            self.skipTest('T16 fixture missing')
        with netcdf_file(case,mmap=False) as nc:
            # Synthetic composite, deliberately preserve all real MMM T16 data.
            variables=dict(nc.variables)
            shape=np.asarray(nc.variables['XB'][:]).shape
            variables['XKEGLF23']=FakeVar(np.ones(shape)*1234,('TIME3','XB'),'CM**2/SEC')
            variables['XKIGLF23']=FakeVar(np.ones(shape)*250,('TIME3','XB'),'CM**2/SEC')
            composite=types.SimpleNamespace(variables=variables)
            models=TestMMMWriter()._write(composite)
            self.assertIn('glf23',models)
            self.assertTrue(all(spec['name'] in models for spec in MODEL_SPECS))
            self.assertEqual(sum(name.startswith('mmm07_') for name in models),5)
            self.assertAlmostEqual(models['glf23'].profiles_1d[0].electrons.energy.d.values[0], .1234)
            self.assertEqual(models['mmm07_total'].code.name,'TRANSP/MMM07')
            self.assertEqual(models['glf23'].code.name,'TRANSP/GLF23')

    def test_glf23_invalid_grid_and_units_fail(self):
        cdf=FakeCDF()
        cdf.variables['XKEGLF23']=FakeVar(np.ones((2,3)),('TIME3','XB'),'EV')
        with self.assertRaisesRegex(ValueError,'unexpected units'):
            build_glf23_plan(cdf)
        cdf.variables['XKEGLF23'].units='CM**2/SEC'
        cdf.variables['XKEGLF23'].dimensions=('TIME3','X')
        with self.assertRaisesRegex(ValueError,'expected'):
            build_glf23_plan(cdf)


def verify_case(nc, case):
    plans = build_theory_plans(nc)
    mmm = plans['mmm07']
    report = mmm_summary(mmm)
    report['case'] = str(case)
    if mmm['active']:
        for name, result in mmm['closure'].items():
            if result['status'] == 'FAIL':
                raise AssertionError(f'{case}: MMM {name} failed closure: '
                                     f'{result["max_normalized_difference"]:.4g}')
        if not report['models']:
            raise AssertionError(f'{case}: active MMM, but no supported model')
    glf = plans['glf23']
    report['glf23'] = {'active': glf['active'],
                       'source_channels': list(glf['fields']) if glf['active'] else []}
    report['tglf'] = plans['tglf']['status']
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='run fast synthetic/writer tests only')
    parser.add_argument('--family', choices=('all','mmm07','glf23'), default='all')
    parser.add_argument('--data-root', type=Path, default=ROOT / 'examples/data')
    parser.add_argument('--report-dir', type=Path, default=ROOT / 'build/mmm-transport')
    ns = parser.parse_args(argv)
    tests = unittest.TestSuite([
        unittest.defaultTestLoader.loadTestsFromTestCase(TestMMMContracts),
        unittest.defaultTestLoader.loadTestsFromTestCase(TestMMMWriter),
    ])
    tests.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestGenericTheory))
    result = unittest.TextTestRunner(verbosity=1).run(tests)
    if not result.wasSuccessful():
        return 1
    if ns.quick:
        return 0
    paths = time_dependent_cdfs(ns.data_root)
    if not paths:
        parser.error(f'no CDF examples in {ns.data_root}')
    results = []
    failed = []
    for path in paths:
        try:
            with netcdf_file(path, mmap=False) as nc:
                item = verify_case(nc, path.relative_to(ns.data_root))
            results.append(item)
            status = 'ACTIVE' if item['active'] else 'INACTIVE'
            print(f'{path.relative_to(ns.data_root)}: PASS (MMM={status}, '
                  f'GLF23={"ACTIVE" if item["glf23"]["active"] else "inactive"}, '
                  f'MMM components={len(item["models"])})')
        except Exception as exc:
            failed.append(f'{path}: {type(exc).__name__}: {exc}')
            print('FAIL', failed[-1])
    ns.report_dir.mkdir(parents=True, exist_ok=True)
    summary = {'cases_checked': len(paths), 'passed': len(results),
               'failed': len(failed), 'cases': results, 'errors': failed}
    (ns.report_dir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    lines = ['# Unified theory anomalous transport source verification', '',
             f'Checked {len(paths)} CDFs; passed {len(results)}, failed {len(failed)}.', '',
             '| CDF | MMM result | Models | GLF23 | Electron closure | Ion closure |',
             '|---|---|---:|---|---|---|']
    for item in results:
        closure = item['closure']
        def show(name):
            c = closure.get(name, {})
            return (f'{c["status"]} ({c["max_normalized_difference"]:.2g})'
                    if 'max_normalized_difference' in c else 'N/A')
        lines.append(f'| {item["case"]} | {"active" if item["active"] else "inactive"} | '
                     f'{len(item["models"])} | {"active" if item["glf23"]["active"] else "inactive"} | {show("electron_heat")} | {show("ion_heat")} |')
    if failed:
        lines.extend(['', '## Errors', ''] + ['- ' + x for x in failed])
    (ns.report_dir / 'summary.md').write_text('\n'.join(lines) + '\n')
    return int(bool(failed))


if __name__ == '__main__':
    raise SystemExit(main())
