#!/usr/bin/python
"""Verify PT_SOLVER multi-region provenance and core_transport associations.

Runs without IMAS-Python.  All CDF examples are inventoried; only explicitly
linked PT_TEMPLATE files are treated as confirmed active model selections.
"""

import argparse
import json
from pathlib import Path
import sys
import unittest
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from input_provenance import InputInventory
from pt_solver_regions import (
    append_pt_regions_xml, annotate_theory_model_xml, append_pt_regions_reference,
    build_pt_region_plan, region_plan_xml, _pt_controls, _trdat_states, _xanom,
)


class IDSString0D:
    def __init__(self, value):
        self.value = value


class TestRegionMapping(unittest.TestCase):
    def test_boundaries_are_per_equation_and_follow_updates(self):
        source = '''LPREDICTIVE_MODE=3
XIMIN_CONF=0.1
XIMAX_CONF=0.8
XIBOUND=0.75
XNBOUND=0.9
XPHIBOUND=0.65
LPREDICT_TE=0
~update_time=1.5
LPREDICT_TE=1
XIBOUND=0.6
'''
        snaps = _trdat_states(source)
        self.assertEqual(len(snaps), 2)
        self.assertEqual(snaps[0]['values']['LPREDICT_TE'], 0)
        self.assertEqual(snaps[1]['values']['LPREDICT_TE'], 1)
        self.assertEqual(snaps[1]['start_time'], 1.5)
        self.assertEqual(snaps[0]['values']['XIBOUND'], .75)
        self.assertEqual(snaps[1]['values']['XIBOUND'], .6)
        plan = {'trdat_present': True, 'trdat_file': '/tmp/sampleTR.DAT',
                'predictive_mode': 3, 'boundary_states': snaps, 'templates': [],
                'model_source': 'no_PT_TEMPLATE', 'trusted_model_assignment': False}
        root = region_plan_xml(plan)
        regions = root.findall('.//radial_domains/region')
        self.assertEqual(len(regions), 18)
        def get(step, model, eq):
            el = root.findall('./boundary_and_equation_timeline/state')[step]
            return next(x for x in el.findall('./radial_domains/region')
                        if x.get('name') == 'pt_' + model and x.get('quantity') == eq)
        # Distinct X*BOUND limits are essential for interpreting region extent.
        self.assertEqual(get(0, 'confinement', 'temperature').get('xi_max'), '0.75')
        self.assertEqual(get(0, 'edge', 'temperature').get('has_extent'), 'false')
        self.assertEqual(get(0, 'edge', 'density').get('xi_max'), '0.9')
        self.assertEqual(get(0, 'edge', 'toroidal_momentum').get('has_extent'), 'false')
        self.assertEqual(get(1, 'confinement', 'temperature').get('xi_max'), '0.6')
        self.assertEqual(get(1, 'edge', 'temperature').get('has_extent'), 'false')

    def test_pt_template_multiple_scheduled_configs(self):
        text = '''2
&pt_updates
pt_times=-1.0, 0.5
/
&pt_controls
pt_update_time=-1.
pt_confinement%glf23%active=.true.
pt_confinement%glf23%xanom=6*.5
pt_confinement%exb%version=3
/
&pt_controls
pt_update_time=0.5
pt_confinement%tglf%active=.true.
pt_confinement%glf23%active=.false.
/
'''
        blocks = _pt_controls(text)
        self.assertEqual(len(blocks), 2)
        self.assertEqual(blocks[0]['start_time'], -1)
        self.assertEqual(blocks[1]['start_time'], .5)
        self.assertEqual(_xanom(blocks[0]['regions']['confinement']['glf23']['xanom']), [.5] * 6)
        self.assertEqual(blocks[1]['regions']['confinement']['glf23']['active'], '.false.')
        plan = {'trdat_present': True, 'trdat_file': 'caseTR.DAT',
                'predictive_mode': 3, 'boundary_states': _trdat_states('LPREDICTIVE_MODE=3'),
                'model_source': 'explicit_PT_TEMPLATE', 'trusted_model_assignment': True,
                'templates': [{'path': 'pt.nml', 'sha256': 'abcdef', 'source_role': 'TEMPLATE',
                               'blocks': blocks}]}
        root = region_plan_xml(plan)
        self.assertEqual(len(root.findall('.//configuration')), 2)
        self.assertEqual(root.find('.//option[@key="xanom"]').get('values'), '0.5,0.5,0.5,0.5,0.5,0.5')
        self.assertEqual(root.find('.//option[@key="version"]').get('scheme'), 'TGYROEXB')
        out = annotate_theory_model_xml('<original_output/>', plan, 'glf23')
        selection = ET.fromstring(out).findall('./pt_solver_region_application/selection')
        self.assertEqual([(s.get('active'), s.get('effective_from')) for s in selection],
                         [('true', 'initial'), ('false', '0.5')])

    def test_wrapped_imas_xml_preserves_source_provenance(self):
        plan = {'trdat_present': True, 'trdat_file': 'test', 'predictive_mode': 3,
                'boundary_states': _trdat_states('LPREDICTIVE_MODE=3'),
                'templates': [], 'model_source': 'no_PT_TEMPLATE',
                'trusted_model_assignment': False}
        wrapped = IDSString0D('<transp2imas_parameters><input_files id="original"/></transp2imas_parameters>')
        text = append_pt_regions_xml(wrapped, plan)
        root = ET.fromstring(text)
        self.assertIsNotNone(root.find('input_files[@id="original"]'))
        self.assertIsNotNone(root.find('pt_solver_regions'))
        self.assertEqual(len(ET.fromstring(append_pt_regions_xml(text, plan)).findall('pt_solver_regions')), 1)
        # A missing/no-region file must not erase input provenance either.
        self.assertIsNotNone(ET.fromstring(append_pt_regions_xml(wrapped, None)).find('input_files'))
        self.assertEqual(ET.fromstring(append_pt_regions_reference('', plan)).find(
            'pt_solver_region_configuration_reference').get('source'), 'no_PT_TEMPLATE')

    def test_end_to_end_readback_detects_changed_region_input(self):
        import types
        import verify_end_to_end as e2e
        trdat = ROOT / 'examples/data/iter/38530P94TR.DAT'
        if not trdat.is_file():
            self.skipTest('ITER TR.DAT unavailable')
        plan = build_pt_region_plan(InputInventory(trdat, '38530P94'))
        serialized = append_pt_regions_xml('<transp2imas_parameters/>', plan)
        ct = types.SimpleNamespace(code=types.SimpleNamespace(
            parameters=IDSString0D(serialized)))
        self.assertEqual(e2e.check_pt_solver_region_readback(ct, plan)[
            'pt_solver_region_configuration'], 'PASS')
        changed = serialized.replace('xi="0.75"', 'xi="0.76"', 1)
        self.assertNotEqual(changed, serialized)
        ct.code.parameters.value = changed
        with self.assertRaisesRegex(AssertionError, 'changed'):
            e2e.check_pt_solver_region_readback(ct, plan)

    def test_input_provenance_on_actual_three_cases(self):
        for path, run in [
            ('iter/38530P94TR.DAT', '38530P94'),
            ('nstx/141711Y05/141711Y05TR.DAT', '141711Y05'),
            ('d3d/147634T16/147634T16TR.DAT', '147634T16'),
        ]:
            tr = ROOT / 'examples/data' / path
            if not tr.parent.is_dir():
                self.skipTest('bundled PT_SOLVER inputs are not available')
            plan = build_pt_region_plan(InputInventory(tr, run))
            self.assertIsNotNone(plan)
            root = region_plan_xml(plan)
            if run == '38530P94':
                self.assertTrue(plan['trusted_model_assignment'])
                self.assertEqual(len(plan['boundary_states']), 3)
                self.assertEqual(plan['boundary_states'][2]['values']['LPREDICT_TE'], 1)
                first = root.find('./boundary_and_equation_timeline/state/radial_domains')
                extent = {x.get('name'): (x.get('xi_min'), x.get('xi_max'))
                          for x in first.findall('region[@quantity="temperature"]')}
                self.assertEqual(extent, {
                    'pt_axial': ('0', '0.1'),
                    'pt_confinement': ('0.1', '0.75'),
                    'pt_edge': ('0.75', '1'),
                })
                self.assertEqual(len(root.findall('.//configuration/model[@name="GLF23"]')), 3)
                self.assertEqual(len(root.findall('.//configuration/model[@name="NEOCH"]')), 3)
                out = annotate_theory_model_xml(IDSString0D('<glf23/>'), plan, 'glf23')
                self.assertEqual(len(ET.fromstring(out).findall('.//selection[@active="true"]')), 3)
            elif run == '141711Y05':
                self.assertTrue(plan['trusted_model_assignment'])
                first = root.find('./boundary_and_equation_timeline/state/radial_domains')
                self.assertEqual(sum(x.get('has_extent') == 'true'
                                     for x in first.findall('region[@quantity="temperature"]')), 1)
                self.assertEqual(len(root.findall('.//configuration/model[@name="CDBM"]')), 1)
            else:
                # The bundled T16 fixture may include its actual TR.DAT.
                # Presence must be detected, not hardcoded from an older fixture set.
                self.assertEqual(plan['trdat_present'], tr.is_file())
                if plan['trdat_present']:
                    self.assertTrue(plan['boundary_states'])
                    self.assertNotEqual(root.find('boundary_and_equation_timeline').get('status'),
                                        'unavailable_without_TR.DAT')
                    if plan['trusted_model_assignment']:
                        self.assertEqual(plan['model_source'], 'explicit_PT_TEMPLATE')
                    else:
                        self.assertNotEqual(plan['model_source'], 'explicit_PT_TEMPLATE')
                        self.assertNotIn('pt_solver_region_application',
                                         annotate_theory_model_xml('<mmm/>', plan, 'mmm07_total'))
                else:
                    self.assertFalse(plan['trusted_model_assignment'])
                    self.assertEqual(plan['model_source'], 'adjacent_unreferenced_candidate')
                    self.assertEqual(root.find('boundary_and_equation_timeline').get('status'),
                                     'unavailable_without_TR.DAT')
                    self.assertNotIn('pt_solver_region_application',
                                     annotate_theory_model_xml('<mmm/>', plan, 'mmm07_total'))

    def test_integration_with_actual_transport_writer(self):
        # A small mock CDF avoids opening all large physical CDFs in quick CI.
        from transport_test_support import FakeConverter, write_core_transport
        from verify_theory_transport import FakeCDF, FakeVar
        nc = FakeCDF()
        import numpy as np
        nc.variables['DVOL'] = FakeVar(np.ones((2, 3)), ('TIME3', 'X'), 'CM**3')
        nc.variables['DAREA'] = FakeVar(np.ones((2, 3)), ('TIME3', 'X'), 'CM**2')
        nc.variables['SURF'] = FakeVar(np.ones((2, 3)), ('TIME3', 'XB'), 'CM**2')
        nc.variables['XKEGLF23'] = FakeVar(np.ones((2, 3)), ('TIME3', 'XB'), 'CM**2/SEC')
        nc.variables['EPTR_MOD'] = FakeVar(np.ones((2, 3)), ('TIME3', 'X'), 'N/CM3/SEC')
        nc.variables['GFLNC_E'] = FakeVar(np.ones((2, 3)), ('TIME3', 'X'), 'N/CM3/SEC')
        nc.variables['QFLNC_E'] = FakeVar(np.ones((2, 3)), ('TIME3', 'X'), 'WATTS/CM3')
        trdat = ROOT / 'examples/data/iter/38530P94TR.DAT'
        if not trdat.is_file():
            self.skipTest('ITER fixture missing')
        inventory = InputInventory(trdat, '38530P94')
        original = FakeConverter
        method = getattr(original, 'input_inventory', None)
        original.input_inventory = lambda self: inventory
        try:
            ct = write_core_transport(nc)
        finally:
            if method is None:
                delattr(original, 'input_inventory')
            else:
                original.input_inventory = method
        root = ET.fromstring(ct.code.parameters)
        self.assertIsNotNone(root.find('pt_solver_regions'))
        # Five MMM component outputs remain five global model entries, not
        # 15 fake region models or three separate IDS occurrences.
        models = {model.identifier.name: model for model in ct.model._items}
        self.assertEqual(len([name for name in models if name.startswith('mmm07')]), 5)
        self.assertIsNone(ET.fromstring(models['mmm07_total'].code.parameters).find(
            'pt_solver_region_application'))  # ITER config specifies GLF23, not MMM.
        self.assertEqual(len(ET.fromstring(models['glf23'].code.parameters).findall(
            './pt_solver_region_application/selection')), 3)
        self.assertEqual(len(ET.fromstring(models['neoclassical'].code.parameters).findall(
            './pt_solver_region_application/selection')), 3)
        self.assertIsNotNone(ET.fromstring(models['transport_solver'].code.parameters).find(
            'pt_solver_region_configuration_reference'))


def inventory_all(root, report_dir=None):
    cases = sorted(root.rglob('*.CDF'))
    if not cases:
        raise AssertionError(f'no bundled CDF files in {root}')
    records = []
    for cdf in cases:
        trdat = cdf.with_name(cdf.stem + 'TR.DAT')
        inventory = InputInventory(trdat, cdf.stem)
        plan = build_pt_region_plan(inventory)
        obj = {'cdf': str(cdf.relative_to(root)),
               'status': ('NO_PT_SOLVER_CONFIGURATION' if plan is None else
                          'EXPLICIT_INPUT' if plan['trusted_model_assignment'] else
                          'CANDIDATE_INPUT' if plan['templates'] else
                          'TRDAT_ONLY'),
               'model_source': plan['model_source'] if plan else None,
               'boundary_time_states': len(plan['boundary_states']) if plan else 0,
               'model_time_states': sum(len(r['blocks']) for r in plan['templates']) if plan else 0}
        records.append(obj)
        print(f"{obj['cdf']}: {obj['status']}; boundary states={obj['boundary_time_states']}; "
              f"model states={obj['model_time_states']}")
        if plan:
            ET.fromstring(ET.tostring(region_plan_xml(plan)))
    if report_dir:
        path = Path(report_dir)
        path.mkdir(parents=True, exist_ok=True)
        (path / 'pt-solver-regions.json').write_text(json.dumps(records, indent=2) + '\n')
    print(f'PT_SOLVER region input coverage: PASS ({len(records)}/{len(cases)} CDF examples inventoried)')
    return records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--data-root', default=str(ROOT / 'examples/data'))
    ap.add_argument('--report-dir', default=None)
    args = ap.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestRegionMapping)
    result = unittest.TextTestRunner(verbosity=1).run(suite)
    if not result.wasSuccessful():
        return 1
    if not args.quick:
        inventory_all(Path(args.data_root), args.report_dir)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
