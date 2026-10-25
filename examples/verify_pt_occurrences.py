#!/usr/bin/python
"""Source-level occurrence contract tests; no installed IMAS required."""
import sys
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from input_provenance import InputInventory
from pt_solver_regions import (MODEL_EXECUTION_SEMANTICS, _pt_controls, _trdat_states, build_pt_region_plan,
                               occurrence_region_models, region_declarations, region_occurrence_xml,
                               region_selection_timeline)

class TestOccurrences(unittest.TestCase):
    def plan(self, path, name):
        return build_pt_region_plan(InputInventory(str(ROOT/path),name))
    def test_three_region_iter(self):
        plan=self.plan(Path('examples/data/iter/38530P94TR.DAT'),'38530P94')
        regions=occurrence_region_models(plan)
        self.assertEqual([r['region'] for r in regions],['axial','confinement','edge'])
        for i,r in enumerate(regions,1):
            root=ET.fromstring(region_occurrence_xml(plan,r,i))
            self.assertEqual(root.get('occurrence'),str(i))
            self.assertEqual(len(root.findall('model_selection')),2)
            self.assertTrue(root.findall('.//input_file/namelist'))
            self.assertTrue(root.findall('.//boundary_and_equation_timeline/state'))
    def test_single_region_nstx(self):
        plan=self.plan(Path('examples/data/nstx/141711Y05/141711Y05TR.DAT'),'141711Y05')
        self.assertEqual([r['region'] for r in occurrence_region_models(plan)],['confinement'])


    SIMULTANEOUS = '''1
&pt_controls
pt_update_time=-1.
pt_confinement%glf23%active=.true.
pt_confinement%tglf%active=.true.
pt_confinement%neoch%active=.true.
/
&pt_controls
pt_update_time=0.5
pt_confinement%glf23%active=.false.
/
'''

    def simultaneous_plan(self):
        blocks = _pt_controls(self.SIMULTANEOUS)
        return {'trdat_present': True, 'trdat_file': 'caseTR.DAT', 'predictive_mode': 3,
                'boundary_states': _trdat_states('LPREDICTIVE_MODE=3'),
                'model_source': 'explicit_PT_TEMPLATE', 'trusted_model_assignment': True,
                'templates': [{'path': 'pt.nml', 'sha256': 'abcdef', 'source_role': 'TEMPLATE', 'blocks': blocks}]}

    def test_the_execution_semantics_are_declared_unresolved_not_asserted(self):
        self.assertEqual(MODEL_EXECUTION_SEMANTICS['status'], 'unresolved')
        self.assertEqual(MODEL_EXECUTION_SEMANTICS['evidence'], 'none supplied')
        plan = self.plan(Path('examples/data/iter/38530P94TR.DAT'), '38530P94')
        for i, item in enumerate(occurrence_region_models(plan), 1):
            root = ET.fromstring(region_occurrence_xml(plan, item, i))
            self.assertEqual(root.get('execution_semantics'), 'unresolved')
            self.assertEqual(root.get('execution_semantics_evidence'), 'none supplied')
            self.assertEqual(root.get('execution_semantics_transp_version'), 'not established')
            self.assertIn('replaces the previous one', root.get('execution_semantics_assumed_rule'))
            self.assertIsNone(root.get('selection_semantics'))                  # the assumption is no longer stated as fact
            self.assertTrue(root.findall('declarations/declaration'))
            for selection in root.findall('model_selection'):
                self.assertEqual(selection.get('interpretation'), 'replacement_assumed_unverified')

    def test_declarations_keep_simultaneous_models_that_the_replacement_interpretation_drops(self):
        plan = self.simultaneous_plan()
        declared = region_declarations(plan, 'confinement')
        self.assertEqual([(d['time'], d['category'], d['model'], d['active']) for d in declared],
                         [(None, 'anomalous', 'glf23', True), (None, 'anomalous', 'tglf', True),
                          (None, 'neoclassical', 'neoch', True), (0.5, 'anomalous', 'glf23', False)])
        # the interpretation (replacement) keeps one anomalous model; the declarations keep what was written
        anomalous = region_selection_timeline(plan)['confinement']['anomalous']
        self.assertEqual([entry['model'] for entry in anomalous], ['tglf'])
        item = occurrence_region_models(plan)[0]
        root = ET.fromstring(region_occurrence_xml(plan, item, 1))
        written = [(d.get('time'), d.get('category'), d.get('model'), d.get('active'))
                   for d in root.findall('declarations/declaration')]
        self.assertEqual(written, [('initial', 'anomalous', 'GLF23', 'true'), ('initial', 'anomalous', 'TGLF', 'true'),
                                   ('initial', 'neoclassical', 'NEOCH', 'true'), ('0.5', 'anomalous', 'GLF23', 'false')])

    def test_declarations_do_not_invent_regions_or_models_the_template_never_mentions(self):
        plan = self.simultaneous_plan()
        self.assertEqual(region_declarations(plan, 'axial'), [])
        self.assertEqual(region_declarations(plan, 'edge'), [])

    def test_combined_aos_allocated_once_before_population(self):
        import ast
        source = (ROOT / 'imas_core_transport.py').read_text()
        module = ast.parse(source)
        cls = next(x for x in module.body if isinstance(x, ast.ClassDef)
                   and x.name == 'imas_core_transport')
        create = next(x for x in cls.body if isinstance(x, ast.FunctionDef)
                      and x.name == 'build')
        text = ast.get_source_segment(source, create)
        self.assertEqual(text.count('ct.model.resize('), 1)
        self.assertIn('ct.model.resize(expected_models)', text)
        self.assertLess(text.index('ct.model.resize(expected_models)'),
                        text.index('self._fill_combined_model(ct.model[0])'))
        self.assertIn('len(ct.model[0].profiles_1d) != self.nt', text)

    def test_regional_slots_are_per_selected_model_and_never_relabel_assembled_output(self):
        source = (ROOT / 'imas_core_transport.py').read_text()
        # The slot is named for the output it holds (the assembled neoclassical result, or the
        # model's own native channels), never renamed to a selected model it does not come from.
        self.assertNotIn('slot.identifier.name = selected.upper()', source)
        self.assertIn("identity = 'neoclassical' if assembled else names[0].upper()", source)
        self.assertIn('slot.profiles_1d.resize(0)', source)

if __name__=='__main__':unittest.main()
