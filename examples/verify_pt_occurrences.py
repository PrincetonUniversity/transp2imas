#!/usr/bin/python
"""Source-level occurrence contract tests; no installed IMAS required."""
import sys
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from input_provenance import InputInventory
from pt_solver_regions import build_pt_region_plan, occurrence_region_models, region_occurrence_xml

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


    def test_combined_aos_allocated_once_before_population(self):
        import ast
        source = (ROOT / 'imas_core_transport.py').read_text()
        module = ast.parse(source)
        cls = next(x for x in module.body if isinstance(x, ast.ClassDef)
                   and x.name == 'imas_core_transport')
        create = next(x for x in cls.body if isinstance(x, ast.FunctionDef)
                      and x.name == 'create_ids')
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
