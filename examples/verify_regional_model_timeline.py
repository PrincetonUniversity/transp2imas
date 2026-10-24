#!/usr/bin/python
# registry ids verified here (mapping_registry.yaml): pt_regional_models
"""Verify time-resolved regional model attribution for PT_SOLVER occurrences.

The converter used to keep only the model selected *last* in each category of a region, so a
GLF23 -> TGLF switch became a single TGLF identifier over the whole record, a model disabled
before the end disappeared (even its region could vanish), and every selected neoclassical
model name was attached to the same assembled neoclassical CDF result.  Now:

* ``region_selection_timeline`` / ``occurrence_region_models`` keep every selected model with its
  intervals (``initial`` ... ``end``), with ``final_models`` as the end-of-run view only,
* one model slot per selected anomalous model; TRANSP publishes one assembled neoclassical
  result, which keeps the name ``neoclassical`` (never a selected model's name) and lists all
  selected neoclassical models with their intervals as metadata,
* the radial validity states carry explicit ``effective_until``.

1. unit tests on synthetic plans (no IMAS),
2. (not with ``--quick``) ITER converted with a GLF23 -> TGLF / NEOCH -> NCLASS switch added to its
   PT_TEMPLATE and read back from HDF5.

    python examples/verify_regional_model_timeline.py [--quick]
"""

import argparse
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pt_solver_regions import (interval_attributes, occurrence_region_models, region_occurrence_xml,
                               region_selection_timeline, region_validity_xml, regional_model_slots)

try:
    import imas
except ImportError:      # pragma: no cover - environment dependent
    imas = None

ITER_DIR = ROOT / 'examples' / 'data' / 'iter'
RUN = {'quick': False}


def block(start, **regions):
    """A template block: ``block(552, confinement={'glf23': 'F', 'tglf': 'T'})``."""
    return {'start_time': start, 'options': {},
            'regions': {region: {model: {'active': flag} for model, flag in models.items()}
                        for region, models in regions.items()}}


def plan_of(*blocks, states=((None, 0.2, 0.8), (552.0, 0.2, 0.8))):
    """States are ``(start, XIMIN_CONF, XIMAX_CONF)``: 0.2/0.8 gives axial, confinement and edge an extent."""
    boundary_states = []
    for t, lo, hi in states:
        values = {'XIMIN_CONF': lo, 'XIMAX_CONF': hi, 'XIBOUND': 1.0, 'XNBOUND': 1.0, 'XPHIBOUND': 1.0}
        boundary_states.append({'start_time': t, 'values': values, 'origins': {k: 'TR.DAT' for k in values}})
    return {'trusted_model_assignment': True, 'boundary_states': boundary_states, 'trdat_present': True,
            'trdat_file': 'TR.DAT', 'predictive_mode': 3, 'model_source': 'explicit_PT_TEMPLATE',
            'templates': [{'path': 'pt.nml', 'sha256': 'x', 'text': '&pt_controls /', 'blocks': list(blocks)}]}


def selection(item, category, model):
    return next(s for s in item['selections'] if s['category'] == category and s['model'] == model)


class TimelineTests(unittest.TestCase):
    def test_a_model_switch_keeps_both_models_with_their_intervals(self):
        plan = plan_of(block(-1, confinement={'glf23': '.true.'}),
                       block(552.0, confinement={'glf23': '.false.', 'tglf': '.true.'}))
        item = occurrence_region_models(plan)[0]
        self.assertEqual(item['region'], 'confinement')
        self.assertEqual(selection(item, 'anomalous', 'glf23')['intervals'], [(None, 552.0)])
        self.assertEqual(selection(item, 'anomalous', 'tglf')['intervals'], [(552.0, None)])
        self.assertEqual(item['final_models'], {'anomalous': 'tglf', 'neoclassical': None})

    def test_a_replacement_without_an_explicit_disable_still_ends_the_old_interval(self):
        plan = plan_of(block(-1, confinement={'glf23': 'T'}), block(552.0, confinement={'tglf': 'T'}))
        item = occurrence_region_models(plan)[0]
        self.assertEqual(selection(item, 'anomalous', 'glf23')['intervals'], [(None, 552.0)])
        self.assertEqual(selection(item, 'anomalous', 'tglf')['intervals'], [(552.0, None)])

    def test_a_model_disabled_before_the_end_is_not_lost_and_keeps_its_region(self):
        plan = plan_of(block(-1, edge={'glf23': 'T'}), block(552.0, edge={'glf23': 'F'}))
        items = occurrence_region_models(plan)
        self.assertEqual([i['region'] for i in items], ['edge'])           # the old code returned []
        self.assertEqual(selection(items[0], 'anomalous', 'glf23')['intervals'], [(None, 552.0)])
        self.assertEqual(items[0]['final_models'], {'anomalous': None, 'neoclassical': None})

    def test_re_enabling_adds_a_second_interval(self):
        plan = plan_of(block(-1, axial={'mmm': 'T'}), block(2.0, axial={'mmm': 'F'}), block(5.0, axial={'mmm': 'T'}))
        self.assertEqual(selection(occurrence_region_models(plan)[0], 'anomalous', 'mmm')['intervals'],
                         [(None, 2.0), (5.0, None)])

    def test_models_enabled_together_at_one_instant_leave_no_empty_interval(self):
        plan = plan_of(block(-1, confinement={'glf23': 'T'}),
                       block(3.0, confinement={'tglf': 'T', 'mmm': 'T'}))
        item = occurrence_region_models(plan)[0]
        self.assertEqual(selection(item, 'anomalous', 'glf23')['intervals'], [(None, 3.0)])
        # tglf was replaced by mmm at the instant it was enabled: no zero-length interval, and
        # it is not listed as a model that was ever really selected
        self.assertEqual([s['model'] for s in item['selections'] if s['category'] == 'anomalous'],
                         ['glf23', 'mmm'])
        self.assertEqual(selection(item, 'anomalous', 'mmm')['intervals'], [(3.0, None)])
        self.assertEqual(item['final_models']['anomalous'], 'mmm')

    def test_blocks_are_ordered_chronologically_not_by_file_position(self):
        plan = plan_of(block(552.0, confinement={'glf23': 'F', 'tglf': 'T'}), block(-1, confinement={'glf23': 'T'}))
        item = occurrence_region_models(plan)[0]
        self.assertEqual(selection(item, 'anomalous', 'glf23')['intervals'], [(None, 552.0)])

    def test_regions_and_categories_are_independent_and_unreadable_flags_are_ignored(self):
        plan = plan_of(block(-1, axial={'glf23': 'T'}, edge={'tglf': 'T'}, confinement={'neoch': 'maybe'}),
                       block(552.0, axial={'neoch': 'T'}))
        timeline = region_selection_timeline(plan)
        self.assertEqual(timeline['axial']['anomalous'], [{'model': 'glf23', 'intervals': [(None, None)]}])
        self.assertEqual(timeline['axial']['neoclassical'], [{'model': 'neoch', 'intervals': [(552.0, None)]}])
        self.assertEqual(timeline['edge']['anomalous'][0]['model'], 'tglf')
        self.assertEqual(timeline['confinement']['neoclassical'], [])         # 'maybe' is never guessed
        self.assertEqual([i['region'] for i in occurrence_region_models(plan)], ['axial', 'edge'])

    def test_a_region_without_radial_extent_gets_no_occurrence(self):
        plan = plan_of(block(-1, axial={'glf23': 'T'}), states=((None, 0.0, 1.0),))   # XIMIN_CONF=0: no axial extent
        self.assertEqual(occurrence_region_models(plan), [])

    def test_xml_carries_every_selection_with_explicit_intervals(self):
        plan = plan_of(block(-1, confinement={'glf23': 'T'}), block(552.0, confinement={'tglf': 'T'}))
        item = occurrence_region_models(plan)[0]
        root = ET.fromstring(region_occurrence_xml(plan, item, 1))
        stored = {(n.get('category'), n.get('name')): [(i.get('effective_from'), i.get('effective_until'))
                                                        for i in n.findall('interval')]
                  for n in root.findall('model_selection')}
        self.assertEqual(stored, {('anomalous', 'GLF23'): [('initial', '552')],
                                  ('anomalous', 'TGLF'): [('552', 'end')]})
        self.assertEqual(interval_attributes(None, None), {'effective_from': 'initial', 'effective_until': 'end'})

    def test_validity_states_have_explicit_end_times(self):
        plan = plan_of(block(-1, axial={'glf23': 'T'}), states=((None, 0.1, 0.9), (2.0, 0.2, 0.8), (7.5, 0.3, 0.7)))
        states = region_validity_xml(plan, 'confinement').findall('state')
        self.assertEqual([(s.get('effective_from'), s.get('effective_until')) for s in states],
                         [('initial', '2'), ('2', '7.5'), ('7.5', 'end')])


class SlotPolicyTests(unittest.TestCase):
    NC, NO_NC = True, False        # whether the CDF has neoclassical channels

    def item(self, *blocks):
        return occurrence_region_models(plan_of(*blocks))[0]

    def test_anomalous_models_each_get_a_slot_and_there_are_no_placeholders(self):
        item = self.item(block(-1, confinement={'glf23': 'T'}), block(552.0, confinement={'tglf': 'T'}))
        slots = regional_model_slots(item, self.NC)
        self.assertEqual([(s['category'], s['selected'][0]['model']) for s in slots],
                         [('anomalous', 'glf23'), ('anomalous', 'tglf')])      # no neoclassical selected: no slot

    def test_one_assembled_neoclassical_slot_lists_every_selected_model(self):
        item = self.item(block(-1, confinement={'neoch': 'T'}), block(552.0, confinement={'nclass': 'T'}))
        slots = regional_model_slots(item, self.NC)
        self.assertEqual(len(slots), 1)
        self.assertTrue(slots[0]['assembled'])
        self.assertEqual([s['model'] for s in slots[0]['selected']], ['neoch', 'nclass'])

    def test_without_neoclassical_channels_each_model_gets_a_settings_only_slot(self):
        item = self.item(block(-1, confinement={'neoch': 'T'}), block(552.0, confinement={'nclass': 'T'}))
        slots = regional_model_slots(item, self.NO_NC)
        self.assertEqual([s['selected'][0]['model'] for s in slots], ['neoch', 'nclass'])
        self.assertFalse(any(s.get('assembled') for s in slots))


def read_occurrences(cdf, outdir, dd_version='4.1.1'):
    result = subprocess.run(
        [sys.executable, str(ROOT / 'convert.py'), '-i', str(cdf), '-b', 'HDF5', '-u', str(outdir),
         '--dd-version', dd_version, '-new', '--ids', 'core_transport'],
        cwd=ROOT, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stdout[-1500:] + result.stderr[-1500:])
    path = Path(re.search(r'ids stored in (\S+)', result.stdout).group(1))
    database, _major, shot, run = path.parts[-4:]
    entry = imas.DBEntry(imas.ids_defs.HDF5_BACKEND, database, int(shot), int(run), str(path.parents[3]),
                         dd_version=dd_version)
    entry.open()
    occurrences = {}
    try:
        for occurrence in range(1, 6):
            try:
                ids = entry.get('core_transport', occurrence=occurrence)
            except Exception:                                  # first empty occurrence ends the list
                break
            occurrences[occurrence] = ids
    finally:
        entry.close()
    return occurrences


def slot_view(model):
    xml = ET.fromstring(str(model.code.parameters))
    scope = xml.find('pt_solver_occurrence')
    intervals = {n.get('name'): [(i.get('effective_from'), i.get('effective_until')) for i in n.findall('interval')]
                 for n in scope.findall('selected_model')}
    return str(model.identifier.name), scope, intervals


SWITCH_BLOCK = '''
&pt_controls
 pt_update_time = 552.
 pt_confinement%glf23%active = .false.
 pt_confinement%tglf%active = .true.
 pt_confinement%neoch%active = .false.
 pt_confinement%nclass%active = .true.
/
'''


@unittest.skipIf(imas is None, 'IMAS-Python is required')
class ConvertedSwitchTests(unittest.TestCase):
    def test_iter_with_a_glf23_to_tglf_and_neoch_to_nclass_switch(self):
        if RUN['quick'] or not ITER_DIR.is_dir():
            self.skipTest('end-to-end readback skipped (--quick or ITER example not present)')
        with tempfile.TemporaryDirectory() as tmp:
            case = Path(tmp) / 'case'
            shutil.copytree(ITER_DIR, case)
            template = case / 'p94_pt.nml'
            template.write_text(template.read_text() + SWITCH_BLOCK)
            occurrences = read_occurrences(case / '38530P94.CDF', Path(tmp) / 'db')
        regions = {str(ET.fromstring(str(ids.code.parameters)).get('region')): ids for ids in occurrences.values()}
        self.assertEqual(sorted(regions), ['axial', 'confinement', 'edge'])
        # confinement: both anomalous models and the assembled neoclassical slot
        names = [slot_view(m)[0] for m in regions['confinement'].model]
        self.assertEqual(names, ['GLF23', 'TGLF', 'neoclassical'])
        glf23, tglf, neoclassical = [slot_view(m) for m in regions['confinement'].model]
        self.assertEqual(glf23[2], {'GLF23': [('initial', '552')]})
        self.assertEqual(tglf[2], {'TGLF': [('552', 'end')]})
        self.assertEqual(glf23[1].get('output_attribution'), 'native_model_channels')
        self.assertEqual(tglf[1].get('cdf_model_output_available'), 'false')       # no verified TGLF channels
        self.assertEqual(tglf[1].get('output_scope'), 'settings_only')
        self.assertEqual(len(regions['confinement'].model[1].profiles_1d), 0)
        self.assertGreater(len(regions['confinement'].model[0].profiles_1d), 0)
        # neoclassical: one assembled slot, never renamed to NEOCH/NCLASS, both selections listed
        self.assertEqual(neoclassical[0], 'neoclassical')
        self.assertEqual(neoclassical[1].get('output_attribution'), 'assembled_result_not_attributed_to_a_selected_model')
        self.assertEqual(neoclassical[2], {'NEOCH': [('initial', '552')], 'NCLASS': [('552', 'end')]})
        for slot in regions['confinement'].model:
            self.assertNotIn(str(slot.identifier.name), ('NEOCH', 'NCLASS'))
        # the other regions did not switch
        for region in ('axial', 'edge'):
            self.assertEqual([slot_view(m)[0] for m in regions[region].model], ['GLF23', 'neoclassical'])
            self.assertEqual(slot_view(regions[region].model[0])[2], {'GLF23': [('initial', 'end')]})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='skip the HDF5 end-to-end readback')
    args, rest = parser.parse_known_args(argv)
    RUN['quick'] = args.quick
    result = unittest.main(argv=[sys.argv[0]] + rest, exit=False, verbosity=1).result
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
