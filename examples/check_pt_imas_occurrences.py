#!/usr/bin/python
"""Read back actual core_transport occurrence 0 and active PT_SOLVER regions.

Run from checkout: python3 examples/check_pt_imas_occurrences.py --machine ITER \
    --shot 38530 --run 1694 --source-run 38530P94 \
    --trdat examples/data/iter/38530P94TR.DAT --user pankin
Requires installed IMAS-Python and the HDF5 database entry from convert.py.
"""
import argparse
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from input_provenance import InputInventory
from pt_solver_regions import build_pt_region_plan, interval_attributes, occurrence_region_models


def string(value):
    value = getattr(value, 'value', value)
    if isinstance(value, bytes):
        return value.decode('utf-8')
    return str(value)


def inspect_occurrence(ct, index, expected=None):
    params = string(ct.code.parameters)
    root = ET.fromstring(params)
    print(f'Occurrence {index}: time={len(ct.time)}, models={len(ct.model)}, xml={root.tag}')
    assert len(ct.time) > 0, f'occurrence {index}: missing time axis'
    if index == 0:
        assert root.find('.//pt_solver_regions') is not None, 'occurrence 0: missing PT_SOLVER XML'
        return
    assert root.tag == 'pt_solver_region_occurrence', (index, root.tag)
    assert root.get('occurrence') == str(index)
    assert root.get('region') == expected['region'], (root.get('region'),expected)
    selections = root.findall('model_selection')
    expected_pairs = {(x['category'], x['model'].upper()) for x in expected['selections']}
    assert {(x.get('category'), x.get('name')) for x in selections} == expected_pairs, \
        ([x.attrib for x in selections], expected_pairs)
    assert root.find('./pt_solver_regions/boundary_and_equation_timeline/state') is not None
    assert root.find('.//input_file/namelist') is not None, 'referenced PT_TEMPLATE namelist XML absent'
    seen = set()
    for model in ct.model:
        model_xml = ET.fromstring(string(model.code.parameters))
        scope = model_xml.find('pt_solver_occurrence')
        assert scope is not None
        assert scope.get('region') == expected['region']
        category = scope.get('category')
        attribution = scope.get('output_attribution')
        identity = string(model.identifier.name)
        listed = scope.get('selected_models').split(',')
        if attribution == 'assembled_result_not_attributed_to_a_selected_model':
            # the assembled neoclassical result is never renamed to a selected model
            assert identity == 'neoclassical' and scope.get('model') == 'neoclassical'
        else:
            assert identity == scope.get('model') == listed[0]
        for node in scope.findall('selected_model'):
            expected_item = next(x for x in expected['selections']
                                 if x['category'] == category and x['model'].upper() == node.get('name'))
            stored = [(i.get('effective_from'), i.get('effective_until')) for i in node.findall('interval')]
            want = [interval_attributes(a, b) for a, b in expected_item['intervals']]
            assert stored == [(w['effective_from'], w['effective_until']) for w in want], (stored, want)
            seen.add((category, node.get('name')))
        state = scope.get('cdf_model_output_available')
        assert state in ('true', 'false')
        if state == 'false':
            assert scope.get('output_scope') == 'settings_only'
        print(f'  {category}: slot {identity!r} selected={scope.get("selected_models")}, '
              f'CDF output={state}, profile slices={len(model.profiles_1d)}')
    assert seen == expected_pairs, (seen, expected_pairs)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--machine',required=True)
    ap.add_argument('--shot',required=True,type=int)
    ap.add_argument('--run',required=True,type=int)
    ap.add_argument('--user',default=None)
    ap.add_argument('--dd-version',default='4.1.1')
    ap.add_argument('--source-run',required=True,help='Full TRANSP run id, e.g. 38530P94')
    ap.add_argument('--trdat',required=True,type=Path)
    a = ap.parse_args()
    try:
        import imas
    except ImportError as exc:
        raise SystemExit('IMAS-Python not installed in current environment') from exc
    plan = build_pt_region_plan(InputInventory(a.trdat.resolve(), a.source_run))
    assert plan and plan['trusted_model_assignment'], 'No authoritative PT_TEMPLATE found'
    regions = occurrence_region_models(plan)
    print('Expected regional occurrences:',[(i+1,r['region'],[(x['category'],x['model'],x['intervals']) for x in r['selections']]) for i,r in enumerate(regions)])
    backend = imas.ids_defs.HDF5_BACKEND
    db = imas.DBEntry(backend, a.machine, a.shot, a.run,
                      a.user or __import__('getpass').getuser(),
                      dd_version=a.dd_version)
    db.open()
    try:
        for index in range(len(regions)+1):
            try:ct = db.get('core_transport', occurrence=index)
            except Exception as exc:
                raise AssertionError(f'missing required occurrence {index}: {exc}') from exc
            inspect_occurrence(ct,index,None if index==0 else regions[index-1])
        # Explicitly probe the first forbidden occurrence. A missing occurrence
        # may manifest as IMAS exception or empty IDS depending on backend.
        extra=len(regions)+1
        try:
            unwanted=db.get('core_transport',occurrence=extra)
        except Exception:
            unwanted=None
        if unwanted is not None:
            try:
                present=bool(len(unwanted.time) or len(unwanted.model))
            except Exception:
                present=True
            assert not present, f'unexpected unused occurrence {extra}'
    finally:
        db.close()
    print('PASS: all expected regional occurrences read back and no extra occurrence found')

if __name__=='__main__':main()
