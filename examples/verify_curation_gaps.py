#!/usr/bin/python
"""Permanent regression tests for the six curation gaps found in the t2i4 review (P0.1 to P0.6).

Each class reproduces the reviewed behaviour on a synthetic fixture and asserts the corrected one:

* P0.1 a v1 sidecar runs through load -> engine -> report (and the CLI), with an explicit migration record;
* P0.2 a resolved control that changes inside one declared interval is a mismatch on that part;
* P0.3 evidence is covered in time only where its time applicability meets the assessed native support;
* P0.4 shared-resource aliases are found throughout the evidence ancestry, selectors keep containers apart;
* P0.5 an event never rewrites a raw integrity finding; exclusion is a separate policy result;
* P0.6 an update that cannot be resolved invalidates the array elements it may affect.

No IMAS, no real fixture, nothing here is experimental data.

    python examples/verify_curation_gaps.py
"""

import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import unittest.mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from curation import MANIFEST_SCHEMA, MANIFEST_SCHEMA_V1
from curation.engine import CurationEngine
from curation import intervals
from curation.inputs import declared_timeline
from curation.manifest import load_manifest, validate_manifest
from curation.reporting import write_report
from curation.source import NativeSource
import verify_curation as base
from verify_curation import Variable, run_screen, variables


def timeline(body):
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'TR.DAT'
        path.write_text(body)
        return declared_timeline(path)


def v2(**sections):
    return {'schema_version': MANIFEST_SCHEMA, **sections}


class IntervalHelperTests(unittest.TestCase):
    """The one interval helper behind control, evidence and reference scopes."""

    def test_partition_cuts_at_every_boundary_of_every_list(self):
        pieces = intervals.partition([[(None, 0.5), (0.5, None)], [(None, 0.2), (0.2, None)]], 0.0, 1.0)
        self.assertEqual(pieces, [(0.0, 0.2), (0.2, 0.5), (0.5, 1.0)])

    def test_boundaries_outside_the_support_do_not_cut_it(self):
        self.assertEqual(intervals.partition([[(None, 5.0), (5.0, None)]], 0.0, 1.0), [(0.0, 1.0)])
        self.assertEqual(intervals.partition([[(None, -1.0), (-1.0, None)]], 0.0, 1.0), [(0.0, 1.0)])

    def test_a_boundary_on_the_support_end_leaves_the_single_instant_piece(self):
        self.assertEqual(intervals.partition([[(None, 1.0), (1.0, None)]], 0.0, 1.0), [(0.0, 1.0), (1.0, 1.0)])

    def test_an_unknown_support_side_stays_unbounded(self):
        self.assertEqual(intervals.partition([[(None, 0.5), (0.5, None)]], None, None), [(None, 0.5), (0.5, None)])

    def test_covering_is_half_open(self):
        declared = [(None, 0.5), (0.5, None)]
        self.assertEqual(intervals.covering(declared, (0.0, 0.5)), 0)
        self.assertEqual(intervals.covering(declared, (0.5, 1.0)), 1)
        self.assertIsNone(intervals.covering([(0.6, None)], (0.0, 0.6)))             # a gap is not covered
        self.assertEqual(intervals.covering([(None, 0.3)], (None, 0.3)), 0)

    def test_support_intersects_window_and_native_spans(self):
        self.assertEqual(intervals.support({'time_min': 0.2}, [(0.0, 1.0), (0.0, 2.0)]), (0.2, 2.0))
        self.assertEqual(intervals.support({'time_max': 0.5}, [(0.0, 1.0)]), (0.0, 0.5))
        self.assertEqual(intervals.support({}, []), (None, None))

    def test_subtract_reports_the_uncovered_parts(self):
        self.assertEqual(intervals.subtract((0.0, 1.0), [(0.2, 0.4), (0.6, 0.7)]), [[0.0, 0.2], [0.4, 0.6], [0.7, 1.0]])
        self.assertEqual(intervals.subtract((0.0, 1.0), [(0.0, 1.0)]), [])

    def test_overlap_does_not_count_touching_endpoints(self):
        self.assertFalse(intervals.overlap((0.0, 1.0), (1.0, 2.0)))
        self.assertTrue(intervals.overlap((0.0, 1.0), (0.9, 2.0)))
        self.assertTrue(intervals.overlap((None, None), (5.0, 6.0)))


class TimelineConsumptionTests(unittest.TestCase):
    """The rules consume indexed and array values from the namelist adapter, and report the parser used."""

    def screen(self, body, manifest=None):
        return run_screen(manifest=manifest, timeline=timeline(body))

    def test_array_assignments_feed_the_heating_fraction_check_elementwise(self):
        _r, found = self.screen('FFULLA=0.5,0.6\nFHALFA=0.2,0.5\n')
        finding = found['HEATING.FRACTIONS']
        self.assertEqual(finding.status, 'fail')                                   # 0.6 + 0.5 > 1 for the second source
        declared = {d['full_key']: d['invalid'] for d in finding.evidence['declarations']}
        self.assertEqual(declared, {'FFULLA(1)': False, 'FFULLA(2)': True})

    def test_repeat_notation_resolves_each_element_for_reference_comparison(self):
        manifest = v2(input_expectations=[{'name': 'EINJA(2)', 'expected': 80.0, 'reference': 'r', 'units': 'keV'}])
        _r, found = self.screen('EINJA=3*80.0\n', manifest)
        self.assertEqual(found['HEATING.REFERENCE.0'].status, 'pass')
        _r, found = self.screen('EINJA=3*80.0\nEINJA(2)=70.0\n', manifest)        # a later element update is element-wise
        self.assertEqual(found['HEATING.REFERENCE.0'].status, 'warning')

    def test_resolved_controls_can_name_an_indexed_element(self):
        controls = [{'name': 'EINJA(3)', 'resolved': True, 'segments': [{'start': None, 'end': None, 'value': 80.0}]}]
        _r, found = self.screen('EINJA=3*80.0\n', v2(resolved_inputs={'controls': controls}))
        self.assertEqual(found['INPUT.RESOLVED.CONTROL.EINJA(3)'].status, 'pass')
        controls[0]['segments'][0]['value'] = 75.0
        _r, found = self.screen('EINJA=3*80.0\n', v2(resolved_inputs={'controls': controls}))
        self.assertEqual(found['INPUT.RESOLVED.CONTROL.EINJA(3)'].status, 'warning')

    def test_a_sparse_array_keeps_its_actual_indices(self):
        manifest = v2(input_expectations=[{'name': 'EINJA(3)', 'expected': 80.0, 'reference': 'r', 'units': 'keV'}])
        _r, found = self.screen('EINJA(3)=80.0\n', manifest)
        self.assertEqual(found['HEATING.REFERENCE.0'].status, 'pass')
        manifest = v2(input_expectations=[{'name': 'EINJA(1)', 'expected': 80.0, 'reference': 'r', 'units': 'keV'}])
        _r, found = self.screen('EINJA(3)=80.0\n', manifest)
        self.assertEqual(found['HEATING.REFERENCE.0'].status, 'not_assessed')      # element 1 was never declared

    def test_the_timeline_finding_reports_the_parser_and_what_it_leaves_unresolved(self):
        _r, found = self.screen('A(1)=1.0\n')
        evidence = found['INPUT.TIMELINE'].evidence
        self.assertEqual(evidence['parser']['kind'], 'f90nml')
        self.assertTrue(evidence['parser']['version'])
        self.assertIn('not the resolved TRANSP configuration', found['INPUT.TIMELINE'].message)
        _r, found = self.screen('A=1.0,2.0\n')                                     # an assumed lower bound stays reported
        self.assertEqual(found['INPUT.TIMELINE'].status, 'warning')
        self.assertTrue(found['INPUT.TIMELINE'].evidence['unresolved'])

    def test_the_limited_fallback_is_reported_and_never_presented_as_the_full_parser(self):
        import builtins
        real_import = builtins.__import__

        def refuse(name, *args, **kwargs):
            if name == 'f90nml':
                raise ImportError('mocked: f90nml unavailable')
            return real_import(name, *args, **kwargs)
        with unittest.mock.patch.object(builtins, '__import__', refuse):
            _r, found = self.screen('A(1)=1.0\n')
        finding = found['INPUT.TIMELINE']
        self.assertEqual(finding.evidence['parser']['kind'], 'scalar_fallback')
        self.assertEqual(finding.status, 'warning')
        self.assertIn('limited scalar reader', finding.message)


class V1MigrationTests(unittest.TestCase):
    """P0.1"""

    def write(self, tmp, document):
        path = Path(tmp) / 'evidence-v1.json'
        path.write_text(json.dumps(document))
        return path

    def test_a_minimal_v1_file_runs_through_load_engine_and_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(tmp, {'schema_version': MANIFEST_SCHEMA_V1})
            manifest = load_manifest(path)
            report, found = run_screen(manifest=manifest)
            self.assertEqual(found['RUN.IDENTITY'].status, 'not_assessed')
            paths, document = write_report(report, Path(tmp) / 'out')
            self.assertEqual(document['effective_manifest']['migration']['from'], MANIFEST_SCHEMA_V1)

    def test_a_populated_v1_file_keeps_its_content_through_the_engine(self):
        raw = {'schema_version': MANIFEST_SCHEMA_V1,
               'evidence': {'rotation': {'origin': 'measured', 'source': 'ref', 'roles': ['validation'],
                                         'uncertainty': 'provided'}},
               'input_expectations': [{'name': 'EINJA(1)', 'expected': 80.0, 'reference': 'r', 'units': 'keV'}],
               'thresholds': {'newton_iterations_warning': 10}}
        with tempfile.TemporaryDirectory() as tmp:
            manifest = load_manifest(self.write(tmp, raw))
        _report, found = run_screen(manifest=manifest, timeline=timeline('EINJA(1)=80.0\n'))
        self.assertEqual(found['HEATING.REFERENCE.0'].status, 'pass')
        self.assertEqual(found['EVIDENCE.rotation'].status, 'pass')
        self.assertEqual(manifest['evidence']['rotation']['usage'], ['independent_validation'])

    def test_the_migration_record_keeps_the_original_schema_and_raw_hash(self):
        raw = {'schema_version': MANIFEST_SCHEMA_V1}
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(tmp, raw)
            manifest = load_manifest(path)
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
        migration = manifest['migration']
        self.assertEqual((manifest['schema_version'], migration['from']), (MANIFEST_SCHEMA, MANIFEST_SCHEMA_V1))
        self.assertEqual(migration['raw_sha256'], digest)
        self.assertTrue(migration['steps'])

    def test_normalisation_is_idempotent_and_never_mutates_the_caller(self):
        raw = {'schema_version': MANIFEST_SCHEMA_V1, 'evidence': {'rotation': {'origin': 'assumed', 'roles': ['input']}}}
        before = copy.deepcopy(raw)
        once = validate_manifest(raw)
        self.assertEqual(raw, before)
        twice = validate_manifest(once)
        self.assertEqual(once, twice)
        self.assertEqual(validate_manifest(twice), twice)
        once['evidence']['rotation']['origin'] = 'x'                # the result is independent of the input
        self.assertEqual(raw, before)

    def test_v2_fields_in_a_raw_v1_document_are_rejected_by_policy(self):
        for key, value in (('events', []), ('run_identity', {}), ('resolved_inputs', {}), ('evidence_sources', [])):
            with self.assertRaises(ValueError, msg=key) as caught:
                validate_manifest({'schema_version': MANIFEST_SCHEMA_V1, key: value})
            self.assertIn('inputs.v1', str(caught.exception))

    def test_a_v2_document_still_works_and_has_no_migration_from_v1(self):
        manifest = validate_manifest(v2())
        self.assertEqual(manifest['migration']['from'], MANIFEST_SCHEMA)
        self.assertEqual(manifest['migration']['steps'], [])

    def test_the_real_cli_accepts_a_v1_sidecar_on_a_generated_cdf(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            base.ReportingAndCliTests().write_cdf(tmp / 'run.CDF')
            sidecar = self.write(tmp, {'schema_version': MANIFEST_SCHEMA_V1,
                                       'evidence': {'rotation': {'origin': 'assumed', 'source': 's'}}})
            before = hashlib.sha256((tmp / 'run.CDF').read_bytes()).hexdigest()
            done = subprocess.run([sys.executable, str(ROOT / 'curate.py'), '-i', str(tmp / 'run.CDF'),
                                   '--evidence', str(sidecar), '--output-dir', str(tmp / 'out')],
                                  capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stderr + done.stdout)
            document = json.loads((tmp / 'out' / 'curation.json').read_text())
            self.assertEqual(document['effective_manifest']['migration']['from'], MANIFEST_SCHEMA_V1)
            self.assertEqual(document['source']['evidence_file']['sha256'],
                             document['effective_manifest']['migration']['raw_sha256'])
            self.assertEqual(hashlib.sha256((tmp / 'run.CDF').read_bytes()).hexdigest(), before)


class ResolvedControlIntervalTests(unittest.TestCase):
    """P0.2"""

    def control(self, segments, name='MODEQ'):
        return v2(resolved_inputs={'controls': [{'name': name, 'resolved': True, 'segments': segments}]})

    def finding(self, manifest, body='MODEQ=1\n', data=None):
        _r, found = run_screen(data, manifest, timeline(body))
        return found['INPUT.RESOLVED.CONTROL.MODEQ']

    def test_a_runtime_change_inside_one_declared_interval_is_a_mismatch(self):
        finding = self.finding(self.control([{'start': None, 'end': 0.5, 'value': 1}, {'start': 0.5, 'end': None, 'value': 2}]))
        self.assertEqual(finding.status, 'warning')
        rows = finding.evidence['intervals']
        self.assertEqual([(r['start_time'], r['end_time'], r['agrees']) for r in rows],
                         [(0.0, 0.5, True), (0.5, 1.0, False)])
        self.assertEqual([r['resolved'] for r in rows], [1, 2])

    def test_agreement_over_the_whole_support_passes(self):
        finding = self.finding(self.control([{'start': None, 'end': None, 'value': 1}]))
        self.assertEqual(finding.status, 'pass')
        self.assertEqual(len(finding.evidence['intervals']), 1)

    def test_a_gap_in_the_resolved_segments_is_not_a_pass(self):
        finding = self.finding(self.control([{'start': None, 'end': 0.3, 'value': 1}, {'start': 0.6, 'end': None, 'value': 1}]))
        self.assertEqual(finding.status, 'not_assessed')
        gaps = [r for r in finding.evidence['intervals'] if r['resolved'] is None]
        self.assertEqual([(g['start_time'], g['end_time']) for g in gaps], [(0.3, 0.6)])

    def test_a_short_overlap_at_the_end_is_detected(self):
        finding = self.finding(self.control([{'start': None, 'end': 0.99, 'value': 1}, {'start': 0.99, 'end': None, 'value': 2}]))
        self.assertEqual(finding.status, 'warning')
        bad = [r for r in finding.evidence['intervals'] if not r['agrees']]
        self.assertEqual([(b['start_time'], b['end_time']) for b in bad], [(0.99, 1.0)])

    def test_a_change_exactly_on_a_declared_boundary_is_half_open(self):
        body = 'MODEQ=1\n~update_time=0.5\nMODEQ=2\n'
        finding = self.finding(self.control([{'start': None, 'end': 0.5, 'value': 1}, {'start': 0.5, 'end': None, 'value': 2}]), body)
        self.assertEqual(finding.status, 'pass')
        shifted = self.finding(self.control([{'start': None, 'end': 0.6, 'value': 1}, {'start': 0.6, 'end': None, 'value': 2}]), body)
        self.assertEqual(shifted.status, 'warning')                         # [0.5, 0.6) is declared 2 but ran 1

    def test_changes_outside_the_assessed_support_do_not_contaminate_it(self):
        manifest = self.control([{'start': None, 'end': 0.5, 'value': 1}, {'start': 0.5, 'end': None, 'value': 2}])
        manifest['assessment_window'] = {'time_min': 0.0, 'time_max': 0.4}
        self.assertEqual(self.finding(manifest).status, 'pass')
        beyond = self.control([{'start': None, 'end': 5.0, 'value': 1}, {'start': 5.0, 'end': None, 'value': 9}])
        self.assertEqual(self.finding(beyond).status, 'pass')               # the source ends at 1.0 s

    def test_an_unknown_initial_state_stays_unknown(self):
        body = '~update_time=0.5\nMODEQ=1\n'                                  # nothing declared before 0.5
        finding = self.finding(self.control([{'start': None, 'end': None, 'value': 1}]), body)
        self.assertEqual(finding.status, 'not_assessed')
        self.assertEqual(finding.evidence['intervals'][0]['declared'], None)

    def test_value_types_are_not_coerced_inside_the_partition(self):
        finding = self.finding(self.control([{'start': None, 'end': None, 'value': True}]), 'MODEQ=1\n')
        self.assertEqual(finding.status, 'warning')


class EvidenceTimeCoverageTests(unittest.TestCase):
    """P0.3"""

    def evidence(self, applicability=None, **extra):
        record = {'origin': 'measured', 'source': 'ref', 'uncertainty': 'provided', **extra}
        if applicability is not None:
            record['applicability'] = applicability
        return v2(evidence={'electron_temperature': record})

    def test_disjoint_evidence_supplies_no_coverage_and_is_preserved(self):
        _r, found = run_screen(manifest=self.evidence({'time_min': 10, 'time_max': 20}))
        for name in ('EVIDENCE.electron_temperature', 'UNCERTAINTY.electron_temperature'):
            self.assertEqual(found[name].status, 'not_assessed', name)
        coverage = found['EVIDENCE.electron_temperature'].evidence['coverage']
        self.assertEqual((coverage['state'], coverage['record_interval']), ('disjoint', [10.0, 20.0]))
        self.assertEqual(found['EVIDENCE.electron_temperature'].evidence['origin'], 'measured')   # record kept

    def test_partial_evidence_leaves_an_explicit_unassessed_interval(self):
        _r, found = run_screen(manifest=self.evidence({'time_min': 0.0, 'time_max': 0.4}))
        finding = found['EVIDENCE.electron_temperature']
        self.assertEqual(finding.status, 'warning')
        coverage = finding.evidence['coverage']
        self.assertEqual(coverage['state'], 'partial')
        self.assertEqual(coverage['unassessed_intervals'], [[0.4, 1.0]])
        self.assertEqual(found['UNCERTAINTY.electron_temperature'].evidence['coverage']['state'], 'partial')

    def test_full_coverage_and_no_time_bound_are_unchanged(self):
        for applicability in (None, {'time_min': -1.0, 'time_max': 5.0}):
            _r, found = run_screen(manifest=self.evidence(applicability))
            self.assertEqual(found['EVIDENCE.electron_temperature'].status, 'pass')
            self.assertEqual(found['EVIDENCE.electron_temperature'].evidence['coverage']['state'], 'full')

    def test_the_coverage_is_intersected_with_the_assessment_window(self):
        manifest = self.evidence({'time_min': 0.0, 'time_max': 0.4})
        manifest['assessment_window'] = {'time_min': 0.0, 'time_max': 0.3}
        _r, found = run_screen(manifest=manifest)
        self.assertEqual(found['EVIDENCE.electron_temperature'].evidence['coverage']['state'], 'full')
        manifest['assessment_window'] = {'time_min': 5.0, 'time_max': 6.0}
        _r, found = run_screen(manifest=manifest)
        self.assertEqual(found['EVIDENCE.electron_temperature'].evidence['coverage']['state'], 'unknown')

    def test_time_and_time3_are_not_interchangeable(self):
        data = variables()
        data['TIME3'] = Variable([0.0, 2.0], ('TIME3',), 'SECONDS')              # the profile clock reaches 2 s, TIME only 1 s
        manifest = self.evidence({'time_min': 0.0, 'time_max': 1.0})
        _r, found = run_screen(data, manifest)
        self.assertEqual(found['EVIDENCE.electron_temperature'].evidence['coverage']['state'], 'partial')
        self.assertEqual(found['EVIDENCE.electron_temperature'].evidence['coverage']['clock'], 'TIME3')

    def test_a_conflicting_declared_identity_cannot_authorize_a_reference(self):
        manifest = v2(run_identity={'device': 'D3D'}, input_expectations=[
            {'name': 'EINJA(1)', 'expected': 80.0, 'reference': 'r', 'units': 'keV', 'applicability': {'device': 'D3D'}}])
        _r, found = run_screen(manifest=manifest, timeline=timeline('EINJA(1)=80.0\n'), device='NSTX')
        self.assertEqual(found['RUN.IDENTITY'].status, 'fail')
        self.assertEqual(found['HEATING.REFERENCE.0'].status, 'not_assessed')   # neither side is trusted


class AliasLineageTests(unittest.TestCase):
    """P0.4"""

    def manifest(self, sources, constraint_parents, validator_parents, **extra):
        return v2(evidence_sources=sources, evidence={
            'electron_temperature': {'id': 'fit', 'origin': 'fitted', 'source': 'fit-a', 'parents': constraint_parents,
                                     'usage': ['constraint']},
            'ion_temperature': {'id': 'check', 'origin': 'measured', 'source': 'check-b', 'parents': validator_parents,
                                'usage': ['independent_validation'], **extra}})

    def test_same_id_ancestry_is_still_found(self):
        _r, found = run_screen(manifest=self.manifest([{'id': 'raw', 'description': 'd'}], ['raw'], ['raw']))
        shared = found['INDEPENDENCE.ion_temperature'].evidence['shared']
        self.assertEqual((shared[0]['kind'], shared[0]['via']), ('same_node', ['raw']))

    def test_differently_named_ancestors_with_one_hash_are_possible_shared_provenance(self):
        sources = [{'id': 'raw-a', 'description': 'd', 'sha256': 'a' * 64}, {'id': 'raw-b', 'description': 'd', 'sha256': 'a' * 64}]
        _r, found = run_screen(manifest=self.manifest(sources, ['raw-a'], ['raw-b']))
        finding = found['INDEPENDENCE.ion_temperature']
        self.assertEqual(finding.status, 'warning')
        kinds = {item['kind'] for item in finding.evidence['shared']}
        self.assertEqual(kinds, {'same_resource_selector_unknown'})
        self.assertIn('raw-a', json.dumps(finding.evidence['shared']))
        self.assertIn('raw-b', json.dumps(finding.evidence['shared']))

    def test_equal_selectors_on_aliased_hashes_are_the_same_resource(self):
        selector = {'dataset': 'TE', 'channel': 3}
        sources = [{'id': 'raw-a', 'description': 'd', 'sha256': 'a' * 64, 'selector': selector},
                   {'id': 'raw-b', 'description': 'd', 'sha256': 'a' * 64, 'selector': dict(selector)}]
        _r, found = run_screen(manifest=self.manifest(sources, ['raw-a'], ['raw-b']))
        self.assertEqual({i['kind'] for i in found['INDEPENDENCE.ion_temperature'].evidence['shared']}, {'same_resource'})

    def test_a_shared_container_with_different_channels_keeps_the_distinction(self):
        sources = [{'id': 'raw-a', 'description': 'd', 'sha256': 'a' * 64, 'selector': {'channel': 3}},
                   {'id': 'raw-b', 'description': 'd', 'sha256': 'a' * 64, 'selector': {'channel': 4}}]
        _r, found = run_screen(manifest=self.manifest(sources, ['raw-a'], ['raw-b']))
        finding = found['INDEPENDENCE.ion_temperature']
        self.assertEqual(finding.status, 'pass')
        self.assertEqual({i['kind'] for i in finding.evidence['shared']}, {'shared_container_different_selector'})
        self.assertIn('not proof of independence', finding.message)

    def test_undeclared_dependence_is_not_reported_as_independence(self):
        sources = [{'id': 'raw-a', 'description': 'd', 'sha256': 'a' * 64}, {'id': 'raw-b', 'description': 'd', 'sha256': 'b' * 64}]
        _r, found = run_screen(manifest=self.manifest(sources, ['raw-a'], ['raw-b']))
        self.assertNotIn('INDEPENDENCE.ion_temperature', found)

    def test_a_hash_alias_deep_in_the_ancestry_is_found(self):
        sources = [{'id': 'raw-a', 'description': 'd', 'sha256': 'a' * 64},
                   {'id': 'raw-b', 'description': 'd', 'sha256': 'a' * 64},
                   {'id': 'mid-a', 'description': 'd', 'parents': ['raw-a']},
                   {'id': 'mid-b', 'description': 'd', 'parents': ['raw-b']}]
        _r, found = run_screen(manifest=self.manifest(sources, ['mid-a'], ['mid-b']))
        self.assertEqual(found['INDEPENDENCE.ion_temperature'].status, 'warning')

    def test_usage_is_scoped_in_time_but_ancestry_is_kept(self):
        sources = [{'id': 'raw', 'description': 'd'}]
        manifest = self.manifest(sources, ['raw'], ['raw'])
        manifest['evidence']['electron_temperature']['applicability'] = {'time_min': 50.0, 'time_max': 60.0}
        report, found = run_screen(manifest=manifest)
        self.assertNotIn('INDEPENDENCE.ion_temperature', found)              # the constraint does not act in the window
        self.assertEqual(report.evidence_graph['raw']['children'], ['check', 'fit'])   # but the ancestry is retained

    def test_missing_references_and_cycles_stay_invalid(self):
        with self.assertRaises(ValueError):
            validate_manifest(v2(evidence_sources=[{'id': 'a', 'description': 'd', 'parents': ['zzz']}]))
        with self.assertRaises(ValueError):
            validate_manifest(v2(evidence_sources=[{'id': 'a', 'description': 'd', 'parents': ['b']},
                                                   {'id': 'b', 'description': 'd', 'parents': ['a']}]))


class EventPolicySeparationTests(unittest.TestCase):
    """P0.5"""

    def negative_te(self):
        data = variables()
        te = np.array(data['TE'].values, dtype=float)
        te[1, 1] = -5.0
        data['TE'] = Variable(te, ('TIME3', 'X'), 'EV')
        return data

    def event(self, rules, **extra):
        return {'id': 'pellet', 'type': 'pellet', 'start': 0.9, 'end': 1.1, 'rules': rules, **extra}

    def test_the_raw_finding_keeps_the_negative_sample(self):
        report, found = run_screen(self.negative_te(), v2(events=[self.event(['PROFILE.TE.NEGATIVE'])]))
        raw = found['PROFILE.TE.NEGATIVE']
        self.assertEqual(raw.status, 'fail')
        self.assertEqual(raw.scope['flagged_count'], 1)
        self.assertEqual(raw.scope['evaluated_count'], 6)
        self.assertEqual(report.masks['PROFILE_TE_NEGATIVE__flagged'].sum(), 1)

    def test_the_transport_policy_can_still_exclude_it_as_a_separate_result(self):
        report, found = run_screen(self.negative_te(), v2(events=[self.event(['PROFILE.TE.NEGATIVE'])]))
        record = report.policy_exclusions['PROFILE.TE.NEGATIVE']
        self.assertEqual((record['excluded_count'], record['flagged_excluded_count'], record['events']), (3, 1, ['pellet']))
        self.assertTrue(record['covers_all_flagged'])
        self.assertEqual(record['eligible_count'], 3)
        self.assertEqual(report.masks['PROFILE_TE_NEGATIVE__excluded'].sum(), 3)
        summary = report.as_dict()['summary']['applications']
        self.assertEqual(summary['heat_transport'], 'review_required')                 # no longer blocked: only unrelated warnings
        self.assertEqual(found['PROFILE.TE.NEGATIVE'].status, 'fail')                  # while the finding itself stays
        self.assertEqual(report.as_dict()['summary']['excused_findings']['heat_transport'], ['PROFILE.TE.NEGATIVE'])

    def test_without_an_event_the_application_is_blocked(self):
        report, _found = run_screen(self.negative_te())
        self.assertEqual(report.as_dict()['summary']['applications']['heat_transport'], 'blocked_by_screening_findings')

    def test_an_event_that_does_not_cover_every_flagged_sample_does_not_excuse_the_finding(self):
        data = self.negative_te()
        te = np.array(data['TE'].values, dtype=float)
        te[0, 2] = -1.0                                                          # a second negative at t = 0, outside the event
        data['TE'] = Variable(te, ('TIME3', 'X'), 'EV')
        report, _found = run_screen(data, v2(events=[self.event(['PROFILE.TE.NEGATIVE'])]))
        self.assertFalse(report.policy_exclusions['PROFILE.TE.NEGATIVE']['covers_all_flagged'])
        self.assertEqual(report.as_dict()['summary']['applications']['heat_transport'], 'blocked_by_screening_findings')

    def test_when_every_sample_is_excluded_the_policy_is_unassessed_not_the_raw_finding(self):
        event = self.event(['PROFILE.TE.NEGATIVE'], start=-1.0, end=5.0)
        report, found = run_screen(self.negative_te(), v2(events=[event]))
        self.assertEqual(found['PROFILE.TE.NEGATIVE'].status, 'fail')
        record = report.policy_exclusions['PROFILE.TE.NEGATIVE']
        self.assertEqual((record['eligible_count'], record['policy_status']), (0, 'not_assessed'))

    def test_mandatory_source_integrity_rules_cannot_be_waived_by_a_wildcard(self):
        data = variables()
        data['TIME3'] = Variable([1.0, 0.0], ('TIME3',), 'SECONDS')              # a reversed clock
        report, found = run_screen(data, v2(events=[self.event(['*'], start=-5.0, end=5.0)]))
        self.assertEqual(found['GRID.TIME3'].status, 'fail')
        self.assertNotIn('GRID.TIME3', report.policy_exclusions)
        inventory = found['EVENTS.INVENTORY'].evidence
        self.assertIn('GRID.TIME3', inventory['events'][0]['mandatory_rules_not_waivable'])

    def test_missingness_and_event_membership_stay_distinct(self):
        data = variables()
        te = np.ma.array(np.array(data['TE'].values, dtype=float), mask=np.zeros((2, 3), bool))
        te.mask[1, 0] = True
        data['TE'] = Variable(te, ('TIME3', 'X'), 'EV')
        report, found = run_screen(data, v2(events=[self.event(['PROFILE.TE.FINITE'])]))
        self.assertEqual(found['PROFILE.TE.FINITE'].status, 'warning')           # the masked sample is still reported
        record = report.policy_exclusions['PROFILE.TE.FINITE']
        self.assertEqual((record['flagged_excluded_count'], record['excluded_count']), (1, 3))

    def test_legacy_event_rules_are_recorded_as_a_legacy_analyst_policy(self):
        manifest = validate_manifest(v2(events=[self.event(['PROFILE.TE.*'])]))
        event = manifest['events'][0]
        self.assertEqual(event['policy'], {'kind': 'legacy_exclusion', 'rules': ['PROFILE.TE.*']})
        self.assertEqual(event['origin'], 'analyst_declaration')


class ArrayInvalidationTests(unittest.TestCase):
    """P0.6"""

    def values(self, body):
        return timeline(body)['states'][-1]['values']

    def test_an_unsupported_whole_array_update_invalidates_every_prior_element(self):
        values = self.values('A(1)=1\nA(2)=2\nA=$unsupported\n')
        self.assertNotIn('A(1)', values)
        self.assertNotIn('A(2)', values)

    def test_a_malformed_update_invalidates_the_array_and_records_the_affected_keys(self):
        tl = timeline('A(1)=1\nA(2)=2\nB(1)=5\nA=$unsupported\n')
        gap = tl['gaps'][0]
        self.assertEqual(sorted(gap['invalidated']), ['A(1)', 'A(2)'])
        self.assertEqual(tl['states'][-1]['values'], {'B(1)': 5})                   # unaffected arrays are kept

    def test_an_unsupported_slice_with_known_bounds_leaves_other_elements(self):
        values = self.values('A(1)=1\nA(2)=2\nA(3)=3\nA(2:3)=$unsupported\n')
        self.assertEqual(values, {'A(1)': 1})

    def test_an_unsupported_slice_with_unknown_bounds_invalidates_the_base_array(self):
        values = self.values('A(1)=1\nA(2)=2\nA(:)=$unsupported\n')
        self.assertEqual(values, {})

    def test_an_unsupported_element_update_invalidates_only_that_element(self):
        values = self.values('A(1)=1\nA(2)=2\nA(2)=$unsupported\n')
        self.assertEqual(values, {'A(1)': 1})

    def test_a_later_valid_assignment_after_invalidation_is_current(self):
        values = self.values('A(1)=1\nA=$unsupported\nA(1)=7\n')
        self.assertEqual(values, {'A(1)': 7})

    def test_an_explicit_zero_and_an_invalidated_element_are_distinct(self):
        values = self.values('A(1)=0\nA(2)=2\nA(2)=$unsupported\n')
        self.assertEqual(values, {'A(1)': 0})
        self.assertNotIn('A(2)', values)

    def test_a_stale_value_is_not_used_by_the_heating_fraction_check(self):
        body = 'FFULLA(1)=0.9\nFHALFA(1)=0.9\nFFULLA=$unsupported\nFHALFA(1)=0.05\n'
        _r, found = run_screen(timeline=timeline(body))
        self.assertEqual(found['HEATING.FRACTIONS'].status, 'not_assessed')      # the 0.9 full fraction is gone


if __name__ == '__main__':
    unittest.main(verbosity=1)
