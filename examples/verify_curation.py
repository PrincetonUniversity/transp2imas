#!/usr/bin/python
"""Pure and classic-CDF integration tests for the independent curation framework.

No IMAS, netCDF4, or external fixture is required. SciPy is used only to generate
small classic CDFs for CLI tests. Run: python3 examples/verify_curation.py
"""

import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from curation import MANIFEST_SCHEMA, MANIFEST_SCHEMA_V1
from curation.engine import CurationEngine
from curation.inputs import declared_timeline, literal
from curation.manifest import load_manifest, validate_manifest
from curation.pilot import PILOT_SCHEMA, build_manifest, load_pilot
from curation.reporting import write_report
from curation.source import NativeSource, density_scale


class Variable:
    def __init__(self, values, dimensions, units='', **attrs):
        self.values = np.ma.array(values, copy=True) if np.ma.isMaskedArray(values) else np.asarray(values)
        self.shape = self.values.shape
        self.dimensions = dimensions
        self.units = units
        self.long_name = 'Synthetic fixture - not measured data'
        for key, value in attrs.items():
            setattr(self, key, value)

    def __getitem__(self, index):
        return self.values[index]


def variables():
    x = np.tile([0.1, 0.5, 0.9], (2, 1))
    return {'TIME': Variable([0., 0.4, 1.], ('TIME',), 'SECONDS'),
            'TIME3': Variable([0., 1.], ('TIME3',), 'SECONDS'),
            'X': Variable(x, ('TIME3', 'X')),
            'NE': Variable(np.full_like(x, 1e13), ('TIME3', 'X'), 'N/CM**3'),
            'TE': Variable(1000 - 500 * x, ('TIME3', 'X'), 'EV'),
            'TI': Variable(1500 - 700 * x, ('TIME3', 'X'), 'EV'),
            'NPTRITER': Variable([0, 4, 2], ('TIME',))}


def run_screen(data=None, manifest=None, timeline=None, **attributes):
    source = NativeSource(types.SimpleNamespace(variables=data if data is not None else variables(), **attributes))
    manifest = manifest or {'schema_version': MANIFEST_SCHEMA}
    report = CurationEngine(source, manifest, timeline).run()
    return report, {f.id: f for f in report.findings}


class NativeSourceTests(unittest.TestCase):
    def test_masked_and_fill_values_are_missing_not_physical_numbers(self):
        for var in (Variable(np.ma.array([1, -999], mask=[False, True]), ('T',)),
                    Variable([1, -999], ('T',), _FillValue=-999)):
            source = NativeSource(types.SimpleNamespace(variables={'A': var}))
            values = source.array('A')
            self.assertEqual(values[0], 1)
            self.assertTrue(np.isnan(values[1]))

    def test_cache_returns_read_only_values_and_reads_once(self):
        source = NativeSource(types.SimpleNamespace(variables={'A': Variable([1., 2.], ('T',))}))
        first = source.array('A')
        self.assertIs(source.array('A'), first)
        self.assertEqual(source.read_counts, {'A': 1})
        with self.assertRaises(ValueError):
            first[0] = -1

    def test_lru_respects_byte_budget(self):
        source = NativeSource(types.SimpleNamespace(variables={k: Variable([1., 2.], ('T',)) for k in 'AB'}),
                              cache_bytes=16)
        source.array('A'); source.array('B'); source.array('A')
        self.assertEqual(source.read_counts['A'], 2)

    def test_large_variable_is_not_read(self):
        source = NativeSource(types.SimpleNamespace(variables={'A': Variable([1., 2.], ('T',))}), max_variable_bytes=8)
        with self.assertRaises(ValueError):
            source.array('A')
        self.assertEqual(source.read_counts, {})

    def test_packed_fill_is_masked_before_scaling(self):
        var = Variable([10, -99], ('T',), _FillValue=-99, scale_factor=.1, add_offset=5)
        source = NativeSource(types.SimpleNamespace(variables={'A': var}), packed=True)
        np.testing.assert_allclose(source.array('A'), [6, np.nan], equal_nan=True)

    def test_valid_range_is_not_exposed_as_valid_data(self):
        var = Variable([-1, 2, 11], ('T',), valid_range=[0, 10])
        source = NativeSource(types.SimpleNamespace(variables={'A': var}))
        np.testing.assert_allclose(source.array('A'), [np.nan, 2, np.nan], equal_nan=True)

    def test_scalar_missing_values(self):
        source = NativeSource(types.SimpleNamespace(variables={'A': Variable(-9., (), missing_value=-9.)}))
        self.assertTrue(np.isnan(source.array('A')))

    def test_units_are_not_inferred_from_magnitude(self):
        self.assertEqual(density_scale('N/CM**3'), 1e6)
        self.assertEqual(density_scale('m^-3'), 1)
        with self.assertRaises(ValueError):
            density_scale('')


class ManifestTests(unittest.TestCase):
    def test_duplicate_json_keys_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'evidence.json'
            path.write_text('{"schema_version":"' + MANIFEST_SCHEMA + '","evidence":{},"evidence":{}}')
            with self.assertRaisesRegex(ValueError, 'duplicate JSON key'):
                load_manifest(path)

    def test_unknown_keys_rejected(self):
        with self.assertRaisesRegex(ValueError, 'unknown keys'):
            validate_manifest({'schema_version': MANIFEST_SCHEMA, 'evidnce': {}})

    def test_boolean_numbers_and_nonfinite_thresholds_rejected(self):
        for value in (True, float('nan'), -1):
            with self.assertRaises(ValueError):
                validate_manifest({'schema_version': MANIFEST_SCHEMA, 'thresholds': {'weak_gradient_relative': value}})

    def test_bad_window_is_not_silently_fixed(self):
        with self.assertRaises(ValueError):
            validate_manifest({'schema_version': MANIFEST_SCHEMA, 'assessment_window': {'time_min': 2, 'time_max': 1}})

    def test_duplicate_species_and_electron_reuse_rejected(self):
        for species in ([{'density': 'ND', 'charge': 1}] * 2, [{'density': 'NE', 'charge': 1}]):
            with self.assertRaises(ValueError):
                validate_manifest({'schema_version': MANIFEST_SCHEMA, 'composition': {'species': species}})

    def test_residual_threshold_requires_normalization_rationale(self):
        with self.assertRaises(ValueError):
            validate_manifest({'schema_version': MANIFEST_SCHEMA, 'thresholds': {'residual_limits': {'RESL2TE': {'limit': 1e-3}}}})

    def test_point_uncertainties_not_accepted_as_derivative_uncertainties(self):
        with self.assertRaises(ValueError):
            validate_manifest({'schema_version': MANIFEST_SCHEMA, 'gradient_uncertainties': {
                'TE': {'sigma_variable': 'TE_SIGMA', 'coordinate': 'X', 'method': 'pointwise_error', 'source': 'test'}}})

    def test_input_manifest_is_not_mutated(self):
        original = {'schema_version': MANIFEST_SCHEMA}
        validate_manifest(original)
        self.assertEqual(original, {'schema_version': MANIFEST_SCHEMA})


class ManifestV2Tests(unittest.TestCase):
    """Manifest v2: strict new sections, v1 compatibility, and a validated evidence graph."""

    def v2(self, **sections):
        return {'schema_version': MANIFEST_SCHEMA, **sections}

    def rejects(self, document, fragment):
        with self.assertRaises(ValueError) as caught:
            validate_manifest(document)
        self.assertIn(fragment, str(caught.exception))

    def test_a_v1_document_stays_valid_and_gains_empty_v2_sections(self):
        result = validate_manifest({'schema_version': MANIFEST_SCHEMA_V1,
                                    'evidence': {'rotation': {'origin': 'measured', 'roles': ['validation']}}})
        self.assertEqual(result['events'], [])
        self.assertEqual(result['evidence']['rotation']['usage'], ['independent_validation'])

    def test_a_v1_document_may_not_use_v2_keys(self):
        self.rejects({'schema_version': MANIFEST_SCHEMA_V1, 'events': []}, 'need schema_version')
        self.rejects({'schema_version': MANIFEST_SCHEMA_V1,
                      'evidence': {'rotation': {'origin': 'measured', 'id': 'x'}}}, 'need schema_version')

    def test_legacy_roles_map_to_one_usage_vocabulary(self):
        result = validate_manifest(self.v2(evidence={'heating': {
            'origin': 'measured', 'roles': ['input', 'calibration'], 'independent_validation': True}}))
        self.assertEqual(result['evidence']['heating']['usage'],
                         ['calibration_target', 'constraint', 'independent_validation'])

    def test_unknown_usage_hash_and_identifier_are_rejected(self):
        self.rejects(self.v2(evidence={'rotation': {'origin': 'measured', 'usage': ['proof']}}), 'usage')
        self.rejects(self.v2(evidence={'rotation': {'origin': 'measured', 'sha256': 'abc'}}), 'SHA-256')
        self.rejects(self.v2(evidence={'rotation': {'origin': 'measured', 'id': 'has space'}}), 'identifier')

    def test_parents_must_exist_and_may_not_form_a_cycle(self):
        self.rejects(self.v2(evidence_sources=[{'id': 'a', 'description': 'x', 'parents': ['ghost']}]), 'ghost')
        self.rejects(self.v2(evidence_sources=[{'id': 'a', 'description': 'x', 'parents': ['b']},
                                               {'id': 'b', 'description': 'y', 'parents': ['a']}]), 'cycle')
        self.rejects(self.v2(evidence_sources=[{'id': 'a', 'description': 'x'}, {'id': 'a', 'description': 'y'}]),
                     'duplicate')

    def test_an_override_needs_value_reason_and_reviewer(self):
        for override in ({'value': 1.0, 'reason': 'r'}, {'value': 1.0, 'reviewer': 'p'}, {'reason': 'r', 'reviewer': 'p'},
                         {'value': float('nan'), 'reason': 'r', 'reviewer': 'p'}):
            self.rejects(self.v2(input_expectations=[{'name': 'X', 'expected': 1.0, 'override': override}]), 'override')
        validate_manifest(self.v2(input_expectations=[{'name': 'X', 'expected': 1.0, 'reference': 'r', 'units': '1',
                                                       'override': {'value': 2.0, 'reason': 'r', 'reviewer': 'p'}}]))

    def test_applicability_ranges_must_be_ordered(self):
        self.rejects(self.v2(input_expectations=[{'name': 'X', 'expected': 1.0,
                                                  'applicability': {'shot_min': 5, 'shot_max': 1}}]), 'shot_min')

    def test_events_need_a_type_a_window_and_rules(self):
        base = {'id': 'e1', 'type': 'sawtooth', 'start': 1.0, 'end': 1.1, 'rules': ['GRADIENT.*']}
        validate_manifest(self.v2(events=[base]))
        self.rejects(self.v2(events=[{**base, 'type': 'quake'}]), 'type')
        self.rejects(self.v2(events=[{**base, 'start': 2.0}]), 'exceeds')
        self.rejects(self.v2(events=[{**base, 'rules': []}]), 'rules')
        self.rejects(self.v2(events=[base, base]), 'duplicate')

    def test_controls_must_be_ordered_nonoverlapping_scalars_that_state_their_kind(self):
        def control(segments, **extra):
            return self.v2(resolved_inputs={'controls': [{'name': 'MODEQ', 'resolved': True,
                                                          'segments': segments, **extra}]})
        validate_manifest(control([{'start': 0, 'end': 1, 'value': 1}, {'start': 1, 'end': None, 'value': 2}]))
        self.rejects(control([{'start': 0, 'end': 2, 'value': 1}, {'start': 1, 'end': 3, 'value': 2}]), 'overlap')
        self.rejects(control([{'start': 0, 'end': 1, 'value': [1, 2]}]), 'scalar')
        self.rejects(control([{'start': 1, 'end': 0, 'value': 1}]), 'before')
        self.rejects(self.v2(resolved_inputs={'controls': [{'name': 'M', 'segments': [{'value': 1}]}]}), 'resolved')

    def test_duplicate_controls_are_rejected_after_normalisation(self):
        segments = [{'start': 0, 'end': 1, 'value': 1}]
        self.rejects(self.v2(resolved_inputs={'controls': [
            {'name': 'MODEQ', 'resolved': True, 'segments': segments},
            {'name': 'mod eq', 'resolved': True, 'segments': segments}]}), 'duplicate')


class EngineV2Tests(unittest.TestCase):
    """Events, run identity, evidence lineage, applicability, overrides and resolved inputs in the engine."""

    def v2(self, **sections):
        return {'schema_version': MANIFEST_SCHEMA, **sections}

    def timeline(self, body):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'TR.DAT'
            path.write_text(body)
            return declared_timeline(path)

    def negative_te(self):
        data = variables()
        te = np.array(data['TE'].values, dtype=float)
        te[1, 1] = -5.0
        data['TE'] = Variable(te, ('TIME3', 'X'), 'EV')
        return data

    # events
    def test_a_negative_sample_inside_a_declared_event_is_excluded_and_accounted_for(self):
        _r, found = run_screen(self.negative_te())
        self.assertEqual(found['PROFILE.TE.NEGATIVE'].status, 'fail')
        event = {'id': 'saw1', 'type': 'sawtooth', 'start': 0.9, 'end': 1.1, 'rules': ['PROFILE.TE.NEGATIVE']}
        report, found = run_screen(self.negative_te(), self.v2(events=[event]))
        finding = found['PROFILE.TE.NEGATIVE']
        self.assertEqual(finding.status, 'pass')
        self.assertEqual(finding.scope['excluded_by_events'], 3)                  # the whole TIME3=1.0 row
        self.assertEqual(finding.scope['excluding_events'], ['saw1'])
        self.assertEqual(report.event_exclusions['PROFILE.TE.NEGATIVE']['excluded_count'], 3)
        self.assertEqual(found['PROFILE.TE.ZERO'].scope.get('excluded_by_events', 0), 0)        # other rules untouched
        self.assertEqual(found['EVENTS.INVENTORY'].status, 'pass')

    def test_excluded_samples_are_not_counted_as_evaluated_or_unflagged(self):
        event = {'id': 'e', 'type': 'elm', 'start': 0.9, 'end': 1.1, 'rules': ['PROFILE.TE.*']}
        report, found = run_screen(self.negative_te(), self.v2(events=[event]))
        scope = found['PROFILE.TE.FINITE'].scope
        self.assertEqual((scope['evaluated_count'], scope['excluded_by_events']), (3, 3))
        self.assertTrue(report.masks['PROFILE_TE_FINITE__excluded'][1].all())
        self.assertFalse(report.masks['PROFILE_TE_FINITE__evaluated'][1].any())

    def test_when_every_sample_is_excluded_the_rule_is_not_assessed(self):
        event = {'id': 'all', 'type': 'diagnostic_outage', 'start': -1.0, 'end': 5.0, 'rules': ['PROFILE.TE.NEGATIVE']}
        _r, found = run_screen(self.negative_te(), self.v2(events=[event]))
        self.assertEqual(found['PROFILE.TE.NEGATIVE'].status, 'not_assessed')
        self.assertIn('inside declared events', found['PROFILE.TE.NEGATIVE'].message)

    def test_an_event_matching_no_rule_or_outside_the_source_is_reported(self):
        events = [{'id': 'typo', 'type': 'pellet', 'start': 0.0, 'end': 1.0, 'rules': ['NO.SUCH.RULE']},
                  {'id': 'late', 'type': 'pellet', 'start': 50.0, 'end': 51.0, 'rules': ['PROFILE.TE.*']}]
        report, found = run_screen(manifest=self.v2(events=events))
        self.assertEqual(found['EVENTS.INVENTORY'].status, 'warning')
        self.assertEqual(sorted(found['EVENTS.INVENTORY'].evidence['needs_review']), ['late', 'typo'])

    def test_no_events_means_unassessed_not_none_found(self):
        _r, found = run_screen()
        self.assertEqual(found['EVENTS.INVENTORY'].status, 'not_assessed')

    # identity
    def test_run_identity_is_checked_against_the_source_attributes(self):
        manifest = self.v2(run_identity={'device': 'nstx', 'shot': 141711})
        _r, found = run_screen(manifest=manifest, device='NSTX', shot=np.int32(141711))
        self.assertEqual(found['RUN.IDENTITY'].status, 'pass')
        _r, found = run_screen(manifest=manifest, device='NSTX', shot=np.int32(141712))
        self.assertEqual(found['RUN.IDENTITY'].status, 'fail')
        _r, found = run_screen(manifest=manifest)
        self.assertEqual(found['RUN.IDENTITY'].status, 'not_assessed')
        _r, found = run_screen()
        self.assertEqual(found['RUN.IDENTITY'].status, 'not_assessed')

    # lineage
    def evidence(self, **extra):
        return {'origin': 'measured', 'source': 'ref', **extra}

    def test_a_validator_sharing_an_ancestor_with_a_constraint_is_not_independent(self):
        manifest = self.v2(
            evidence_sources=[{'id': 'efit_raw', 'description': 'EFIT reconstruction'}],
            evidence={'electron_temperature': self.evidence(id='te_input', usage=['constraint'], parents=['efit_raw'],
                                                           source='a'),
                      'ion_temperature': self.evidence(id='ti_check', usage=['independent_validation'],
                                                       parents=['efit_raw'], source='b')})
        _r, found = run_screen(manifest=manifest)
        finding = found['INDEPENDENCE.ion_temperature']
        self.assertEqual(finding.status, 'warning')
        self.assertEqual(finding.evidence['shared'][0]['via'], ['efit_raw'])

    def test_a_validator_with_separate_lineage_raises_no_warning(self):
        manifest = self.v2(
            evidence_sources=[{'id': 'a_raw', 'description': 'x'}, {'id': 'b_raw', 'description': 'y'}],
            evidence={'electron_temperature': self.evidence(id='te', usage=['constraint'], parents=['a_raw'], source='a'),
                      'ion_temperature': self.evidence(id='ti', usage=['independent_validation'], parents=['b_raw'],
                                                       source='b')})
        report, found = run_screen(manifest=manifest)
        self.assertNotIn('INDEPENDENCE.ion_temperature', found)
        self.assertEqual(report.evidence_graph['a_raw']['children'], ['te'])

    def test_a_shared_hash_defeats_a_claim_of_independence(self):
        digest = 'ab' * 32
        manifest = self.v2(evidence={'electron_temperature': self.evidence(usage=['constraint'], sha256=digest, source='a'),
                                     'ion_temperature': self.evidence(usage=['independent_validation'], sha256=digest,
                                                                      source='b')})
        _r, found = run_screen(manifest=manifest)
        self.assertEqual(found['INDEPENDENCE.ion_temperature'].status, 'warning')

    # applicability and overrides
    def expectation(self, **extra):
        return {'name': 'EINJA(1)', 'expected': 80.0, 'reference': 'sheet', 'units': 'keV', **extra}

    def test_a_reference_for_another_shot_range_is_not_applicable(self):
        body = 'EINJA(1)=70.0\n'
        manifest = self.v2(run_identity={'shot': 141711}, input_expectations=[
            self.expectation(applicability={'shot_min': 120000, 'shot_max': 130000})])
        _r, found = run_screen(manifest=manifest, timeline=self.timeline(body))
        self.assertEqual(found['HEATING.REFERENCE.0'].status, 'not_applicable')
        manifest['input_expectations'][0]['applicability'] = {'shot_min': 141000, 'shot_max': 142000}
        _r, found = run_screen(manifest=manifest, timeline=self.timeline(body))
        self.assertEqual(found['HEATING.REFERENCE.0'].status, 'warning')               # applies, and 70 != 80

    def test_an_unknown_shot_makes_applicability_unknown_not_true(self):
        manifest = self.v2(input_expectations=[self.expectation(applicability={'shot_min': 1})])
        _r, found = run_screen(manifest=manifest, timeline=self.timeline('EINJA(1)=80.0\n'))
        self.assertEqual(found['HEATING.REFERENCE.0'].status, 'not_assessed')

    def test_a_reviewed_override_replaces_the_value_and_keeps_the_original(self):
        override = {'value': 70.0, 'reason': 'campaign used 70 keV', 'reviewer': 'physicist'}
        manifest = self.v2(input_expectations=[self.expectation(override=override)])
        _r, found = run_screen(manifest=manifest, timeline=self.timeline('EINJA(1)=70.0\n'))
        finding = found['HEATING.REFERENCE.0']
        self.assertEqual(finding.status, 'pass')
        self.assertTrue(finding.evidence['override_applied'])
        self.assertEqual(finding.evidence['expectation']['expected'], 80.0)
        self.assertIn('physicist', finding.message)

    def test_a_reference_time_window_restricts_the_compared_intervals(self):
        body = 'EINJA(1)=80.0\n~update_time=1.0\nEINJA(1)=10.0\n'
        manifest = self.v2(input_expectations=[self.expectation(applicability={'time_max': 0.5})])
        _r, found = run_screen(manifest=manifest, timeline=self.timeline(body))
        self.assertEqual(found['HEATING.REFERENCE.0'].status, 'pass')                  # the 10 keV interval is not compared

    def test_evidence_for_another_device_does_not_apply(self):
        manifest = self.v2(run_identity={'device': 'NSTX'}, evidence={'rotation': self.evidence(
            applicability={'device': 'DIII-D'})})
        _r, found = run_screen(manifest=manifest)
        self.assertEqual(found['EVIDENCE.rotation'].status, 'not_applicable')

    # resolved inputs
    def test_resolved_controls_are_compared_with_the_declarations(self):
        controls = [{'name': 'NLMDIF', 'resolved': True, 'segments': [{'start': None, 'end': None, 'value': False}]}]
        manifest = self.v2(resolved_inputs={'controls': controls})
        _r, found = run_screen(manifest=manifest, timeline=self.timeline('NLMDIF=.TRUE.\n'))
        self.assertEqual(found['INPUT.RESOLVED.CONTROL.NLMDIF'].status, 'warning')
        _r, found = run_screen(manifest=manifest, timeline=self.timeline('NLMDIF=.FALSE.\n'))
        self.assertEqual(found['INPUT.RESOLVED.CONTROL.NLMDIF'].status, 'pass')

    def test_a_declared_only_control_cannot_confirm_the_run(self):
        controls = [{'name': 'NLMDIF', 'resolved': False, 'segments': [{'value': True}]}]
        _r, found = run_screen(manifest=self.v2(resolved_inputs={'controls': controls}),
                               timeline=self.timeline('NLMDIF=.TRUE.\n'))
        self.assertEqual(found['INPUT.RESOLVED.CONTROL.NLMDIF'].status, 'not_assessed')

    def test_a_value_of_another_type_is_not_coerced_to_agree(self):
        controls = [{'name': 'LPREDICT_NE', 'resolved': True, 'segments': [{'value': True}]}]
        _r, found = run_screen(manifest=self.v2(resolved_inputs={'controls': controls}),
                               timeline=self.timeline('LPREDICT_NE=1\n'))
        self.assertEqual(found['INPUT.RESOLVED.CONTROL.LPREDICT_NE'].status, 'warning')

    def test_input_file_hashes_are_verified(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'TR.DAT'
            path.write_text('A=1\n')
            good = hashlib.sha256(path.read_bytes()).hexdigest()
            for digest, expected in ((good, 'pass'), ('00' * 32, 'fail'), (None, 'not_assessed')):
                item = {'role': 'trdat', 'path': str(path)}
                if digest:
                    item['sha256'] = digest
                _r, found = run_screen(manifest=self.v2(resolved_inputs={'inputs': [item]}))
                self.assertEqual(found['INPUT.RESOLVED.FILES'].status, expected, digest)
            gone = {'role': 'trdat', 'path': str(Path(tmp) / 'absent'), 'sha256': good}
            _r, found = run_screen(manifest=self.v2(resolved_inputs={'inputs': [gone]}))
            self.assertEqual(found['INPUT.RESOLVED.FILES'].status, 'error')

    def test_the_code_version_is_compared_with_the_source(self):
        manifest = self.v2(resolved_inputs={'code': {'name': 'TRANSP', 'version': '21.X.X'}})
        _r, found = run_screen(manifest=manifest, TRANSP_version='21.X.X')
        self.assertEqual(found['INPUT.RESOLVED.CODE'].status, 'pass')
        _r, found = run_screen(manifest=manifest, TRANSP_version='20.1')
        self.assertEqual(found['INPUT.RESOLVED.CODE'].status, 'warning')
        _r, found = run_screen(manifest=manifest)
        self.assertEqual(found['INPUT.RESOLVED.CODE'].status, 'not_assessed')

    def test_the_source_arrays_are_not_modified_by_the_v2_checks(self):
        data = self.negative_te()
        before = {k: np.array(v.values, copy=True) for k, v in data.items()}
        event = {'id': 'e', 'type': 'elm', 'start': 0.9, 'end': 1.1, 'rules': ['PROFILE.*']}
        run_screen(data, self.v2(events=[event]))
        for name, values in before.items():
            np.testing.assert_array_equal(data[name].values, values, name)


class DeclaredInputTests(unittest.TestCase):
    def parse(self, body):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'TR.DAT'
            path.write_text(body)
            return declared_timeline(path)

    def test_d_exponents_and_quoted_bang(self):
        self.assertEqual(literal('1.0D-3'), .001)
        plan = self.parse("A=1D-3\nFILE='dir!name.nml' ! comment\n")
        self.assertEqual(plan['states'][0]['values'], {'A': .001, 'FILE': 'dir!name.nml'})

    def test_timeline_preserves_both_model_states(self):
        plan = self.parse('LPREDICT_NE=0\nNLMDIF=.TRUE.\n~update_time=0.5\nLPREDICT_NE=1\n')
        self.assertEqual([s['values']['LPREDICT_NE'] for s in plan['states']], [0, 1])
        self.assertEqual(plan['states'][0]['end_time'], .5)

    def test_unsupported_new_value_does_not_reuse_old_value(self):
        plan = self.parse('A=1\n~update_time=.5\nA=1+2\n')
        self.assertNotIn('A', plan['states'][1]['values'])
        self.assertTrue(plan['gaps'])

    def test_includes_and_multiple_assignments_are_reported(self):
        plan = self.parse("include 'other.nml'\nA=1, B=2\n")
        self.assertEqual(len(plan['gaps']), 2)
        self.assertNotIn('A', plan['states'][0]['values'])

    def test_nonincreasing_updates_disable_temporal_rules(self):
        plan = self.parse('~update_time=2\n~update_time=1\n')
        self.assertFalse(plan['temporal_valid'])

    def test_a_group_header_may_share_its_line_with_an_assignment(self):
        plan = self.parse('&TRIN LPREDICT_NE=1\nLPREDICT_NMAIN=0\n')
        self.assertEqual(plan['states'][0]['values'], {'LPREDICT_NE': 1.0, 'LPREDICT_NMAIN': 0.0})
        self.assertEqual(plan['gaps'], [])
        plan = self.parse('&INNML\nA=1\n/\n')                         # a bare header and terminator are still skipped
        self.assertEqual(plan['states'][0]['values'], {'A': 1.0})
        self.assertEqual(plan['gaps'], [])

    def test_missing_input_is_not_all_false_controls(self):
        self.assertFalse(declared_timeline()['available'])
        self.assertEqual(declared_timeline()['states'], [])


class ScreeningTests(unittest.TestCase):
    def test_native_time_grids_and_unknown_evidence_remain_distinct(self):
        report, findings = run_screen()
        self.assertEqual(findings['NUMERICS.ITERATION_DOMAIN'].status, 'pass')
        self.assertEqual(findings['NUMERICS.CONVERGENCE'].status, 'not_assessed')
        self.assertEqual(findings['UNCERTAINTY.rotation'].status, 'not_assessed')
        self.assertEqual(len(report.coordinates['coordinate_TIME']), 3)
        self.assertEqual(len(report.coordinates['coordinate_TIME3']), 2)
        self.assertFalse(report.as_dict()['summary']['scientifically_certified'])

    def test_present_but_unreadable_variable_is_an_error_not_missing_evidence(self):
        class BrokenVariable(Variable):
            def __getitem__(self, index):
                raise OSError('synthetic read failure')
        data = variables(); data['TE'] = BrokenVariable(np.ones((2, 3)), ('TIME3', 'X'), 'EV')
        _, findings = run_screen(data)
        self.assertEqual(findings['PROFILE.TE.ASSESSMENT'].status, 'error')

    def test_explicit_numeric_one_in_a_logical_switch_is_recognized(self):
        timeline = {'available': True, 'states': [{'start_time': None, 'end_time': None,
                    'values': {'LPREDICT_NE': 1, 'NLMDIF': 1}, 'source_lines': {}}], 'gaps': []}
        _, findings = run_screen(timeline=timeline)
        self.assertEqual(findings['MODEL.DENSITY_CURRENT_COUPLING'].status, 'warning')

    def test_nonmonotonic_clock_fails(self):
        data = variables(); data['TIME'] = Variable([0, 1, .5], ('TIME',), 'SECONDS')
        _, findings = run_screen(data)
        self.assertEqual(findings['GRID.TIME'].status, 'fail')

    def test_negative_profile_is_not_clipped(self):
        data = variables(); data['TE'].values[1, 1] = -10
        report, findings = run_screen(data)
        self.assertEqual(findings['PROFILE.TE.NEGATIVE'].status, 'fail')
        self.assertEqual(data['TE'].values[1, 1], -10)
        self.assertTrue(report.masks['PROFILE_TE_NEGATIVE__flagged'][1, 1])

    def test_nonfinite_values_are_unavailable_not_negative_density(self):
        data = variables(); data['NE'].values[0, 1] = -999; data['NE']._FillValue = -999
        _, findings = run_screen(data)
        self.assertEqual(findings['PROFILE.NE.FINITE'].status, 'warning')
        self.assertEqual(findings['PROFILE.NE.NEGATIVE'].status, 'pass')

    def test_zero_density_is_an_interval_review_not_a_fabricated_correction(self):
        data = variables(); data['NE'].values[0] = 0
        _, findings = run_screen(data)
        self.assertEqual(findings['PROFILE.NE.ZERO'].status, 'warning')

    def test_radial_shape_mismatch_detected(self):
        data = variables(); data['NE'] = Variable(np.ones((2, 2)), ('TIME3', 'X'), 'm^-3')
        _, findings = run_screen(data)
        self.assertEqual(findings['PROFILE.NE.GRID_MISMATCH'].status, 'fail')

    def test_equal_lengths_do_not_authorize_time_coordinate_relabeling(self):
        data = variables(); data['TE'].dimensions = ('TIME', 'X')
        _, findings = run_screen(data)
        self.assertEqual(findings['PROFILE.TE.ASSESSMENT'].status, 'not_assessed')

    def test_window_limits_exact_native_sample_masks(self):
        manifest = {'schema_version': MANIFEST_SCHEMA,
                    'assessment_window': {'time_min': .5, 'rho_min': .3, 'rho_max': .6}}
        report, _ = run_screen(manifest=manifest)
        expected = np.array([[False, False, False], [False, True, False]])
        np.testing.assert_array_equal(report.masks['PROFILE_TE_FINITE__evaluated'], expected)

    def test_empty_window_is_not_pass(self):
        _, findings = run_screen(manifest={'schema_version': MANIFEST_SCHEMA,
                                          'assessment_window': {'time_min': 99}})
        self.assertEqual(findings['PROFILE.TE.FINITE'].status, 'not_assessed')

    def test_flat_profile_flags_coefficient_identifiability_not_bad_profile(self):
        data = variables(); data['TE'].values[:] = 1000
        _, findings = run_screen(data)
        self.assertEqual(findings['FIT.TE.WEAK_GRADIENT'].status, 'warning')
        self.assertEqual(findings['PROFILE.TE.NEGATIVE'].status, 'pass')

    def test_supplied_derivative_uncertainty_is_used_without_inventing_covariance(self):
        data = variables(); data['DTE_SIGMA'] = Variable(np.full((2, 3), 600), ('TIME3', 'X'), 'EV')
        manifest = {'schema_version': MANIFEST_SCHEMA, 'gradient_uncertainties': {
            'TE': {'sigma_variable': 'DTE_SIGMA', 'coordinate': 'X',
                   'method': 'numpy_gradient_edge_order_1', 'source': 'synthetic derivative ensemble'}}}
        _, findings = run_screen(data, manifest)
        self.assertEqual(findings['FIT.TE.GRADIENT_UNCERTAINTY'].status, 'warning')

    def test_measured_zero_rotation_is_not_automatically_rejected(self):
        data = variables(); data['OMEGA'] = Variable(np.zeros((2, 3)), ('TIME3', 'X'), 'RAD/SEC')
        manifest = {'schema_version': MANIFEST_SCHEMA, 'evidence': {'rotation': {
            'origin': 'measured', 'source': 'synthetic diagnostic', 'uncertainty': 'provided', 'reliability': 'reviewed'}}}
        _, findings = run_screen(data, manifest)
        self.assertEqual(findings['ROTATION.AVAILABILITY'].status, 'pass')
        self.assertEqual(findings['EVIDENCE.rotation'].status, 'pass')

    def test_reused_calibration_source_cannot_claim_independence(self):
        manifest = {'schema_version': MANIFEST_SCHEMA, 'evidence': {'heating': {
            'origin': 'manual', 'source': 'synthetic neutron observation', 'roles': ['calibration', 'validation'],
            'independent_validation': True}}}
        _, findings = run_screen(manifest=manifest)
        self.assertEqual(findings['INDEPENDENCE.heating'].status, 'warning')

    def test_density_and_magnetic_diffusion_trigger_review_not_prohibition(self):
        timeline = {'available': True, 'states': [{'start_time': None, 'end_time': None,
                    'values': {'LPREDICT_NE': 1, 'NLMDIF': True}, 'source_lines': {}}], 'gaps': []}
        _, findings = run_screen(timeline=timeline)
        self.assertEqual(findings['MODEL.DENSITY_CURRENT_COUPLING'].status, 'warning')
        self.assertNotEqual(findings['MODEL.DENSITY_CURRENT_COUPLING'].status, 'fail')

    def test_current_fraction_sum_above_one_fails(self):
        timeline = {'available': True, 'states': [{'start_time': None, 'end_time': None,
                    'values': {'FFULLA(1)': .7, 'FHALFA(1)': .5}, 'source_lines': {}}], 'gaps': []}
        _, findings = run_screen(timeline=timeline)
        self.assertEqual(findings['HEATING.FRACTIONS'].status, 'fail')

    def test_configuration_reference_deviation_warns_without_modification(self):
        timeline = {'available': True, 'states': [{'start_time': None, 'end_time': None,
                    'values': {'EINJA(1)': 90000}, 'source_lines': {'EINJA(1)': 2}}], 'gaps': []}
        manifest = {'schema_version': MANIFEST_SCHEMA, 'input_expectations': [{'name': 'EINJA(1)',
                     'expected': 80000, 'units': 'eV', 'reference': 'synthetic reference', 'atol': 1000}]}
        _, findings = run_screen(manifest=manifest, timeline=timeline)
        self.assertEqual(findings['HEATING.REFERENCE.0'].status, 'warning')
        self.assertEqual(timeline['states'][0]['values']['EINJA(1)'], 90000)

    def test_residual_threshold_needs_explicit_rationale(self):
        data = variables(); data['RESL2TE'] = Variable([.1, .01, .001], ('TIME',))
        manifest = {'schema_version': MANIFEST_SCHEMA, 'thresholds': {'residual_limits': {
            'RESL2TE': {'limit': .02, 'rationale': 'synthetic regression tolerance, not a TRANSP default'}}}}
        _, findings = run_screen(data, manifest)
        self.assertEqual(findings['NUMERICS.RESL2TE.THRESHOLD'].status, 'warning')
        self.assertEqual(findings['NUMERICS.CONVERGENCE'].status, 'not_assessed')


class DeclaredStateReportingTests(unittest.TestCase):
    def timeline(self, body):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'TR.DAT'
            path.write_text(body)
            return declared_timeline(path)

    def test_a_reference_name_with_spaces_matches_the_declaration(self):
        manifest = {'schema_version': MANIFEST_SCHEMA,
                    'input_expectations': [{'name': 'EINJA ( 1 )', 'expected': 80.0, 'reference': 'campaign sheet',
                                            'units': 'keV'}]}
        _report, found = run_screen(manifest=manifest, timeline=self.timeline('EINJA(1)=80.0\n'))
        self.assertEqual(found['HEATING.REFERENCE.0'].status, 'pass')
        manifest['input_expectations'][0]['expected'] = 70.0
        _report, found = run_screen(manifest=manifest, timeline=self.timeline('EINJA(1)=80.0\n'))
        self.assertEqual(found['HEATING.REFERENCE.0'].status, 'warning')

    def test_unusable_update_times_say_the_state_rules_were_not_evaluated(self):
        body = 'LPREDICT_NE=1\nNLMDIF=.TRUE.\n'
        _report, found = run_screen(timeline=self.timeline(body))
        self.assertEqual(found['MODEL.DENSITY_CURRENT_COUPLING'].status, 'warning')       # usable timeline: coupling seen
        _report, found = run_screen(timeline=self.timeline('~update_time=2\n~update_time=1\n' + body))
        for identifier in ('MODEL.DENSITY_CURRENT_COUPLING', 'MODEL.ROTATION_CONTROLS', 'HEATING.FRACTIONS'):
            self.assertEqual(found[identifier].status, 'not_assessed', identifier)
            self.assertIn('Not evaluated', found[identifier].message, identifier)
            self.assertIn('INPUT.TIMELINE', found[identifier].message, identifier)
        self.assertEqual(found['INPUT.TIMELINE'].status, 'warning')


class CorruptSourceTests(unittest.TestCase):
    def test_undecodable_metadata_is_an_error_not_missing_evidence(self):
        data = variables()
        data['TE'] = Variable(data['TE'].values, ('TIME3', 'X'), 'EV', valid_range='not-a-range')
        _report, found = run_screen(data=data)
        self.assertEqual(found['PROFILE.TE.ASSESSMENT'].status, 'error')
        self.assertIn('valid_range', found['PROFILE.TE.ASSESSMENT'].message)

    def test_unsupported_units_stay_not_assessed(self):
        data = variables()
        data['NE'] = Variable(data['NE'].values, ('TIME3', 'X'), 'FURLONGS')
        _report, found = run_screen(data=data)
        statuses = {i: f.status for i, f in found.items() if i.startswith('PROFILE.NE')}
        self.assertNotIn('error', statuses.values(), statuses)


class CompositionTests(unittest.TestCase):
    def example(self):
        data = variables()
        data['ND'] = Variable(np.full((2, 3), .8e19), ('TIME3', 'X'), 'm^-3')
        data['NIMP_TEST'] = Variable(np.full((2, 3), .1e13), ('TIME3', 'X'), 'N/CM**3')
        manifest = {'schema_version': MANIFEST_SCHEMA, 'composition': {
            'complete': True, 'includes_fast_ions': True, 'source': 'synthetic closed inventory',
            'species': [{'density': 'ND', 'charge': 1}, {'density': 'NIMP_TEST', 'charge': 2}]}}
        return data, manifest

    def test_charge_closure_in_mixed_density_units(self):
        data, manifest = self.example()
        _, findings = run_screen(data, manifest)
        self.assertEqual(findings['COMPOSITION.QUASINEUTRALITY'].status, 'pass')

    def test_charge_mismatch_warns(self):
        data, manifest = self.example(); data['ND'].values[:] *= .5
        _, findings = run_screen(data, manifest)
        self.assertEqual(findings['COMPOSITION.QUASINEUTRALITY'].status, 'warning')

    def test_partial_inventory_never_claims_closure(self):
        data, manifest = self.example(); manifest['composition']['includes_fast_ions'] = False
        _, findings = run_screen(data, manifest)
        self.assertEqual(findings['COMPOSITION.QUASINEUTRALITY'].status, 'not_assessed')

    def test_charge_profile_keeps_native_dimensions(self):
        data, manifest = self.example()
        data['Z_TEST'] = Variable(np.full((2, 3), 2.), ('TIME3', 'X'))
        manifest['composition']['species'][1] = {'density': 'NIMP_TEST', 'charge_variable': 'Z_TEST'}
        _, findings = run_screen(data, manifest)
        self.assertEqual(findings['COMPOSITION.QUASINEUTRALITY'].status, 'pass')

    def test_unknown_density_units_are_not_guessed(self):
        data, manifest = self.example(); data['ND'].units = ''
        _, findings = run_screen(data, manifest)
        self.assertEqual(findings['COMPOSITION.QUASINEUTRALITY'].status, 'not_assessed')

    def test_negative_species_not_included_as_a_valid_charge_term(self):
        data, manifest = self.example(); data['ND'].values[:] = -1
        _, findings = run_screen(data, manifest)
        self.assertEqual(findings['COMPOSITION.QUASINEUTRALITY.INPUTS'].status, 'warning')
        self.assertEqual(findings['COMPOSITION.QUASINEUTRALITY'].status, 'not_assessed')


class ReportingAndCliTests(unittest.TestCase):
    def write_cdf(self, path):
        from scipy.io import netcdf_file
        with netcdf_file(path, 'w') as nc:
            nc.createDimension('TIME', 3); nc.createDimension('TIME3', 2); nc.createDimension('X', 3)
            for name, var in variables().items():
                output = nc.createVariable(name, 'd', var.dimensions)
                output[:] = var.values
                output.units = var.units

    def test_report_serializes_without_nonstandard_nan_json(self):
        report, _ = run_screen()
        with tempfile.TemporaryDirectory() as tmp:
            paths, data = write_report(report, tmp)
            json.loads(paths['json'].read_text(), parse_constant=lambda s: self.fail(s))
            with np.load(paths['npz'], allow_pickle=False) as arrays:
                self.assertIn('coordinate_TIME', arrays)
                self.assertEqual(arrays['PROFILE_TE_FINITE__evaluated'].dtype, bool)
            self.assertIn('not a scientific certification', paths['md'].read_text())

    def test_reports_get_the_default_file_permissions(self):
        report, _ = run_screen()
        previous = os.umask(0o022)
        try:
            with tempfile.TemporaryDirectory() as tmp:
                paths, _ = write_report(report, tmp)
                for name, path in paths.items():
                    self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o644, name)
        finally:
            os.umask(previous)

    def test_a_rendering_failure_leaves_the_previous_report_untouched(self):
        report, _ = run_screen()
        with tempfile.TemporaryDirectory() as tmp:
            paths, _ = write_report(report, tmp)
            before = {name: path.read_bytes() for name, path in paths.items()}
            with mock.patch('curation.reporting.markdown_report', side_effect=KeyError('boom')):
                with self.assertRaises(KeyError):
                    write_report(report, tmp)
            self.assertEqual({name: path.read_bytes() for name, path in paths.items()}, before)

    def test_an_unexpected_data_error_exits_2_not_1(self):
        import curate
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'test.CDF'
            self.write_cdf(path)
            for error in (KeyError('X'), IndexError('x'), AttributeError('x')):
                with mock.patch('curate.CurationEngine') as engine:
                    engine.return_value.run.side_effect = error
                    self.assertEqual(curate.main(['-i', str(path), '--output-dir', str(Path(tmp) / 'out')]), 2,
                                     type(error).__name__)

    def test_report_cannot_replace_input_sidecar(self):
        report, _ = run_screen()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'curation.json'; path.write_text('important input')
            with self.assertRaises(ValueError):
                write_report(report, tmp, protected_paths=[path])
            self.assertEqual(path.read_text(), 'important input')

    def test_classic_cdf_cli_and_source_hash_preservation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'test.CDF'; self.write_cdf(path)
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            output = Path(tmp) / 'report'
            result = subprocess.run([sys.executable, str(ROOT / 'curate.py'), '-i', str(path),
                                     '--output-dir', str(output)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), digest)
            self.assertEqual(json.loads((output / 'curation.json').read_text())['source']['sha256'], digest)

    def test_require_complete_changes_exit_status_not_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'test.CDF'; self.write_cdf(path)
            result = subprocess.run([sys.executable, str(ROOT / 'curate.py'), '-i', str(path), '--require-complete',
                                     '--output-dir', str(Path(tmp) / 'report')], capture_output=True, text=True)
            self.assertEqual(result.returncode, 1, result.stderr)

    def test_bad_manifest_is_an_execution_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'test.CDF'; self.write_cdf(path)
            evidence = Path(tmp) / 'bad.json'; evidence.write_text('{"schema_version":"wrong"}')
            result = subprocess.run([sys.executable, str(ROOT / 'curate.py'), '-i', str(path),
                                     '--evidence', str(evidence), '--output-dir', str(Path(tmp) / 'report')],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertFalse((Path(tmp) / 'report').exists())


class PilotAdapterTests(unittest.TestCase):
    """The pilot adapter computes hashes, passes supplied facts through, and refuses an unfilled template."""

    def facts(self, **extra):
        return {'schema_version': PILOT_SCHEMA, 'run_identity': {'device': 'X', 'shot': 1}, **extra}

    def test_the_template_is_refused_until_it_is_filled(self):
        with self.assertRaises(ValueError) as caught:
            load_pilot(ROOT / 'examples' / 'curation' / 'pilot-template.json')
        self.assertIn('TO_BE_SUPPLIED', str(caught.exception))
        self.assertIn('facts.run_identity.device', str(caught.exception))

    def test_file_hashes_are_computed_from_the_files_not_supplied(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / 'TR.DAT').write_text('A=1\n')
            manifest = build_manifest(self.facts(files=[{'role': 'trdat', 'path': 'TR.DAT'}]), tmp)
            item = manifest['resolved_inputs']['inputs'][0]
            self.assertEqual(item['sha256'], hashlib.sha256(b'A=1\n').hexdigest())
            self.assertEqual(manifest['schema_version'], MANIFEST_SCHEMA)
            with self.assertRaises(ValueError):                     # a supplied hash would be a claim, not a computation
                build_manifest(self.facts(files=[{'role': 'trdat', 'path': 'TR.DAT', 'sha256': '0' * 64}]), tmp)
            with self.assertRaises(OSError):
                build_manifest(self.facts(files=[{'role': 'trdat', 'path': 'absent'}]), tmp)

    def test_unknown_sections_and_wrong_schema_are_rejected(self):
        with self.assertRaises(ValueError):
            build_manifest(self.facts(surprise=1), '.')
        with self.assertRaises(ValueError):
            build_manifest({'schema_version': 'other'}, '.')

    def test_the_result_is_validated_like_any_manifest(self):
        with self.assertRaises(ValueError):
            build_manifest(self.facts(events=[{'id': 'e', 'type': 'bogus', 'start': 0, 'end': 1, 'rules': ['*']}]), '.')

    def test_the_cli_runs_a_pilot_and_leaves_the_inputs_untouched(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            ReportingAndCliTests().write_cdf(tmp / 'run.CDF')
            (tmp / 'TR.DAT').write_text('NLMDIF=.TRUE.\n')
            facts = self.facts(files=[{'role': 'trdat', 'path': 'TR.DAT'}],
                               events=[{'id': 'e', 'type': 'elm', 'start': 0.0, 'end': 0.1, 'rules': ['PROFILE.TE.*']}])
            (tmp / 'pilot.json').write_text(json.dumps(facts))
            before = {n: hashlib.sha256((tmp / n).read_bytes()).hexdigest() for n in ('run.CDF', 'TR.DAT', 'pilot.json')}
            done = subprocess.run([sys.executable, str(ROOT / 'curate.py'), '-i', str(tmp / 'run.CDF'),
                                   '--trdat', str(tmp / 'TR.DAT'), '--pilot', str(tmp / 'pilot.json'),
                                   '--output-dir', str(tmp / 'out')], capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stderr + done.stdout)
            self.assertEqual(before, {n: hashlib.sha256((tmp / n).read_bytes()).hexdigest() for n in before})
            document = json.loads((tmp / 'out' / 'curation.json').read_text())
            self.assertEqual(document['source']['evidence_file']['kind'], 'pilot_facts')
            self.assertEqual(document['events'][0]['id'], 'e')
            self.assertIn('INPUT.RESOLVED.FILES', {f['id'] for f in document['findings']})
            done = subprocess.run([sys.executable, str(ROOT / 'curate.py'), '-i', str(tmp / 'run.CDF'), '--pilot',
                                   str(ROOT / 'examples' / 'curation' / 'pilot-template.json'),
                                   '--output-dir', str(tmp / 'out2')], capture_output=True, text=True)
            self.assertEqual(done.returncode, 2)
            self.assertFalse((tmp / 'out2').exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
