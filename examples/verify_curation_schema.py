#!/usr/bin/python
"""Tests of the offline JSON Schemas and the validation pipeline of curation documents.

Stages under test (each one distinct, and none needs IMAS or the network):

    strict JSON decoding -> declared-version structure (JSON Schema) -> migration and canonical normalization
    -> resource/reference validation -> semantic validation

    python examples/verify_curation_schema.py
"""

import copy
import hashlib
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import jsonschema

from curation import (ACCEPTED_MANIFEST_SCHEMAS, MANIFEST_SCHEMA, MANIFEST_SCHEMA_V1, MANIFEST_SCHEMA_V2,
                      REPORT_SCHEMA)
from curation import schema
from curation.manifest import load_manifest, referenced_paths, validate_manifest, validate_resources
from curation.pilot import PILOT_SCHEMA, build_manifest, load_pilot
from curation.policy import load_policy
from curation.reporting import write_report
import verify_curation as base
from verify_curation import run_screen

V1, V2, V3 = MANIFEST_SCHEMA_V1, MANIFEST_SCHEMA_V2, MANIFEST_SCHEMA


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class SchemaFileTests(unittest.TestCase):
    def test_every_packaged_schema_is_valid_draft_2020_12(self):
        for path in sorted(schema.SCHEMA_DIR.glob('*.schema.json')):
            document = json.loads(path.read_text())
            jsonschema.Draft202012Validator.check_schema(document)
            self.assertEqual(document['$schema'], 'https://json-schema.org/draft/2020-12/schema', path.name)
            self.assertTrue(document['$id'].startswith(schema.URN), path.name)

    def test_the_expected_schemas_ship(self):
        names = {p.name.replace('.schema.json', '') for p in schema.SCHEMA_DIR.glob('*.schema.json')}
        self.assertEqual(names, {'common', *schema.SCHEMA_NAMES})

    def test_references_are_local_and_validation_needs_no_network(self):
        for path in schema.SCHEMA_DIR.glob('*.schema.json'):
            for ref in [m for m in json.dumps(json.loads(path.read_text())).split('"$ref": "')[1:]]:
                self.assertTrue(ref.startswith(schema.URN), (path.name, ref[:60]))

        def refuse(*args, **kwargs):
            raise AssertionError('network access attempted during curation validation')
        with mock.patch.object(socket, 'socket', refuse), mock.patch.object(socket, 'getaddrinfo', refuse):
            schema._cache.clear()
            validate_manifest({'schema_version': V3})
            schema._cache.clear()

    def test_schema_validation_does_not_import_imas(self):
        done = subprocess.run([sys.executable, '-c', 'import sys; sys.path.insert(0, %r); import curation.manifest, '
                               'curation.reporting, curate; print("imas" in sys.modules)' % str(ROOT)],
                              capture_output=True, text=True)
        self.assertEqual(done.stdout.strip(), 'False', done.stderr)


class StrictDecodingTests(unittest.TestCase):
    def test_duplicate_keys_nan_infinity_and_overflow_are_rejected_at_decoding(self):
        for text, fragment in (('{"a": 1, "a": 2}', 'duplicate JSON key'), ('{"a": NaN}', 'nonfinite JSON constant'),
                               ('{"a": Infinity}', 'nonfinite'), ('{"a": -Infinity}', 'nonfinite'),
                               ('{"a": 1e999}', 'numeric overflow'), ('{"a": -1e999}', 'numeric overflow')):
            with self.assertRaises(ValueError, msg=text) as caught:
                schema.loads_strict(text)
            self.assertIn(fragment, str(caught.exception))

    def test_a_finite_number_and_a_large_exponent_that_fits_are_accepted(self):
        self.assertEqual(schema.loads_strict('{"a": 1e300, "b": 3}'), {'a': 1e300, 'b': 3})

    def test_programmatically_supplied_nonfinite_scalars_are_rejected(self):
        for bad in (float('nan'), float('inf'), float('-inf')):
            with self.assertRaises(ValueError, msg=repr(bad)) as caught:
                validate_manifest({'schema_version': V3, 'thresholds': {'weak_gradient_relative': bad}})
            self.assertIn('nonfinite', str(caught.exception))
            with self.assertRaises(ValueError):
                validate_manifest({'schema_version': V3, 'input_expectations': [
                    {'name': 'X', 'expected': bad, 'reference': 'r', 'units': '1'}]})

    def test_a_file_with_a_duplicate_key_or_overflow_is_refused_by_the_loader(self):
        with tempfile.TemporaryDirectory() as tmp:
            for text in ('{"schema_version": "%s", "notes": "a", "notes": "b"}' % V3,
                         '{"schema_version": "%s", "thresholds": {"weak_gradient_relative": 1e999}}' % V3):
                path = Path(tmp) / 'e.json'
                path.write_text(text)
                with self.assertRaises(ValueError):
                    load_manifest(path)


class StructuralValidationTests(unittest.TestCase):
    def reject(self, document, *fragments):
        with self.assertRaises(ValueError) as caught:
            validate_manifest(document)
        for fragment in fragments:
            self.assertIn(fragment, str(caught.exception))

    def test_an_error_names_schema_path_value_type_and_condition(self):
        self.reject({'schema_version': V3, 'events': [{'id': 'e', 'type': 'sawtooth', 'start': 'soon', 'end': 1.0,
                                                       'rules': ['*']}]},
                    '[schema evidence.v3]', 'events[0].start', "'soon'", 'type str', "failed 'type'")

    def test_unknown_scientific_keys_are_rejected_everywhere_they_can_occur(self):
        for document in ({'schema_version': V3, 'evidnce': {}},
                         {'schema_version': V3, 'thresholds': {'weak_gradient_relativ': 1.0}},
                         {'schema_version': V3, 'evidence': {'rotation': {'origin': 'measured', 'sorce': 'x'}}},
                         {'schema_version': V3, 'composition': {'complet': True}},
                         {'schema_version': V3, 'evidence': {'not_a_role': {'origin': 'measured'}}},
                         {'schema_version': V3, 'events': [{'id': 'e', 'type': 'elm', 'start': 0, 'end': 1, 'rules': ['*'],
                                                           'bogus': 1}]}):
            with self.assertRaises(ValueError, msg=str(document)):
                validate_manifest(document)

    def test_shot_identity_is_an_integer_and_values_are_not_coerced(self):
        validate_manifest({'schema_version': V3, 'run_identity': {'shot': 141711}})
        for shot in (141711.5, '141711', True, -1):
            self.reject({'schema_version': V3, 'run_identity': {'shot': shot}}, 'run_identity.shot')
        self.reject({'schema_version': V3, 'thresholds': {'weak_gradient_relative': '0.1'}}, 'thresholds.weak_gradient_relative')
        self.reject({'schema_version': V3, 'thresholds': {'weak_gradient_relative': True}}, 'thresholds.weak_gradient_relative')
        self.reject({'schema_version': V3, 'thresholds': {'weak_gradient_relative': -1.0}}, 'thresholds.weak_gradient_relative')

    def test_the_extension_namespace_is_deliberate(self):
        manifest = validate_manifest({'schema_version': V3, 'extensions': {'vendor.my_rule': {'limit': 3}}})
        self.assertEqual(manifest['extensions'], {'vendor.my_rule': {'limit': 3}})
        for bad in ({'plain': {}}, {'vendor.rule': 3}, {'Vendor.Rule': {}}):
            self.reject({'schema_version': V3, 'extensions': bad}, 'extensions')

    def test_event_origins_and_detector_settings_are_validated(self):
        event = {'id': 'e', 'type': 'unclassified_transient', 'start': 0.1, 'end': 0.2, 'rules': ['PROFILE.*'],
                 'origin': 'automatic_candidate', 'review_state': 'unreviewed', 'timing_uncertainty': 0.01,
                 'detector': {'method': 'abrupt_drop', 'version': '1.0', 'settings': {'threshold': 5.0}}}
        self.assertEqual(validate_manifest({'schema_version': V3, 'events': [event]})['events'][0]['origin'],
                         'automatic_candidate')
        self.reject({'schema_version': V3, 'events': [{**event, 'origin': 'oracle'}]}, 'events[0].origin')
        self.reject({'schema_version': V3, 'events': [{**event, 'timing_uncertainty': -1.0}]}, 'events[0].timing_uncertainty')
        self.reject({'schema_version': V3, 'events': [{**event, 'review_state': 'maybe'}]}, 'events[0].review_state')
        detector = {'id': 'd1', 'method': 'abrupt_drop', 'diagnostic': 'TE0', 'clock': 'TIME', 'mode': 'causal',
                    'expected_units': ['EV'], 'max_sample_interval': 0.01}
        validate_manifest({'schema_version': V3, 'event_detectors': [detector]})
        self.reject({'schema_version': V3, 'event_detectors': [{**detector, 'mode': 'whenever'}]}, 'event_detectors[0].mode')
        self.reject({'schema_version': V3, 'event_detectors': [detector, detector]}, 'duplicate detector id')
        self.reject({'schema_version': V3, 'event_detectors': [{**detector, 'parameters': {'window': 4, 'min_baseline': 8}}]},
                    'min_baseline must be smaller than window')

    def test_the_resource_selector_structure_is_validated(self):
        record = {'origin': 'measured', 'selector': {'dataset': 'TE', 'channel': 3}}
        validate_manifest({'schema_version': V3, 'evidence': {'electron_temperature': record}})
        for selector in ({}, {'dataset': ['TE']}, {'bad key': 1}, 'TE'):
            self.reject({'schema_version': V3, 'evidence': {'electron_temperature': {**record, 'selector': selector}}},
                        'selector')

    def test_semantic_errors_are_not_reported_as_schema_failures(self):
        for document, fragment in (
                ({'schema_version': V3, 'evidence_sources': [{'id': 'a', 'description': 'd', 'parents': ['b']},
                                                             {'id': 'b', 'description': 'd', 'parents': ['a']}]}, 'cycle'),
                ({'schema_version': V3, 'evidence_sources': [{'id': 'a', 'description': 'd', 'parents': ['zzz']}]}, 'zzz'),
                ({'schema_version': V3, 'resolved_inputs': {'controls': [{'name': 'M', 'resolved': True, 'segments': [
                    {'start': 0, 'end': 2, 'value': 1}, {'start': 1, 'end': 3, 'value': 2}]}]}}, 'overlap'),
                ({'schema_version': V3, 'composition': {'species': [{'density': 'ND', 'charge': 1}, {'density': 'ND', 'charge': 1}]}}, 'unique'),
                ({'schema_version': V3, 'assessment_window': {'time_min': 2, 'time_max': 1}}, 'exceeds')):
            with self.assertRaises(ValueError) as caught:
                validate_manifest(document)
            self.assertIn(fragment, str(caught.exception))
            self.assertNotIn('[schema', str(caught.exception))             # structurally valid, semantically not


class MigrationAndNormalizationTests(unittest.TestCase):
    def documents(self):
        return {V1: {'schema_version': V1, 'evidence': {'rotation': {'origin': 'measured', 'source': 's', 'roles': ['validation']}},
                     'input_expectations': [{'name': 'EINJA(1)', 'expected': 80.0, 'reference': 'r', 'units': 'keV'}]},
                V2: {'schema_version': V2, 'run_identity': {'device': 'NSTX', 'shot': 141711},
                     'events': [{'id': 'e', 'type': 'pellet', 'start': 0.1, 'end': 0.2, 'rules': ['PROFILE.*']}],
                     'evidence': {'rotation': {'origin': 'measured', 'source': 's', 'usage': ['independent_validation']}}},
                V3: {'schema_version': V3, 'sample_support': {'TIME3': {'averaging_interval': 0.02}},
                     'evidence': {'rotation': {'origin': 'measured', 'source': 's'}}}}

    def test_every_supported_raw_version_normalizes_to_a_canonical_v3_that_validates(self):
        for version, raw in self.documents().items():
            canonical = validate_manifest(raw)
            self.assertEqual(canonical['schema_version'], V3, version)
            self.assertEqual(canonical['migration']['from'], version)
            schema.validate_structure(canonical, 'evidence.v3')
            for section in ('sample_support', 'event_detectors', 'uncertainty', 'extensions', 'run_identity', 'events'):
                self.assertIn(section, canonical, (version, section))

    def test_normalization_is_idempotent_for_every_version_and_never_mutates_input(self):
        for version, raw in self.documents().items():
            before = copy.deepcopy(raw)
            once = validate_manifest(raw)
            self.assertEqual(raw, before, version)
            self.assertEqual(validate_manifest(once), once, version)
            self.assertEqual(validate_manifest(validate_manifest(once)), once, version)

    def test_migration_steps_are_listed_per_version_hop(self):
        steps = {v: validate_manifest(raw)['migration']['steps'] for v, raw in self.documents().items()}
        self.assertEqual(len(steps[V1]), 4)
        self.assertEqual(len(steps[V2]), 2)
        self.assertEqual(steps[V3], [])
        self.assertTrue(all('->' in s for s in steps[V2][:1]))

    def test_older_documents_may_not_carry_newer_fields_or_a_migration_record(self):
        for document in ({'schema_version': V1, 'events': []}, {'schema_version': V1, 'uncertainty': {}},
                         {'schema_version': V2, 'uncertainty': {'covariances': []}},
                         {'schema_version': V2, 'extensions': {}}, {'schema_version': V2, 'sample_support': {}},
                         {'schema_version': V2, 'migration': {'from': V1, 'raw_sha256': '0' * 64, 'steps': []}}):
            with self.assertRaises(ValueError, msg=str(document)) as caught:
                validate_manifest(document)
            self.assertIn('are not part of schema', str(caught.exception))

    def test_the_evidence_identity_does_not_depend_on_the_declared_version(self):
        raw_v1 = {'schema_version': V1, 'evidence': {'rotation': {'origin': 'assumed', 'source': 's'}}}
        raw_v3 = {'schema_version': V3, 'evidence': {'rotation': {'origin': 'assumed', 'source': 's', 'usage': []}}}
        a = base.CurationEngine(base.NativeSource(base.types.SimpleNamespace(variables=base.variables())), raw_v1).run()
        b = base.CurationEngine(base.NativeSource(base.types.SimpleNamespace(variables=base.variables())), raw_v3).run()
        self.assertEqual(a.fingerprints['evidence_sha256'], b.fingerprints['evidence_sha256'])
        self.assertEqual(a.fingerprints['assessment_sha256'], b.fingerprints['assessment_sha256'])

    def test_all_three_versions_run_through_the_real_cli_on_a_generated_cdf(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            base.ReportingAndCliTests().write_cdf(tmp / 'run.CDF')
            for version, raw in self.documents().items():
                (tmp / 'e.json').write_text(json.dumps(raw))
                done = subprocess.run([sys.executable, str(ROOT / 'curate.py'), '-i', str(tmp / 'run.CDF'),
                                       '--evidence', str(tmp / 'e.json'), '--output-dir', str(tmp / 'out')],
                                      capture_output=True, text=True)
                self.assertEqual(done.returncode, 0, (version, done.stderr))
                document = json.loads((tmp / 'out' / 'curation.json').read_text())
                self.assertEqual(document['effective_manifest']['migration']['from'], version)
                self.assertEqual(document['schema_version'], REPORT_SCHEMA)


class DocumentSchemaTests(unittest.TestCase):
    def test_the_report_structure_is_validated_before_it_is_written(self):
        report, _ = run_screen()
        with tempfile.TemporaryDirectory() as tmp:
            _paths, document = write_report(report, tmp)
        schema.validate_structure(document, 'report.v3')
        broken = copy.deepcopy(document)
        broken['summary']['scientifically_certified'] = True
        with self.assertRaises(ValueError) as caught:
            schema.validate_structure(broken, 'report.v3')
        self.assertIn('summary.scientifically_certified', str(caught.exception))
        with tempfile.TemporaryDirectory() as tmp:
            with mock.patch.object(type(report), 'as_dict', return_value=broken):
                with self.assertRaises(ValueError):
                    write_report(report, Path(tmp) / 'out')
                self.assertFalse((Path(tmp) / 'out').exists())             # nothing was touched

    def test_a_policy_file_is_checked_against_the_schema_with_a_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'policy.json'
            for body, fragment in (
                    ({'schema_version': 'transp2imas.curation.policy.v1', 'parameters': {'weak_gradient_relative': -1}},
                     'parameters.weak_gradient_relative'),
                    ({'schema_version': 'transp2imas.curation.policy.v1', 'rules': {'disable': {'numerics.runtime': ''}}},
                     'rules.disable.numerics.runtime'),
                    ({'schema_version': 'transp2imas.curation.policy.v1', 'applications': ['all']}, 'applications[0]'),
                    ({'schema_version': 'transp2imas.curation.policy.v1', 'surprise': 1}, 'surprise')):
                path.write_text(json.dumps(body))
                with self.assertRaises(ValueError) as caught:
                    load_policy(path)
                self.assertIn(fragment, str(caught.exception))
                self.assertIn('[schema policy.v1]', str(caught.exception))
            path.write_text(json.dumps({'schema_version': 'transp2imas.curation.policy.v1', 'name': 'ok'}))
            self.assertEqual(load_policy(path).name, 'ok')

    def test_pilot_facts_are_structurally_validated_and_placeholders_still_refused_first(self):
        facts = {'schema_version': PILOT_SCHEMA, 'run_identity': {'device': 'X', 'shot': 1}}
        self.assertEqual(build_manifest(facts, '.')['schema_version'], V3)
        with self.assertRaises(ValueError) as caught:
            build_manifest({**facts, 'files': [{'role': 'trdat', 'path': 'x', 'sha256': '0' * 64}]}, '.')
        self.assertIn('[schema pilot.v1]', str(caught.exception))                # a supplied hash is structurally refused
        with self.assertRaises(ValueError) as caught:
            build_manifest({**facts, 'run_identity': {'device': 'TO_BE_SUPPLIED', 'shot': 1}}, '.')
        self.assertIn('TO_BE_SUPPLIED', str(caught.exception))
        self.assertNotIn('[schema', str(caught.exception))


class ReferenceResolutionTests(unittest.TestCase):
    def covariance_resource(self, path, digest):
        return {'id': 'cov_te', 'quantity': 'TE', 'kind': 'covariance', 'path': path, 'sha256': digest, 'units': 'EV',
                'coordinate': 'X', 'ordering': 'ascending_coordinate', 'estimator': 'numpy_gradient_edge_order_1',
                'provenance': 'synthetic fixture, not diagnostic data'}

    def write(self, directory, path='cov.npz'):
        np.savez(Path(directory) / path, time=np.array([0.0, 1.0]), covariance=np.zeros((2, 3, 3)))
        document = {'schema_version': V3, 'uncertainty': {'covariances': [
            self.covariance_resource(path, sha256(Path(directory) / path))]}}
        sidecar = Path(directory) / 'e.json'
        sidecar.write_text(json.dumps(document))
        return sidecar

    def test_relative_references_resolve_against_the_sidecar_not_the_working_directory(self):
        with tempfile.TemporaryDirectory() as tmp, tempfile.TemporaryDirectory() as elsewhere:
            sidecar = self.write(tmp)
            previous = os.getcwd()
            os.chdir(elsewhere)
            try:
                manifest = load_manifest(sidecar)
            finally:
                os.chdir(previous)
            spec = manifest['uncertainty']['covariances'][0]
            self.assertEqual(spec['path'], 'cov.npz')                          # the declared path is kept
            self.assertEqual(spec['resolved_path'], str(Path(tmp).resolve() / 'cov.npz'))
            self.assertEqual(referenced_paths(manifest), [Path(tmp).resolve() / 'cov.npz'])
            self.assertEqual(validate_manifest(manifest), manifest)           # idempotent with resolved paths

    def test_a_missing_or_substituted_resource_fails_clearly(self):
        with tempfile.TemporaryDirectory() as tmp:
            sidecar = self.write(tmp)
            np.savez(Path(tmp) / 'cov.npz', time=np.array([0.0]), covariance=np.ones((1, 3, 3)))   # substituted
            with self.assertRaises(ValueError) as caught:
                load_manifest(sidecar)
            self.assertIn('does not match the declared', str(caught.exception))
            (Path(tmp) / 'cov.npz').unlink()
            with self.assertRaises(ValueError) as caught:
                load_manifest(sidecar)
            self.assertIn('is not a file', str(caught.exception))

    def test_the_report_cannot_overwrite_a_referenced_resource(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            out = tmp / 'out'
            out.mkdir()
            np.savez(tmp / 'tmp_cov.npz', time=np.array([0.0, 1.0]), covariance=np.zeros((2, 3, 3)))
            (out / 'curation.json').write_bytes(b'{}')                        # a resource that collides with a report file
            document = {'schema_version': V3, 'uncertainty': {'covariances': [
                self.covariance_resource('out/curation.json', sha256(out / 'curation.json'))]}}
            (tmp / 'e.json').write_text(json.dumps(document))
            manifest = load_manifest(tmp / 'e.json')
            report, _ = run_screen()
            with self.assertRaises(ValueError) as caught:
                write_report(report, out, protected_paths=referenced_paths(manifest))
            self.assertIn('protected input', str(caught.exception))
            self.assertEqual((out / 'curation.json').read_bytes(), b'{}')

    def test_validate_resources_is_a_no_op_without_covariances(self):
        validate_resources(validate_manifest({'schema_version': V3}))


class ValidateOnlyCliTests(unittest.TestCase):
    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(ROOT / 'curate.py'), *args], capture_output=True, text=True)

    def test_validation_needs_no_cdf_and_says_what_it_did_not_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            sidecar = Path(tmp) / 'e.json'
            sidecar.write_text(json.dumps({'schema_version': V1}))
            done = self.run_cli('--validate-only', '--evidence', str(sidecar))
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertIn(f'declared {V1}, canonical {V3}', done.stdout)
            self.assertIn('Not checked without a source', done.stdout)

    def test_a_structural_error_exits_2_and_names_the_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            sidecar = Path(tmp) / 'e.json'
            sidecar.write_text(json.dumps({'schema_version': V3, 'thresholds': {'weak_gradient_relative': 'x'}}))
            done = self.run_cli('--validate-only', '--evidence', str(sidecar))
            self.assertEqual(done.returncode, 2)
            self.assertIn('thresholds.weak_gradient_relative', done.stderr)

    def test_input_is_required_without_validate_only(self):
        done = self.run_cli('--evidence', 'x.json')
        self.assertEqual(done.returncode, 2)
        self.assertIn('--input', done.stderr)

    def test_a_policy_is_validated_too(self):
        with tempfile.TemporaryDirectory() as tmp:
            policy = Path(tmp) / 'p.json'
            policy.write_text(json.dumps({'schema_version': 'transp2imas.curation.policy.v1', 'rules': {'disable': {'typo.rule': 'r'}}}))
            self.assertEqual(self.run_cli('--validate-only', '--policy', str(policy)).returncode, 0)   # names are checked when planning
            policy.write_text(json.dumps({'schema_version': 'transp2imas.curation.policy.v1', 'surprise': 1}))
            done = self.run_cli('--validate-only', '--policy', str(policy))
            self.assertEqual(done.returncode, 2)
            self.assertIn('[schema policy.v1]', done.stderr)


if __name__ == '__main__':
    unittest.main(verbosity=1)
