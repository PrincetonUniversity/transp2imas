#!/usr/bin/python
"""Tests of the rule registry, policy, execution ledger, plugins and provenance of curation.

Source-free and IMAS-free: synthetic data only.  Real-case behaviour of the extracted rules is protected by
``verify_curation.py`` and ``verify_curation_gaps.py``; this file covers the architecture around them:

* duplicate ids, unknown requested rules, unknown prerequisites and cycles are errors;
* planning is deterministic, dependencies first, ``after_all`` rules last;
* a policy changes what is asked, never what a finding says; required rules cannot be disabled;
* every rule has a ledger state, and an exception is an execution error, never ``not_assessed``;
* a check that did not run cannot let an application look screened;
* plugins are explicit opt-in, recorded with version and hash, and cannot come from evidence JSON;
* provenance hashes subpackages and is cached; the assessment fingerprint moves with the policy.

    python examples/verify_curation_rules.py
"""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from curation import MANIFEST_SCHEMA
from curation.context import CurationContext
from curation.engine import CurationEngine
from curation.manifest import validate_manifest
from curation.policy import Policy, is_mandatory
from curation.provenance import code_identity
from curation.registry import Rule, RuleRegistry
from curation.results import Recorder, Report
from curation.rules import coordinates, default_registry
from curation.source import NativeSource
import verify_curation as base
from verify_curation import run_screen, types, variables


def noop(ctx, rec):
    return None


def engine(policy=None, registry=None, manifest=None, data=None):
    source = NativeSource(types.SimpleNamespace(variables=data if data is not None else variables()))
    return CurationEngine(source, manifest or {'schema_version': MANIFEST_SCHEMA}, None, policy, registry)


class RegistryTests(unittest.TestCase):
    def test_a_duplicate_rule_id_is_an_error(self):
        registry = RuleRegistry()
        registry.register(Rule('a.one', '1', 'a', noop))
        with self.assertRaises(ValueError) as caught:
            registry.register(Rule('a.one', '2', 'a', noop))
        self.assertIn('duplicate rule id', str(caught.exception))

    def test_an_unknown_prerequisite_is_an_error(self):
        registry = RuleRegistry()
        registry.register(Rule('a.one', '1', 'a', noop, requires=('a.ghost',)))
        with self.assertRaises(ValueError) as caught:
            registry.validate()
        self.assertIn('a.ghost', str(caught.exception))

    def test_a_cycle_is_an_error_and_names_the_path(self):
        registry = RuleRegistry()
        registry.register(Rule('a.one', '1', 'a', noop, requires=('a.two',)))
        registry.register(Rule('a.two', '1', 'a', noop, requires=('a.one',)))
        with self.assertRaises(ValueError) as caught:
            registry.validate()
        self.assertIn('cycle', str(caught.exception))
        self.assertIn('a.one', str(caught.exception))

    def test_a_mandatory_rule_may_not_depend_on_an_optional_one(self):
        registry = RuleRegistry()
        registry.register(Rule('a.optional', '1', 'a', noop))
        registry.register(Rule('a.required', '1', 'a', noop, requires=('a.optional',), mandatory=True))
        with self.assertRaises(ValueError):
            registry.validate()

    def test_the_plan_is_deterministic_with_dependencies_first_and_after_all_last(self):
        registry = RuleRegistry()
        registry.register(Rule('z.last', '1', 'z', noop, after_all=True))
        registry.register(Rule('b.second', '1', 'b', noop, requires=('c.first',)))
        registry.register(Rule('c.first', '1', 'c', noop))
        registry.register(Rule('a.free', '1', 'a', noop))
        order = [item.rule.id for item in registry.plan(Policy.default())]
        self.assertEqual(order, ['c.first', 'b.second', 'a.free', 'z.last'])
        self.assertEqual(order, [item.rule.id for item in registry.plan(Policy.default())])

    def test_a_policy_naming_an_unknown_rule_is_an_error(self):
        registry = default_registry()
        for policy in (Policy(disabled={'nope.rule': 'because'}), Policy(enabled_only=['nope.rule'])):
            with self.assertRaises(ValueError) as caught:
                registry.plan(policy)
            self.assertIn('nope.rule', str(caught.exception))

    def test_the_default_registry_has_stable_ids_and_versions(self):
        ids = [r.id for r in default_registry().rules]
        self.assertEqual(ids, ['coordinates.grid', 'identity.run', 'profiles.screen', 'uncertainty.gradient_covariance',
                               'evidence.provenance',
                               'evidence.rotation', 'composition.quasineutrality', 'model.declared_inputs',
                               'heating.reference', 'inputs.resolved', 'numerics.runtime', 'events.detect',
                               'events.inventory'])
        self.assertTrue(all(r.version == '1' for r in default_registry().rules))


class PolicyTests(unittest.TestCase):
    def ledger(self, report):
        return {e['rule']: e for e in report.ledger}

    def test_a_disabled_optional_rule_stays_visible_in_the_ledger(self):
        report = engine(Policy(disabled={'numerics.runtime': 'not needed here'})).run()
        entry = self.ledger(report)['numerics.runtime']
        self.assertEqual((entry['state'], entry['note']), ('policy_disabled', 'not needed here'))
        self.assertFalse([f for f in report.findings if f.id.startswith('NUMERICS.')])

    def test_a_required_rule_cannot_be_disabled(self):
        report = engine(Policy(disabled={'coordinates.grid': 'skip the grids'})).run()
        entry = self.ledger(report)['coordinates.grid']
        self.assertEqual(entry['state'], 'executed')
        self.assertIn('ignored', entry['note'])
        self.assertTrue([f for f in report.findings if f.id.startswith('GRID.')])

    def test_coverage_counts_required_and_optional_rules_separately(self):
        report = engine(Policy(enabled_only=['profiles.screen'])).run()
        coverage = report.rule_coverage()
        self.assertEqual(coverage['required']['executed'], coverage['required']['total'])      # required always run
        self.assertEqual(coverage['optional']['executed'], 1)
        self.assertIn('numerics.runtime', coverage['optional']['not_executed'])
        self.assertEqual(coverage['optional']['not_executed']['numerics.runtime'], 'policy_disabled')

    def test_a_check_that_did_not_run_does_not_let_an_application_look_screened(self):
        full = engine().run()
        subset = engine(Policy(enabled_only=['profiles.screen'])).run()
        self.assertEqual(set(subset.ledger[i]['state'] for i in range(len(subset.ledger))) - {'executed'}, {'policy_disabled'})
        for application, status in subset.as_dict()['summary']['applications'].items():
            self.assertNotEqual(status, 'screening_only_not_certified', application)
        self.assertIn('rule_coverage', full.as_dict()['summary'])

    def test_a_policy_never_changes_a_finding(self):
        full = engine().run()
        subset = engine(Policy(enabled_only=['coordinates.grid', 'profiles.screen'])).run()
        wanted = {f.id: f for f in subset.findings}
        for finding in full.findings:
            if finding.id in wanted:
                self.assertEqual((finding.status, finding.message), (wanted[finding.id].status, wanted[finding.id].message))
        self.assertEqual(set(wanted) - {f.id for f in full.findings}, set())

    def test_a_missing_prerequisite_means_the_rule_does_not_run(self):
        report = engine(Policy(disabled={'model.declared_inputs': 'off'})).run()
        states = self.ledger(report)
        self.assertEqual(states['heating.reference']['state'], 'missing_prerequisite')
        self.assertEqual(states['heating.reference']['missing'], ['model.declared_inputs'])
        self.assertFalse([f for f in report.findings if f.id.startswith('HEATING.')])

    def test_policy_parameters_override_sidecar_thresholds_and_are_validated(self):
        with self.assertRaises(ValueError):
            Policy(parameters={'not_a_parameter': 1.0})
        with self.assertRaises(ValueError):
            Policy(parameters={'weak_gradient_relative': float('nan')})
        with self.assertRaises(ValueError):
            Policy(parameters={'weak_gradient_relative': -1.0})
        eng = engine(Policy(parameters={'weak_gradient_relative': 0.5}),
                     manifest={'schema_version': MANIFEST_SCHEMA, 'thresholds': {'weak_gradient_relative': 0.01}})
        self.assertEqual(eng.context.thresholds['weak_gradient_relative'], 0.5)
        self.assertEqual(eng.manifest['thresholds']['weak_gradient_relative'], 0.01)         # the evidence is unchanged

    def test_a_disable_reason_is_required(self):
        with self.assertRaises(ValueError):
            Policy(disabled={'numerics.runtime': ''})

    def test_policy_documents_are_strict(self):
        policy = Policy.from_dict({'schema_version': 'transp2imas.curation.policy.v1', 'name': 'n',
                                   'rules': {'disable': {'numerics.runtime': 'r'}}})
        self.assertEqual(policy.to_dict()['rules']['disable'], {'numerics.runtime': 'r'})
        for bad in ({'schema_version': 'x'}, {'schema_version': 'transp2imas.curation.policy.v1', 'typo': 1},
                    {'schema_version': 'transp2imas.curation.policy.v1', 'rules': {'typo': []}}):
            with self.assertRaises(ValueError):
                Policy.from_dict(bad)

    def test_mandatory_patterns_cover_source_integrity_and_not_physics_screens(self):
        for identifier in ('GRID.TIME3', 'GRID.X', 'PROFILE.TE.GRID_MISMATCH', 'PROFILE.NE.UNITS', 'NUMERICS.RESL2TE.READ',
                           'RUN.IDENTITY', 'PROFILE.TE.ASSESSMENT'):
            self.assertTrue(is_mandatory(identifier), identifier)
        for identifier in ('PROFILE.TE.NEGATIVE', 'PROFILE.TE.FINITE', 'MODEL.DENSITY_CURRENT_COUPLING', 'FIT.TE.WEAK_GRADIENT'):
            self.assertFalse(is_mandatory(identifier), identifier)


class ExecutionLedgerTests(unittest.TestCase):
    def failing_registry(self, exception):
        registry = default_registry()

        def boom(ctx, rec):
            raise exception
        registry.register(Rule('test.boom', '1', 'test', boom, 'raises'))
        return registry

    def test_an_exception_is_an_execution_error_not_a_missing_evidence_finding(self):
        report = engine(registry=self.failing_registry(KeyError('programming bug'))).run()
        entry = {e['rule']: e for e in report.ledger}['test.boom']
        self.assertEqual(entry['state'], 'execution_error')
        self.assertEqual(entry['error']['type'], 'KeyError')
        finding = {f.id: f for f in report.findings}['RULE.test.boom.EXECUTION']
        self.assertEqual(finding.status, 'error')
        self.assertTrue(report.has_execution_errors())
        self.assertEqual(report.as_dict()['summary']['applications']['heat_transport'], 'blocked_by_screening_findings')

    def test_other_rules_still_run_after_one_raises(self):
        report = engine(registry=self.failing_registry(RuntimeError('x'))).run()
        self.assertTrue([f for f in report.findings if f.id.startswith('GRID.')])
        self.assertEqual({e['rule']: e['state'] for e in report.ledger}['events.inventory'], 'executed')

    def test_a_rule_whose_prerequisite_raised_is_not_run(self):
        registry = default_registry()

        def boom(ctx, rec):
            raise ValueError('data-shaped programming error')
        registry.register(Rule('test.first', '1', 'test', boom))
        registry.register(Rule('test.second', '1', 'test', noop, requires=('test.first',)))
        report = engine(registry=registry).run()
        states = {e['rule']: e['state'] for e in report.ledger}
        self.assertEqual((states['test.first'], states['test.second']), ('execution_error', 'missing_prerequisite'))

    def test_findings_carry_the_rule_and_version_that_produced_them(self):
        report = engine().run()
        found = {f.id: f for f in report.findings}
        self.assertEqual((found['GRID.TIME'].rule, found['GRID.TIME'].rule_version), ('coordinates.grid', '1'))
        self.assertEqual(found['PROFILE.TE.NEGATIVE'].rule, 'profiles.screen')
        self.assertEqual(found['EVENTS.INVENTORY'].rule, 'events.inventory')

    def test_ledger_counts_match_the_findings_of_each_rule(self):
        report = engine().run()
        for entry in report.ledger:
            mine = [f for f in report.findings if f.rule == entry['rule']]
            self.assertEqual(sum(entry['findings'].values()), len(mine), entry['rule'])

    def test_the_cli_exits_2_after_writing_the_report_when_a_rule_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            base.ReportingAndCliTests().write_cdf(tmp / 'run.CDF')
            (tmp / 'boom_plugin.py').write_text(textwrap.dedent('''\
                from curation.registry import Rule
                __version__ = '0.0.1-test'

                def explode(ctx, rec):
                    raise AttributeError('plugin bug')

                def register(registry):
                    registry.register(Rule('plugin.explode', '1', 'plugin', explode, 'raises'))
                '''))
            done = subprocess.run([sys.executable, str(ROOT / 'curate.py'), '-i', str(tmp / 'run.CDF'),
                                   '--plugin', 'boom_plugin', '--output-dir', str(tmp / 'out')],
                                  capture_output=True, text=True, env={'PYTHONPATH': str(tmp), 'PATH': '/usr/bin:/bin'})
            self.assertEqual(done.returncode, 2, done.stderr)
            self.assertIn('plugin.explode', done.stderr)
            document = json.loads((tmp / 'out' / 'curation.json').read_text())
            ledger = {e['rule']: e for e in document['execution_ledger']}
            self.assertEqual(ledger['plugin.explode']['state'], 'execution_error')
            plugin = document['registry']['plugins'][0]
            self.assertEqual((plugin['module'], plugin['version'], plugin['rules']), ('boom_plugin', '0.0.1-test', ['plugin.explode']))
            self.assertEqual(len(plugin['sha256']), 64)


class PluginTests(unittest.TestCase):
    def test_plugins_load_only_by_explicit_registration_and_are_recorded(self):
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, 'my_rules.py').write_text(textwrap.dedent('''\
                from curation.registry import Rule
                __version__ = '2.0'

                def screen(ctx, rec):
                    rec.report.add('PLUGIN.SEEN', 'pass', 'plugin rule ran')

                def register(registry):
                    registry.register(Rule('plugin.screen', '1', 'plugin', screen, 'x'))
                '''))
            sys.path.insert(0, tmp)
            try:
                registry = default_registry()
                registry.load_plugin('my_rules')
            finally:
                sys.path.remove(tmp)
                sys.modules.pop('my_rules', None)
            report = engine(registry=registry).run()
        self.assertIn('PLUGIN.SEEN', {f.id for f in report.findings})
        self.assertEqual(report.registry['plugins'][0]['rules'], ['plugin.screen'])
        self.assertEqual(report.registry['plugins'][0]['version'], '2.0')

    def test_a_module_without_register_is_refused(self):
        with self.assertRaises(ValueError):
            default_registry().load_plugin('json')

    def test_evidence_json_cannot_name_a_plugin(self):
        for key in ('plugins', 'rules', 'plugin', 'module'):
            with self.assertRaises(ValueError, msg=key):
                validate_manifest({'schema_version': MANIFEST_SCHEMA, key: ['os.system']})


class IndependentRuleTests(unittest.TestCase):
    """A rule runs on its own against a context and recorder: no engine, no IMAS, no experimental data."""

    def test_the_grid_rule_runs_standalone(self):
        source = NativeSource(types.SimpleNamespace(variables=variables()))
        manifest = validate_manifest({'schema_version': MANIFEST_SCHEMA})
        ctx = CurationContext(source, manifest)
        report = Report(source.provenance(), manifest, ctx.timeline)
        coordinates.check_coordinates(ctx, Recorder(ctx, report))
        found = {f.id: f.status for f in report.findings}
        self.assertEqual((found['GRID.TIME'], found['GRID.TIME3'], found['GRID.X']), ('pass', 'pass', 'pass'))

    def test_the_context_does_not_modify_the_source_arrays(self):
        data = variables()
        before = {k: v.values.copy() for k, v in data.items()}
        engine(data=data).run()
        for name, values in before.items():
            self.assertTrue((data[name].values == values).all(), name)


class ProvenanceAndFingerprintTests(unittest.TestCase):
    def test_provenance_hashes_subpackages_shared_modules_and_dependencies(self):
        identity = code_identity()
        self.assertIn('curation/rules/profiles.py', identity['modules'])
        self.assertIn('curation/engine.py', identity['modules'])
        self.assertIn('fortran_literals.py', identity['modules'])
        self.assertIn('numpy', identity['dependencies'])
        self.assertTrue(all(len(h) == 64 for h in identity['modules'].values()))

    def test_provenance_is_computed_once_per_process(self):
        self.assertIs(code_identity(), code_identity())

    def test_the_assessment_fingerprint_moves_with_the_policy_but_the_evidence_fingerprint_does_not(self):
        base = engine().run().fingerprints
        tightened = engine(Policy(parameters={'weak_gradient_relative': 0.5})).run().fingerprints
        subset = engine(Policy(enabled_only=['coordinates.grid'])).run().fingerprints
        self.assertEqual(base['evidence_sha256'], tightened['evidence_sha256'])
        self.assertEqual(base['evidence_sha256'], subset['evidence_sha256'])
        self.assertEqual(len({base['assessment_sha256'], tightened['assessment_sha256'], subset['assessment_sha256']}), 3)

    def test_the_same_inputs_give_the_same_fingerprints_and_ledger(self):
        a, b = engine().run(), engine().run()
        self.assertEqual(a.fingerprints, b.fingerprints)
        self.assertEqual([e['rule'] for e in a.ledger], [e['rule'] for e in b.ledger])

    def test_the_report_carries_policy_registry_ledger_and_coverage(self):
        document = engine().run().as_dict()
        for key in ('execution_ledger', 'policy', 'registry', 'fingerprints'):
            self.assertIn(key, document)
        self.assertEqual(document['policy']['name'], 'default')
        self.assertEqual(len(document['registry']['rules']), 13)
        json.dumps(document, allow_nan=False)


class CliPolicyTests(unittest.TestCase):
    def test_the_cli_applies_a_policy_file_and_records_its_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            base.ReportingAndCliTests().write_cdf(tmp / 'run.CDF')
            policy = tmp / 'policy.json'
            policy.write_text(json.dumps({'schema_version': 'transp2imas.curation.policy.v1', 'name': 'grids-only',
                                          'rules': {'enable_only': ['coordinates.grid']}}))
            done = subprocess.run([sys.executable, str(ROOT / 'curate.py'), '-i', str(tmp / 'run.CDF'), '--policy', str(policy),
                                   '--output-dir', str(tmp / 'out')], capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stderr)
            document = json.loads((tmp / 'out' / 'curation.json').read_text())
            self.assertEqual(document['policy']['name'], 'grids-only')
            self.assertEqual(document['source']['policy_file']['path'], str(policy.resolve()))
            states = {e['rule']: e['state'] for e in document['execution_ledger']}
            self.assertEqual(states['profiles.screen'], 'policy_disabled')
            self.assertEqual(states['identity.run'], 'executed')                     # required
            self.assertIn('curation/engine.py', document['source']['code_identity']['modules'])
            blocked = subprocess.run([sys.executable, str(ROOT / 'curate.py'), '-i', str(tmp / 'run.CDF'), '--policy', str(policy),
                                      '--require-complete', '--output-dir', str(tmp / 'out2')], capture_output=True, text=True)
            self.assertEqual(blocked.returncode, 1)                                  # rules did not all execute

    def test_the_cli_refuses_a_policy_that_names_an_unknown_rule(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            base.ReportingAndCliTests().write_cdf(tmp / 'run.CDF')
            policy = tmp / 'policy.json'
            policy.write_text(json.dumps({'schema_version': 'transp2imas.curation.policy.v1',
                                          'rules': {'disable': {'typo.rule': 'r'}}}))
            done = subprocess.run([sys.executable, str(ROOT / 'curate.py'), '-i', str(tmp / 'run.CDF'), '--policy', str(policy),
                                   '--output-dir', str(tmp / 'out')], capture_output=True, text=True)
            self.assertEqual(done.returncode, 2)
            self.assertIn('typo.rule', done.stderr)
            self.assertFalse((tmp / 'out').exists())


if __name__ == '__main__':
    unittest.main(verbosity=1)
