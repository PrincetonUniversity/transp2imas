#!/usr/bin/python
"""Evidence-aware screening orchestration; findings are not a scientific certification.

The engine only plans and executes rules.  Each rule (``curation.rules``) reads the immutable ``CurationContext`` and
writes through the ``Recorder``; the engine records in the execution ledger what happened to every rule (executed,
disabled by policy, missing a prerequisite, or raised) so that a check that did not run is never mistaken for one
that passed.  An exception in a rule is an execution error, never a ``not_assessed`` finding.  No quantity is
repaired here.
"""

import hashlib
import json

from curation.context import CurationContext
from curation.manifest import validate_manifest
from curation.policy import Policy
from curation.results import APPLICATIONS, Finding, Recorder, Report
from curation.rules import default_registry

__all__ = ['CurationEngine', 'Finding', 'Report', 'APPLICATIONS']


def _digest(document):
    return hashlib.sha256(json.dumps(document, sort_keys=True, allow_nan=False).encode()).hexdigest()


class CurationEngine:
    """A standalone, injectable engine; input evidence remains analyst-declared."""

    def __init__(self, source, manifest, timeline=None, policy=None, registry=None):
        self.source = source
        self.manifest = validate_manifest(manifest)
        self.policy = policy or Policy.default()
        self.registry = registry or default_registry()
        self.context = CurationContext(source, self.manifest, timeline, self.policy.parameters)
        self.timeline = self.context.timeline
        self.report = Report(source.provenance(), self.manifest, self.timeline)
        self.recorder = Recorder(self.context, self.report)
        self.window = self.context.window
        self.thresholds = self.context.thresholds

    def run(self):
        report = self.report
        executed = set()
        for item in self.registry.plan(self.policy):
            rule = item.rule
            entry = {'rule': rule.id, 'version': rule.version, 'family': rule.family, 'required': rule.mandatory,
                     'state': 'executed', 'note': item.note, 'findings': {}}
            report.ledger.append(entry)
            if item.state == 'policy_disabled':
                entry.update(state='policy_disabled', note=item.note)
                continue
            missing = [r for r in rule.requires if r not in executed]
            if missing:
                entry.update(state='missing_prerequisite', missing=missing)
                continue
            report.current_rule, report.current_version = rule.id, rule.version
            first = len(report.findings)
            try:
                rule.func(self.context, self.recorder)
                executed.add(rule.id)
            except Exception as exc:          # a programming error in a rule is an execution error, not missing evidence
                entry.update(state='execution_error', error={'type': type(exc).__name__, 'message': str(exc)})
                report.add(f'RULE.{rule.id}.EXECUTION', 'error',
                           f'Rule {rule.id} raised {type(exc).__name__}: {exc}', applications=list(APPLICATIONS))
            for finding in report.findings[first:]:
                entry['findings'][finding.status] = entry['findings'].get(finding.status, 0) + 1
        report.current_rule = report.current_version = ''
        report.read_counts = dict(self.source.read_counts)
        report.policy = self.policy.to_dict()
        report.registry = self.registry.describe()
        report.fingerprints = self.fingerprints()
        return report

    def fingerprints(self):
        """Evidence identity and effective assessment identity, kept apart.

        The evidence fingerprint covers the scientific content of the sidecar only (not its migration record).  The
        assessment fingerprint adds everything that changes what was asked of it: the policy, the effective screening
        parameters and the version of every rule, so a threshold change alters it though the evidence is unchanged.
        """
        evidence = {k: v for k, v in self.manifest.items() if k != 'migration'}
        evidence_sha = _digest(evidence)
        assessment = {'evidence': evidence_sha, 'policy': self.policy.to_dict(),
                      'parameters': self.context.thresholds,
                      'rules': [[r.id, r.version] for r in self.registry.rules],
                      'plugins': self.registry.plugins}
        return {'evidence_sha256': evidence_sha, 'policy_sha256': self.policy.fingerprint(),
                'assessment_sha256': _digest(assessment)}
