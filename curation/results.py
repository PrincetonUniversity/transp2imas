#!/usr/bin/python
"""Findings, the report, the rule-execution ledger and the recorder rules write through.

A *finding* is an observation with a scientific status.  The *ledger* is separate: it records what happened to each
rule (executed, disabled by policy, missing a prerequisite, or raised), so that a check that did not run can never be
mistaken for one that passed.  Neither is changed by an application policy; policy only decides which findings block
which application (see ``Report.application_status``).
"""

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import hashlib
import json
import re

import numpy as np

from curation import REPORT_SCHEMA, __version__
from curation import applicability, events as event_rules, policy as policy_rules
from curation.source import unit_key

APPLICATIONS = ('heat_transport', 'momentum_transport', 'current_evolution', 'surrogate_training')
ROLE_APPS = {
    'electron_temperature': ('heat_transport', 'current_evolution', 'surrogate_training'),
    'ion_temperature': ('heat_transport', 'momentum_transport', 'surrogate_training'),
    'electron_density': APPLICATIONS, 'rotation': ('momentum_transport', 'heat_transport', 'surrogate_training'),
    'composition': APPLICATIONS, 'heating': APPLICATIONS,
}
STATUSES = ('pass', 'warning', 'fail', 'not_assessed', 'not_applicable', 'error')
# what happened to a rule, as distinct from what it found
LEDGER_STATES = ('executed', 'policy_disabled', 'missing_prerequisite', 'execution_error')


@dataclass
class Finding:
    id: str
    status: str
    message: str
    variables: list = field(default_factory=list)
    applications: list = field(default_factory=lambda: list(APPLICATIONS))
    evidence: dict = field(default_factory=dict)
    metrics: dict = field(default_factory=dict)
    scope: dict = field(default_factory=dict)
    recommendation: str = ''
    rule: str = ''
    rule_version: str = ''


@dataclass
class Report:
    source: dict
    manifest: dict
    timeline: dict
    findings: list = field(default_factory=list)
    masks: dict = field(default_factory=dict, repr=False)
    coordinates: dict = field(default_factory=dict, repr=False)
    channel_inventory: dict = field(default_factory=dict)
    read_counts: dict = field(default_factory=dict)
    run_identity: dict = field(default_factory=dict)
    evidence_graph: dict = field(default_factory=dict)
    policy_exclusions: dict = field(default_factory=dict)
    events: list = field(default_factory=list)
    detected_events: list = field(default_factory=list)
    uncertainty: dict = field(default_factory=dict)
    ledger: list = field(default_factory=list)
    policy: dict = field(default_factory=dict)
    registry: dict = field(default_factory=dict)
    fingerprints: dict = field(default_factory=dict)
    current_rule: str = ''
    current_version: str = ''

    def add(self, identifier, status, message, **kwargs):
        if status not in STATUSES or any(f.id == identifier for f in self.findings):
            raise ValueError(f'invalid or duplicate finding: {identifier}, {status}')
        finding = Finding(identifier, status, message, rule=self.current_rule, rule_version=self.current_version,
                          **kwargs)
        self.findings.append(finding)
        return finding

    # -- the ledger: what happened to the rules
    def required_not_executed(self):
        return [e['rule'] for e in self.ledger if e['required'] and e['state'] != 'executed']

    def not_executed(self):
        """Every rule that did not run to completion, required or optional."""
        return [e['rule'] for e in self.ledger if e['state'] != 'executed']

    def has_execution_errors(self):
        return any(e['state'] == 'execution_error' for e in self.ledger)

    def rule_coverage(self):
        """Required (source-integrity) and optional rules counted separately; a disabled rule stays visible."""
        coverage = {}
        for label, required in (('required', True), ('optional', False)):
            rows = [e for e in self.ledger if e['required'] == required]
            coverage[label] = {'total': len(rows), 'executed': sum(e['state'] == 'executed' for e in rows),
                               'not_executed': {e['rule']: e['state'] for e in rows if e['state'] != 'executed'}}
        return coverage

    # -- application policy over raw findings
    def _excused(self, finding, application):
        """Whether a declared policy exclusion covers every flagged sample of ``finding`` for ``application``.

        The raw finding is not changed; this only decides whether it blocks that application.  Source-integrity
        rules never get a record, so they are never excused.
        """
        record = self.policy_exclusions.get(finding.id)
        return bool(record and record['covers_all_flagged'] and application in record['applications'])

    def application_status(self, application):
        findings = [f for f in self.findings if application in f.applications]
        live = [f for f in findings if not self._excused(f, application)]
        if any(f.status in ('fail', 'error') for f in live):
            return 'blocked_by_screening_findings'
        if any(f.status == 'warning' for f in live):
            return 'review_required'
        unassessed = [f for f in live if f.status == 'not_assessed' or self._policy_unassessed(f, application)]
        if unassessed or self.not_executed():                 # a check that did not run is not a check that passed
            return 'insufficient_evidence'
        return 'screening_only_not_certified'

    def _policy_unassessed(self, finding, application):
        record = self.policy_exclusions.get(finding.id)
        return bool(record and record['policy_status'] == 'not_assessed' and application in record['applications'])

    def excused_findings(self):
        return {a: [f.id for f in self.findings if a in f.applications and f.status in ('warning', 'fail')
                    and self._excused(f, a)] for a in APPLICATIONS}

    def as_dict(self):
        counts = {s: sum(f.status == s for f in self.findings) for s in STATUSES}
        applications = self.policy.get('applications') or list(APPLICATIONS)
        return {'schema_version': REPORT_SCHEMA, 'curation_version': __version__,
                'generated_at_utc': datetime.now(timezone.utc).isoformat(),
                'source': self.source,
                'manifest_sha256': hashlib.sha256(json.dumps(self.manifest, sort_keys=True,
                                                             allow_nan=False).encode()).hexdigest(),
                'effective_manifest': self.manifest,
                'fingerprints': self.fingerprints,
                'declared_input_timeline': self.timeline,
                'summary': {'counts': counts, 'scientifically_certified': False,
                            'excused_findings': self.excused_findings(),
                            'rule_coverage': self.rule_coverage(),
                            'scope': 'read-only source screening; no IDS backend validation',
                            'applications': {a: self.application_status(a) for a in applications}},
                'findings': [asdict(f) for f in self.findings],
                'execution_ledger': self.ledger, 'policy': self.policy, 'registry': self.registry,
                'run_identity': self.run_identity, 'evidence_graph': self.evidence_graph,
                'events': self.events, 'detected_events': self.detected_events,
                'policy_exclusions': self.policy_exclusions, 'uncertainty': self.uncertainty,
                'channel_inventory': self.channel_inventory, 'source_read_counts': self.read_counts,
                'mask_storage': {'format': 'numpy_npz',
                                 'keys': sorted(set(self.masks) | set(self.coordinates)),
                                 'semantics': 'evaluated masks distinguish unassessed samples from unflagged samples; '
                                              '__excluded masks mark evaluated samples a declared event makes ineligible '
                                              'under a policy; the raw finding and its masks are never changed by them'}}


class Recorder:
    """The write side rules use: findings, masks and coordinates go to the report; the context is only read."""

    def __init__(self, ctx, report):
        self.ctx = ctx
        self.report = report

    def selection(self, name):
        """Select native samples; never infer a grid from equal array lengths."""
        ctx = self.ctx
        data = ctx.source.array(name)
        meta = ctx.source.meta(name)
        dims = tuple(meta['dimensions'])
        mask = np.ones(data.shape, dtype=bool)
        for axis, dim in enumerate(dims):
            if dim not in ('TIME', 'TIME3'):
                continue
            if not ctx.source.has(dim):
                raise ValueError(f'{name}: missing native time coordinate {dim}')
            times = ctx.source.array(dim)
            if (times.ndim != 1 or len(times) != data.shape[axis] or
                    not np.all(np.isfinite(times)) or not np.all(np.diff(times) > 0)):
                raise ValueError(f'{name}: invalid native coordinate {dim}')
            if unit_key(ctx.source.meta(dim)['units']) not in ('S', 'SEC', 'SECS', 'SECOND', 'SECONDS'):
                raise ValueError(f'{dim}: native time units must be seconds')
            selected = np.ones(times.shape, dtype=bool)
            if 'time_min' in ctx.window:
                selected &= times >= ctx.window['time_min']
            if 'time_max' in ctx.window:
                selected &= times <= ctx.window['time_max']
            reshape = [1] * data.ndim
            reshape[axis] = len(times)
            mask &= selected.reshape(reshape)
            self.report.coordinates['coordinate_' + dim] = times
        # Save an explicit native radial coordinate even without a radial cut.
        for radial in ('X', 'XB'):
            if radial in dims and ctx.source.has(radial):
                radius = ctx.source.array(radial)
                if radius.shape == data.shape and tuple(ctx.source.meta(radial)['dimensions']) == dims:
                    self.report.coordinates['coordinate_' + radial] = radius
        if any(key.startswith('rho_') for key in ctx.window) and data.ndim >= 2:
            radial = next((d for d in dims if d in ('X', 'XB')), None)
            if radial is None or not ctx.source.has(radial):
                raise ValueError(f'{name}: requested radial window has no explicit X/XB coordinate')
            radius = ctx.source.array(radial)
            if radius.shape != data.shape or tuple(ctx.source.meta(radial)['dimensions']) != dims:
                raise ValueError(f'{name}: radial window requires the matching native X/XB grid')
            if not ctx.normalized_radius(radial):
                raise ValueError(f'{radial}: normalized-radius unit is not established')
            if 'rho_min' in ctx.window:
                mask &= radius >= ctx.window['rho_min']
            if 'rho_max' in ctx.window:
                mask &= radius <= ctx.window['rho_max']
            mask &= np.isfinite(radius)
            self.report.coordinates['coordinate_' + radial] = radius
        return mask

    def _scope(self, name, identifier, evaluated, flagged):
        key = re.sub('[^A-Za-z0-9_]', '_', identifier)
        self.report.masks[key + '__evaluated'] = np.asarray(evaluated, dtype=bool)
        self.report.masks[key + '__flagged'] = np.asarray(flagged, dtype=bool)
        points = np.argwhere(flagged)
        return {'source_dimensions': self.ctx.source.meta(name)['dimensions'],
                'shape': list(evaluated.shape), 'evaluated_count': int(evaluated.sum()),
                'flagged_count': int(flagged.sum()), 'mask_prefix': key,
                'flagged_index_examples': points[:12].tolist(), 'index_examples_truncated': len(points) > 12,
                'assessment_window': self.ctx.window}

    def _event_membership(self, name, identifier, shape):
        """Samples of ``name`` inside a declared event that applies to rule ``identifier``.

        Returns ``(mask, {event id: mask}, {clock: basis})``.  Membership is a fact about the clock; it is neither
        missingness nor an integrity result.  Source-integrity rules have no membership: no event may touch them.
        """
        ctx = self.ctx
        mask = np.zeros(shape, dtype=bool)
        each, bases = {}, {}
        if policy_rules.is_mandatory(identifier) or not any(
                event_rules.acts_as_exclusion(e) and applicability.rule_matches(e['rules'], identifier)
                for e in ctx.events):
            return mask, each, bases
        for axis, dim in enumerate(ctx.source.meta(name)['dimensions']):
            if dim not in ('TIME', 'TIME3') or not ctx.source.has(dim) or axis >= len(shape):
                continue
            times = ctx.source.array(dim)
            if times.ndim != 1 or len(times) != shape[axis]:
                continue
            inside, per_event, basis = event_rules.membership(ctx.events, identifier, times,
                                                              ctx.manifest['sample_support'].get(dim))
            reshape = [1] * len(shape)
            reshape[axis] = len(times)
            mask |= inside.reshape(reshape)
            for event_id, hit in per_event.items():
                full = np.broadcast_to(hit.reshape(reshape), shape)
                each[event_id] = each.get(event_id, np.zeros(shape, dtype=bool)) | full
            bases[dim] = basis
        return mask, each, bases

    def _record_policy(self, identifier, name, evaluated, flagged, member, each, bases):
        """Record which evaluated samples are ineligible for which applications, leaving the finding as observed."""
        excluded = member & evaluated
        event_ids = [e for e, hit in each.items() if (hit & evaluated).any()]
        if not event_ids or not excluded.any():
            return
        key = re.sub('[^A-Za-z0-9_]', '_', identifier)
        self.report.masks[key + '__excluded'] = excluded
        eligible = int((evaluated & ~excluded).sum())
        by_id = {e['id']: e for e in self.ctx.events}
        applications = sorted({a for e in event_ids for a in by_id[e].get('applications', APPLICATIONS)})
        self.report.policy_exclusions[identifier] = {
            'variable': name, 'events': list(event_ids), 'applications': applications,
            'origins': {e: by_id[e]['origin'] for e in event_ids},
            'policy': 'legacy_exclusion', 'support_basis': bases,
            'total_count': int(evaluated.size), 'unassessed_count': int((~evaluated).sum()),
            'excluded_count': int(excluded.sum()),
            'excluded_by_event': {e: int((each[e] & evaluated).sum()) for e in event_ids},
            'flagged_excluded_count': int((flagged & excluded).sum()), 'eligible_count': eligible,
            'covers_all_flagged': bool(flagged.any() and not (flagged & ~excluded).any()),
            'policy_status': 'not_assessed' if eligible == 0 else 'eligible'}

    def samples(self, identifier, name, evaluated, flagged, message, *, severity='warning',
                applications=APPLICATIONS, metrics=None, recommendation=''):
        evaluated, flagged = np.asarray(evaluated, dtype=bool), np.asarray(flagged, dtype=bool)
        n = int(evaluated.sum())
        status = 'not_assessed' if n == 0 else severity if np.any(flagged) else 'pass'
        finding = self.report.add(identifier, status, message if n else 'No assessable samples in the selected window.',
                                  variables=[name], applications=list(applications), metrics=metrics or {},
                                  scope=self._scope(name, identifier, evaluated, flagged), recommendation=recommendation)
        member, each, bases = self._event_membership(name, identifier, evaluated.shape)
        self._record_policy(identifier, name, evaluated, flagged, member, each, bases)
        return finding

    def unavailable(self, identifier, message, names=(), apps=APPLICATIONS):
        return self.report.add(identifier, 'not_assessed', message, variables=list(names), applications=list(apps))

    def source_problem(self, identifier, exc, names=(), apps=APPLICATIONS):
        # A backend/read failure is not merely absent scientific evidence.
        from curation.source import SourceDataError
        if isinstance(exc, (OSError, RuntimeError, TypeError, SourceDataError)):
            return self.report.add(identifier, 'error', str(exc), variables=list(names), applications=list(apps))
        return self.unavailable(identifier, str(exc), names, apps)
