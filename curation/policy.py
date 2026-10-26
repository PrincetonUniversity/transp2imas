#!/usr/bin/python
"""Application policy: what a declared exclusion may and may not do to a raw finding.

A *raw finding* records an observation about the source (a value is negative, a sample is missing, a clock is
reversed).  A *policy decision* says whether a sample is eligible for an application.  Events and waivers act on
eligibility only: they never rewrite the finding.

Source-integrity rules are *mandatory*: no event, wildcard or waiver excludes samples from them, because a bad
clock, grid, unit or read error invalidates the assessment itself rather than a transient within it.
"""

from curation.applicability import rule_matches

MANDATORY_RULES = ('GRID.*', 'PROFILE.*.GRID_MISMATCH', 'PROFILE.*.UNITS', 'PROFILE.*.ASSESSMENT', '*.READ',
                   'RUN.IDENTITY')


def is_mandatory(identifier):
    """Whether ``identifier`` is a source-integrity rule that no declared exclusion may touch."""
    return rule_matches(MANDATORY_RULES, identifier)


POLICY_SCHEMA = 'transp2imas.curation.policy.v1'
# screening parameters a policy may override (positive finite numbers); everything else is evidence, not policy
POLICY_PARAMETERS = ('quasineutrality_relative', 'weak_gradient_relative', 'gradient_sigma_multiplier',
                     'newton_iterations_warning')
ALL_APPLICATIONS = ('heat_transport', 'momentum_transport', 'current_evolution', 'surrogate_training')


class Policy:
    """Which rules run, with which screening parameters, and which applications are reported.

    A policy changes *what is asked*, never what a finding says.  Rules it disables are recorded in the execution
    ledger as ``policy_disabled`` (not silently absent); required source-integrity rules cannot be disabled, and a
    request to do so is kept in ``ignored_disable_requests`` rather than honoured.
    """

    def __init__(self, name='default', enabled_only=None, disabled=None, parameters=None, applications=None):
        self.name = name
        self.enabled_only = None if enabled_only is None else tuple(enabled_only)
        self.disabled = dict(disabled or {})
        self.parameters = dict(parameters or {})
        self.applications = tuple(applications) if applications else ALL_APPLICATIONS
        for key, value in self.parameters.items():
            if key not in POLICY_PARAMETERS:
                raise ValueError(f'policy.parameters.{key}: not a policy parameter (allowed: {POLICY_PARAMETERS})')
            if isinstance(value, bool) or not isinstance(value, (int, float)) or value != value or abs(value) == float('inf') or value <= 0:
                raise ValueError(f'policy.parameters.{key}: expected a positive finite number')
        for application in self.applications:
            if application not in ALL_APPLICATIONS:
                raise ValueError(f'policy.applications: unknown application {application!r}')
        for rule_id, reason in self.disabled.items():
            if not isinstance(reason, str) or not reason.strip():
                raise ValueError(f'policy.rules.disable.{rule_id}: a reason is required')

    @classmethod
    def default(cls):
        return cls()

    @classmethod
    def from_dict(cls, document):
        """Build a policy from a validated ``transp2imas.curation.policy.v1`` document."""
        if document.get('schema_version') != POLICY_SCHEMA:
            raise ValueError(f'policy.schema_version must be {POLICY_SCHEMA!r}')
        unknown = set(document) - {'$schema', 'schema_version', 'name', 'rules', 'parameters', 'applications', 'notes'}
        if unknown:
            raise ValueError(f'policy: unknown keys: {", ".join(sorted(unknown))}')
        rules = document.get('rules', {})
        unknown = set(rules) - {'enable_only', 'disable'}
        if unknown:
            raise ValueError(f'policy.rules: unknown keys: {", ".join(sorted(unknown))}')
        return cls(document.get('name', 'unnamed'), rules.get('enable_only'), rules.get('disable'),
                   document.get('parameters'), document.get('applications'))

    def to_dict(self):
        """The canonical, JSON-ready form; its hash is the policy fingerprint."""
        return {'schema_version': POLICY_SCHEMA, 'name': self.name,
                'rules': {'enable_only': None if self.enabled_only is None else list(self.enabled_only),
                          'disable': dict(sorted(self.disabled.items()))},
                'parameters': dict(sorted(self.parameters.items())), 'applications': list(self.applications)}

    def fingerprint(self):
        import hashlib
        import json
        return hashlib.sha256(json.dumps(self.to_dict(), sort_keys=True, allow_nan=False).encode()).hexdigest()


def load_policy(path):
    """Read a policy file: strict JSON, then the packaged schema, then the semantic checks of ``Policy``."""
    from pathlib import Path
    from curation.schema import loads_strict, validate_structure
    document = loads_strict(Path(path).read_text(encoding='utf-8'))
    validate_structure(document, 'policy.v1')
    return Policy.from_dict(document)
