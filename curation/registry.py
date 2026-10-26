#!/usr/bin/python
"""Explicit rule registration, dependencies and deterministic planning.

A rule is a function ``rule(ctx, rec)`` with a stable id, a version, a family and declared prerequisites.  Rules are
registered explicitly; nothing is discovered from evidence JSON, and an external rule is trusted code that the
operator opts into (``load_plugin``), recorded with its module version and file hash.

Planning is deterministic: registration order, dependencies first, ``after_all`` rules last.  Duplicate ids, unknown
requested rules, unknown dependencies and cycles are errors, not silent omissions.  A rule whose prerequisite did not
execute is *not run* (state ``missing_prerequisite``): unknown never becomes a pass by omission.
"""

from dataclasses import dataclass, field
import hashlib
import importlib


@dataclass(frozen=True)
class Rule:
    id: str
    version: str
    family: str
    func: object = field(repr=False)
    title: str = ''
    requires: tuple = ()            # ids of rules that must have executed first
    mandatory: bool = False         # source-integrity rule: cannot be disabled by policy
    after_all: bool = False         # runs after every other rule (it reads their findings)
    parameters: tuple = ()          # names of the screening parameters the rule reads


@dataclass(frozen=True)
class PlanItem:
    rule: Rule
    state: str                      # 'run' or 'policy_disabled'
    note: str = ''


class RuleRegistry:
    def __init__(self):
        self._rules = {}
        self.plugins = []

    def register(self, rule):
        if rule.id in self._rules:
            raise ValueError(f'duplicate rule id: {rule.id}')
        self._rules[rule.id] = rule
        return rule

    def __contains__(self, rule_id):
        return rule_id in self._rules

    def __getitem__(self, rule_id):
        return self._rules[rule_id]

    @property
    def rules(self):
        return list(self._rules.values())

    def validate(self):
        """Raise on unknown dependencies, cycles, and a mandatory rule that depends on an optional one."""
        for rule in self._rules.values():
            for dependency in rule.requires:
                if dependency not in self._rules:
                    raise ValueError(f'rule {rule.id}: unknown prerequisite {dependency!r}')
                if rule.mandatory and not self._rules[dependency].mandatory:
                    raise ValueError(f'mandatory rule {rule.id} depends on optional rule {dependency}')
        self._order()                                         # raises on a cycle

    def _order(self):
        ordinary = [r for r in self._rules.values() if not r.after_all]
        last = [r for r in self._rules.values() if r.after_all]
        done, order, active = set(), [], []

        def visit(rule):
            if rule.id in done:
                return
            if rule.id in active:
                raise ValueError('rule prerequisites form a cycle: ' + ' -> '.join(active + [rule.id]))
            active.append(rule.id)
            for dependency in rule.requires:
                if dependency in self._rules:
                    visit(self._rules[dependency])
            active.pop()
            done.add(rule.id)
            order.append(rule)

        for rule in ordinary + last:
            visit(rule)
        # a prerequisite pulled in from the ``after_all`` group must still come before its dependents; the visit
        # order guarantees that, so the only remaining requirement is that after_all rules end up last
        return [r for r in order if not r.after_all] + [r for r in order if r.after_all]

    def plan(self, policy):
        """``[PlanItem]`` in execution order for ``policy``; raises on a request naming an unknown rule."""
        self.validate()
        for rule_id in list(policy.disabled) + list(policy.enabled_only or ()):
            if rule_id not in self._rules:
                raise ValueError(f'policy names unknown rule {rule_id!r}; known: {", ".join(sorted(self._rules))}')
        items = []
        for rule in self._order():
            note = ''
            wanted = policy.enabled_only is None or rule.id in policy.enabled_only
            disabled = rule.id in policy.disabled or not wanted
            if disabled and rule.mandatory:
                note = 'a request to disable this required rule was ignored'
                disabled = False
            if disabled:
                items.append(PlanItem(rule, 'policy_disabled', policy.disabled.get(rule.id, 'not in enable_only')))
            else:
                items.append(PlanItem(rule, 'run', note))
        return items

    def describe(self):
        """JSON-ready identity of the registry: every rule's id, version and contract, and the plugins loaded."""
        return {'rules': [{'id': r.id, 'version': r.version, 'family': r.family, 'title': r.title,
                           'requires': list(r.requires), 'mandatory': r.mandatory, 'after_all': r.after_all,
                           'parameters': list(r.parameters)} for r in self._rules.values()],
                'plugins': list(self.plugins)}

    def load_plugin(self, module_name):
        """Opt-in registration of trusted external rules.

        ``module_name`` comes from the operator (a command-line flag), never from evidence JSON.  The module must
        define ``register(registry)``; its ``__version__`` and the SHA-256 of its source file are recorded.
        """
        module = importlib.import_module(module_name)
        if not callable(getattr(module, 'register', None)):
            raise ValueError(f'plugin {module_name}: no register(registry) function')
        before = set(self._rules)
        module.register(self)
        path = getattr(module, '__file__', None)
        digest = hashlib.sha256(open(path, 'rb').read()).hexdigest() if path else None
        self.plugins.append({'module': module_name, 'version': getattr(module, '__version__', None),
                             'sha256': digest, 'rules': sorted(set(self._rules) - before)})
        self.validate()
