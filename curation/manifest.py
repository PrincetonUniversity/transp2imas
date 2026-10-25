#!/usr/bin/python
"""Strict validation for a small, explicit evidence sidecar; no assumed defaults.

Evidence strings are analyst declarations, not independently verified facts.
"""

import copy
import json
import math
from pathlib import Path

from curation import MANIFEST_SCHEMA

ROLES = ('electron_temperature', 'ion_temperature', 'electron_density',
         'rotation', 'composition', 'heating')
ORIGINS = ('measured', 'fitted', 'manual', 'automatic', 'assumed', 'inferred',
           'predicted', 'missing', 'unknown')
DEFAULT_THRESHOLDS = {'quasineutrality_relative': 0.02,
                      'weak_gradient_relative': 0.001, 'gradient_sigma_multiplier': 2.0}


def _keys(obj, allowed, where):
    if not isinstance(obj, dict):
        raise ValueError(f'{where}: expected an object')
    unknown = set(obj) - set(allowed)
    if unknown:
        raise ValueError(f'{where}: unknown keys: {", ".join(sorted(unknown))}')


def _number(value, where, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f'{where}: expected a finite number')
    if positive and value <= 0:
        raise ValueError(f'{where}: must be positive')


def validate_manifest(document):
    """Return an independent normalized copy; reject typos and ambiguous contracts."""
    _keys(document, ('schema_version', 'assessment_window', 'evidence', 'thresholds',
                     'composition', 'input_expectations', 'gradient_uncertainties', 'notes'), 'manifest')
    if document.get('schema_version') != MANIFEST_SCHEMA:
        raise ValueError(f'manifest.schema_version must be {MANIFEST_SCHEMA!r}')
    result = copy.deepcopy(document)
    window = result.setdefault('assessment_window', {})
    _keys(window, ('time_min', 'time_max', 'rho_min', 'rho_max'), 'assessment_window')
    for key, value in window.items():
        _number(value, f'assessment_window.{key}')
    for prefix in ('time', 'rho'):
        lo, hi = window.get(prefix + '_min'), window.get(prefix + '_max')
        if lo is not None and hi is not None and lo > hi:
            raise ValueError(f'assessment_window: {prefix}_min exceeds {prefix}_max')
    for key in ('rho_min', 'rho_max'):
        if key in window and not 0 <= window[key] <= 1:
            raise ValueError(f'assessment_window.{key}: expected normalized radius in [0, 1]')
    evidence = result.setdefault('evidence', {})
    _keys(evidence, ROLES, 'evidence')
    for role, item in evidence.items():
        _keys(item, ('origin', 'source', 'uncertainty', 'roles', 'reliability',
                     'independent_validation', 'notes'), 'evidence.' + role)
        if item.get('origin', 'unknown') not in ORIGINS:
            raise ValueError(f'evidence.{role}: unsupported origin')
        if item.get('uncertainty', 'unknown') not in ('provided', 'unknown', 'not_applicable'):
            raise ValueError(f'evidence.{role}: unsupported uncertainty status')
        if item.get('reliability', 'unknown') not in ('reviewed', 'unreliable', 'unknown'):
            raise ValueError(f'evidence.{role}: unsupported reliability status')
        for key in ('source', 'notes'):
            if key in item and not isinstance(item[key], str):
                raise ValueError(f'evidence.{role}.{key}: expected text')
        roles = item.get('roles', [])
        if not isinstance(roles, list) or any(r not in ('input', 'calibration', 'validation',
                                                     'initial_condition', 'boundary_condition') for r in roles):
            raise ValueError(f'evidence.{role}.roles: unsupported role')
        if 'independent_validation' in item and not isinstance(item['independent_validation'], bool):
            raise ValueError(f'evidence.{role}.independent_validation: expected boolean')
    thresholds = result.setdefault('thresholds', {})
    _keys(thresholds, (*DEFAULT_THRESHOLDS, 'newton_iterations_warning', 'residual_limits'), 'thresholds')
    for key, default in DEFAULT_THRESHOLDS.items():
        thresholds.setdefault(key, default)
    for key, value in thresholds.items():
        if key != 'residual_limits':
            _number(value, 'thresholds.' + key, positive=True)
    limits = thresholds.setdefault('residual_limits', {})
    _keys(limits, ('RESL2TE', 'RESL2TI', 'RESL2OMG'), 'thresholds.residual_limits')
    for name, spec in limits.items():
        _keys(spec, ('limit', 'rationale'), 'residual_limits.' + name)
        _number(spec.get('limit'), name + '.limit', positive=True)
        if not isinstance(spec.get('rationale'), str) or not spec['rationale'].strip():
            raise ValueError(f'{name}: a residual threshold requires its normalization/rationale')
    comp = result.get('composition')
    if comp is not None:
        _keys(comp, ('electron_density', 'species', 'complete', 'includes_fast_ions', 'source'), 'composition')
        for flag in ('complete', 'includes_fast_ions'):
            if flag in comp and not isinstance(comp[flag], bool):
                raise ValueError(f'composition.{flag}: expected boolean')
        if not isinstance(comp.get('electron_density', 'NE'), str):
            raise ValueError('composition.electron_density: expected a CDF variable name')
        if not isinstance(comp.get('source', ''), str):
            raise ValueError('composition.source: expected text')
        species = comp.get('species', [])
        if not isinstance(species, list):
            raise ValueError('composition.species: expected a list')
        seen = set()
        for item in species:
            _keys(item, ('density', 'charge', 'charge_variable'), 'composition.species')
            name = item.get('density')
            if not isinstance(name, str) or not name or name in seen:
                raise ValueError('composition.species: density names must be nonempty and unique')
            if name == comp.get('electron_density', 'NE'):
                raise ValueError('electron density cannot also be a positive-ion species')
            seen.add(name)
            if ('charge' in item) == ('charge_variable' in item):
                raise ValueError('each species needs exactly one of charge or charge_variable')
            if 'charge' in item:
                _number(item['charge'], name + '.charge', positive=True)
            elif not isinstance(item['charge_variable'], str) or not item['charge_variable']:
                raise ValueError('charge_variable must be a nonempty CDF name')
    expectations = result.setdefault('input_expectations', [])
    if not isinstance(expectations, list):
        raise ValueError('input_expectations: expected a list')
    for item in expectations:
        _keys(item, ('name', 'expected', 'atol', 'rtol', 'reference', 'units'), 'input_expectations')
        if not all(isinstance(item.get(k), str) and item[k].strip() for k in ('name', 'reference', 'units')):
            raise ValueError('input_expectations need name, reference and native units (or "1")')
        _number(item.get('expected'), 'input_expectations.expected')
        for key in ('atol', 'rtol'):
            item.setdefault(key, 0.0)
            _number(item[key], key)
            if item[key] < 0:
                raise ValueError(f'{key}: expected a nonnegative tolerance')
    gradients = result.setdefault('gradient_uncertainties', {})
    _keys(gradients, ('TE', 'TI'), 'gradient_uncertainties')
    for name, spec in gradients.items():
        _keys(spec, ('sigma_variable', 'coordinate', 'method', 'source'), 'gradient_uncertainties.' + name)
        if spec.get('coordinate') != 'X' or spec.get('method') != 'numpy_gradient_edge_order_1':
            raise ValueError('gradient uncertainty must refer to d(profile)/dX using numpy_gradient_edge_order_1')
        if not all(isinstance(spec.get(k), str) and spec[k] for k in ('sigma_variable', 'source')):
            raise ValueError('gradient uncertainty needs sigma_variable and source')
    if 'notes' in result and not isinstance(result['notes'], str):
        raise ValueError('notes must be text')
    return result


def load_manifest(path=None):
    document = {'schema_version': MANIFEST_SCHEMA}
    if path is not None:
        def reject_constant(value):
            raise ValueError(f'nonfinite JSON constant: {value}')
        def unique_object(pairs):
            result = {}
            for key, value in pairs:
                if key in result:
                    raise ValueError(f'duplicate JSON key: {key}')
                result[key] = value
            return result
        document = json.loads(Path(path).read_text(encoding='utf-8'), parse_constant=reject_constant,
                              object_pairs_hook=unique_object)
    return validate_manifest(document)
