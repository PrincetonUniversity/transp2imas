#!/usr/bin/python
"""Strict validation for a small, explicit evidence sidecar; no assumed defaults.

Evidence strings are analyst declarations, not independently verified facts.
"""

import copy
import json
import math
from pathlib import Path

import hashlib
import re

from curation import ACCEPTED_MANIFEST_SCHEMAS, MANIFEST_SCHEMA, MANIFEST_SCHEMA_V1, MANIFEST_SCHEMA_V2
from curation.schema import loads_strict, reject_nonfinite, validate_structure

ROLES = ('electron_temperature', 'ion_temperature', 'electron_density',
         'rotation', 'composition', 'heating')
ORIGINS = ('measured', 'fitted', 'manual', 'automatic', 'assumed', 'inferred',
           'predicted', 'missing', 'unknown')
# What evidence did for the run (v2).  The legacy v1 ``roles`` map onto these: input -> constraint,
# calibration -> calibration_target, validation -> independent_validation.
USAGE = ('constraint', 'initial_condition', 'boundary_condition', 'calibration_target', 'independent_validation')
LEGACY_USAGE = {'input': 'constraint', 'calibration': 'calibration_target', 'validation': 'independent_validation',
                'initial_condition': 'initial_condition', 'boundary_condition': 'boundary_condition'}
CONSTRAINING = ('constraint', 'initial_condition', 'boundary_condition', 'calibration_target')
EVENT_TYPES = ('startup', 'termination', 'sawtooth', 'elm', 'pellet', 'heating_transition', 'diagnostic_outage',
               'unclassified_transient')
BOUNDARIES = ('injected', 'launched', 'absorbed', 'other')
INPUT_ROLES = ('trdat', 'namelist', 'pt_template', 'ufile', 'other')
V2_TOP_LEVEL = ('run_identity', 'evidence_sources', 'resolved_inputs', 'events')
V3_TOP_LEVEL = ('sample_support', 'event_detectors', 'uncertainty', 'extensions')
SCHEMA_OF = {MANIFEST_SCHEMA_V1: 'evidence.v1', MANIFEST_SCHEMA_V2: 'evidence.v2', MANIFEST_SCHEMA: 'evidence.v3'}
MIGRATION_KEY = 'migration'       # written by this module; a raw v1 document may not carry it
EVENT_ORIGINS = ('experimental_observation', 'actuator_command', 'modeled_transp', 'analyst_declaration',
                 'automatic_candidate')
V2_EVIDENCE = ('id', 'parents', 'usage', 'applicability', 'revision', 'sha256', 'units', 'power_accounting_boundary',
               'override', 'selector')
V2_EXPECTATION = ('applicability', 'evidence_id', 'revision', 'sha256', 'override', 'power_accounting_boundary')
ID_PATTERN = re.compile(r'[A-Za-z0-9_.:-]+\Z')
SHA256 = re.compile(r'[0-9a-fA-F]{64}\Z')
DEFAULT_THRESHOLDS = {'quasineutrality_relative': 0.02,
                      'weak_gradient_relative': 0.001, 'gradient_sigma_multiplier': 2.0}


def _number(value, where, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f'{where}: expected a finite number')
    if positive and value <= 0:
        raise ValueError(f'{where}: must be positive')


def _text(value, where, *, allow_empty=False):
    if not isinstance(value, str) or (not allow_empty and not value.strip()):
        raise ValueError(f'{where}: expected nonempty text')


def _identifier(value, where):
    if not isinstance(value, str) or not ID_PATTERN.match(value):
        raise ValueError(f'{where}: expected an identifier of letters, digits and _ . : -')


def _hash(value, where):
    if not isinstance(value, str) or not SHA256.match(value):
        raise ValueError(f'{where}: expected a 64-digit hexadecimal SHA-256')


def _validate_applicability(spec, where):
    """Where a configuration or reference applies: device, campaign, a shot range and a time window (all optional)."""
    for key in ('device', 'campaign'):
        if key in spec:
            _text(spec[key], f'{where}.{key}')
    for key in ('shot_min', 'shot_max', 'time_min', 'time_max'):
        if key in spec:
            _number(spec[key], f'{where}.{key}')
    for lo, hi in (('shot_min', 'shot_max'), ('time_min', 'time_max')):
        if lo in spec and hi in spec and spec[lo] > spec[hi]:
            raise ValueError(f'{where}: {lo} exceeds {hi}')


def _validate_override(spec, where):
    """A reviewed, deliberate deviation: the value, why, and who reviewed it.  None of the three may be left out."""
    if 'value' not in spec or isinstance(spec['value'], bool) or not isinstance(spec['value'], (int, float, str)):
        raise ValueError(f'{where}.value: expected a number or text')
    if isinstance(spec['value'], float) and not math.isfinite(spec['value']):
        raise ValueError(f'{where}.value: expected a finite number')
    _text(spec.get('reason'), f'{where}.reason')
    _text(spec.get('reviewer'), f'{where}.reviewer')


def _validate_reference_v2(item, where):
    """Fields shared by an evidence record and an input expectation (v2)."""
    if 'applicability' in item:
        _validate_applicability(item['applicability'], where + '.applicability')
    if 'override' in item:
        _validate_override(item['override'], where + '.override')
    if 'revision' in item:
        _text(item['revision'], where + '.revision')
    if 'sha256' in item:
        _hash(item['sha256'], where + '.sha256')
    if 'evidence_id' in item:
        _identifier(item['evidence_id'], where + '.evidence_id')
    if 'power_accounting_boundary' in item and item['power_accounting_boundary'] not in BOUNDARIES:
        raise ValueError(f'{where}.power_accounting_boundary: expected one of {BOUNDARIES}')


def _selector(spec, where):
    """Resource selector: which dataset/channel/fit/calibration inside a container the record refers to."""
    if not isinstance(spec, dict) or not spec:
        raise ValueError(f'{where}: expected a nonempty object naming the dataset, channel, fit or calibration')
    for key, value in spec.items():
        if not ID_PATTERN.match(key) or isinstance(value, bool) or not isinstance(value, (str, int, float)):
            raise ValueError(f'{where}.{key}: expected a text or number value')
        if isinstance(value, float) and not math.isfinite(value):
            raise ValueError(f'{where}.{key}: expected a finite number')


def _validate_evidence_v2(role, item):
    where = 'evidence.' + role
    _validate_reference_v2(item, where)
    if 'id' in item:
        _identifier(item['id'], where + '.id')
    if 'parents' in item:
        if not isinstance(item['parents'], list):
            raise ValueError(f'{where}.parents: expected a list of evidence ids')
        for parent in item['parents']:
            _identifier(parent, where + '.parents')
    if 'units' in item:
        _text(item['units'], where + '.units')
    if 'selector' in item:
        _selector(item['selector'], where + '.selector')
    if 'usage' in item:
        if not isinstance(item['usage'], list) or any(u not in USAGE for u in item['usage']):
            raise ValueError(f'{where}.usage: expected a list drawn from {USAGE}')


def _validate_v2_sections(result):
    """Run identity, evidence sources, resolved inputs and events (v2), and the evidence graph's references."""
    identity = result.setdefault('run_identity', {})
    for key in ('device', 'campaign'):
        if key in identity:
            _text(identity[key], 'run_identity.' + key)
    if 'shot' in identity:
        _number(identity['shot'], 'run_identity.shot')

    sources = result.setdefault('evidence_sources', [])
    if not isinstance(sources, list):
        raise ValueError('evidence_sources: expected a list')
    ids = {}
    for index, node in enumerate(sources):
        where = f'evidence_sources[{index}]'
        _identifier(node.get('id'), where + '.id')
        _text(node.get('description'), where + '.description')
        if 'source' in node:
            _text(node['source'], where + '.source')
        if 'sha256' in node:
            _hash(node['sha256'], where + '.sha256')
        if 'selector' in node:
            _selector(node['selector'], where + '.selector')
        node.setdefault('parents', [])
        if not isinstance(node['parents'], list):
            raise ValueError(f'{where}.parents: expected a list of evidence ids')
        for parent in node['parents']:
            _identifier(parent, where + '.parents')
        if any(u not in USAGE for u in node.setdefault('usage', [])):
            raise ValueError(f'{where}.usage: expected a list drawn from {USAGE}')
        if node['id'] in ids:
            raise ValueError(f'{where}.id: duplicate evidence id {node["id"]!r}')
        ids[node['id']] = where
    for role, item in result.get('evidence', {}).items():
        record_id = item.get('id', 'evidence:' + role)
        if record_id in ids:
            raise ValueError(f'evidence.{role}.id: duplicate evidence id {record_id!r}')
        ids[record_id] = 'evidence.' + role
    graph = {}
    for node in sources:
        graph[node['id']] = node['parents']
    for role, item in result.get('evidence', {}).items():
        if 'parents' in item:
            graph[item.get('id', 'evidence:' + role)] = item['parents']
    for node_id, parents in graph.items():
        for parent in parents:
            if parent not in ids:
                raise ValueError(f'{ids[node_id]}: parent {parent!r} is not a declared evidence id')
    _reject_cycles(graph)
    for expectation in result.get('input_expectations', []):
        if 'evidence_id' in expectation and expectation['evidence_id'] not in ids:
            raise ValueError(f'input_expectations: evidence_id {expectation["evidence_id"]!r} is not a declared evidence id')

    resolved = result.setdefault('resolved_inputs', {})
    if 'code' in resolved:
        for key in ('name', 'version'):
            _text(resolved['code'].get(key), 'resolved_inputs.code.' + key)
        if 'executable_sha256' in resolved['code']:
            _hash(resolved['code']['executable_sha256'], 'resolved_inputs.code.executable_sha256')
        if 'source' in resolved['code']:
            _text(resolved['code']['source'], 'resolved_inputs.code.source')
    inputs = resolved.setdefault('inputs', [])
    if not isinstance(inputs, list):
        raise ValueError('resolved_inputs.inputs: expected a list')
    for index, item in enumerate(inputs):
        where = f'resolved_inputs.inputs[{index}]'
        if item.get('role') not in INPUT_ROLES:
            raise ValueError(f'{where}.role: expected one of {INPUT_ROLES}')
        _text(item.get('path'), where + '.path')
        for key in ('revision', 'source'):
            if key in item:
                _text(item[key], f'{where}.{key}')
        if 'sha256' in item:
            _hash(item['sha256'], where + '.sha256')
    controls = resolved.setdefault('controls', [])
    if not isinstance(controls, list):
        raise ValueError('resolved_inputs.controls: expected a list')
    seen = set()
    for index, control in enumerate(controls):
        where = f'resolved_inputs.controls[{index}]'
        _text(control.get('name'), where + '.name')
        key = re.sub(r'\s+', '', control['name']).upper()
        if key in seen:
            raise ValueError(f'{where}.name: duplicate control {control["name"]!r}')
        seen.add(key)
        if not isinstance(control.get('resolved'), bool):
            raise ValueError(f'{where}.resolved: state whether this is the resolved runtime value (true) or a declaration (false)')
        segments = control.get('segments')
        if not isinstance(segments, list) or not segments:
            raise ValueError(f'{where}.segments: expected a nonempty list of time segments')
        previous_end = None
        for number, segment in enumerate(segments):
            at = f'{where}.segments[{number}]'
            for edge in ('start', 'end'):
                if segment.get(edge) is not None:
                    _number(segment[edge], f'{at}.{edge}')
            if 'value' not in segment or isinstance(segment['value'], (list, dict)) or segment['value'] is None:
                raise ValueError(f'{at}.value: expected one scalar value')
            lo, hi = segment.get('start'), segment.get('end')
            if lo is not None and hi is not None and lo >= hi:
                raise ValueError(f'{at}: start must be before end')
            if previous_end is not None and (lo is None or lo < previous_end):
                raise ValueError(f'{at}: segments must be in time order and must not overlap')
            previous_end = hi if hi is not None else math.inf

    events = result.setdefault('events', [])
    if not isinstance(events, list):
        raise ValueError('events: expected a list')
    event_ids = set()
    for index, event in enumerate(events):
        where = f'events[{index}]'
        _identifier(event.get('id'), where + '.id')
        if event['id'] in event_ids:
            raise ValueError(f'{where}.id: duplicate event id {event["id"]!r}')
        event_ids.add(event['id'])
        if event.get('type') not in EVENT_TYPES:
            raise ValueError(f'{where}.type: expected one of {EVENT_TYPES}')
        _number(event.get('start'), where + '.start')
        _number(event.get('end'), where + '.end')
        if event['start'] > event['end']:
            raise ValueError(f'{where}: start exceeds end')
        rules = event.get('rules')
        if not isinstance(rules, list) or not rules or any(not isinstance(r, str) or not r for r in rules):
            raise ValueError(f'{where}.rules: name the rules (shell-style patterns over rule ids) the event applies to; '
                             'an event that applies to nothing is not a declaration')
        for key in ('source', 'note'):
            if key in event:
                _text(event[key], f'{where}.{key}')
        if event.setdefault('origin', 'analyst_declaration') not in EVENT_ORIGINS:
            raise ValueError(f'{where}.origin: expected one of {EVENT_ORIGINS}')
        # The legacy ``rules`` list is an analyst exclusion *policy*, not an assertion about event physics.
        legacy = {'kind': 'legacy_exclusion', 'rules': list(rules)}
        if event.setdefault('policy', legacy) != legacy:
            raise ValueError(f'{where}.policy: only the legacy exclusion policy over the declared rules is supported')
        if 'applications' in event and (not isinstance(event['applications'], list) or not event['applications']
                                        or any(not isinstance(a, str) for a in event['applications'])):
            raise ValueError(f'{where}.applications: expected a nonempty list of application names')


    # --- v3: sample support, detectors, covariance resources, extensions (structure is the schema's job)
    result.setdefault('sample_support', {})
    result.setdefault('extensions', {})
    detector_ids = set()
    for index, detector in enumerate(result.setdefault('event_detectors', [])):
        if detector['id'] in detector_ids:
            raise ValueError(f'event_detectors[{index}].id: duplicate detector id {detector["id"]!r}')
        detector_ids.add(detector['id'])
        params = detector.get('parameters', {})
        if 'window' in params and 'min_baseline' in params and params['min_baseline'] >= params['window']:
            raise ValueError(f'event_detectors[{index}].parameters: min_baseline must be smaller than window')
    covariances = result.setdefault('uncertainty', {}).setdefault('covariances', [])
    seen_ids, seen_quantities = set(), set()
    for index, spec in enumerate(covariances):
        if spec['id'] in seen_ids:
            raise ValueError(f'uncertainty.covariances[{index}].id: duplicate id {spec["id"]!r}')
        if spec['quantity'] in seen_quantities:
            raise ValueError(f'uncertainty.covariances[{index}].quantity: {spec["quantity"]} already has a covariance; '
                             'two covariances for one quantity would be ambiguous (merge them upstream)')
        seen_ids.add(spec['id'])
        seen_quantities.add(spec['quantity'])


def _reject_cycles(graph):
    """Raise on a cycle among evidence ids (a source cannot be its own ancestor)."""
    state = {}

    def visit(node, path):
        if state.get(node) == 'done':
            return
        if state.get(node) == 'active':
            raise ValueError('evidence parents form a cycle: ' + ' -> '.join(path + [node]))
        state[node] = 'active'
        for parent in graph.get(node, ()):
            visit(parent, path + [node])
        state[node] = 'done'

    for node in graph:
        visit(node, [])


def _canonical_sha256(document):
    return hashlib.sha256(json.dumps(document, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _migration_record(document, raw_sha256):
    """The migration history of a canonical (v3) manifest.

    A raw v1 or v2 document is migrated by an explicit, listed set of steps; a canonical document that already
    carries a record keeps it unchanged, which makes normalization idempotent.
    """
    if MIGRATION_KEY in document:
        return copy.deepcopy(document[MIGRATION_KEY])             # structure already checked by the v3 schema
    steps = []
    if document['schema_version'] == MANIFEST_SCHEMA_V1:
        steps += [f'{MANIFEST_SCHEMA_V1} -> {MANIFEST_SCHEMA_V2}: added empty run_identity, evidence_sources, '
                  'resolved_inputs and events sections',
                  'evidence usage derived from the v1 roles and independent_validation fields']
    if document['schema_version'] in (MANIFEST_SCHEMA_V1, MANIFEST_SCHEMA_V2):
        steps += [f'{MANIFEST_SCHEMA_V2} -> {MANIFEST_SCHEMA}: added empty sample_support, event_detectors, uncertainty '
                  'and extensions sections',
                  'declared events keep their legacy exclusion policy and take origin analyst_declaration']
    return {'from': document['schema_version'], 'raw_sha256': raw_sha256 or _canonical_sha256(document), 'steps': steps}


def validate_raw(document):
    """Stages 1-2, without touching ``document``: finite numbers, declared version, structure and compatibility.

    A document is judged against *its own* declared schema.  A v1 document may not use v2 fields and a v2 document
    may not use v3 fields (or a migration record): it is migrated by this module, not extended in place.  Structure
    comes from the packaged JSON Schema; the cross-field checks that follow are semantic, not structural.
    """
    if not isinstance(document, dict):
        raise ValueError('manifest: expected an object')
    reject_nonfinite(document, 'manifest')
    version = document.get('schema_version')
    if version not in ACCEPTED_MANIFEST_SCHEMAS:
        raise ValueError(f'manifest.schema_version must be one of {ACCEPTED_MANIFEST_SCHEMAS}')
    newer = {MANIFEST_SCHEMA_V1: V2_TOP_LEVEL + V3_TOP_LEVEL + (MIGRATION_KEY,),
             MANIFEST_SCHEMA_V2: V3_TOP_LEVEL + (MIGRATION_KEY,), MANIFEST_SCHEMA: ()}[version]
    used = [key for key in newer if key in document]
    if used:
        raise ValueError(f'manifest: {", ".join(used)} are not part of schema {version}; declare {MANIFEST_SCHEMA} '
                         '(an older document is migrated by this module, not by adding newer fields to it)')
    validate_structure(document, SCHEMA_OF[version])


def _resolve_references(manifest, base_dir):
    """Add ``resolved_path`` next to each declared file reference; the declared ``path`` is kept as written.

    Relative paths resolve against ``base_dir`` (the sidecar's directory), never the working directory, so a sidecar
    means the same thing wherever it is run.  An already-resolved reference is left alone (idempotent).
    """
    if base_dir is None:
        return
    references = list(manifest['resolved_inputs'].get('inputs', [])) + list(manifest['uncertainty']['covariances'])
    for item in references:
        if 'resolved_path' not in item:
            path = Path(item['path'])
            item['resolved_path'] = str(path if path.is_absolute() else Path(base_dir) / path)


def referenced_paths(manifest):
    """Every file the evidence refers to (resolved where the sidecar location was known): protected inputs."""
    paths = []
    for item in list(manifest['resolved_inputs'].get('inputs', [])) + list(manifest['uncertainty']['covariances']):
        paths.append(Path(item.get('resolved_path', item['path'])))
    return paths


def validate_resources(manifest):
    """Stage 4, resource validation: covariance resources must exist and match their declared SHA-256.

    A missing or substituted file is an error here, before any assessment, not a quiet ``not_assessed``.  Input files
    declared under ``resolved_inputs`` are compared by the rules instead (their disagreement is a finding).
    """
    for spec in manifest['uncertainty']['covariances']:
        path = Path(spec.get('resolved_path', spec['path']))
        if not path.is_file():
            raise ValueError(f'uncertainty.covariances[{spec["id"]}].path: {spec["path"]!r} is not a file '
                             f'(resolved to {str(path)!r})')
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != spec['sha256'].lower():
            raise ValueError(f'uncertainty.covariances[{spec["id"]}]: sha256 {digest} of {str(path)!r} does not match '
                             f'the declared {spec["sha256"]}')


def validate_manifest(document, raw_sha256=None, base_dir=None):
    """Return an independent canonical (v3) copy; reject typos and ambiguous contracts.

    Pipeline: ``validate_raw`` (finite numbers, declared version, JSON Schema) -> migration to the canonical v3
    representation with a migration record -> the semantic checks below (ordering, graph, units and numerical
    constraints no schema can express).  Applying this to its own result changes nothing, and the caller's document
    is never modified.  ``raw_sha256`` is the hash of the file the document was read from, when there is one;
    ``base_dir`` is the directory file references resolve against (the sidecar's own).
    """
    validate_raw(document)
    result = copy.deepcopy(document)
    result.pop('$schema', None)                       # an editor hint, not evidence
    window = result.setdefault('assessment_window', {})
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
    for role, item in evidence.items():
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
        _validate_evidence_v2(role, item)
        usage = {LEGACY_USAGE[r] for r in roles} | set(item.get('usage', ()))
        if item.get('independent_validation'):
            usage.add('independent_validation')
        item['usage'] = sorted(usage)                       # one vocabulary downstream, whichever spelling was used
    thresholds = result.setdefault('thresholds', {})
    for key, default in DEFAULT_THRESHOLDS.items():
        thresholds.setdefault(key, default)
    for key, value in thresholds.items():
        if key != 'residual_limits':
            _number(value, 'thresholds.' + key, positive=True)
    limits = thresholds.setdefault('residual_limits', {})
    for name, spec in limits.items():
        _number(spec.get('limit'), name + '.limit', positive=True)
        if not isinstance(spec.get('rationale'), str) or not spec['rationale'].strip():
            raise ValueError(f'{name}: a residual threshold requires its normalization/rationale')
    comp = result.get('composition')
    if comp is not None:
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
        _validate_reference_v2(item, 'input_expectations')
        if not all(isinstance(item.get(k), str) and item[k].strip() for k in ('name', 'reference', 'units')):
            raise ValueError('input_expectations need name, reference and native units (or "1")')
        _number(item.get('expected'), 'input_expectations.expected')
        for key in ('atol', 'rtol'):
            item.setdefault(key, 0.0)
            _number(item[key], key)
            if item[key] < 0:
                raise ValueError(f'{key}: expected a nonnegative tolerance')
    gradients = result.setdefault('gradient_uncertainties', {})
    for name, spec in gradients.items():
        if spec.get('coordinate') != 'X' or spec.get('method') != 'numpy_gradient_edge_order_1':
            raise ValueError('gradient uncertainty must refer to d(profile)/dX using numpy_gradient_edge_order_1')
        if not all(isinstance(spec.get(k), str) and spec[k] for k in ('sigma_variable', 'source')):
            raise ValueError('gradient uncertainty needs sigma_variable and source')
    if 'notes' in result and not isinstance(result['notes'], str):
        raise ValueError('notes must be text')
    _validate_v2_sections(result)
    _resolve_references(result, base_dir)
    result[MIGRATION_KEY] = _migration_record(document, raw_sha256)
    result['schema_version'] = MANIFEST_SCHEMA
    return result


def read_json(path):
    """Strict JSON: duplicate keys, NaN/Infinity and numeric overflow (``1e999``) are rejected."""
    return loads_strict(Path(path).read_text(encoding='utf-8'))


def load_manifest(path=None):
    """Read, validate and migrate an evidence sidecar; file references resolve against its directory."""
    if path is None:
        return validate_manifest({'schema_version': MANIFEST_SCHEMA})
    path = Path(path)
    manifest = validate_manifest(read_json(path), raw_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                 base_dir=path.resolve().parent)
    validate_resources(manifest)
    return manifest
