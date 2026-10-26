#!/usr/bin/python
"""Offline JSON Schema validation (Draft 2020-12) of curation documents, and strict JSON decoding.

Every schema ships in ``curation/schemas`` and is resolved through a local ``referencing`` registry: nothing is
fetched from the network.  Validation is *structural* only: a schema-valid sidecar is not authenticated evidence,
and cross-field, unit, coordinate and graph checks happen in later stages (see ``curation.manifest``).

Decoding rejects what a schema type alone would let through in Python: duplicate keys, ``NaN``/``Infinity``
constants and numeric overflow such as ``1e999``; ``reject_nonfinite`` covers programmatically built documents.
"""

import json
import math
from pathlib import Path

try:
    import jsonschema
    from referencing import Registry, Resource
    from referencing.jsonschema import DRAFT202012
except ImportError as exc:                                          # pragma: no cover - environment dependent
    raise ImportError('curation schema validation needs jsonschema>=4.18 and referencing; install the curation '
                      'requirements: pip install -r requirements-curation.txt') from exc

SCHEMA_DIR = Path(__file__).resolve().parent / 'schemas'
URN = 'urn:transp2imas:curation:'
SCHEMA_NAMES = ('evidence.v1', 'evidence.v2', 'evidence.v3', 'pilot.v1', 'policy.v1', 'report.v3')
_cache = {}


def _load():
    if not _cache:
        documents = {path.stem.replace('.schema', ''): json.loads(path.read_text(encoding='utf-8'))
                     for path in sorted(SCHEMA_DIR.glob('*.schema.json'))}
        _cache['documents'] = documents
        _cache['registry'] = Registry().with_resources(
            (document['$id'], Resource.from_contents(document, default_specification=DRAFT202012))
            for document in documents.values())
    return _cache


def schema_document(name):
    return _load()['documents'][name]


def validator(name):
    cache = _load()
    return jsonschema.Draft202012Validator(cache['documents'][name], registry=cache['registry'])


def _path(error):
    parts = []
    for item in error.absolute_path:
        parts.append(f'[{item}]' if isinstance(item, int) else f'.{item}')
    return ''.join(parts).lstrip('.') or '<document>'


def describe(name, error):
    """One line: schema, JSON path, the offending value and type, and the failed condition."""
    value = error.instance
    shown = repr(value)
    shown = shown if len(shown) <= 80 else shown[:77] + '...'
    return (f'[schema {name}] {_path(error)}: {error.message} '
            f'(value {shown} of type {type(value).__name__}; failed {error.validator!r}: '
            f'{json.dumps(error.validator_value, default=str)[:100]})')


def validate_structure(document, name):
    """Raise ``ValueError`` listing every structural violation of ``document`` against schema ``name``."""
    errors = sorted(validator(name).iter_errors(document), key=lambda e: [str(p) for p in e.absolute_path])
    if errors:
        # report the most specific error of each failing location, not the whole oneOf/anyOf tree
        best = {}
        for error in errors:
            key = tuple(str(p) for p in error.absolute_path)
            best.setdefault(key, jsonschema.exceptions.best_match([error]))
        raise ValueError('; '.join(describe(name, e) for e in list(best.values())[:8])
                         + ('' if len(best) <= 8 else f'; and {len(best) - 8} more'))


def reject_nonfinite(value, where='document'):
    """Raise on a NaN or infinite number anywhere in a programmatically built document."""
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f'{where}: nonfinite number {value!r}')
    if isinstance(value, dict):
        for key, item in value.items():
            reject_nonfinite(item, f'{where}.{key}')
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            reject_nonfinite(item, f'{where}[{index}]')


def _reject_constant(name):
    raise ValueError(f'nonfinite JSON constant: {name}')


def _parse_float(text):
    value = float(text)
    if not math.isfinite(value):
        raise ValueError(f'numeric overflow in JSON number {text}')
    return value


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'duplicate JSON key: {key}')
        result[key] = value
    return result


def loads_strict(text):
    """Strict JSON decoding: duplicate keys, ``NaN``/``Infinity`` and overflowing numbers are errors."""
    return json.loads(text, parse_constant=_reject_constant, parse_float=_parse_float, object_pairs_hook=_unique_object)
