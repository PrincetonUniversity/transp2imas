#!/usr/bin/python
"""Pilot adapter: analyst-supplied facts about one run -> a v2 evidence manifest.

Nothing is invented.  The facts file carries what a person established (run identity, code version, resolved
controls, evidence origins and lineage, declared events, campaign references).  The adapter adds only what can be
computed from files on disk, the SHA-256 of each listed input, and refuses a facts file that still holds a
``TO_BE_SUPPLIED`` placeholder, so a template cannot be run by accident.  The result goes through the same
``validate_manifest`` as any sidecar.
"""

from pathlib import Path

from curation import MANIFEST_SCHEMA
from curation.manifest import read_json, validate_manifest
from curation.resolved import sha256_of

PILOT_SCHEMA = 'transp2imas.curation.pilot.v1'
PLACEHOLDER = 'TO_BE_SUPPLIED'
PASSED_THROUGH = ('assessment_window', 'evidence', 'evidence_sources', 'events', 'input_expectations', 'composition',
                  'thresholds', 'gradient_uncertainties', 'run_identity', 'notes')


def placeholders(value, where='facts'):
    """JSON paths of every placeholder string in a loaded document."""
    if isinstance(value, str):
        return [where] if PLACEHOLDER in value else []
    if isinstance(value, dict):
        return [p for key, item in value.items() for p in placeholders(item, f'{where}.{key}')]
    if isinstance(value, list):
        return [p for number, item in enumerate(value) for p in placeholders(item, f'{where}[{number}]')]
    return []


def build_manifest(facts, base):
    """The v2 manifest for ``facts`` (a loaded document); relative file paths resolve against ``base``."""
    if facts.get('schema_version') != PILOT_SCHEMA:
        raise ValueError(f'pilot facts schema_version must be {PILOT_SCHEMA!r}')
    unknown = set(facts) - set(PASSED_THROUGH) - {'schema_version', 'files', 'resolved_inputs'}
    if unknown:
        raise ValueError(f'pilot facts: unknown keys: {", ".join(sorted(unknown))}')
    open_items = placeholders(facts)
    if open_items:
        raise ValueError(f'pilot facts still contain {PLACEHOLDER}: ' + ', '.join(open_items))
    manifest = {'schema_version': MANIFEST_SCHEMA}
    manifest.update({key: facts[key] for key in PASSED_THROUGH if key in facts})
    resolved = dict(facts.get('resolved_inputs', {}))
    inputs = list(resolved.get('inputs', []))
    for item in facts.get('files', []):
        if not isinstance(item, dict) or set(item) - {'role', 'path', 'revision', 'source'}:
            raise ValueError('pilot facts: files entries take role, path, revision, source only (the hash is computed)')
        path = Path(item['path'])
        path = path if path.is_absolute() else Path(base) / path
        inputs.append({**item, 'path': str(path), 'sha256': sha256_of(path)})      # an absent file is an OSError
    if inputs:
        resolved['inputs'] = inputs
    if resolved:
        manifest['resolved_inputs'] = resolved
    return validate_manifest(manifest)


def load_pilot(path):
    """Load a pilot facts file (strict JSON) and build its manifest."""
    path = Path(path)
    return build_manifest(read_json(path), path.resolve().parent)
