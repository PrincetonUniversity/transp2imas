#!/usr/bin/python
"""Human-readable and machine-readable curation reports; never write source data."""

import json
import os
from pathlib import Path
import tempfile

import numpy as np


def _escape(value):
    return str(value).replace('|', '\\|').replace('\n', ' ').replace('\r', ' ')


def markdown_report(document):
    """Markdown report with evidence gaps and native-coordinate scope."""
    summary = document['summary']
    counts = summary['counts']
    source = document['source']
    lines = ['# TRANSP curation screening report', '',
             '**This report is not a scientific certification or an IMAS write/readback test.**', '',
             f"Source: `{source.get('path') or 'in-memory dataset'}`  ",
             f"Source SHA-256: `{source.get('sha256') or 'not available'}`  ",
             f"Evaluator: `{document['curation_version']}`; schema: `{document['schema_version']}`", '',
             '## Screening summary', '',
             '| Status | Findings |', '|---|---:|']
    lines += [f'| {name} | {number} |' for name, number in counts.items()]
    lines += ['', '## Application-specific preliminary status', '',
              '| Application | Status |', '|---|---|']
    lines += [f'| {name} | {status} |' for name, status in summary['applications'].items()]
    lines += ['', 'Missing uncertainties are not zero uncertainties. A passed accounting check does not validate',
              'a diagnostic, fitted profile, heating configuration, or transport model.', '',
              '## Findings', '', '| Rule | Status | Variables | Finding |', '|---|---|---|---|']
    for finding in document['findings']:
        lines.append(f"| `{finding['id']}` | {finding['status']} | {_escape(', '.join(finding['variables']))} | {_escape(finding['message'])} |")
    lines += ['', '## Details requiring review', '']
    for finding in document['findings']:
        if finding['status'] in ('pass', 'not_applicable'):
            continue
        lines += [f"### {finding['id']} ({finding['status']})", '', finding['message'], '']
        if finding['recommendation']:
            lines += [finding['recommendation'], '']
        for key in ('metrics', 'scope', 'evidence'):
            if finding[key]:
                lines += [f'**{key.capitalize()}**', '', '```json',
                          json.dumps(finding[key], indent=2, sort_keys=True, allow_nan=False), '```', '']
    lines += ['## Runtime diagnostic inventory', '',
              '| Variable | Native dimensions | Units | Finite selected samples | Range |', '|---|---|---|---:|---|']
    for name, item in document['channel_inventory'].items():
        limits = f"{item.get('minimum', 'unavailable')} to {item.get('maximum', 'unavailable')}"
        lines.append(f"| `{name}` | {_escape(', '.join(item['dimensions']))} | {_escape(item['units'])} | {item['finite_samples']} | {limits} |")
    lines += ['', '## Limitations', '',
              'Only documented rules implemented by this evaluator were screened. Raw diagnostic reconstruction,',
              'uncertainty propagation through TRANSP, species-inventory completeness, campaign-reference approval,',
              'event segmentation, and actual solver convergence are not certified. Input provenance is declared',
              'by the sidecar; no external reference has been independently retrieved or authenticated.', '',
              'The JSON report preserves the effective thresholds, input timeline, and evidence declarations.',
              'The NPZ file contains separate evaluated and flagged masks; an unevaluated point is not a pass.', '']
    return '\n'.join(lines)


def _umask():
    mask = os.umask(0)
    os.umask(mask)
    return mask


def _replace_with_default_mode(temporary, path):
    """``mkstemp`` creates 0600; give the report the permissions a normally created file would have."""
    os.chmod(temporary, 0o666 & ~_umask())
    os.replace(temporary, path)


def _atomic_text(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as stream:
            stream.write(value)
        _replace_with_default_mode(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def write_report(report, directory, *, protected_paths=()):
    """Write three sidecars; refuse any destination resolving to an input file.

    All three documents are rendered before any file is touched, so a rendering failure leaves the
    previous report as it was.  Individual files are then replaced atomically with the default file
    permissions; the bundle is not a transactional database update.  Original CDF, TR.DAT and evidence
    files are never opened for write.
    """
    directory = Path(directory)
    targets = {ext: directory / ('curation.' + ext) for ext in ('json', 'md', 'npz')}
    protected = {Path(p).resolve() for p in protected_paths if p is not None}
    if any(p.resolve() in protected for p in targets.values()):
        raise ValueError('A report destination would overwrite a protected input file')
    directory.mkdir(parents=True, exist_ok=True)
    document = report.as_dict()
    document['mask_storage']['path'] = targets['npz'].name
    # Render everything first: JSON validity and the Markdown text can fail, and must not leave a mixed bundle.
    serialized = json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + '\n'
    markdown = markdown_report(document)
    fd, temporary = tempfile.mkstemp(prefix='.curation-masks-', dir=directory)
    try:
        with os.fdopen(fd, 'wb') as stream:
            np.savez_compressed(stream, **report.masks, **report.coordinates)
        _replace_with_default_mode(temporary, targets['npz'])
    finally:
        Path(temporary).unlink(missing_ok=True)
    _atomic_text(targets['json'], serialized)
    _atomic_text(targets['md'], markdown)
    return targets, document
