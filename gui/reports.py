#!/usr/bin/python
"""View models of a curation report and of a conversion report.

Pure functions over decoded JSON: nothing here touches the filesystem or the network.  A view model keeps every
status and message the report carries; it only reshapes (and bounds the size of) what the page has to draw.  A report
written by an older curation version is shown with what it has and a notice, never silently upgraded.
"""

CERTIFICATION_NOTICE = ('A curation report is a screening record. It is not a scientific certification, and a passed '
                        'rule means only that its stated predicate passed.')
STATUS_ORDER = ('fail', 'error', 'warning', 'not_assessed', 'not_applicable', 'pass')
MAX_ASSIGNMENTS_SHOWN = 200
CURATION_PREFIX = 'transp2imas.curation.report.'


def is_curation_report(document):
    return isinstance(document, dict) and str(document.get('schema_version', '')).startswith(CURATION_PREFIX)


def is_conversion_report(document):
    return (isinstance(document, dict) and 'writers' in document and 'format_version' in document
            and 'selected_ids' in document)


def graph_layout(nodes, column=240, row=84):
    """Layered layout of the evidence graph: parents left of children, rows ordered by id.

    Returns ``{'nodes': [...], 'edges': [{'from', 'to'}], 'width', 'height'}``; edges point from a node to what was
    derived from it, as the report's ``children`` say.  A cycle cannot occur in a validated report; if one is
    present the depth search stops rather than looping.
    """
    depth = {}

    def level(identifier, trail=()):
        if identifier in depth:
            return depth[identifier]
        if identifier in trail:
            return 0
        parents = [p for p in nodes[identifier].get('parents', []) if p in nodes]
        depth[identifier] = 1 + max((level(p, trail + (identifier,)) for p in parents), default=-1)
        return depth[identifier]

    for identifier in nodes:
        level(identifier)
    layers = {}
    for identifier in sorted(nodes):
        layers.setdefault(depth[identifier], []).append(identifier)
    placed = []
    for d, identifiers in sorted(layers.items()):
        for index, identifier in enumerate(identifiers):
            node = nodes[identifier]
            placed.append({'id': identifier, 'x': d * column, 'y': index * row, 'kind': node.get('kind'),
                           'role': node.get('role'), 'usage': node.get('usage', []), 'source': node.get('source', ''),
                           'sha256': node.get('sha256'), 'selector': node.get('selector'), 'scope': node.get('scope'),
                           'parents': node.get('parents', []), 'description': node.get('description')})
    edges = [{'from': p, 'to': identifier} for identifier, node in sorted(nodes.items())
             for p in node.get('parents', []) if p in nodes]
    return {'nodes': placed, 'edges': edges, 'width': (max(layers, default=0) + 1) * column,
            'height': max((len(v) for v in layers.values()), default=1) * row}


def _counts(finding):
    scope = finding.get('scope') or {}
    return {key: scope.get(key) for key in ('evaluated_count', 'flagged_count', 'shape') if key in scope}


def curation_view(document, mask_prefixes=()):
    """Everything the curation tabs draw, from one decoded ``curation.json``."""
    schema = document.get('schema_version', '')
    summary = document.get('summary', {})
    findings = document.get('findings', [])
    exclusions = document.get('policy_exclusions', {})
    notices = [CERTIFICATION_NOTICE]
    if schema != CURATION_PREFIX + 'v3':
        notices.append(f'This report is {schema or "of unknown format"}; the rule ledger, policy and fingerprints '
                       'appear only in report.v3, so those tabs may be empty.')
    prefixes = set(mask_prefixes)
    view_findings = []
    for finding in findings:
        prefix = (finding.get('scope') or {}).get('mask_prefix')
        view_findings.append({
            **{k: finding.get(k) for k in ('id', 'status', 'message', 'variables', 'applications', 'recommendation',
                                           'rule', 'rule_version', 'metrics', 'evidence', 'scope')},
            'counts': _counts(finding), 'has_masks': bool(prefix and prefix in prefixes),
            'mask_prefix': prefix, 'policy_exclusion': exclusions.get(finding['id'])})
    timeline = document.get('declared_input_timeline') or {}
    assignments = timeline.get('assignments', [])
    graph = document.get('evidence_graph') or {}
    manifest = document.get('effective_manifest') or {}
    return {
        'kind': 'curation', 'notices': notices,
        'header': {'curation_version': document.get('curation_version'), 'schema_version': schema,
                   'generated_at_utc': document.get('generated_at_utc'), 'source': document.get('source', {}),
                   'policy_name': (document.get('policy') or {}).get('name'),
                   'fingerprints': document.get('fingerprints'), 'manifest_sha256': document.get('manifest_sha256')},
        'summary': {'counts': summary.get('counts', {}), 'applications': summary.get('applications', {}),
                    'rule_coverage': summary.get('rule_coverage'), 'excused_findings': summary.get('excused_findings'),
                    'scientifically_certified': summary.get('scientifically_certified', False)},
        'findings': view_findings,
        'ledger': document.get('execution_ledger', []),
        'policy': document.get('policy'), 'registry': document.get('registry'),
        'run_identity': document.get('run_identity'),
        'migration': manifest.get('migration'),
        'events': {'declared': document.get('events', []), 'detected': document.get('detected_events', []),
                   'policy_exclusions': exclusions, 'sample_support': manifest.get('sample_support', {})},
        'evidence': {'graph': graph_layout(graph) if graph else None,
                     'coverage': [{'id': f['id'], 'status': f['status'], **(f.get('evidence') or {}).get('coverage', {})}
                                  for f in findings if (f.get('evidence') or {}).get('coverage')]},
        'uncertainty': document.get('uncertainty', {}),
        'timeline': {'available': timeline.get('available'), 'parser': timeline.get('parser'),
                     'files': timeline.get('files'), 'gaps': timeline.get('gaps', []),
                     'unresolved': timeline.get('unresolved', []), 'states': len(timeline.get('states', [])),
                     'temporal_valid': timeline.get('temporal_valid'), 'source': timeline.get('source'),
                     'assignments_total': len(assignments), 'assignments_truncated': timeline.get('assignments_truncated'),
                     'assignments': assignments[:MAX_ASSIGNMENTS_SHOWN]},
        'code_identity': (document.get('source') or {}).get('code_identity'),
        'mask_prefixes': sorted(prefixes),
        'status_order': list(STATUS_ORDER)}


def conversion_view(document):
    """The conversion report as the Convert tab draws it: header, publication, per-IDS table, mapping journal."""
    publication = document.get('publication') or {}
    rows, journal = [], []
    for writer in document.get('writers', []):
        mapping = writer.get('mapping') or {}
        reads = writer.get('reads') or {}
        entries = 0
        for kind, record in sorted(mapping.items()):           # kinds: skipped, missing_optional, adjusted
            for field in (record or {}).get('fields', []):
                entries += field.get('count', 1)
                journal.append({'ids': writer.get('ids'), 'kind': kind, 'target': field.get('target'),
                                'reason': field.get('reason'), 'count': field.get('count', 1),
                                'truncated': bool(record.get('truncated'))})
        rows.append({'ids': writer.get('ids'), 'occurrences': writer.get('occurrences'), 'status': writer.get('status'),
                     'seconds': writer.get('seconds'), 'write_count': writer.get('write_count'),
                     'write_seconds': writer.get('write_seconds'), 'consumed': len(writer.get('consumed', [])),
                     'variables_read': reads.get('variables'), 'mapping_entries': entries,
                     'error': writer.get('error')})
    accepted = publication.get('accepted', [])
    return {
        'kind': 'conversion',
        'notices': ['A passing conversion report records what was requested and accepted by the backend; it is not a '
                    'scientific validation of the converted data.'],
        'header': {k: document.get(k) for k in ('status', 'machine', 'shot', 'run', 'occurrence', 'backend', 'dd_version',
                                                'mapping_policy', 'source_cdf', 'output_directory', 'error',
                                                'format_version')},
        'selected_ids': document.get('selected_ids', []), 'written_ids': document.get('written_ids', []),
        'publication': {'status': publication.get('status'), 'accepted': accepted, 'kind': publication.get('kind'),
                        'failed': publication.get('failed'), 'unattempted': publication.get('unattempted', [])},
        'writers': rows, 'mapping_journal': journal,
        'setup_reads': document.get('setup_reads', [])}
