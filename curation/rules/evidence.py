#!/usr/bin/python
"""Run identity, evidence provenance, coverage, uncertainty declarations and lineage."""

from curation import applicability, evidence as evidence_graph, intervals
from curation.results import ROLE_APPS

# the native clock each evidence role's quantity lives on (profiles on TIME3, scalar histories on TIME)
ROLE_CLOCK = {'electron_temperature': 'TIME3', 'ion_temperature': 'TIME3', 'electron_density': 'TIME3',
              'rotation': 'TIME3', 'composition': 'TIME3', 'heating': 'TIME'}

def check_identity(ctx, rec):
    declared = ctx.manifest['run_identity']
    found = {'device': ctx.source.global_attribute('device'), 'shot': ctx.source.global_attribute('shot')}
    rows, status = [], 'pass'
    for key in ('device', 'shot'):
        if key not in declared:
            continue
        have = found[key]
        if have is None:
            rows.append({'item': key, 'declared': declared[key], 'source': None, 'match': None})
            status = status if status == 'fail' else 'not_assessed'
            continue
        match = (str(have).strip().lower() == str(declared[key]).strip().lower() if key == 'device'
                 else int(have) == int(declared[key]))
        rows.append({'item': key, 'declared': declared[key], 'source': have, 'match': match})
        if not match:
            status = 'fail'
    rec.report.run_identity = {'declared': declared, 'source_attributes': found, 'effective': ctx.identity,
                                'conflicts': list(ctx.identity_conflicts),
                                'run_id': ctx.source.global_attribute('Runid'),
                                'transp_version': ctx.source.global_attribute('TRANSP_version')}
    if not rows:
        rec.unavailable('RUN.IDENTITY', 'No run identity was declared; the source is not tied to a named run.')
        return
    message = {'pass': 'The declared run identity matches the source file.',
               'fail': 'The declared run identity disagrees with the source file; the supplied evidence may belong to another run.',
               'not_assessed': 'The source file does not state part of the declared run identity; it is unverified.'}[status]
    rec.report.add('RUN.IDENTITY', status, message, evidence={'checks': rows})


def _coverage(ctx, role, record):
    """Where the record's time applicability meets the assessed native support of the role's quantity.

    ``state`` is ``full``, ``partial``, ``disjoint`` or ``unknown`` (the record has time bounds but the support
    is empty or not established).  A record with no time bounds covers whatever support there is.  TIME and
    TIME3 are separate clocks: the role names the one its quantity lives on.
    """
    spec = record.get('applicability') or {}
    clock = ROLE_CLOCK[role]
    lo, hi = intervals.support(ctx.window, ctx.native_spans((clock,)))
    t0, t1 = spec.get('time_min'), spec.get('time_max')
    result = {'clock': clock, 'support': [lo, hi], 'record_interval': [t0, t1], 'covered_interval': None,
              'unassessed_intervals': []}
    if t0 is None and t1 is None:
        result.update(state='full', covered_interval=[lo, hi])
    elif lo is None or hi is None or lo > hi:
        result['state'] = 'unknown'
    elif (t0 is not None and t0 > hi) or (t1 is not None and t1 < lo):
        result['state'] = 'disjoint'
    else:
        covered = [lo if t0 is None else max(lo, t0), hi if t1 is None else min(hi, t1)]
        missing = intervals.subtract((lo, hi), [tuple(covered)])
        result.update(state='partial' if missing else 'full', covered_interval=covered, unassessed_intervals=missing)
    return result


def _scope_of(role, record, coverage):
    """The node's effective time scope for lineage: the covered part, ``False`` if none, ``None`` if unscoped."""
    if coverage['state'] == 'disjoint':
        return False
    if coverage['state'] == 'unknown':
        spec = record.get('applicability') or {}
        return (spec.get('time_min'), spec.get('time_max'))
    return None if coverage['covered_interval'] is None else tuple(coverage['covered_interval'])


def check_evidence(ctx, rec):
    not_applicable = {}
    for role, record in ctx.manifest['evidence'].items():
        state, reasons = applicability.judge(record.get('applicability'), ctx.identity)
        if state != applicability.APPLIES:
            not_applicable[role] = (state, reasons)
    coverage = {role: _coverage(ctx, role, record) for role, record in ctx.manifest['evidence'].items()
                if role not in not_applicable}
    for role, apps in ROLE_APPS.items():
        record = ctx.manifest['evidence'].get(role, {})
        origin = record.get('origin', 'unknown')
        source = record.get('source', '')
        if role in not_applicable:
            state, reasons = not_applicable[role]
            status = 'not_applicable' if state == applicability.NOT_APPLICABLE else 'not_assessed'
            rec.report.add('EVIDENCE.' + role, status,
                            'The declared evidence does not apply to this run: ' + '; '.join(reasons)
                            if state == applicability.NOT_APPLICABLE else
                            'Whether the declared evidence applies cannot be established: ' + '; '.join(reasons),
                            applications=list(apps), evidence=record)
            rec.report.add('UNCERTAINTY.' + role, 'not_assessed', 'No applicable evidence, hence no uncertainty.',
                            applications=list(apps))
            continue
        if origin in ('unknown', 'missing') or not source:
            status, message = 'not_assessed', 'Origin/reference is missing or unknown; automatic retrieval is not validation.'
        elif record.get('reliability') == 'unreliable':
            status, message = 'warning', 'Evidence is explicitly marked unreliable.'
        elif origin in ('assumed', 'inferred', 'predicted'):
            status, message = 'warning', 'This quantity is assumption/model dependent, not an independent measurement.'
        else:
            status, message = 'pass', 'An origin and reference were declared; their accuracy was not independently established.'
        uncertainty = record.get('uncertainty', 'unknown')
        u_status = 'pass' if uncertainty == 'provided' else 'not_assessed'
        u_message = ('Uncertainty is declared as provided; propagation/calibration remains unverified.'
                     if uncertainty == 'provided' else 'No usable uncertainty specification; unknown is not zero.')
        evidence, u_evidence = dict(record), {}
        if role in coverage:
            cov = coverage[role]
            evidence['coverage'], u_evidence['coverage'] = cov, cov
            if cov['state'] in ('disjoint', 'unknown'):
                reason = ('The declared evidence does not overlap the assessed interval' if cov['state'] == 'disjoint'
                          else 'The declared evidence interval cannot be placed on the assessed support')
                status, message = 'not_assessed', f'{reason}; the quantity has no evidence there (the record is kept).'
                u_status, u_message = 'not_assessed', message
            elif cov['state'] == 'partial':
                note = f' It covers only part of the assessed interval; see unassessed_intervals.'
                if status == 'pass':
                    status = 'warning'
                message += note
                if u_status == 'pass':
                    u_status = 'warning'
                u_message += note
        rec.report.add('EVIDENCE.' + role, status, message, applications=list(apps), evidence=evidence)
        rec.report.add('UNCERTAINTY.' + role, u_status, u_message, applications=list(apps), evidence=u_evidence)
    # Independence is judged on declared lineage (parents, resources, selectors); it is never a proof.
    evidence = ctx.manifest['evidence']
    excluded = [evidence[role].get('id', 'evidence:' + role) for role in not_applicable]
    scopes = {evidence[role].get('id', 'evidence:' + role): _scope_of(role, evidence[role], cov)
              for role, cov in coverage.items()}
    nodes = evidence_graph.build_graph(ctx.manifest, excluded, scopes)
    rec.report.evidence_graph = evidence_graph.summary(nodes)
    for role, record in evidence.items():
        identifier = record.get('id', 'evidence:' + role)
        if 'independent_validation' not in record['usage'] or identifier not in nodes or role not in ROLE_APPS:
            continue
        shared = evidence_graph.shared_lineage(nodes, identifier)
        if not shared:
            continue
        dependent = [item for item in shared if item['kind'] in evidence_graph.DEPENDENT]
        apps = list(ROLE_APPS[role])
        details = {'shared_with': sorted({item['with'] for item in shared}), 'shared': shared,
                   'source': record.get('source')}
        if dependent:
            rec.report.add('INDEPENDENCE.' + role, 'warning',
                            'The claimed independent validator shares, or may share, lineage with evidence that constrained the run.',
                            applications=apps, evidence=details)
        else:
            rec.report.add('INDEPENDENCE.' + role, 'pass',
                            'The validator and a constraining record share a container but declare different selectors; '
                            'this is not proof of independence (undeclared dependence is invisible).',
                            applications=apps, evidence=details)


def check_rotation(ctx, rec):
    names = [name for name in ('OMEGA', 'VTOR') if ctx.source.has(name)]
    apps = ('momentum_transport', 'heat_transport', 'surrogate_training')
    record = ctx.manifest['evidence'].get('rotation', {})
    if not names:
        rec.report.add('ROTATION.AVAILABILITY', 'warning',
                        'No recognized rotation output; do not replace missing rotation by zero.',
                        applications=list(apps), recommendation='Provide measurement provenance or explicit rotation scenarios.')
    else:
        rec.report.add('ROTATION.AVAILABILITY', 'pass',
                        'Rotation output exists; its presence does not establish measurement provenance.',
                        variables=names, applications=list(apps), evidence=record)
