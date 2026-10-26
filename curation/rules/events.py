#!/usr/bin/python
"""Inventory of declared events and the policy exclusions they produced."""

import numpy as np

from curation import detectors, events as event_rules



def check_events(ctx, rec):
    spans = {}
    for dim in ('TIME', 'TIME3'):
        if ctx.source.has(dim):
            times = ctx.source.array(dim)
            if times.ndim == 1 and times.size and np.all(np.isfinite(times)):
                spans[dim] = (float(times[0]), float(times[-1]))
    emitted = {f.id for f in rec.report.findings}
    rows = event_rules.inventory(ctx.events, spans, emitted)
    rec.report.events = rows
    if not rows:
        rec.unavailable('EVENTS.INVENTORY', 'No events were declared; transients are not detected automatically.')
        return
    problems = [r['id'] for r in rows if r['outside_source_time'] or r['unused']]
    rec.report.add('EVENTS.INVENTORY', 'warning' if problems else 'pass',
                    'Declared events were applied as analyst exclusions; a declaration outside the source time or matching '
                    'no rule is reported rather than ignored.' if problems else
                    'Declared events were applied as analyst exclusions; they were not detected or verified here.',
                    evidence={'events': rows, 'needs_review': problems,
                              'policy_exclusions': rec.report.policy_exclusions})


def check_event_detectors(ctx, rec):
    """Optional candidate detection on explicitly configured, unit-checked diagnostics.

    Candidates are recorded as ``automatic_candidate`` events that act as nothing until an analyst accepts them.  A
    declared event excludes samples from the detector only if it names ``EVENTS.DETECT.<id>`` in its rules, so the
    exclusion is a recorded policy.  No candidates is never evidence of an event-free plasma.
    """
    for config in ctx.manifest['event_detectors']:
        identifier = 'EVENTS.DETECT.' + config['id']
        name, clock = config['diagnostic'], config['clock']
        settings = dict(config.get('parameters', {}))
        common = dict(diagnostic=name, mode=config['mode'], expected_units=set(config['expected_units']),
                      max_sample_interval=config['max_sample_interval'], **settings)
        try:
            if not (ctx.source.has(name) and ctx.source.has(clock)):
                result = detectors.abrupt_drop(None, None, units=config['expected_units'][0], **common)
                excluded = None
            else:
                meta = ctx.source.meta(name)
                values, times = ctx.source.array(name), ctx.source.array(clock)
                dims = tuple(meta['dimensions'])
                if dims == (clock, 'X') and 'radial_index' in config and config['radial_index'] < values.shape[1]:
                    values = values[:, config['radial_index']]
                elif dims != (clock,):
                    raise ValueError(f'{name}: dimensions {list(dims)} are not the declared clock {clock} (a profile needs '
                                     'radial_index); no other axis is borrowed')
                excluded, _each, _basis = event_rules.membership(ctx.events, identifier, times,
                                                                 ctx.manifest['sample_support'].get(clock))
                result = detectors.abrupt_drop(times, values, units=meta['units'], excluded=excluded, **common)
        except ValueError as exc:
            rec.report.add(identifier, 'error', f'Detector configuration or data problem: {exc}', variables=[name])
            continue
        document = result.as_dict()
        for candidate in document['candidates']:
            rec.report.detected_events.append({**candidate, 'detector': config['id'], 'diagnostic': name, 'clock': clock})
        if result.status != 'ok':
            status, message = 'not_assessed', f'Detection was not possible ({result.status}); this is not evidence of an event-free plasma.'
        elif document['candidates']:
            status, message = 'warning', (f"{len(document['candidates'])} unclassified transient candidate(s); scores are "
                                          'uncalibrated and the candidates exclude nothing until an analyst accepts them.')
        else:
            status, message = 'pass', 'No candidates above the threshold; this is not evidence of an event-free plasma.'
        rec.report.add(identifier, status, message, variables=[name],
                       evidence={'detection': document, 'excluded_by_declared_events': None if excluded is None else int(excluded.sum())})
