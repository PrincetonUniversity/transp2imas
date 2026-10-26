#!/usr/bin/python
"""Declared-input timeline, model-control review flags and resolved-input comparison."""

import re

from curation import intervals, resolved
from curation.inputs import normalize_name
from curation.results import ROLE_APPS



def check_declared_inputs(ctx, rec):
    if not ctx.timeline.get('available'):
        rec.unavailable('INPUT.TIMELINE', 'No TR.DAT supplied; controls and default values remain unknown.')
    else:
        timeline = ctx.timeline
        gaps = timeline.get('gaps', [])
        parser = timeline.get('parser') or {'kind': 'unspecified'}
        unresolved = timeline.get('unresolved', [])
        limited = parser.get('kind') != 'f90nml'
        concerns = bool(gaps or unresolved or limited)
        message = ('Explicit declarations were parsed (' + parser['kind'] + '); this is the declared state in file order, '
                   'not the resolved TRANSP configuration: defaults, omitted bounds, preprocessing and model-execution '
                   'semantics are not applied.')
        if limited:
            message += ' The limited scalar reader was used, so arrays, repeats and file references are unresolved.'
        if unresolved:
            message += f' {len(unresolved)} construct(s) remain unresolved (see evidence.unresolved).'
        rec.report.add('INPUT.TIMELINE', 'warning' if concerns else 'pass', message,
                       evidence={'parser': parser, 'parser_gaps': gaps, 'unresolved': unresolved,
                                 'files': timeline.get('files'), 'source': timeline.get('source'),
                                 'assignments_recorded': len(timeline.get('assignments', [])),
                                 'assignments_truncated': timeline.get('assignments_truncated', False)})
    temporal = ctx.timeline.get('temporal_valid', True)
    states = list(ctx.states_in_window()) if temporal else []
    coupled = []
    conflicts = []
    fractions = []
    for state in states:
        vals = state['values']
        density = [k for k in ('LPREDICT_NE', 'LPREDICT_NMAIN', 'LPREDICT_NIMP')
                   if isinstance(vals.get(k), (int, float)) and not isinstance(vals[k], bool) and vals[k] > 0]
        diffusion = [k for k in ('NLMDIF', 'NLISODIF') if vals.get(k) is True or (isinstance(vals.get(k), (int, float)) and vals[k] == 1)]
        if density and diffusion:
            coupled.append({'start_time': state['start_time'], 'end_time': state['end_time'],
                            'density_controls': density, 'diffusion_controls': diffusion})
        if vals.get('NLVWNC') in (True, 1) and vals.get('NLIVPO') in (True, 1):
            conflicts.append({'start_time': state['start_time'], 'end_time': state['end_time']})
        for name, value in vals.items():
            if re.fullmatch(r'FFULLA(?:\(\d+\))?', name):
                other = name.replace('FFULLA', 'FHALFA')
                half = vals.get(other)
                if (isinstance(value, (int, float)) and not isinstance(value, bool)
                        and isinstance(half, (int, float)) and not isinstance(half, bool)):
                    fractions.append({'full_key': name, 'half_key': other, 'full': value, 'half': half,
                                      'start_time': state['start_time'], 'end_time': state['end_time'],
                                      'invalid': not (0 <= value <= 1 and 0 <= half <= 1 and value + half <= 1 + 1e-9)})
    if not temporal:
        # no declared state could be placed in time: say so instead of reporting a clean "nothing found"
        reason = ('Not evaluated: the TR.DAT update times are unusable, so no declared state could be placed '
                  'in time (see INPUT.TIMELINE).')
        for identifier in ('MODEL.DENSITY_CURRENT_COUPLING', 'MODEL.ROTATION_CONTROLS', 'HEATING.FRACTIONS'):
            rec.unavailable(identifier, reason)
    else:
        _report_declared_models(ctx, rec, coupled, conflicts, fractions)


def _report_declared_models(ctx, rec, coupled, conflicts, fractions):
    rec.report.add('MODEL.DENSITY_CURRENT_COUPLING', 'warning' if coupled else 'not_assessed',
                    'Density and magnetic diffusion are both declared active; inspect composition, resistivity, boundaries and update coupling.'
                    if coupled else 'No simultaneous activation established from explicit supported controls; consistency is not certified.',
                    evidence={'intervals': coupled}, recommendation='This combination is not prohibited; document the closure.')
    rec.report.add('MODEL.ROTATION_CONTROLS', 'warning' if conflicts else 'not_assessed',
                    'Both radial-potential choices are declared; verify compatibility for the actual TRANSP version.'
                    if conflicts else 'No explicit double selection found; omitted/default controls are not resolved.',
                    evidence={'intervals': conflicts}, applications=list(ROLE_APPS['rotation']))
    rec.report.add('HEATING.FRACTIONS', 'fail' if any(x['invalid'] for x in fractions) else 'pass' if fractions else 'not_assessed',
                    'FFULLA/FHALFA are current fractions; the implied third fraction must be nonnegative. No power-fraction equivalence is assumed.',
                    evidence={'declarations': fractions})


def check_resolved_inputs(ctx, rec):
    spec = ctx.manifest['resolved_inputs']
    if not (spec.get('code') or spec['inputs'] or spec['controls']):
        rec.unavailable('INPUT.RESOLVED', 'No resolved-input evidence was supplied; controls and files are unverified.')
        return
    code = spec.get('code')
    if code:
        version = ctx.source.global_attribute('TRANSP_version')
        if version is None:
            status, message = 'not_assessed', 'The source file does not state its TRANSP version.'
        elif str(version).strip() != code['version'].strip():
            status, message = 'warning', 'The declared code version differs from the version recorded in the source file.'
        else:
            status, message = 'pass', 'The declared code version matches the version recorded in the source file.'
        rec.report.add('INPUT.RESOLVED.CODE', status, message,
                        evidence={'declared': code, 'source_version': version})
    else:
        rec.unavailable('INPUT.RESOLVED.CODE', 'The code version used for the run was not declared.')
    rows = resolved.check_files(spec['inputs'])
    if not rows:
        rec.unavailable('INPUT.RESOLVED.FILES', 'No input files were declared with hashes.')
    else:
        states = {row['status'] for row in rows}
        status = ('fail' if 'mismatch' in states else 'error' if 'unreadable' in states
                  else 'not_assessed' if 'unhashed' in states else 'pass')
        rec.report.add('INPUT.RESOLVED.FILES', status,
                        'Declared input files compared by SHA-256; a match proves only that this file is the declared one.',
                        evidence={'files': rows})
    _check_controls(ctx, rec, spec['controls'])


def _check_controls(ctx, rec, controls):
    temporal = ctx.timeline.get('temporal_valid', True)
    states = ctx.timeline.get('states', []) if temporal and ctx.timeline.get('available') else []
    lo, hi = intervals.support(ctx.window, ctx.native_spans())
    declared = [(state['start_time'], state['end_time']) for state in states]
    for control in controls:
        key = normalize_name(control['name'])
        identifier = 'INPUT.RESOLVED.CONTROL.' + key
        if not control['resolved']:
            rec.report.add(identifier, 'not_assessed',
                            'Supplied as a declaration, not as a resolved runtime value; it cannot confirm the run.',
                            variables=[key], evidence={'control': control})
            continue
        if not states:
            rec.report.add(identifier, 'not_assessed',
                            'No declared TR.DAT state exists to compare with the resolved control.',
                            variables=[key], evidence={'control': control, 'intervals': [], 'unknown_intervals': 0})
            continue
        runtime = [(seg.get('start'), seg.get('end')) for seg in control['segments']]
        records, unknown = [], 0
        for piece in intervals.partition([declared, runtime], lo, hi):
            state_index = intervals.covering(declared, piece)
            segment = intervals.covering(runtime, piece)
            declared_value = None if state_index is None else states[state_index]['values'].get(key)
            resolved_value = None if segment is None else control['segments'][segment]['value']
            known = declared_value is not None and segment is not None
            records.append({'start_time': piece[0], 'end_time': piece[1], 'declared': declared_value,
                            'resolved': resolved_value, 'segment': segment,
                            'agrees': resolved.same_value(declared_value, resolved_value) if known else None})
            unknown += not known
        if any(r['agrees'] is False for r in records):
            status, message = 'warning', 'The TR.DAT declaration differs from the resolved runtime value on part of the support; the declaration is not what ran there.'
        elif not records or unknown:
            status, message = 'not_assessed', 'Part of the assessed support has no declared value or no resolved runtime value.'
        else:
            status, message = 'pass', 'The declaration agrees with the resolved runtime value over the whole assessed support.'
        rec.report.add(identifier, status, message, variables=[key],
                        evidence={'control': control, 'intervals': records, 'unknown_intervals': unknown,
                                  'support': [lo, hi]})
