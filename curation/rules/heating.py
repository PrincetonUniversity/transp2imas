#!/usr/bin/python
"""Declared heating configuration compared with analyst-supplied campaign references."""

from curation import applicability
from curation.inputs import normalize_name



def check_heating_reference(ctx, rec):
    temporal = ctx.timeline.get('temporal_valid', True)
    states = list(ctx.states_in_window()) if temporal else []
    _report_heating_reference(ctx, rec, states)


def _report_heating_reference(ctx, rec, states):
    expectations = ctx.manifest['input_expectations']
    if not expectations:
        rec.unavailable('HEATING.REFERENCE', 'No campaign/configuration reference expectations were supplied.')
    for index, item in enumerate(expectations):
        key = normalize_name(item['name'])
        identifier = f'HEATING.REFERENCE.{index}'
        records = []
        missing = 0
        state_of_applicability, reasons = applicability.judge(item.get('applicability'), ctx.identity)
        if state_of_applicability != applicability.APPLIES:
            not_applicable = state_of_applicability == applicability.NOT_APPLICABLE
            rec.report.add(identifier, 'not_applicable' if not_applicable else 'not_assessed',
                            ('This reference does not apply to this run: ' if not_applicable else
                             'Whether this reference applies cannot be established: ') + '; '.join(reasons),
                            variables=[key], evidence={'expectation': item})
            continue
        override = item.get('override')
        if override is not None and not isinstance(override['value'], (int, float)):
            rec.report.add(identifier, 'error', 'The override value is not numeric and cannot replace a numeric reference.',
                            variables=[key], evidence={'expectation': item})
            continue
        expected = override['value'] if override else item['expected']
        time_spec = item.get('applicability')
        for state in states:
            if not applicability.overlaps(time_spec, state['start_time'], state['end_time']):
                continue
            value = state['values'].get(key)
            if not isinstance(value, (int, float)) or isinstance(value, bool):
                missing += 1
                continue
            matches = abs(value - expected) <= item['atol'] + item['rtol'] * abs(expected)
            records.append({'start_time': state['start_time'], 'end_time': state['end_time'],
                            'value': value, 'matches_reference': matches,
                            'source_line': state['source_lines'].get(key)})
        status = ('warning' if any(not r['matches_reference'] for r in records)
                  else 'not_assessed' if not records or missing else 'pass')
        note = (' A reviewed override replaced the reference value (' + override['reason'] + '; reviewer '
                + override['reviewer'] + '); the original is kept in the evidence.') if override else ''
        rec.report.add(identifier, status,
                        'Explicit native-unit declarations compared to the supplied reference; deviations are not automatically corrected.' + note,
                        variables=[key], evidence={'expectation': item, 'effective_expected': expected,
                                                   'override_applied': override is not None,
                                                   'intervals': records, 'unknown_intervals': missing})
