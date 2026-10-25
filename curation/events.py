#!/usr/bin/python
"""Analyst-declared events (sawteeth, ELMs, pellets, heating transitions, outages) and the samples they exclude.

Events are never detected here.  Each declared event names the rule ids it applies to; samples inside its closed
time window are *excluded* from those rules' evaluations and are accounted for as excluded, which is different
from being unflagged (they were not assessed) and different from being outside the assessment window.
"""

import numpy as np

from curation.applicability import rule_matches


def exclusion(events, identifier, times):
    """``(mask, event_ids)``: samples of the 1-D ``times`` inside any event that applies to rule ``identifier``."""
    mask = np.zeros(np.shape(times), dtype=bool)
    used = []
    for event in events:
        if not rule_matches(event['rules'], identifier):
            continue
        inside = (times >= event['start']) & (times <= event['end'])
        if inside.any():
            used.append(event['id'])
        mask |= inside
    return mask, used


def inventory(events, spans, emitted):
    """Per-event facts for the report.

    ``spans`` maps a time dimension to ``(first, last)`` of its native coordinate; ``emitted`` is the set of finding
    ids the run produced.  An event is *outside the source* when its window misses every native span, and
    *unused* when none of its rules matched any emitted finding.
    """
    rows = []
    for event in events:
        outside = all(event['end'] < lo or event['start'] > hi for lo, hi in spans.values()) if spans else True
        matched = sorted(i for i in emitted if rule_matches(event['rules'], i))
        rows.append({'id': event['id'], 'type': event['type'], 'start': event['start'], 'end': event['end'],
                     'rules': event['rules'], 'outside_source_time': outside, 'matched_findings': matched,
                     'unused': not matched, 'source': event.get('source')})
    return rows
