#!/usr/bin/python
"""Events (sawteeth, ELMs, pellets, heating transitions, outages, unclassified transients) and the samples they cover.

An event is a catalog entry, not a detection: it has an *origin* (experimental observation, actuator command,
modeled TRANSP event, analyst declaration, or automatic candidate), an interval, a timing uncertainty and a review
state.  A commanded pellet is not a delivered pellet and a modeled sawtooth is not independent experimental
evidence, so origin is kept and reported with every use.

Membership is a fact about the clock.  A sample is a *member* of an event when its support interval meets the
event's closed window ``[start - pad, end + pad]`` (``pad`` = timing uncertainty + declared padding); endpoints are
inclusive and simultaneous events union.  The support interval is the declared sample averaging interval of that
clock (``sample_support``), else the sample is treated as a point, which is labelled an approximation.  TIME and
TIME3 are separate clocks and are never borrowed from each other.

Membership makes samples ineligible under the event's declared (legacy exclusion) policy, but only for events that
may act as exclusions (see ``acts_as_exclusion``), and never for source-integrity rules.  It does not change the raw
finding, and it is different from missingness and from being outside the assessment window.
"""

import numpy as np

from curation.applicability import rule_matches
from curation.policy import is_mandatory

POINT = 'point_sample_approximation'
DECLARED = 'declared_averaging_interval'


def acts_as_exclusion(event):
    """Whether an event may make samples ineligible.

    An automatic candidate is only a candidate until an analyst accepts it, and a rejected event excludes nothing.
    """
    if event.get('review_state') == 'reviewed_rejected':
        return False
    if event['origin'] == 'automatic_candidate':
        return event.get('review_state') == 'reviewed_accepted'
    return True


def padding(event):
    return event.get('timing_uncertainty', 0.0) + event.get('padding', 0.0)


def sample_window(times, support):
    """``(lo, hi, basis)``: the interval each sample represents (the point itself when no support is declared)."""
    times = np.asarray(times, dtype=float)
    if not support:
        return times, times, POINT
    width = support['averaging_interval']
    if support.get('alignment', 'center') == 'trailing':
        return times - width, times, DECLARED
    return times - width / 2.0, times + width / 2.0, DECLARED


def membership(events, identifier, times, support=None):
    """``(mask, {event id: mask}, basis)`` of the samples of the 1-D ``times`` inside events acting on ``identifier``."""
    times = np.asarray(times, dtype=float)
    lo, hi, basis = sample_window(times, support)
    union = np.zeros(times.shape, dtype=bool)
    each = {}
    for event in events:
        if not acts_as_exclusion(event) or not rule_matches(event['rules'], identifier):
            continue
        pad = padding(event)
        inside = (hi >= event['start'] - pad) & (lo <= event['end'] + pad)
        each[event['id']] = inside
        union |= inside
    return union, each, basis


def exclusion(events, identifier, times):
    """Compatibility form of ``membership`` for point samples: ``(mask, ids of events that hit a sample)``."""
    mask, each, _ = membership(events, identifier, times)
    return mask, [event_id for event_id, inside in each.items() if inside.any()]


def inventory(events, spans, emitted):
    """Per-event facts for the report.

    ``spans`` maps a time dimension to ``(first, last)`` of its native coordinate; ``emitted`` is the set of finding
    ids the run produced.  An event is *outside the source* when its padded window misses every native span, and
    *unused* when none of its rules matched a finding that an exclusion may apply to (or the event cannot act as an
    exclusion).  Source-integrity findings that a rule pattern reached are listed apart: an event can never waive them.
    """
    rows = []
    for event in events:
        pad = padding(event)
        outside = (all(event['end'] + pad < lo or event['start'] - pad > hi for lo, hi in spans.values())
                   if spans else True)
        reached = sorted(i for i in emitted if rule_matches(event['rules'], i))
        acts = acts_as_exclusion(event)
        matched = [i for i in reached if not is_mandatory(i)] if acts else []
        rows.append({'id': event['id'], 'type': event['type'], 'origin': event['origin'], 'start': event['start'],
                     'end': event['end'], 'padding': pad, 'rules': event['rules'], 'policy': event['policy'],
                     'review_state': event.get('review_state'), 'acts_as_exclusion': acts,
                     'independent_experimental_evidence': event['origin'] == 'experimental_observation',
                     'reference': event.get('reference'), 'detector': event.get('detector'),
                     'applications': event.get('applications'), 'outside_source_time': outside,
                     'matched_findings': matched,
                     'mandatory_rules_not_waivable': [i for i in reached if is_mandatory(i)],
                     'unused': not matched, 'source': event.get('source')})
    return rows
