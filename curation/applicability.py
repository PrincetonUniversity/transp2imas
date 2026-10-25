#!/usr/bin/python
"""Whether a declared reference or configuration applies to this run.

The answer is one of ``applies``, ``not_applicable`` or ``unknown``.  ``unknown`` means a criterion was declared
but the run does not establish the fact it needs (no shot number, no device); it is never read as "applies".
"""

import fnmatch

APPLIES, NOT_APPLICABLE, UNKNOWN = 'applies', 'not_applicable', 'unknown'


def judge(spec, identity):
    """``(state, reasons)`` for an applicability ``spec`` against ``identity`` (device, shot, campaign).

    Time bounds are not judged here: they restrict which declared intervals a reference is compared against
    (see :func:`overlaps`).
    """
    if not spec:
        return APPLIES, []
    unknown, mismatch = [], []
    for key in ('device', 'campaign'):
        if key in spec:
            have = identity.get(key)
            if have is None:
                unknown.append(f'{key} of the run is not established')
            elif str(have).strip().lower() != str(spec[key]).strip().lower():
                mismatch.append(f'{key} {have!r} is not {spec[key]!r}')
    if 'shot_min' in spec or 'shot_max' in spec:
        shot = identity.get('shot')
        if shot is None:
            unknown.append('shot of the run is not established')
        else:
            if 'shot_min' in spec and shot < spec['shot_min']:
                mismatch.append(f'shot {shot} is before {spec["shot_min"]}')
            if 'shot_max' in spec and shot > spec['shot_max']:
                mismatch.append(f'shot {shot} is after {spec["shot_max"]}')
    if mismatch:
        return NOT_APPLICABLE, mismatch
    if unknown:
        return UNKNOWN, unknown
    return APPLIES, []


def overlaps(spec, start, end):
    """Whether the interval ``[start, end)`` (``None`` = unbounded) meets the time bounds of ``spec``."""
    lo, hi = (spec or {}).get('time_min'), (spec or {}).get('time_max')
    if lo is not None and end is not None and end <= lo:
        return False
    if hi is not None and start is not None and start > hi:
        return False
    return True


def rule_matches(patterns, identifier):
    """Whether any shell-style pattern in ``patterns`` matches the rule/finding id (case-sensitive)."""
    return any(fnmatch.fnmatchcase(identifier, pattern) for pattern in patterns)
