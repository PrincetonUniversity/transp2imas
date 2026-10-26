#!/usr/bin/python
"""Half-open time intervals: one helper for every scope comparison in curation.

Declared TR.DAT states and resolved-control segments are half-open ``[start, end)``; ``None`` is an unbounded side
(unknown, never a finite guess).  An *assessed support* is the closed span ``[lo, hi]`` of the quantity's native
clock intersected with the assessment window; a side that is not established is ``None``.

``partition`` cuts that support at every boundary of every supplied interval list, so that each piece lies inside
exactly one interval of each list (or inside none: a gap).  A boundary exactly at ``hi`` leaves the single-instant
piece ``[hi, hi]``, because the sample at ``hi`` belongs to the state that starts there.
"""

import math


def bound(value, side):
    """``None`` -> the infinite bound of that side, so that comparisons need no special cases."""
    if value is None:
        return -math.inf if side == 'start' else math.inf
    return value


def support(window, spans):
    """``(lo, hi)`` of the assessed support: window limits intersected with the native spans (``None`` = unknown).

    ``spans`` is a list of ``(first, last)`` native clock spans; with several (TIME and TIME3) the support is their
    union span, with none the native side is unknown and only the window limits apply.
    """
    lo = min(s[0] for s in spans) if spans else None
    hi = max(s[1] for s in spans) if spans else None
    if window.get('time_min') is not None:
        lo = window['time_min'] if lo is None else max(lo, window['time_min'])
    if window.get('time_max') is not None:
        hi = window['time_max'] if hi is None else min(hi, window['time_max'])
    return lo, hi


def partition(lists, lo=None, hi=None):
    """Atomic ``(start, end)`` pieces of ``[lo, hi]`` cut at every finite boundary in ``lists``.

    ``lists`` is an iterable of iterables of ``(start, end)`` intervals.  A piece's ``None`` side is only possible
    when the corresponding support side is unknown.  Boundaries outside ``(lo, hi]`` do not cut the support.
    """
    cuts = set()
    for intervals in lists:
        for start, end in intervals:
            cuts.update(edge for edge in (start, end) if edge is not None)
    low, high = bound(lo, 'start'), bound(hi, 'end')
    inner = sorted(c for c in cuts if low < c < high)
    edges = [lo] + inner + [hi]
    pieces = [(edges[k], edges[k + 1]) for k in range(len(edges) - 1)]
    if hi is not None and hi in cuts and (lo is None or hi > lo):
        pieces.append((hi, hi))
    if lo is not None and hi is not None and lo > hi:
        return []
    if lo is not None and hi is not None and lo == hi:
        return [(lo, hi)]
    return pieces


def covering(intervals, piece):
    """Index of the interval of ``intervals`` that contains the whole atomic ``piece``, else ``None``."""
    start = bound(piece[0], 'start')
    for index, (lo, hi) in enumerate(intervals):
        if bound(lo, 'start') <= start < bound(hi, 'end'):
            return index
    return None


def overlap(a, b):
    """Whether the half-open intervals ``a`` and ``b`` share a point (touching at an endpoint is not overlap)."""
    return bound(a[0], 'start') < bound(b[1], 'end') and bound(b[0], 'start') < bound(a[1], 'end')


def subtract(whole, parts):
    """The pieces of the closed interval ``whole = (lo, hi)`` not covered by the closed intervals ``parts``.

    Used for "what part of the assessed support has no evidence"; the result is a sorted list of ``[start, end]``.
    """
    lo, hi = whole
    pieces = [[lo, hi]]
    for plo, phi in sorted(parts):
        nxt = []
        for a, b in pieces:
            if phi <= a or plo >= b:
                nxt.append([a, b])
                continue
            if plo > a:
                nxt.append([a, plo])
            if phi < b:
                nxt.append([phi, b])
        pieces = nxt
    return [p for p in pieces if p[1] > p[0]]
