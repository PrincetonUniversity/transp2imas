#!/usr/bin/python
"""Optional event-CANDIDATE detector ``abrupt_drop``: a sudden decrease of one scalar time history.

A candidate is an interval in which a diagnostic trace fell abruptly relative to its own recent noise.  It is never
classified (type is always ``unclassified_transient``; an ELM, sawtooth or pellet is not inferred), its origin is
``automatic_candidate``, and ``score`` is an uncalibrated ratio, not a probability.  The detector cannot establish
that a plasma was event-free: no candidates with missing, excluded, gapped or too slowly sampled data says nothing
about events in those intervals (see the ``notes`` and ``coverage`` of the result).

Method (version ``METHOD_VERSION``), per segment of valid samples (finite, unmasked, not policy-excluded, and not
separated by a time gap > ``max_gap_factor`` x the median spacing; windows never cross a segment boundary and
excluded samples are neither used nor bridged): step k has increment d_k = y_k - y_{k-1}; its score is
-d_k / (1.4826 x MAD of the other increments in the baseline window).  ``causal`` uses the ``window`` increments
before k; ``noncausal`` the increments within +-window//2 of k (k excluded), truncated at the segment ends, and must
not feed forecasting.  The local level is the median of up to ``window`` samples before the step (past only, both
modes) and a relative drop needs a positive level.  A step is flagged when score >= ``threshold`` and
-d_k/level >= ``min_relative_drop``; consecutive flagged steps merge.  Timing is only known to the sample spacing.
"""

from dataclasses import dataclass, field
import math

import numpy as np

from curation.source import unit_key

METHOD = 'abrupt_drop'
METHOD_VERSION = '1.0'
EVENT_TYPE = 'unclassified_transient'
ORIGIN = 'automatic_candidate'
MODES = ('causal', 'noncausal')
STATUSES = ('ok', 'missing_diagnostic', 'no_data', 'insufficient_cadence', 'insufficient_baseline')
NOT_EVIDENCE = ('absence of candidates is not evidence of an event-free plasma: unassessed intervals, '
                'slow sampling, missing or excluded samples and detector insensitivity are not assessed as quiet')
MAD_SCALE = 1.4826
REASONS = ('missing', 'excluded', 'gap', 'insufficient_baseline')


def _num(value):
    """JSON-safe float (None for NaN/inf)."""
    value = float(value)
    return value if math.isfinite(value) else None


@dataclass
class DetectionResult:
    diagnostic: str
    units: str
    mode: str
    status: str
    settings: dict
    candidates: list = field(default_factory=list)
    coverage: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)
    method: str = METHOD
    method_version: str = METHOD_VERSION

    def as_dict(self):
        """JSON-serializable summary (no NaN/inf)."""
        return {'diagnostic': self.diagnostic, 'units': self.units, 'mode': self.mode, 'method': self.method,
                'method_version': self.method_version, 'status': self.status, 'settings': dict(self.settings),
                'candidates': [dict(c) for c in self.candidates], 'coverage': dict(self.coverage),
                'notes': list(self.notes)}


def _check_text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{name} must be an explicit non-empty string, got {value!r}')


def _validate(times, units, expected_units, diagnostic, mode, time_units, max_sample_interval, threshold,
              min_relative_drop, window, min_baseline, max_gap_factor, timing_padding):
    _check_text(diagnostic, 'diagnostic')
    if mode not in MODES:
        raise ValueError(f'mode must be one of {MODES}, got {mode!r}')
    _check_text(units, 'units')
    if not expected_units or isinstance(expected_units, str):
        raise ValueError('expected_units must be a non-empty collection of unit strings')
    if unit_key(units) not in {unit_key(u) for u in expected_units}:
        raise ValueError(f'units {units!r} of {diagnostic!r} not in expected {sorted(expected_units)}')
    if not isinstance(time_units, str) or unit_key(time_units) != 'S':
        raise ValueError(f"time units must be declared 's', got {time_units!r}")
    for name, value in (('max_sample_interval', max_sample_interval), ('threshold', threshold)):
        if not (math.isfinite(value) and value > 0):
            raise ValueError(f'{name} must be finite and positive')
    if not (math.isfinite(min_relative_drop) and min_relative_drop >= 0):
        raise ValueError('min_relative_drop must be finite and >= 0')
    if not (math.isfinite(max_gap_factor) and max_gap_factor > 1):
        raise ValueError('max_gap_factor must be finite and > 1')
    if not (math.isfinite(timing_padding) and timing_padding >= 0):
        raise ValueError('timing_padding must be finite and >= 0')
    if int(min_baseline) != min_baseline or min_baseline < 3 or int(window) != window or window < min_baseline:
        raise ValueError('need integer window >= min_baseline >= 3')


def _times(times):
    t = np.array(times, dtype=float)
    if t.ndim != 1 or not np.all(np.isfinite(t)) or (t.size > 1 and not np.all(np.diff(t) > 0)):
        raise ValueError('times must be a 1-D finite strictly increasing array')
    return t


def _merge(spans):
    out = []
    for lo, hi in spans:
        if out and lo <= out[-1][1]:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return out


def abrupt_drop(times, values, units, diagnostic, mode, expected_units, max_sample_interval, time_units='s',
                excluded=None, threshold=5.0, min_relative_drop=0.05, window=20, min_baseline=8,
                max_gap_factor=3.0, timing_padding=0.0):
    """Detect abrupt-decrease candidates; ``values=None`` means the diagnostic is absent.  Inputs are not modified.

    ``max_sample_interval`` (s) is the slowest median spacing the caller accepts for the physics of interest.
    Raises ValueError for undeclared/mismatched units, missing identifier, unknown mode or bad times/parameters.
    """
    _validate(times, units, expected_units, diagnostic, mode, time_units, max_sample_interval, threshold,
              min_relative_drop, window, min_baseline, max_gap_factor, timing_padding)
    window, min_baseline = int(window), int(min_baseline)
    settings = {'threshold': float(threshold), 'min_relative_drop': float(min_relative_drop), 'window': window,
                'min_baseline': min_baseline, 'max_gap_factor': float(max_gap_factor),
                'max_sample_interval': float(max_sample_interval), 'timing_padding': float(timing_padding),
                'time_units': 's', 'expected_units': sorted(str(u) for u in expected_units),
                'mad_scale': MAD_SCALE, 'candidate_type': EVENT_TYPE, 'origin': ORIGIN}
    result = DetectionResult(diagnostic, units, mode, 'ok', settings, notes=[NOT_EVIDENCE])
    if mode == 'noncausal':
        result.notes.append('noncausal baseline uses future samples: not for use as a forecasting feature')
    empty = {'total_steps': 0, 'assessed': 0, 'unassessed': 0, 'unassessed_by_reason': dict.fromkeys(REASONS, 0),
             'unassessed_intervals': [], 'total_samples': 0}

    def finish(status, note):
        result.status, result.coverage = status, dict(empty)
        result.notes.append(note)
        return result

    if values is None or times is None:
        return finish('missing_diagnostic', f'diagnostic {diagnostic!r} is absent; nothing was assessed')
    t = _times(times)
    y = np.ma.filled(np.ma.asarray(values, dtype=float), np.nan).astype(float)
    if y.shape != t.shape:
        raise ValueError('times and values must be 1-D arrays of equal length')
    n = t.size
    exc = np.zeros(n, dtype=bool) if excluded is None else np.array(excluded, dtype=bool)
    if exc.shape != t.shape:
        raise ValueError('excluded must match times in shape')
    finite = np.isfinite(y)
    valid = finite & ~exc
    empty['total_samples'] = n
    if n < 2 or not valid.any():
        return finish('no_data', 'fewer than two samples or no valid samples')
    dt = np.diff(t)
    median_dt = float(np.median(dt))
    result.settings['median_sample_interval'] = median_dt
    if median_dt > max_sample_interval:
        return finish('insufficient_cadence', f'median spacing {median_dt:g} s exceeds max_sample_interval '
                      f'{max_sample_interval:g} s; no assessment made')
    gap = dt > max_gap_factor * median_dt
    steps = n - 1                                   # step k (1..n-1) joins samples k-1 and k; index k-1 below
    reason = [None] * steps
    for i in range(steps):
        if not (finite[i] and finite[i + 1]):
            reason[i] = 'missing'
        elif exc[i] or exc[i + 1]:
            reason[i] = 'excluded'
        elif gap[i]:
            reason[i] = 'gap'
    inseg = np.array([r is None for r in reason])
    d = np.where(inseg, np.diff(np.where(finite, y, 0.0)), np.nan)
    half = window // 2
    flagged = np.zeros(steps, dtype=bool)
    score = np.full(steps, np.nan)
    rel = np.full(steps, np.nan)
    for i in range(steps):
        if not inseg[i]:
            continue
        lo_seg = i
        while lo_seg > 0 and inseg[lo_seg - 1]:
            lo_seg -= 1
        hi_seg = i
        while hi_seg < steps - 1 and inseg[hi_seg + 1]:
            hi_seg += 1
        if mode == 'causal':
            lo, hi = max(lo_seg, i - window), i - 1
        else:
            lo, hi = max(lo_seg, i - half), min(hi_seg, i + half)
        base = np.array([d[j] for j in range(lo, hi + 1) if j != i], dtype=float)
        if base.size < min_baseline:
            reason[i] = 'insufficient_baseline'
            continue
        sigma = MAD_SCALE * np.median(np.abs(base - np.median(base)))
        if not sigma > 0:
            reason[i] = 'insufficient_baseline'     # zero spread: the score would be unbounded, so not assessed
            continue
        score[i] = -d[i] / sigma
        seg_first_sample = lo_seg
        level = float(np.median(y[max(seg_first_sample, i - window + 1):i + 1]))
        if level > 0:
            rel[i] = -d[i] / level
        flagged[i] = bool(score[i] >= threshold and rel[i] >= min_relative_drop)
    counts = dict.fromkeys(REASONS, 0)
    for r in reason:
        if r is not None:
            counts[r] += 1
    spans = _merge([[float(t[i]), float(t[i + 1])] for i in range(steps) if reason[i] is not None])
    unassessed = sum(counts.values())
    result.coverage = {'total_steps': steps, 'assessed': steps - unassessed, 'unassessed': unassessed,
                       'unassessed_by_reason': counts, 'unassessed_intervals': spans, 'total_samples': n,
                       'excluded_samples': int(exc.sum()), 'missing_samples': int((~finite).sum())}
    if unassessed == steps:
        result.status = 'insufficient_baseline' if counts['insufficient_baseline'] else 'no_data'
        result.notes.append('no step had a usable baseline; nothing was assessed')
        return result
    i = 0
    while i < steps:
        if not flagged[i]:
            i += 1
            continue
        j = i
        while j + 1 < steps and flagged[j + 1]:
            j += 1
        t0, t1 = float(t[i]), float(t[j + 1])
        result.candidates.append({
            'type': EVENT_TYPE, 'origin': ORIGIN, 'diagnostic': diagnostic, 'units': units, 'method': METHOD,
            'method_version': METHOD_VERSION, 'mode': mode,
            'start': t0 - timing_padding, 'end': t1 + timing_padding, 'interval_closed': True,
            'unpadded_interval': [t0, t1], 'timing_uncertainty': t1 - t0, 'timing_padding': float(timing_padding),
            'score': _num(np.max(score[i:j + 1])), 'score_is_probability': False,
            'relative_drop': _num(np.max(rel[i:j + 1])),
            'sample_before_index': i, 'sample_after_index': j + 1, 'first_step_index': i + 1,
            'last_step_index': j + 1, 'review_state': 'unreviewed'})
        i = j + 1
    return result
