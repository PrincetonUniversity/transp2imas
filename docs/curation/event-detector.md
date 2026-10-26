# Optional event-candidate detector: `abrupt_drop`

Module `curation/detectors.py`, test `examples/verify_curation_detectors.py`. Method name `abrupt_drop`, version `1.0`.

## Purpose and non-purpose

Purpose: flag intervals where one explicitly identified, unit-checked scalar time history (for example a central
electron temperature trace) decreased abruptly relative to its own recent noise. These are **candidates** for an
analyst to review.

Not a purpose: classification. The type is always `unclassified_transient` and the origin always
`automatic_candidate`. A sudden decrease of a trace is not thereby a sawtooth crash, ELM or pellet; diagnostic
glitches, noise bursts and heating changes produce the same signature. The detector never establishes that a plasma
was event-free, and its output never acts as an exclusion until an analyst reviews and promotes it
(`review_state` is `unreviewed`).

## Inputs

`abrupt_drop(times, values, units, diagnostic, mode, expected_units, max_sample_interval, ...)`; `times` and
`values` are plain 1-D arrays (`values` may be a masked array or contain NaN; `values=None` means the diagnostic is
absent). Inputs are never modified.

* `times`: native clock, seconds (`time_units='s'` must be declared), finite and strictly increasing, else ValueError.
* `units` must be given and must match, after upper-casing and blank removal, one of `expected_units` (e.g.
  `{'EV','KEV'}`), else ValueError. Units are never guessed from magnitude.
* `diagnostic`: required explicit variable/channel name. `mode`: required, `'causal'` or `'noncausal'`, no default.
* `max_sample_interval` (s): required; chosen by the caller for the physics of interest.
* `excluded`: optional boolean mask of samples a policy excludes.

## Algorithm (exact)

1. A sample is valid if finite, unmasked and not excluded. Segments are maximal runs of valid samples with no time
   step larger than `max_gap_factor` x the median spacing of `times`. Windows never cross a segment boundary;
   excluded samples are neither used for baseline or level nor bridged.
2. For each pair of consecutive samples in one segment, step k has increment `d_k = y_k - y_{k-1}`.
3. Baseline increments (step k itself excluded): `causal`, the `window` increments before k; `noncausal`, the
   increments within `window//2` either side of k. Both are truncated at the segment ends. Fewer than
   `min_baseline` increments, or zero spread, leaves the step unassessed (`insufficient_baseline`).
4. `score_k = -d_k / (1.4826 x MAD)`, MAD the median absolute deviation of the baseline increments about their median.
5. Local level: median of up to `window` samples of the segment ending at `y_{k-1}` (past only, both modes).
   `relative_drop = -d_k / level`, defined only for a positive level.
6. Step k is flagged when `score >= threshold` and `relative_drop >= min_relative_drop`. Consecutive flagged steps
   merge into one candidate.
7. Candidate interval: closed `[t_{k0-1}, t_{k1}]` (the drop happened between two samples), widened by
   `timing_padding` on both sides. `timing_uncertainty` is the unpadded interval width (one sample spacing, more
   for merged steps); timing error up to that width is inherent. Also reported: `score` (maximum step score),
   `relative_drop` (maximum), and `sample_before_index`, `sample_after_index`, `first_step_index`, `last_step_index`.

Defaults: `threshold=5.0`, `min_relative_drop=0.05`, `window=20`, `min_baseline=8`, `max_gap_factor=3.0`,
`timing_padding=0.0`. `max_sample_interval`, `expected_units`, `units`, `diagnostic`, `mode` have none.

## Causal and noncausal

`causal` reads only samples at or before the evaluated step for the noise baseline (testable: truncating the
future leaves earlier results unchanged), so it can be used when preparing forecasting data. `noncausal` uses a
centered window, may react to changes after the step, and must not be used for forecasting features; its result
carries a note saying so. The two agree for stationary noise and differ when the variance changes around the step.

## Status, coverage and what "no candidates" means

Status: `ok`; `missing_diagnostic` (no data given); `no_data` (fewer than two samples, no valid sample, or no step
assessed for another reason); `insufficient_cadence` (median spacing > `max_sample_interval`, nothing assessed);
`insufficient_baseline` (no step had a usable baseline). Every result carries the note that absence of candidates
is not evidence of an event-free plasma.

`coverage` reports `total_steps`, `assessed` (steps with a valid baseline), `unassessed` (missing, excluded, gap,
insufficient baseline, with counts per reason) and `unassessed_intervals` as `[start, end]` lists; always
`assessed + unassessed == total_steps`. No candidates with status `ok` means only that no assessed step met the
criteria. It says nothing about unassessed intervals, drops smaller than `min_relative_drop`, events on other
diagnostics, or events too fast or slow for this method and cadence.

## Failure modes

* False positives: noise bursts (non-Gaussian or heavier-tailed noise than the MAD scale assumes), diagnostic
  glitches, spikes followed by recovery, a quiet baseline directly before a legitimate noise increase (causal mode).
* False negatives: drops under `min_relative_drop` even with a huge score; drops spread over many samples; drops
  in a noisy baseline; drops inside missing, excluded or gapped intervals (deliberately not bridged); the first
  `min_baseline` steps of every segment.
* Timing error up to the sample spacing; slow cadence is rejected rather than interpolated.
* Traces with non-positive level cannot produce a relative drop and so cannot be flagged.

## Score is not a probability

`score` is a ratio of a drop to a robust noise scale. Its relation to the chance of a real event is not calibrated;
no calibration data is used. Keep it as a score (`score_is_probability: false`); do not threshold it as if it were
a probability or combine it as one.

## Synthetic example (all numbers made up)

```python
import numpy as np
from curation.detectors import abrupt_drop

rng = np.random.default_rng(1)                       # SYNTHETIC trace, not a measurement
t = np.arange(400) * 1.0e-3
te = 2.0 + 0.01 * rng.standard_normal(400)           # keV
te[250:] *= 0.7                                      # synthetic abrupt drop between samples 249 and 250
r = abrupt_drop(t, te, 'keV', 'SYNTH_TE0', 'causal', {'EV', 'KEV'}, max_sample_interval=0.01)
print(r.status, [(c['start'], c['end'], round(c['score'])) for c in r.candidates])
print(r.coverage['assessed'] + r.coverage['unassessed'] == r.coverage['total_steps'])
```

The candidate interval is `[0.249, 0.250]`; the score varies with the seed. Synthetic traces test the algorithm's
documented behaviour only, never the physical validity of a candidate.
