# Events: origins, membership, policy and the optional detector

An event is a **catalog entry**, not a detection. Curation uses events to say which samples an application should not
use, never to say that a physical event happened. Declared events and the declarations they replace are described in
[evidence-format.md](evidence-format.md); the detector algorithm is in [event-detector.md](event-detector.md).

## What an event records

| Field | Meaning |
|---|---|
| `id`, `type` | identity; `startup`, `termination`, `sawtooth`, `elm`, `pellet`, `heating_transition`, `diagnostic_outage` or `unclassified_transient` |
| `start`, `end` | the closed interval in seconds |
| `origin` | `experimental_observation`, `actuator_command`, `modeled_transp`, `analyst_declaration` (default), `automatic_candidate` |
| `reference`, `source`, `diagnostic` | where the event comes from (record, diagnostic, log) |
| `detector` | `{method, version, settings}` when an algorithm produced it |
| `timing_uncertainty`, `padding` | seconds added on both sides of the interval before membership is computed |
| `review_state` | `unreviewed`, `reviewed_accepted`, `reviewed_rejected` |
| `rules`, `applications` | the rule ids (shell-style patterns) and applications the exclusion applies to |
| `policy` | `{"kind": "legacy_exclusion", "rules": [...]}`, derived from `rules` |

Origin matters. A *commanded* pellet is not a delivered pellet; a *modeled* sawtooth is TRANSP's own assumption and not
independent experimental evidence; an *automatic candidate* is a score above a threshold. The inventory reports
`independent_experimental_evidence` as true only for `experimental_observation`, and every exclusion record lists the
origin of the events behind it. Migrating legacy `events[].rules` does **not** assert that the event physics was
verified: they become `legacy_exclusion` policies of an `analyst_declaration`.

## Membership

A sample is a *member* of an event when its support interval meets the closed window
`[start - pad, end + pad]`, with `pad = timing_uncertainty + padding`.

* The support interval comes from `sample_support[clock]` (`averaging_interval`, `alignment` `center` or `trailing`) when
  declared. With none declared the sample is treated as a point, which is **labelled an approximation**
  (`support_basis: point_sample_approximation`, otherwise `declared_averaging_interval`).
* Endpoints are inclusive; simultaneous events union their members; each event's own members are also kept so counts
  can be attributed (`excluded_by_event`).
* `TIME` and `TIME3` are separate clocks. Membership is computed on the clock named by the variable's own dimension and
  is never borrowed across, even when the two arrays have equal lengths.

An event excludes samples only if it **may act as an exclusion**: not a `reviewed_rejected` event, and not an
`automatic_candidate` until an analyst marks it `reviewed_accepted`. It never acts on source-integrity rules
(`GRID.*`, grid mismatch, units, read errors, run identity), including with `rules: ["*"]`; those rules reached by a
pattern are listed in the inventory as `mandatory_rules_not_waivable`.

## What is recorded, and what is not changed

For each rule an exclusion applies to, `report.policy_exclusions[rule id]` records `total_count`, `unassessed_count`
(outside the assessment window or otherwise not evaluated), `excluded_count`, `eligible_count`, the counts per event,
how many *flagged* samples are excluded, whether the exclusion covers every flagged sample, and `policy_status`
(`not_assessed` when nothing eligible remains). The identities are exact:
`excluded + eligible = evaluated` and `total = evaluated + unassessed`. The `*__excluded` mask is in `curation.npz`.

The raw finding, its status and its `__evaluated`/`__flagged` masks are **unchanged**. Event membership is also kept
distinct from missingness: a masked sample inside an event is still reported as missing.

`EVENTS.INVENTORY` warns about an event whose padded window misses the source time, or that matched no rule it could act
on, instead of ignoring it.

## The optional detector

```json
"event_detectors": [{
  "id": "te0_drop", "method": "abrupt_drop", "diagnostic": "TE0", "clock": "TIME", "mode": "causal",
  "expected_units": ["EV"], "max_sample_interval": 0.005,
  "parameters": {"threshold": 5.0, "min_relative_drop": 0.05}
}]
```

The diagnostic is named explicitly, its units are checked against `expected_units` (a mismatch is an `error`
finding), a 2-D profile needs a `radial_index`, and no other axis is borrowed. Candidates are recorded in
`report.detected_events` as `unclassified_transient`, origin `automatic_candidate`, `review_state: unreviewed`, with an
uncalibrated `score` (never a probability) and a timing uncertainty equal to the sample spacing. They exclude nothing
until an analyst copies one into `events` as `reviewed_accepted`.

* **No candidates is not evidence of an event-free plasma.** The finding says so; a missing diagnostic, too slow a
  cadence, or too little baseline is `not_assessed`, never a pass.
* Samples can be removed from the detector only by a declared event that names `EVENTS.DETECT.<id>` in its rules, so
  the exclusion is a recorded policy. The detector never bridges missing samples or excluded samples.
* `mode` is explicit: `causal` uses only samples at or before each step (usable when preparing forecasting data);
  `noncausal` may center its noise baseline and must not feed forecasting features.

See [event-detector.md](event-detector.md) for the algorithm, defaults, false-positive/negative behaviour and the
limits of testing it on synthetic traces.
