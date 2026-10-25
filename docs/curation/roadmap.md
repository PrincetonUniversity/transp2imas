# Curation roadmap and verification plan

## Stage 1: implemented foundation

Independent read-only CDF access; native grids; declared input/evidence records; local evaluated/flagged
masks; deterministic rules; preliminary application summaries; JSON/Markdown/NPZ reports; strict
sidecar validation; synthetic corruption and CLI tests. This is a screening tool, not a data-repair
or scientific-certification system.

Before operational adoption, review rule tolerances and missing-data policies with the relevant
TRANSP/diagnostic experts. Run the package in the site's normal environment and retain generated
reports with the input/model revision.

## Stage 1 next changes: source integrity and validated model state

1. Fix and independently test the NBI time-axis/source-identity issues identified in the t2i2 review.
2. Harden the converter's own namelist scalar parsing and masked-data handling. The curation reader
   handles these safely but does not repair the converter.
3. Define a shared, stable source/provenance protocol between conversion and curation after those
   contracts are tested. Avoid building a second independent physics-mapping implementation.
4. Ingest **resolved** run inputs, code version and referenced PT model namelists. Preserve raw
   assignments, uncertainty about omitted defaults and time-dependent validity. Verify additive
   versus replacing model semantics against TRANSP before turning them into hard rules.
5. Add authoritative source/campaign geometry and waveform records, approved override diffs, and
   power-accounting boundaries. Never force a transport match by correcting inputs automatically.
6. Add event masks and explicit multi-diagnostic fit/provenance adapters.

**Status of items 1 and 2 (converter side).** Done, each verified against the bundled cases:

- *NBI clocks and identities.* The time axis was zeroed for every beam without half- and third-energy
  components (both ITER beams, all 14 DIII-D 163303U09 beams, two of three NSTX beams: all-zero times in the
  written IDS). Native `TIME` is now written unchanged and an off beam is expressed by its power. Units keep
  their own `PINJ<n>` numbers, a single species total applies to every unit, several species are not matched by
  list position, and the species-only fallback uses the channel's own suffix (`verify_nbi_mapping.py`,
  `nbi_native_time`).
- *Fortran numeric literals.* `D` exponents (`1.0D-3`, `2.5d0`, `1D3`) were read as 0 after a printed error; the
  bundled NSTX `TR.DAT` files contain `TBONA(1) = 9.205d-2` and `XLBTNA(3) = 1.142d3`. They are read, a `!` inside
  quotes is not a comment, and a value that is not a number is a logged unknown instead of 0.
- *Masked source data.* A masked (fill-value) sample is a named broken mapping instead of a fill value passed on as
  a physical number. No bundled CDF has masked samples, so this closes a latent gap; the missing-data policy for the
  IDS itself is still to be decided.
- *Profile shape.* `check_profile_shape` accepted a wrong time length and any other dimensionality; it now requires
  exactly `(nt, nx)` for a present profile.

Still open from the list: item 3 (a shared source/provenance protocol), 4 to 6 (resolved inputs, authoritative
geometry and events), and the NBI torque total when `TQRPL` is nonzero and `FTQRIP` is unavailable.

Acceptance: representative real runs plus deliberately incorrect source identities/units/times,
valid missing-data cases, and independent expert-reviewed expected outcomes. A bare parser/test pass
is not a valid scientific acceptance criterion.

## Stage 2: uncertainty-aware transport products

Implement correlated profile/gradient ensembles, systematic diagnostic shifts, equilibrium/radial
registration variations, joint charge-constrained composition ensembles, rotation hypotheses and
heating/source perturbations. Use consistent TRANSP reruns for effects that change the solution;
do not alter an impurity array in an already computed result and call it propagated uncertainty.

Separate flux validity from coefficient identifiability. Include covariance/shared-source effects
in diagnostic/model comparisons. Treat Monte Carlo seed variation separately from physical input
variation. Record ensemble identity, member status, failed members and parent data.

Acceptance: synthetic known-truth uncertainty calibration; coverage of weak and sign-changing
gradients; reproducible same-seed ensembles; no train/test leakage among related realizations.

## Stage 3: application-specific curated datasets

Calibrate acceptance policies for power balance, heat transport, momentum studies, current-profile
evolution, model validation and surrogate training. Add independent diagnostic comparisons and
sensitivity-based root-cause ranking. Support reviewed decisions and waivers without erasing findings.
Report precision as well as agreement; inflating uncertainty must not produce a high-quality badge.

Integrate appropriate typed IMAS validity/fitting/numerics fields after DD3/DD4 semantic tests.
Keep richer uncertainty/decision provenance in versioned records. Regional occurrences remain model
views; uncertainty ensemble members use separate dataset/run identities.

Acceptance: site IMAS write/readback on supported DD versions, exact source/record lineage, independent
expected-IDS inventory, and reviewer assessment of representative datasets. No mock test establishes
scientific correctness or serialization by itself.

## Test layers

- Unit: literal parsing, masks, unit conversions, source identity, strict manifest validation,
  thresholds and rule boundaries.
- Mutation: wrong beam energy/fractions, reordered times, negative profile, unknown rotation,
  composition omissions/double counts, masked fill values, unsuitable uncertainty arrays.
- Integration: real classic CDF and netCDF4/HDF5 inputs, JSON/Markdown/NPZ roundtrip, CLI exits,
  protected input paths, native TIME versus TIME3, user-selected windows.
- Scientific: expert labels, conservation/accounting identities, genuinely independent diagnostics,
  event-aware policies and known-truth uncertainty tests.
- Backend: later IMAS export validation and numerical readback, with occurrence/model-aware expectations.

Do not silently skip every real fixture and call the run verified. Source-only archives can pass
`make verify-curation`; full-repository checks must keep their fixture requirements explicit.
