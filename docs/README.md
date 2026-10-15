# transp2imas documentation

This directory contains the maintained technical and developer documentation for
`transp2imas`. The top-level README is intentionally focused on installation and
basic use; the documents here explain the design, validation strategy, mapping
policy, and release process in more detail.

## Start here

- [All-example verification](examples-verification.md) -- automatic CDF discovery,
  all-run reports, PF/waves coverage, and known input/current warnings.
- [Plotting converted IDSs](plotting.md) -- real `plot_ids.py` commands and options.
- [Input provenance](input-provenance.md) -- XML, source references, and model code libraries.
- [Neoclassical theory transport](neoclassical-transport.md) -- GFLNC/QFLNC divergences and SI flux reconstruction.
- [Neutral profiles](neutral-profiles.md) -- independent wall/volume/recombination density, temperature and rotation.
- [Sawtooth crash history](sawteeth.md) -- C-Mod `TLASTSAW` events, typed q=1 radius, native `TIME` and model/event limitations.
- [Unified theory transport](theory-transport.md) -- MMM07 and GLF23 outputs, TGLF-ready registry, source coefficient closures, and model-specific output XML.

- [Architecture](architecture.md) -- module layout, control flow, shared helper
  classes, and where to add new conversion logic.
- [Development guide](development.md) -- local environment, coding conventions,
  Python header policy, and contribution workflow.
- [Passive conductors](pf-passive.md) -- PS_n current conventions, authority-bound
  geometry JSON, DD3/DD4 outline rules, and verification procedure.
- [Testing and validation](testing.md) -- what each `make` target checks and which
  checks require an installed IMAS Data Dictionary.
- [Mapping workflow](mapping-workflow.md) -- the group-by-group procedure used to
  extend TRANSP variable coverage without guessing IMAS semantics.
  This document also describes `examples/audit_mapping_coverage.py` and the
  generated Markdown/CSV/JSON source-coverage reports used to select the next group.
- [Release process](release-process.md) -- version files, Conventional Commits,
  Release Please, GitHub settings, and how the initial `v0.9.0` baseline is
  established.
- [Technical reference](technical-reference.md) -- detailed CDF-to-IDS mappings,
  units, sign conventions, compatibility notes, provenance, and historical
  development audits.

## Documentation policy

Documentation should be updated in the same pull request as a behavior or
mapping change. In particular:

1. New or changed CDF-to-IMAS mappings belong in `technical-reference.md`.
2. New development or validation procedures belong in `development.md` or
   `testing.md`.
3. Changes to versioning or GitHub Actions belong in `release-process.md`.
4. New IDS groups should be added to the status/roadmap section in
   `mapping-workflow.md`.

The pre-0.9 internal 8.x labels retained in the technical reference are
historical development milestones, not public semantic versions.

- [pf_active](pf-active.md) - audited coil/supply channels, KT_ target-current exclusion, machine geometry, and DD3/DD4 circuits.

- [Real ITER/NSTX ISolver PF source validation](isolver-fixtures.md) — coil and passive-current CDF channels, CC/KK effective turns, file layout and large-fixture guidance.

- [TIME/TIME3 mapping and equilibrium regression](time-grids.md)

- [End-to-end IMAS readback](end-to-end.md): convert every example, reopen written IDSs and validate native TIME3.
