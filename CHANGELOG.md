# Changelog

Notable changes to transp2imas are recorded here. Future entries and GitHub
releases are managed automatically by Release Please from Conventional Commit
messages.

## 0.9.0 - 2026-09-28

Initial public pre-1.0 baseline of the IMAS-Python converter.

- Support DD 3.39 and DD 4.1.1 conversion paths with centralized COCOS/sign handling.
- Expand high-confidence TRANSP mappings across core profiles, equilibrium, summary,
  core sources, core transport, radiation, edge sources, actuator IDSs, and NBI
  distribution sources.
- Populate NBI parallel current from `CURB * CONPLJB / BZ`.
- Preserve populated `core_sources` entries with IMAS-Python arrays of structures.
- Add workflow, dataset FAIR metadata, source hashes, and per-IDS provenance.
- Add source-side regression checks for the bundled NSTX cases and static writer/API
  validation.
