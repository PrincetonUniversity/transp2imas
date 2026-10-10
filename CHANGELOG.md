# Changelog

Notable changes to transp2imas are recorded here. Future entries and GitHub
releases are managed automatically by Release Please from Conventional Commit
messages.

## 0.9.0 - 2026-09-28

Initial public pre-1.0 baseline of the IMAS-Python converter.

Code-quality fixes included in the baseline:

- fix duplicate atomic-species dictionary keys for lithium and Hassium/Holmium;
- correct the Dysprosium key and Technetium atomic metadata discovered in the same audit;
- remove process-global `numpy.seterr(all="raise")` from the conversion entry point;
- use explicit imports in `convert.py`;
- replace broad numerical exception handling in the temperature-peaking/shared-array paths with local NumPy error handling and specific expected exceptions;
- add a source-only duplicate-literal-key check to `make check`.

- Support DD 3.39 and DD 4.1.1 conversion paths with centralized COCOS/sign handling.
- Expand high-confidence TRANSP mappings across core profiles, equilibrium, summary,
  core sources, core transport, radiation, edge sources, actuator IDSs, and NBI
  distribution sources.
- Populate NBI parallel current from `CURB * CONPLJB / BZ`.
- Preserve populated `core_sources` entries with IMAS-Python arrays of structures.
- Add workflow, dataset FAIR metadata, source hashes, and per-IDS provenance.
- Add source-side regression checks for the bundled NSTX cases and static writer/API
  validation.
- Validate and extend the `waves` writer with bundled DIII-D TORBEAM/GENRAY EC
  and C-Mod multi-frequency ICRF cases. IC coherent waves now use resolved
  `PIC_Fj` absorbed-power profiles with `RFPWRABj` closure; `QICHA` is only a
  single-wave fallback. `PEICH`, `PIICH`, and `QRFBTOT` remain plasma-source
  quantities rather than being incorrectly summed as wave absorption.
- Add `examples/verify_waves_mapping.py` and `make verify-waves` for portable
  RF source-side regression checks.
- Add `examples/audit_mapping_coverage.py` and `make audit` to classify CDF variables
  as converter-referenced, derived/support, validation-only, deliberately deferred,
  or still unclassified, with Markdown/CSV/JSON reports.
