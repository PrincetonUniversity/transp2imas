# Mapping workflow and roadmap

## Objective

The project goal is to preserve as much useful TRANSP information in IMAS as
possible while keeping every mapping physically defensible and testable. Work
is organized by related IDS/physics groups rather than by isolated variable
names.

## Procedure for each group

For each candidate group:

1. **Inventory the source.** Identify relevant CDF variables, dimensions,
   attributes, namelist controls, and optional companion files.
2. **Establish semantics.** Determine what each TRANSP quantity actually means,
   not just what its variable name suggests.
3. **Record source units.** Define the exact conversion into the IMAS unit.
4. **Identify the grid.** Determine time base, zone/zone-boundary placement,
   radial coordinate, species dimension, and any interpolation required.
5. **Resolve signs.** Apply COCOS/device-direction conventions centrally and
   document them.
6. **Choose the IMAS destination.** Confirm that the node has the same physical
   meaning and representation.
7. **Implement conservatively.** Missing optional inputs should normally skip
   the optional field/group rather than fabricate data.
8. **Validate independently.** Add source-side closure/integral/identity checks
   before considering the mapping complete.
9. **Check DD compatibility.** Verify both DD3 and DD4 paths when the project
   claims support for both.
10. **Document the mapping.** Update the technical reference and this roadmap.


## Mapping coverage audit

Use `examples/audit_mapping_coverage.py` to turn the group-by-group objective into
a reproducible source inventory. The audit scans every variable in a TRANSP CDF
and compares it with static references in the converter and source-side
verification programs. It classifies each variable as:

- `mapped/consumed` -- referenced by converter code;
- `support/derived-input` -- a grid, geometry, normalization, or integration
  input used to construct another mapped quantity;
- `metadata/provenance` -- used only by metadata/provenance code;
- `validation-only` -- known to a verifier but not the converter;
- `deferred` -- reviewed and deliberately left unmapped with a recorded reason;
- `unclassified` -- no static converter/verification/deferral reference was
  found and the quantity remains a candidate for review.

Run the standard report with:

```bash
make audit CASE=141711Y01
```

This writes Markdown, CSV, and JSON inventories below `build/audit/`. A direct
invocation can be used for any external TRANSP CDF:

```bash
python examples/audit_mapping_coverage.py /path/to/run.CDF \
    --markdown build/audit/run.md \
    --csv build/audit/run.csv \
    --json build/audit/run.json
```

The coverage percentage is intentionally conservative and must not be treated
as a physics-completeness metric. Static analysis can miss dynamically composed
CDF names, and a referenced quantity may serve only as a grid or normalization.
The full table is meant to drive the next manual physics audit, not replace it.
When a previously unclassified quantity is investigated, either map it, add a
source-side validation reference, or record a specific deferral reason.

## Current broad status at the 0.9.0 baseline

### Substantially covered

- `core_profiles`
- `equilibrium`
- `summary`
- `core_sources`
- `core_transport`
- `radiation`
- `edge_sources`
- `waves` -- EC and IC absorbed/deposited wave output is now validated against
  dedicated DIII-D and C-Mod regression cases.
- NBI `distribution_sources`
- `nbi`, `ec_launchers`, `ic_antennas`, `lh_antennas`
- `tf`, `pf_active` (source-verified waveforms and optional external geometry), `wall`
- `workflow`, `dataset_fair`, and per-IDS provenance

### Implemented but dependent on richer source data

- `pf_active` -- mapped KK_/VV_/VS_ waveforms, validated per-channel time and unit conversions,
  optional DD3/DD4 coil/circuit device model, and intentionally deferred KT_
  target-current references. No regression CDF contains active PF output; see
  [pf_active](pf-active.md).
- `pf_passive` -- implemented for source-side PS_n currents and explicitly supplied
  authoritative geometry, including DD3/DD4 outline conversion. Actual passive
  current/geometry identity awaits a PS_n regression CDF; see [pf_passive](pf-passive.md).
- `sawteeth` -- requires explicit sawtooth crash-event data.

### Highest-value next group

The preferred next physics target is full fast-ion `distributions` from a
TRANSP/NUBEAM run containing the full FBM phase-space arrays and their energy,
pitch, species, radial, and time coordinates. `distribution_sources` already
provides a natural validated stepping stone but must not be confused with a
full phase-space distribution.

If NUBEAM/FBM extension is postponed, `pf_passive` is implemented with
explicit geometry and source-only synthetic regression, but it still needs a real
`PS_n` CDF with matching authoritative machine description. The next unfilled
physics IDS is `sawteeth`, contingent on a validated crash-event source case.

## Intentionally unmapped or deferred examples

Some quantities are deliberately left unmapped until the correct IMAS semantic
and metric/grid treatment are established. Examples include neutral source-state
temperatures/rotation, boundary energy transport quantities, reduced fast-ion
moments when full phase-space coordinates are absent, and quantities on an
incompatible radial grid.

"Unmapped" is not equivalent to "forgotten". The technical reference should
record important deferred quantities and the reason for deferral.

### Neoclassical and neutral-state extension

`core_transport` neoclassical model sources are reconstructed from native
TRANSP flux divergences, never directly copied as transport flux.
`core_profiles` neutral state is mapped per wall/volume/recombination
population with verified density/temperature units and separate code XML for
rotation on DD4. See [neoclassical-neutrals.md](neoclassical-neutrals.md).
Regression: `make verify-neoclassical-neutrals`, then `make verify-all`.

