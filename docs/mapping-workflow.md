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

## Current broad status at the 0.9.0 baseline

### Substantially covered

- `core_profiles`
- `equilibrium`
- `summary`
- `core_sources`
- `core_transport`
- `radiation`
- `edge_sources`
- NBI `distribution_sources`
- `nbi`, `ec_launchers`, `ic_antennas`, `lh_antennas`
- `tf`, `pf_active`, `wall`
- `workflow`, `dataset_fair`, and per-IDS provenance

### Implemented but dependent on richer source data

- `waves` -- requires spatial RF absorption/deposition data.
- `pf_passive` -- requires passive-structure currents plus trustworthy structure
  geometry/electrical definitions.
- `sawteeth` -- requires explicit sawtooth crash-event data.

### Highest-value next group

The preferred next physics target is full fast-ion `distributions` from a
TRANSP/NUBEAM run containing the full FBM phase-space arrays and their energy,
pitch, species, radial, and time coordinates. `distribution_sources` already
provides a natural validated stepping stone but must not be confused with a
full phase-space distribution.

If an FBM regression case is not available, the next target should be chosen by
which validation-quality source case is available first: `pf_passive`, `waves`,
or `sawteeth`.

## Intentionally unmapped or deferred examples

Some quantities are deliberately left unmapped until the correct IMAS semantic
and metric/grid treatment are established. Examples include neutral source-state
temperatures/rotation, boundary energy transport quantities, reduced fast-ion
moments when full phase-space coordinates are absent, and quantities on an
incompatible radial grid.

"Unmapped" is not equivalent to "forgotten". The technical reference should
record important deferred quantities and the reason for deferral.
