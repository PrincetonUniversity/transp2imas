# Architecture

## Purpose

`transp2imas` converts one TRANSP output CDF, plus its associated input namelist
and referenced model inputs when available, into a set of IMAS Interface Data
Structures (IDSs). The design deliberately separates common NetCDF/IMAS helpers
from one-writer-per-IDS mapping modules.

## Top-level control flow

The normal conversion path is:

1. `convert.py` parses command-line options and resolves the input CDF.
2. `nc2imas.py` opens the TRANSP CDF, parses metadata/namelists, determines the
   target Data Dictionary behavior, and provides shared conversion helpers.
3. `imas_class.py` manages IMAS database-entry creation/opening and common IDS
   write operations.
4. Individual `imas_*.py` modules populate one IDS family at a time and submit
   what they built.
5. `ConversionRunner` (in `conversion_plan.py`) calls each selected writer in
   order, then flushes that writer's IDSs to the database and journals the
   outcome. Common provenance is attached before each IDS is submitted.
6. Source-side verification tools under `examples/` independently check units,
   identities, sign assumptions, and mapping prerequisites.

## Main modules

### Input and shared conversion layer

- `netcdf_class.py` -- TRANSP CDF variable/dimension/attribute access.
- `nc2imas.py` -- shared mapping functions, sign/COCOS handling, version
  metadata, namelist/template capture, code/provenance XML, and helper methods
  used by IDS writers.
- `utilities.py` -- general helper functions.
- `atoms_class.py` -- atomic/species helper logic.

### IMAS access layer

- `imas_class.py` -- backend selection, `DBEntry` lifecycle, IDS factory access
  and the readback manifest of accepted puts.
- `dd_adapter.py` -- Data Dictionary version questions (`at_least`, `before`,
  `is_dd3`, `is_dd4`) and target-version selection; no IMAS dependency.
- `imas_store.py` -- `WriteRequest` and `ImasStore`: a writer's `putdb_imas*`
  call submits a request; with deferral on, the request is held until the runner
  flushes it, so constructing a writer does not itself write to the database.

### Composition status

The review asks for `nc2imas` to be replaced gradually by composed services.
Done so far, behind the unchanged `nc2imas` interface: DD version logic
(`DDAdapter`), write queuing (`ImasStore`, `WriteRequest`), the writer loop
(`ConversionRunner`), the COCOS/direction sign policy (`cocos_signs.resolve_signs`
-> immutable `CocosSigns`) and conversion provenance (`code_provenance`: source-file
records, `ids_properties` text, `code.parameters` XML) and the assignment
primitives (`ids_assignment`: `is_valid`, array arithmetic, `fill_data`, `fill_value`,
which carry the strict/lenient mapping policy). A writer whose build raises submits none of its own output. A backend failure during a multi-IDS
flush can leave earlier accepted IDSs in the entry; this is not rollback/transactional storage. The
journal records accepted writes and the failure, and later writers do not run. Tests: `verify_conversion_runner.py`, `verify_cocos_signs.py`,
`verify_code_provenance.py`, `verify_ids_assignment.py`, `verify_conversion_services.py`.

**Writer-facing context.** Writers no longer take the whole converter. `conversion_context.py`
groups the 75 names they used into six read-only views, each with an explicit list of the names it
forwards: `cdf` (reading the TRANSP file), `namelist`, `dd` (version predicates, `new_ids`),
`signs` (COCOS signs and the vacuum field), `grid` (time bases, unit factors, interpolation) and
`fill` (assigning to IDS nodes, identities, properties and code blocks), plus `ctx.mapping` for the
ledger. `IdsWriter` exposes them as `self.cdf`, `self.signs`, and so on, derived from `self.im`, so a
writer reads `self.signs.ip_sign` rather than reaching into `self.im`. A view forwards to the converter
at call time, so a test double that offers the listed names works as the owner, and a name a view does
not list is an error that names the view. `verify_conversion_context.py` fails if a list names something
`nc2imas` lacks, if a writer uses an unlisted name or a view outside an `IdsWriter`, or if a writer
reaches a view-owned name through `self.im` again.

**What moved out of `nc2imas`.** The implementations behind the views are now separate modules, and
`nc2imas` keeps the historical method names as one-line delegators (it was 909 lines, now about 500):

- `namelist_reader.py`, `name_patterns.py` -- loading and querying the TR.DAT namelist, number parsing and
  the `?` pattern builders, as functions over an explicit `nmlist` dict and logger;
- `grid_services.py` -- reading the `TIME`/`TIME3` bases, time interpolation, moving a 1-D history between
  the two bases, zone-to-enclosed quantities, radial gradients (`radial_interpolation` already held the
  radial mapping);
- `ids_identity.py` -- the DD-version identity setters (species name or label, `z_n` type, EC beam and NBI
  unit identity, EC mode) over a `dd_at_least` predicate, plus `IdentityHelper`, which
  `DistributionsBuilder` takes instead of the converter;
- `conversion_provenance.py` -- `ProvenanceRecorder` owns the input inventory, source-file records and
  `code.parameters` XML caches and builds `ids_properties` and `code`; `reload_namelist` resets it;
- `ids_assignment.py`, `cocos_signs.py`, `code_provenance.py`, `dd_adapter.py`, `imas_store.py` --
  extracted earlier;
- CDF access and species discovery live in `netcdf_class.py` (including `get_ncattr_text`) and the IMAS
  connection in `imas_class.py`; `nc2imas` still inherits both as its two I/O adapters.

What stays in `nc2imas` is state (the namelist dict, time bases, signs, caches' owner), lifecycle, the
delegators, and a few methods that orchestrate five or more of its own services:
`fill_profiles1d_grid`, `vacuum_field_reference`/`vacuum_b0_magnitude`, `fill_vac_tor_field`,
`check_profile_shape`, `fill_data_nc`. Moving those out would mean passing about ten callables into each,
which is the same coupling in another place, so they stay with the state they read. No writer module uses
`self.im` any more (`verify_conversion_context.py` enforces it); a capability a writer needs, such as
`supports_occurrences`, is a named member of the context. Tests of the extracted modules:
`verify_conversion_services.py` (pure) and, for end-to-end equivalence, `compare_conversions.py`.

**Proving a refactor changed nothing.** `examples/compare_conversions.py snapshot DIR` converts four
bundled cases (ITER, DIII-D 147634T16, NSTX 141711Y05, C-Mod 89846; 67 IDSs) and `compare A B`
reads two snapshots back and diffs every IDS node. Two runs of unchanged code differ only in
`ids_properties/creation_date`, which is the only field it ignores, so identical output before and
after a change is a strong check: the move to views was verified this way.

### Writer interface

`writer_base.IdsWriter` is the target interface: construction only stores the converter;
`build()` creates and fills the IDS objects and returns them as a list (or `(ids, occurrence)`
pairs); `run()` submits that list; `run_writer()` is what the runner calls. A build that fails
submits nothing, and what a writer will write is a value that can be inspected before any
submission. All 22 writers use it (`tf`, `wall`, `lh_antennas`, `ic_antennas`, `nbi`,
`ec_launchers`, `sawteeth`, `dataset_fair`, `workflow`, `summary`, `transport_solver_numerics`,
`pf_active`, `pf_passive`, `radiation`, `distribution_sources`, `edge_sources`, `waves`,
`core_sources`, `core_profiles`, `equilibrium`, `distributions`, `core_transport`), and
`verify_writer_interface.py` fails if `convert.py` imports a writer that does not. A writer
reads its grids at the start of `build()`, not in its constructor; the constructor of a writer
with options (`ec_launchers`, `pf_active`, `pf_passive`, `equilibrium`) only stores them.
`core_transport` returns occurrence 0 plus the region view occurrences as `(ids, occurrence)`
pairs, and `equilibrium.build()` times `_build()`. `run_writer` still ignores an object that is
not an `IdsWriter`, so an external writer that works in its constructor keeps working.
To migrate one: subclass `IdsWriter`, rename `create_ids` to `build`, return the IDS list
instead of calling `putdb_imas`, return `[]` where it used to write nothing, and update any test
that built the writer expecting a put.

### IDS writers

Each `imas_<ids>.py` file owns the mapping for its IDS. Writers should avoid
reimplementing unit/sign/grid logic that belongs in `nc2imas.py`.

Current writer modules include:

- machine/actuator IDSs: `wall`, `tf`, `pf_active`, `pf_passive`, `nbi`,
  `ec_launchers`, `ic_antennas`, `lh_antennas`;
- plasma-state IDSs: `core_profiles`, `equilibrium`, `summary`;
- source/transport IDSs: `core_sources`, `core_transport`, `radiation`,
  `edge_sources`, `waves`, `distribution_sources`;
- event/metadata IDSs: `sawteeth`, `workflow`, `dataset_fair`.

## Data Dictionary compatibility

The converter supports both the DD3 and DD4 families where practical. The code
must never assume that a field name or array structure is identical across Data
Dictionary families. Compatibility logic should be centralized when it is used
by more than one writer.

DD4 output follows COCOS 17. The authoritative sign rules and examples are in
`technical-reference.md`.

## Error-handling philosophy

A missing optional TRANSP quantity should normally cause the specific optional
mapping to be skipped, not abort the entire conversion. Conversely, an
ambiguous physical meaning, incompatible radial grid, or unvalidated unit/sign
conversion is a reason not to map the quantity until it is resolved.

The objective is broad coverage with defensible semantics, not maximum node
population at any cost.

## Adding a new writer

A new IDS writer should:

1. live in `imas_<ids_name>.py`;
2. have `#!/usr/bin/python` as line 1 and a module docstring immediately after;
3. subclass `IdsWriter` and use the context views for CDF access (`self.cdf`), units and
   grids (`self.grid`), signs (`self.signs`), DD compatibility (`self.dd`) and IDS filling
   (`self.fill`) rather than reaching into `self.im`; a name no view offers is added to the
   right list in `conversion_context.py`, which `verify_conversion_context.py` checks;
4. avoid constructing values that cannot be validated from source data;
5. be registered in the conversion sequence only after source-side checks are
   available;
6. add or update documentation and regression checks in the same change.


## Resource lifetime

`netcdf_connection` owns the input TRANSP CDF and implements `__enter__`/`__exit__` plus an idempotent `close_ncfile()`. `nc2imas` extends that lifetime rule to both the input CDF and the IMAS `DBEntry`: `convert.py` uses `with nc2imas(...) as nci`, and `nc2imas.close_imas_database()` closes the IMAS entry and the NetCDF dataset. `__del__` remains only a last-resort cleanup path and is not required for correct operation.

The `pf_active` writer delegates CDF decoding, optional machine JSON
validation, and DD-specific element/circuit population to
`pf_active_mapping.py` (source-only, no IMAS dependency). `KT_` is a target
reference and is not written to an actual supply-current field; see
[`pf-active.md`](pf-active.md).

## Where the time goes

`--report-json` records, per writer, `seconds` (total) and `write_seconds` (time in the backend: IMAS validation,
put and compression); the difference is reading, transforming and assigning. Nine of the twelve bundled cases take
1-2 s. Three (C-Mod 89846, NSTX 142202U01 and 120967E03, with many time slices or fast-ion data) take 25-100 s, and
there about 60% of the time is the backend: on C-Mod, 24 of the 31 s inside `put` is IMAS-Python's own
`validate`. That is not reduced by changes here, and the conversion does not validate a second time.

**What was measured, and what it showed.** cProfile and a per-variable read count (now part of the journal as
`reads`, see [mapping status](mapping-status.md)) on the three large cases:

- Reading the CDF is cheap: 7,200 reads in 1.2 s on C-Mod and 10,200 reads in 1.8 s on NSTX 142202U01. A cache of
  source arrays would save milliseconds, so none was added. The repeated reads that exist are one call site
  (`DAREA` 4,764 times on C-Mod) and per-slice reads of whole time series.
- The real redundancy was **loop-invariant work inside per-time-slice loops**: arrays that do not depend on the
  slice recomputed (and re-read) for every slice, which makes the cost grow with the square of the number of
  slices. Each was hoisted so it is computed once; the output is unchanged (`compare_conversions.py`).

| where | what was recomputed per slice | change | effect (own time) |
|---|---|---|---|
| `equilibrium.trace_boundary` | a full matplotlib Figure and Axes just to contour `PSIRZ` | `contour_lines.line_contours` calls the same `contourpy` generator without a Figure | 10.2 ms to 0.6 ms per slice (16x, identical vertices); NSTX 142202U01 writer 6.9 s to 2.1 s |
| `core_transport._fill_grid` | `DAREA` read and two cumulative sums over all times, for every model and slice | `_enclosed_geometry`, cached on the grid | C-Mod 2.3 s to 1.3 s (0.6 s measured alone) |
| `radiation` line radiation | enclosed power and radial interpolation of every species and charge state | computed once before the time loop | C-Mod 2.1 s to 0.4 s |
| `core_sources._fill_channel` | enclosed particles and power of each ion species | computed once before the time loop | small |
| `equilibrium` | `RMCB0`/`YMCB0` whole series | read once | small |

`import matplotlib.pyplot` was also moved from module level to the debug plot that uses it (about 0.5 s of import
time on every conversion). Absolute totals of two end-to-end runs are only comparable on the same machine
state: the backend share, which these changes do not touch, varied by 35% between two runs on a shared
workstation, so compare a writer's own time or use an A/B in one process as the table does.

Not done, on purpose: caching of source arrays (no cost to save, above), the backend's validation and write, and
the FI archives below (already parallel). To look for the next candidate, read `reads.most_repeated` in a journal:
a variable read hundreds of times by one writer sits in a loop.

Reuse implemented earlier: `radial_interpolation.interpolate_radial_profiles` caches the X/XB grid mapping by the
content of the coordinate arrays (bounded cache, original loop for non-monotone or NaN grids), which removed the
repeated per-slice `np.interp` loop that dominated the `radiation` writer. Equivalence test:
`examples/verify_radial_interpolation.py`; for the contour, `verify_equilibrium_boundary_trace.py` compares
`line_contours` with `Axes.contour` vertex for vertex on synthetic and real flux maps.

NUBEAM FI_TAR archives are the other large cost: reaching the equilibrium members means decompressing most of a
gzip stream of hundreds of MB. `nubeam_fi.read_archives_in_order` reads the archives of a run with up to four
`spawn` processes (threads gave 25%, processes 2x on two archives), keeps a bounded window in memory and yields in
the original order. `TRANSP2IMAS_FI_WORKERS=1` reads them sequentially. NSTX 120967E03 `distributions`: 40 s to 21 s.

## Independent scientific curation

`curate.py` and `curation/` screen source data and explicit evidence without constructing `nc2imas`,
importing IMAS or opening an output database. `NativeSource` provides copied read-only arrays,
explicit fill/missing handling and a bounded cache; `CurationEngine` emits typed findings and local
evaluated/flagged masks; the strict evidence parser rejects unknown contracts instead of guessing.
The first-stage sidecars are JSON, Markdown and NPZ. Scientific acceptance, full model-state/default
resolution, uncertainty ensembles and typed IMAS curation export are future stages.

This is not an alternative physics converter. Reuse a common source/provenance protocol after the
converter's missing-data and namelist contracts have been hardened and regression tested. Keep
source-scientific findings separate from failed mappings or incomplete backend serialization.
See [curation design](curation/design.md).
