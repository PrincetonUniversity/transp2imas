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
4. Individual `imas_*.py` modules populate one IDS family at a time.
5. Common provenance is attached before each IDS is written.
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

- `imas_class.py` -- backend selection, `DBEntry` lifecycle, IDS factory access,
  DD-version handling, and common put/write behavior.

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
3. use existing `nc2imas` helpers for CDF access, units, signs, provenance, and
   DD compatibility where possible;
4. avoid constructing values that cannot be validated from source data;
5. be registered in the conversion sequence only after source-side checks are
   available;
6. add or update documentation and regression checks in the same change.


## Resource lifetime

`netcdf_connection` owns the input TRANSP CDF and implements `__enter__`/`__exit__` plus an idempotent `close_ncfile()`. `nc2imas` extends that lifetime rule to both the input CDF and the IMAS `DBEntry`: `convert.py` uses `with nc2imas(...) as nci`, and `nc2imas.close_imas_database()` closes the IMAS entry and the NetCDF dataset. `__del__` remains only a last-resort cleanup path and is not required for correct operation.
