# Converter input contracts

What the converter requires of its inputs and what it does when they are not met. The scientific side of the same
policies is in [curation/scientific-policy.md](curation/scientific-policy.md); the per-variable mapping is in
[mapping-registry.md](mapping-registry.md).

## Namelist values

`fortran_literals.py` is the one scalar-literal parser, shared with curation. It reads integers, reals (with `D` or
`E` exponents), Fortran logicals and quoted strings. There is no `eval`.

- A value that is not one supported scalar (an array, an expression, `nan`, `inf`, a number with trailing text, an
  empty value) is **unknown**: `get_var_` returns `None` and logs a warning. It is never `0`, `False`, or an earlier
  assignment.
- `parse_number` raises `ValueError` for such text; callers that need a number decide what unknown means.
- This is not a Fortran parser: defaults, includes and runtime updates are not resolved.

## Profile shape and coordinates

A present profile must have exactly `(nt, nx)`, with a leading `TIME3` axis and a radial axis from the declared
`X`/`XB` family. `check_profile_shape` rejects a wrong time length, another dimensionality, or a radial dimension that
does not match the grid, and the error names the variable and what it found. The check uses the CDF's declared
dimensions, never array length alone, so equal lengths on `TIME` and `TIME3` are not mistaken for one grid.

## Missing and masked data

netCDF4 masks `_FillValue`, `missing_value`, `valid_range` and the default fill. A masked sample that reaches a
writer is a broken mapping (`MappingError` in strict mode); in lenient mode the variable is skipped and the reason is
recorded in the mapping journal. An absent variable and a variable that is present and zero are different, and are
never merged.

## NBI units and species

- Every unit and its energy trace are written on the native `TIME` axis unchanged. An off beam is expressed by its
  zero power, not by the clock. A trace whose length differs from `TIME`, or a `TIME` that is not finite and strictly
  increasing, is a broken mapping.
- A unit keeps its `PINJ<n>` number and is ordered by that number, not by discovery order.
- A species is attached only when exactly one `PINJ_<sp>` total exists and the injector powers add up to it (relative
  tolerance `SPECIES_CLOSURE_RTOL = 1e-4`). Several species, or a total the injectors do not close on, leave the
  species unset and are reported. Closure is evidence of consistency, not a proof of the species.
- A species total with no injector (`PINJ_D` alone) is an aggregate over unknown injectors; no unit is invented.

## NBI torque

The `core_sources` toroidal-momentum total is the sum of `TQBE`, `TQBI`, `TQBTH`, `TQJXB` and, when present,
`TQABDF`, and is written only when it is complete. A missing required term withholds the total (a zero term does not),
and the `TQJXB` component is still written when it exists. A zero or absent `TQRPL` adds nothing. A nonzero `TQRPL`
contributes `FTQRIP*TQRPL` if `FTQRIP` is a finite scalar or a `TIME3` history, otherwise the total is withheld and
the XML record says `status="partial"`, `total_written="false"`. The bundled runs have a zero or absent `TQRPL` and no
`FTQRIP`, so the weighted path is covered only by synthetic tests.

## Publication

IMAS writes are sequential, not transactional. See the publication section of the scientific policy for the report
and the completion marker.
