# TRANSP CDF TIME and TIME3 coordinates

TRANSP CDF uses **two distinct time coordinates**.  Do not use the length of an
array to infer which coordinate applies: read the NetCDF variable's declared
`dimensions` and the actual values of `TIME` and `TIME3`.

| Native coordinate | Example variables in DIII-D 147634T13 | Typical IMAS destination |
|---|---|---|
| `TIME` | `PCURC`, `RAXIS`, `YAXIS`, `PSI0_TR`, `PVOL`, `TE0`, `TI0`, `Q95`, `PLFLXA`, `RMCB*`, `RMSB*`, `YMCB*`, `YMSB*` | Scalar observables, or interpolated onto a target IDS time grid |
| `TIME3` | `X`, `XB`, `TE`, `TI`, `DVOL`, `PLFLX2PI`, `Q`, `RMC00`, `ELONG`, `TRIANG`, `PSIRZ`, `B_FIELD`, `RGRID`, `ZGRID` | Profile quantities and `equilibrium.time_slice` |

In **147634T13**, `TIME` contains **4 samples**, while `TIME3` contains **7**:

```text
TIME  = [2.543573, 2.564947, 2.585168, 2.600000]
TIME3 = [2.543573, 2.553663, 2.564947, 2.574295,
         2.585168, 2.594486, 2.600000]
```

For `equilibrium`, the native 1-D/2-D profiles and magnetic fields are on
`TIME3`.  Therefore `equilibrium.time` and the seven `time_slice` entries both
use `TIME3`, regardless of the number of scalar `TIME` points.

To fill per-slice global quantities and the LCFS boundary consistently:

1. Read a `TIME`-indexed scalar (e.g., `RAXIS`, `PCURC`, `Q95`) and interpolate
   it onto `TIME3` when placing it inside `equilibrium.time_slice[i]`.
2. Read the `TIME`-indexed Fourier LCFS coefficients (`RMCB*`, `RMSB*`,
   `YMCB*`, `YMSB*`), interpolate each coefficient onto `TIME3`, then
   evaluate the boundary outline for **every** profile-time slice.
3. Read `TIME3`-indexed boundary shaping (`ELONG`, `TRIANG`, `TRIANGU`,
   `TRIANGL`) directly at the last radial boundary, without interpolating it
   backwards onto `TIME`.
4. Retain each variable's actual CDF time-coordinate values.  In particular,
   equal-length `TIME` and `TIME3` arrays may still have different timestamps.
   The shared `nc2imas.get_ncvar_on_time()` helper uses the variable's declared
   NetCDF dimension and will not confuse equal-length grids.

The underlying scalar measurements still retain their original `TIME` values
in the source CDF and input provenance.  Mapping them to equilibrium `TIME3`
only resamples their values for a time-consistent IDS.  `np.interp` uses
endpoint values if target times lie just outside the source range; source
coverage should be reviewed if this is substantial.

## Verification

```bash
make verify-equilibrium
make verify-core-profiles
make verify-all
```

`verify_equilibrium_timegrids.py` checks **all discovered CDFs**, performs
synthetic tests with **equal-length but different-valued time arrays**, and
executes the converter's actual boundary-moment conversion on each CDF using
a source-backed adapter.  These are **source-side checks**, not an actual
IMAS-Python `equilibrium.validate()` and database write.

To reproduce the DD 4.1.1 write test in an IMAS-enabled environment:

```bash
python3 convert.py -i examples/data/d3d/147634T13/147634T13.CDF \
  -b HDF5 -new --dd-version 4.1.1
```

Do not bypass IMAS validation with `IMAS_AL_DISABLE_VALIDATE=1`: that conceals
coordinate errors rather than fixing them.
