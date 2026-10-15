# TRANSP sawteeth mapping

The new `sawteeth` writer uses **native IMAS typed fields** and the **CDF `TIME`**
coordinate, which can be different from the profile coordinate `TIME3`.
`core_profiles` and `equilibrium` continue to use `TIME3` independently.

## Observed CDF quantities and semantics

| TRANSP CDF | Dimension | Meaning | IMAS DD3.39/DD4.1.1 |
|---|---|---|---|
| `TIME` | `TIME` | Scalar-output observation times, s | `sawteeth.time` |
| `TLASTSAW` | `TIME` | Time of most recently reported crash, s | `sawteeth.diagnostics.previous_crash_time` |
| Change in `TLASTSAW` | `TIME` | A **new reported** crash since the previous observation | Identifiable from the typed previous-crash history; not assigned a trigger mechanism |
| Intervals between successive **observed new** crashes | `TIME` | Last *verified* sawtooth period, s | `sawteeth.diagnostics.previous_period` |
| `XIQ1` | `TIME` | Normalized `xi` radius of `q=1` surface | `sawteeth.diagnostics.rho_tor_norm_q1` |
| `XSAWBOUND` | `TIME` | TRANSP "Adhoc sawtooth boundary" | **Not yet mapped** to a mixing or inversion radius |
| `TESAW`, `DNSAW`, `DUSAW`, `XISAW` | scalar | Model settings/auxiliary parameters, not verified post-crash profiles | Preserved via source/input provenance where available; no invented profile |

`IMAS_MISSING_FLOAT = -9e40` is used for invalid/absent `XIQ1` samples
(the source often stores `-0.1` or `-0.2`) and for the undetermined early
`sawteeth.diagnostics.previous_period` samples. These are **missing values**,
not physically zero radii or periods. A period is reported only after two
**observed new** changes to `TLASTSAW`. The initial `TLASTSAW` value could
represent a crash preceding the CDF interval or an initialization value;
we do **not** count it as a new event or use it to infer the first period.

The IDS uses `ids_properties.homogeneous_time = 1`, per the IMAS Data
Dictionary. No separate event-time arrays or duplicate radial coordinates
are embedded in model XML. Existing `code.parameters` continues to preserve
original TRANSP inputs.

**Intentionally not filled:** `crash_trigger` and
`diagnostics.previous_crash_trigger` are *mechanism codes* (ideal kink,
resistive kink, etc.), not Boolean "crash happened" flags. Neither the
available `TLASTSAW` history nor `NLSAW=.T.` identifies the mechanism.
Similarly, a valid `XSAWBOUND` control is not yet proven equal to the
physical mixing or inversion radius, and ordinary `TIME3` plasma profiles
cannot be labeled "post-crash" without precise crash-aligned source data.

## Real examples

C-Mod `89846.CDF` has 397 `TIME` samples and **38 changes** in `TLASTSAW`
beyond its initial baseline; 184 samples contain positive, valid `XIQ1`.
`XSAWBOUND` is zero throughout. Its `89846TR.DAT` enables sawtooth modeling
(`NLSAW=.T.`), but does not by itself identify a crash-trigger mechanism.

ITER `38530P94.CDF` has a constant `TLASTSAW=550 s` for all five output
times; no new crash can be asserted, so the `sawteeth` writer skips it.
Other supplied cases have no `TLASTSAW` and are also skipped even when they
contain `XIQ1` or the scalar `TESAW` placeholder.

## Run the converter

```bash
python3 convert.py -i examples/data/cmod/89846.CDF -b HDF5 -new \
    --dd-version 4.1.1 --ids sawteeth --report-json build/cmod-sawteeth.json
```

Source regression for **every** supplied CDF (including C-Mod, ITER and
all NSTX/DIII-D cases):

```bash
make check
make verify-sawteeth
make verify-all
```

`make verify-sawteeth` creates `build/sawteeth/sawteeth.json` and
`build/sawteeth/sawteeth.csv`. The CSV separately lists the event timestamps taken from `TLASTSAW`
and the subsequent CDF `TIME` sample at which each event was observed.
To test a real IMAS write and coordinate readback:

```bash
python3 examples/verify_end_to_end.py --case 89846 --ids sawteeth \
    --dd-version 4.1.1 --report-dir build/e2e-sawteeth
```

The source regression and mock IDS write/readback do **not** replace a real
DD4.1.1 or DD3.39 write/readback on a machine with IMAS-Python installed.
