# `pf_passive`: TRANSP CDF currents and authoritative device geometry

`pf_passive` represents **axisymmetric, unconnected passive-conductor loops** and
their electrical currents. It is not a replacement for `wall`, and externally
powered/connected coil circuits belong in `pf_active`. A TRANSP CDF may contain
passive-structure current outputs but **does not** reliably contain the machine
geometry, resistivity, effective turns, or a mapping from each current output to
a specific piece of hardware.

The 0.9.0 converter therefore supports two modes:

1. **Current-only legacy mode** (no new command-line options): when `PS_1`,
   `PS_2`, etc. exist, it writes their currents and labels as incomplete passive
   loops, with a high-visibility warning. No geometry is invented. An IDS with
   only currents is **not suitable for electromagnetic field reconstruction**.
2. **Authoritative-geometry mode:** pass a separate validated machine JSON
   description with `--pf-passive-geometry PATH`. The converter writes the
   specified conductor elements/geometry/resistance and attaches each CDF
   current *only* to its explicitly mapped passive loop. Geometry-only loops
   (no associated `PS_n` in that run) are also written without currents.

Neither mode assumes `PS_n` sign is the same as `Ip` or `Bphi`; the sign of the
physical loop current is independent of the COCOS signs. When joining a current
to geometry, the manifest **requires** `cdf_current_sign` to be `1` or `-1`
according to the device's established loop-current orientation. The element
`turns_with_sign` is a separate physical orientation and is also required.

## Obtain the machine geometry

Use an existing verified device `pf_passive` IDS, ISolver/equilibrium-machine
configuration, or an approved machine engineering description. The **identity
between `PS_n` and device loop names** must be confirmed from TRANSP machine
configuration, not guessed from text labels such as `<BVFB>`. Translate that
information into the JSON manifest described below. All dimensional inputs
use **metres**, angles **radians**, resistance **ohms**, resistivity **ohm.metres**.

This repository contains only
[`examples/pf_passive_geometry.synthetic.json`](../examples/pf_passive_geometry.synthetic.json),
a **deliberately invented synthetic fixture** for automated tests. DO NOT use
that file for a tokamak conversion and DO NOT interpret its coordinates or
resistance values as machine data.

## JSON manifest format

```json
{
  "schema_version": 1,
  "machine": "YOUR_MACHINE_NAME",
  "source_reference": "Authoritative machine-configuration record / revision and date",
  "units": "m",
  "loops": [
    {
      "name": "DEVICE_PASSIVE_LOOP_1",
      "description": "Machine-verified description",
      "cdf_variable": "PS_1",
      "cdf_current_sign": 1,
      "resistance_ohm": 0.0001,
      "elements": [
        {
          "name": "DEVICE_PASSIVE_ELEMENT_1",
          "turns_with_sign": 1,
          "geometry": {
            "type": "rectangle",
            "r": 1.2,
            "z": 0.1,
            "width": 0.1,
            "height": 0.1
          }
        }
      ]
    }
  ]
}
```

**The sample values above are arbitrary, for schema illustration only.**
Replace every loop identifier, sign, coordinate and electrical property with
validated machine values. The `machine` field must match the CDF device name
(case-insensitively). `source_reference` is required and identifies the
provenance of the machine geometry.

Loop fields:

| Field | Meaning | Requirement |
| --- | --- | --- |
| `name` | Unique machine passive-loop identifier | Required |
| `description` | Optional machine description | Optional |
| `cdf_variable` | Matching TRANSP output (`PS_1`, `PS_2`, etc.) | Optional; absent means geometry-only |
| `cdf_current_sign` | +1 or -1 applied to CDF current before writing | Required when `cdf_variable` present |
| `resistance_ohm` | Loop resistance (ohms) | Optional, nonnegative |
| `resistivity_ohm_m` | Material resistivity (ohm.metres) | Optional, nonnegative |
| `elements` | One or more physical cross-section elements | Required |

Each element requires `name`, a **signed and nonzero** `turns_with_sign`, and a
`geometry` object. Optional element fields are `description`, `identifier` and
`area_m2`. Zero or missing geometry fields are **not** filled by guesswork.

Supported IMAS geometry shapes (choose one `geometry.type` per element):

- `rectangle`: `r`, `z` (centre); positive `width`, `height`.
- `outline`: arrays `r`, `z` describing a non-degenerate polygon. Supply
  vertices once (either open or explicitly closed); the converter removes the
  duplicate vertex when validating and rewrites **open DD3**, **closed DD4**.
- `oblique`: `r`, `z`, `length_alpha`, `length_beta`, `alpha`, `beta`.
- `arcs_of_circle`: aligned arrays `r`, `z`, `curvature_radii`.
- `annulus`: `r`, `z`, `radius_inner`, `radius_outer`.
- `thick_line`: `first_point:{r,z}`, `second_point:{r,z}`, `thickness`.

The first point of an outline must be repeated at the end **only in DD4**.
Both cases are handled automatically based on the output Data Dictionary.

When a geometry manifest is used, every **present** CDF `PS_n` channel must be
listed exactly once; otherwise the converter raises a `ValueError` instead of
silently omitting/misidentifying a current. Conversely, a manifest may include
additional machine loops whose channels were not saved in this CDF; those
loops remain geometry-only. This permits a full device description even when
no passive currents are recorded.

## Run and verify

```bash
# Normal conversion: currents only if PS_n is present
python3 convert.py -i /path/to/RUN.CDF -b HDF5 -new -u /path/to/output

# With physical conductor data (after you prepare a real machine manifest)
python3 convert.py -i /path/to/RUN.CDF -b HDF5 -new -u /path/to/output \
    --pf-passive-geometry /path/to/machine-pf-passive.json

# Equivalent Makefile override
make convert CDF=/path/to/RUN.CDF PF_PASSIVE_GEOMETRY=/path/to/machine-pf-passive.json

# Synthetic + real-case negative tests; no IMAS required
make verify-pf-passive

# Standard converter tests
make check
make verify-all
```

Check `log.out` for the number of loops written, current traces and populated
geometries. For a genuine passive-current regression, independently compare
`loop[i].current` against the corresponding CDF `PS_n` after applying only the
manifest's `cdf_current_sign` and the documented A/kA unit conversion. Do not
apply `Ip` direction or COCOS sign to a passive circuit current.

### Time coordinates and units

`PS_n` variables must be finite **1D** curves on `TIME` or `TIME3` (or an
unambiguously length-matched legacy time dimension). The writer uses the
native vector in **`loop[i].time`** and writes `pf_passive.time` only as a
common/global reference. If loop times differ, the IDS is marked
`homogeneous_time=0`; otherwise it is `1`.

Accepted unit attributes are `A` (amperes), `kA`, `MA` (megaampere), and `mA` (milliampere); SI prefix case matters. When TRANSP omits
units, the converter preserves the historical assumption that `PS_n` currents
are in amperes; obtaining a real example with a unit attribute will be part
of final physical validation.

### Provenance

`ids_properties.comment` records the manifest name, SHA256 hash and machine
source reference; `code.parameters` also embeds the entire input JSON in XML.
Geometry cannot be silently synthesized from a CDF. Invalid machine names,
nonfinite coordinates, invalid polygons, duplicate mappings, absent turns
or unspecified `cdf_current_sign` cause explicit errors.

## Validation status and remaining work

The older six NSTX/DIII-D/C-Mod CDFs contain **no `PS_n` channels**. We now
also have two real ISolver fixtures: ITER `38530P94` saves 119 passive-current
channels and NSTX `142202U01` saves 82. `make verify-isolver` tests every
source current/time pair, the order of PS indices, and the label extracted
from each CDF `long_name`; the source-side checks pass. Synthetic tests still
verify sign handling, geometry validity and DD3/DD4 polygon conventions.
These tests do **not** establish actual loop geometry, winding orientations,
resistivity or electrical connectivity. We still need the corresponding
verified ISolver passive-structure model for physical magnetic-field modeling. Then add that source case to the physical regression suite and run
`make verify-dd` and an end-to-end conversion under IMAS-Python.

Data Dictionary references:

- [DD4.1.1 pf_passive](https://imas-data-dictionary.readthedocs.io/en/4.1.1/generated/ids/pf_passive.html)
- [DD3 pf_passive](https://imas-data-dictionary.readthedocs.io/en/3.41.0-doc/generated/ids/pf_passive.html)

See [ITER/NSTX ISolver source cases](isolver-fixtures.md) for the new file layout and complete PS_n inventories.
