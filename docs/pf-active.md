# `pf_active`: active poloidal-field coils, supplies and circuits

The `pf_active` writer uses electrical waveforms saved by TRANSP and, **only
when supplied**, an authoritative machine description to populate IMAS
`pf_active.coil`, `pf_active.supply`, and `pf_active.circuit`.

The main goal is correctness rather than guessing a machine description from
CDF variable names. A CDF waveform normally **cannot** determine the coil's
physical cross-section, number of turns, winding direction, supply-to-coil
connectivity, or whether a target current equals the actual current. These
properties must come from a reviewed ISolver machine model, an existing machine
`pf_active` IDS, or another authoritative engineering source.

## Electrical CDF variable mapping

| TRANSP variable | Source meaning | IMAS mapping without machine JSON | With validated machine JSON |
| --- | --- | --- | --- |
| `KK_<coil>` | Actual/simulated coil circuit-current output | `coil[].current.data`, source sign unchanged | `coil[].current`; can additionally supply `supply[].current` and/or `circuit[].current` when the machine description explicitly connects them |
| `VV_<coil>` | Coil terminal-voltage output | `coil[].voltage.data` | `coil[].voltage`, with explicit voltage sign |
| `VS_<supply>` | Power-supply/source voltage | `supply[].voltage.data` | `supply[].voltage`, with explicit voltage sign |
| `KT_<circuit>` | **Target/reference current**, not actual measured current | **Not written to `pf_active`** | **Still not used as actual current**; intended future destination is `pulse_schedule.pf_active` current reference |
| `CC_<coil>` | Equilibrium coil **ampere-turns** (`CC_ ≈ N·KK_` in the ITER/NSTX cases) | **Never directly substituted** | Can be divided by authoritative `cc_effective_turns` and mapped to `coil[].current` with an explicit sign |
| `PF_<coil>` | Legacy auxiliary placeholder/machine-specific variable | Not written | Not written without a dedicated physical audit |

The `KK_/VV_/VS_` meanings and time grids have now been validated against two
real ISolver runs: ITER `38530P94` (11 coils) and NSTX `142202U01` (9 coils).
The earlier six generic NSTX/DIII-D/C-Mod cases lack these channels and remain
useful negative fixtures. See [real ISolver source validation](isolver-fixtures.md).
`KD_*` (desired current) and `KS_*` (JIT sigma) are also recognized and
explicitly deferred; they are not actual coil/supply measurements.

**Important correction:** the old `imas_pf_active.py` wrote `KT_*` into the
`supply.voltage` array and then overwrote it with `VS_*`. Simply moving `KT_*`
into `supply.current` would also be incorrect, because `KT_*` is a target
rather than measured current. The new converter recognizes `KT_*`, logs the
omission, and leaves the actual current unset until a genuine actual-current
channel has been matched. It never turns a target into a measured trace.

## CDF-only mode (no machine model)

The writer scans all `KK_*` and `VV_*` names and **takes their union**, so it
retains a coil having only a current or only a voltage. It independently
collects `VS_*` supply channels. Different coils/supplies are not assumed to
share a suffix. No artificial one-supply-per-coil mapping is constructed.

Each populated dynamic node gets its own native CDF time array:

- `coil[i].current.time`, `coil[i].voltage.time`;
- `supply[i].current.time`, `supply[i].voltage.time`;
- `circuit[i].current.time`, `circuit[i].voltage.time` (when provided).

`pf_active.time` is the generic reference axis. When all written dynamic traces
share the same time coordinate, `homogeneous_time=1`; otherwise it is `0` and
each node retains its exact independent time. The writer never silently
interpolates `TIME` into `TIME3` or vice versa.

CDF-only mode writes no guessed coil geometry, turns, resistance, or circuit
connections. Winding polarity is **not physically verified**: the source
values are preserved as recorded, and the IDS comment and log make that
limitation explicit. Such an IDS is not suitable for magnetic-field
reconstruction until a device model has been added.

If no `KK_/VV_/VS_` waveforms are present and no authoritative geometry is
specified, `pf_active` is **not written**. The old `pf_active.time = [-999]`
sentinel has been removed.

## Authoritative machine JSON (`--pf-active-geometry`)

Only use geometry and topology from verified device engineering files. The
repository includes
[`examples/pf_active_geometry.synthetic.json`](../examples/pf_active_geometry.synthetic.json)
for regression tests. **All values and circuit connections in this sample are
invented and must not be used for a real tokamak.**

A minimal authored manifest has this form:

```json
{
  "schema_version": 1,
  "machine": "REAL_DEVICE_NAME",
  "source_reference": "Engineering device definition and revision/date",
  "units": "m",
  "coils": [
    {
      "name": "REAL_COIL_ID",
      "description": "Machine coil description",
      "cdf_current": "KK_REAL_ID",
      "current_sign": 1,
      "cdf_voltage": "VV_REAL_ID",
      "voltage_sign": 1,
      "resistance_ohm": 0.001,
      "elements": [
        {
          "name": "REAL_ELEMENT_ID",
          "turns_with_sign": 1,
          "geometry": {
            "type": "rectangle", "r": 1.0, "z": 0.0,
            "width": 0.1, "height": 0.1
          }
        }
      ]
    }
  ],
  "supplies": [],
  "circuits": []
}
```

**The numbers and identifiers above are illustrative, not measured machine
data.** Replace them with approved source values; the source file must not be
committed or used as a physical device geometry without validation.

The JSON is deliberately strict:

- `machine` must match the CDF device, and `source_reference` must identify the
  origin and revision of the geometry/circuit model.
- `units` is `m`; geometry angles are in radians, resistances in ohms.
- Every element needs a nonzero, **signed** `turns_with_sign`. If multiple
  physical sub-elements have opposite winding sense, signs may differ.
- Each bound waveform needs an explicit `current_sign` or `voltage_sign` of
  `+1` or `-1`. **This is not the plasma `Ip`/`Bphi`/COCOS sign.**
- A coil current may refer to `KK_` or, after case-specific confirmation, `CC_`;
  a coil voltage must refer to `VV_`. Supply voltage uses `VS_`. Supply
  current may be explicitly tied to a physically equivalent `KK_`/`CC_` trace,
  with its independently verified supply sign.
- A `KT_` target is rejected as an actual-current source.
- Every present `KK_/VV_/VS_` source must be accounted for by a coil/supply
  entry (not just a circuit). Unmapped channels or duplicate coil/supply
  binding names are rejected. Coils in the manifest with no matching waveform
  are retained as valid **geometry-only** coils.
- The six geometry shapes supported by `pf_passive` are reused: `rectangle`,
  `outline`, `oblique`, `arcs_of_circle`, `annulus`, `thick_line`. Element
  outlines are open under DD3 and explicitly closed under DD4.

### Circuit connectivity

Circuits are entirely optional. **Do not infer one circuit per coil.** Only
include a circuit when its actual node connectivity is known. Each circuit
accepts a name, description and one or both of the following:

- `connections_dd3`: integer node-by-side matrix. Its number of columns is
  `2 * (N_supplies + N_coils)`; entries must be `0` or `1`.
- `connections_dd4`: integer node-by-component matrix. Its number of columns
  is `N_supplies + N_coils`; entries must be `-1`, `0`, or `1`.

These are *different* DD semantics, not different encodings of the same
matrix. The converter does **not** infer one from the other. A circuit is
written only in the DD family for which the manifest explicitly supplies a
connection matrix. The source may optionally specify `cdf_current`,
`current_sign`, `cdf_voltage`, and `voltage_sign` when there is an explicitly
verified one-to-one correspondence between that waveform and the circuit.
The resulting `circuit[i].current/voltage` nodes receive their own time axes.

Example: see the *synthetic* JSON test fixture. The converter validates matrix
shape and entry types but cannot independently validate electromagnetic
orientation or circuit topology without an actual machine source.

### Provenance and validation

In manifest mode, the converter appends the source file name, SHA256 and
`source_reference` to `ids_properties.comment`, and embeds the manifest text
inside `code.parameters`. In CDF-only mode it documents explicitly that the
geometry and sign choices were not verified.

```bash
# CDF-only conversion
python3 convert.py -i /path/to/RUN.CDF -b HDF5 -new -u /path/to/imasdb

# With approved machine geometry and circuits
python3 convert.py -i /path/to/RUN.CDF -b HDF5 -new -u /path/to/imasdb \
    --pf-active-geometry /path/to/machine-pf-active.json

# Equivalent Makefile invocation
make convert CDF=/path/to/RUN.CDF \
    PF_ACTIVE_GEOMETRY=/path/to/machine-pf-active.json

# Source-only synthetic test; no IMAS needed
make verify-pf-active

# Separate, actual schema/end-to-end tests on an IMAS-Python installation
make verify-dd
```

Portable tests use a synthetic CDF with `KK_/VV_/KT_/VS_/CC_`, varied grids,
SI unit prefixes, DD3/DD4 element outlines and distinct circuit matrices.
A fake IDS writer checks that no target current is written as actual current,
that all available dynamic values retain their true times, and that the
writer does not insert dummy times. The existing six bundled real CDFs are
also checked for absence of these output families.

**Outstanding:** the attached ITER/NSTX CDFs provide authentic
`KK_/VV_/VS_` waveforms. A matching independently verified device model
is still needed for physical coil geometry and wiring. Check unit
metadata, `KK_` and `CC_` sign conventions, winding orientation, actual supply
connections, and compare converted currents/voltages with independent device
references. Run `make verify-dd` and an end-to-end IMAS write on a supported
site installation before treating PF active geometry/current as physically
validated.

Relevant schema references:

- [IMAS DD4.1.1 `pf_active`](https://imas-data-dictionary.readthedocs.io/en/4.1.1/generated/ids/pf_active.html)
- [IMAS DD3.41 `pf_active`](https://imas-data-dictionary.readthedocs.io/en/3.41.0-doc/generated/ids/pf_active.html)
- [IMAS `pulse_schedule` (planned KT_ references)](https://imas-data-dictionary.readthedocs.io/en/4.1.1/generated/ids/pulse_schedule.html)

## `CC_` ampere-turns: physical-current conversion

In real ISolver CDFs `CC_<coil>` has an equilibrium coil-current label and
`A` units, but differs substantially from the circuit-current `KK_<coil>`:
`CC_CS3U/KK_CS3U = 548` in ITER and `CC_OH/KK_OH = 963` in NSTX. The
factor is constant to CDF float precision across time, consistent with an
effective series-turn count. `CC_` is **not** a drop-in circuit-current trace.

Use `KK_` whenever available. If an authoritative device model specifically
requires the `CC_` representation, the geometry manifest must give the
signed/nonzero effective turns explicitly, for example:

```json
"cdf_current": "CC_REAL_COIL",
"current_sign": 1,
"cc_effective_turns": 548
```

The IDS current is then `current_sign * CC_ / cc_effective_turns`. The
converter rejects a `CC_` coil current without `cc_effective_turns` or any
attempt to use `CC_` as a supply/circuit current. The example number is
observational, **not** a device-model endorsement. Obtain the effective-turns
value and winding convention from the reviewed ISolver hardware description.

For full source verification, run `make verify-isolver` and inspect the active
and passive channel inventories in `build/pf-isolver/`. The ITER/NSTX files
verify currents, voltage time series, `CC_/KK_` closure and CDF channel labels,
but **do not provide the authoritative cross-sections or circuit connectivity**.
