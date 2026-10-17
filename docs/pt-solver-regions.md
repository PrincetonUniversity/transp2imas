# PT_SOLVER multi-region conversion and provenance

## Why axial / confinement / edge are **not separate IMAS occurrences**

A TRANSP PT_SOLVER execution advances **one coupled, assembled plasma-profile
solution**. `pt_axial`, `pt_confinement`, and `pt_edge` are radial *domains in
which different reduced/anomalous and neoclassical transport models are
selected*. They are not distinct independent simulations or competing
alternative solutions. The converter therefore writes **one**
`core_transport` occurrence for the CDF (with the usual output run/occurrence
identification), not three overlapping occurrences or a synthetic
`alternative_solution` hierarchy.

- `core_transport.model[combined]` contains observed/effective transport and
  the appropriate fluxes reconstructed from observed balances.
- `core_transport.model[transport_solver]`, when available, contains the
  **assembled PT_SOLVER modeled fluxes** from TRANSP's `*_MOD` divergences.
- `core_transport.model[glf23]`, `model[mmm07_total]`, individual MMM reduced
  components, and the neoclassical models retain their **actual CDF
  coefficient/flux data** as separate model contributions. The MMM total must
  not be summed again with its components.
- The model arrays and their real `TIME3`, `X`, and `XB` coordinates remain in
  **typed native IMAS grids and output fields**, exactly as before. **No CDF
  coefficient profile is duplicated, artificially clipped, set to zero
  outside a selected region, or relabeled as a region-specific calculation.**

Distinct occurrences are appropriate only if a source actually provides
**independent complete solutions** (e.g., a separately evaluated input or
alternative run). A regional PT_SOLVER option alone does not justify that.
Similarly, there is no evidence in these CDFs of an independent
`alternative_solution` for each region.

## Region coordinate and equation-specific outer boundaries

The PT_SOLVER radial coordinate is

\[\xi=\sqrt{\Phi_{\mathrm{tor}}/\Phi_{\mathrm{tor,bdy}}}.\]

This is the radial convention underlying `rho_tor_norm` stored on the normal
`core_transport` grids; **the complete `X` and `XB` coordinate arrays are
already stored in typed fields**. `XIMIN_CONF` etc. are *input boundary
parameters*, not additional copies of radial grid coordinates.

For each predicted equation, let `B` be its outer predictive boundary:

| Equation(s) | TR.DAT boundary | Default |
|---|---|---:|
| Electron and ion temperature, `Te`, `Ti` | `XIBOUND` | 1.0 |
| Particle/density transport | `XNBOUND` | 1.0 |
| Toroidal angular momentum `Pphi` | `XPHIBOUND` | 1.0 |

The radial domains for a given `B` are

| Region | Intended ξ range | Notes |
|---|---|---|
| `pt_axial` | `[0, min(XIMIN_CONF, B)]` | May have zero extent |
| `pt_confinement` | `[XIMIN_CONF, min(XIMAX_CONF, B)]` | Empty if the upper limit is not greater than the lower |
| `pt_edge` | `[XIMAX_CONF, B]` | May have zero extent, or disappear when `B <= XIMAX_CONF` |

`XIMIN_CONF` defaults to 0.0 and `XIMAX_CONF` to 1.0 when an authoritative
`TR.DAT` is supplied but does not override them. **Defaults are not applied
when the original TR.DAT is absent**, because the executed configuration is
unknown. Input changes at `~update_time` are preserved as separate boundary
and equation-selection states; region limits are recomputed for each state.

These regions end at the *predicted equation's* boundary. They do **not**
characterize the pedestal/boundary-condition region beyond that boundary.
Furthermore, `LPREDICT_TE`, `LPREDICT_TI`, `LPREDICT_NE`, `LPREDICT_PPHI`, etc.
determine which equations are actually advanced. An active GLF23 coefficient
selection for all three regions **does not prove all equations were predicted
for the entire run**: `LPREDICT_*` can change during a run.

## Where we save the settings in IMAS

The XML below holds **model/solver inputs** for which DD 4.1.1 provides no
dedicated typed regional-domain/activation structure. It does *not* serialize
`TIME3`, `X`, `XB`, or any coordinate arrays in lieu of IMAS grids.

| IMAS location | Content |
|---|---|
| `core_transport.time` | Native `TIME3` array (seconds), typed IMAS field |
| `core_transport.model[].profiles_1d[it].time` | Native time for each transport model |
| `...grid_d.rho_tor_norm` | Native `XB` coefficient grid |
| `...grid_flux.rho_tor_norm` | Appropriate existing coefficient/flux grid; for MMM reduced components, native `X` used as a documented secondary coordinate |
| `core_transport.code.parameters / pt_solver_regions` | `XIMIN_CONF`, `XIMAX_CONF`, `XIBOUND`, `XNBOUND`, `XPHIBOUND`, derived per-equation regional limits, `LPREDICT_*`, change times and PT_TEMPLATE model selections |
| `core_transport.model[transport_solver].code.parameters` | Reference to the assembled-solution configuration, where modeled fluxes exist |
| `core_transport.model[glf23] / model[mmm07_*] / model[neoclassical].code.parameters` | **Explicit** region-use selections attached only when a referenced PT_TEMPLATE supports them; reduced-model outputs remain separate |
| `workflow` and `core_transport.code.library[].parameters` | Original `TR.DAT` and discoverable external namelists, unmodified original text, provenance and SHA-256 |

Model selections are listed per `&pt_controls` block and `pt_update_time`,
including **explicit** `%active`, `%xanom` (six transport-channel factors
`Te,Ti,Pphi,Ni,Nz,Ne`), and the `%exb%version` when supplied. An omitted model
switch in a subsequent update is reported as **unspecified**, rather than
silently assumed active or inactive. The complete original namelist remains
available for interpreting more complicated solver-update semantics.

The code is deliberately conservative about *provenance*: a referenced
`PT_TEMPLATE` constitutes evidence of the specified model selection; an
unreferenced adjacent `.nml` file is archived as a **candidate** only. A
nonzero GLF23 or MMM coefficient in a CDF establishes an actual output but
not the regional model selection used to obtain it.

### ITER `38530P94`

`38530P94TR.DAT` explicitly specifies:

```fortran
XIMIN_CONF=0.1
XIMAX_CONF=0.75
XIBOUND=1.00
XNBOUND=1.00
XPHIBOUND=1.00
PT_TEMPLATE='p94_pt.nml'
```

`p94_pt.nml` explicitly activates `GLF23` and `NEOCH` for **all three**
regions. At the start of the run, the temperature prediction switches in
`TR.DAT` are zero; a later `~update_time = 552.0` changes
`LPREDICT_TE`, `LPREDICT_TI`, and `LPREDICT_PPHI` to 1. The converter keeps
both states (along with the intervening RF-only update). The initial radial
ranges, for all three equation categories, are axial `[0,0.1]`, confinement
`[0.1,0.75]`, and edge `[0.75,1.0]`.

The `GLF23` model entry is **not multiplied three times**. Its source CDF
coefficients remain a single complete output array; metadata identifies
which PT_SOLVER regions selected that model.

### NSTX `141711Y05`

The authoritative `TR.DAT` sets `XIMIN_CONF=0`, `XIMAX_CONF=1`, and
`XIBOUND=1`, and references `y05_pt.nml`. Only `pt_confinement` has nonzero
radial extent (the axial and edge regions are empty). The template explicitly
selects CDBM and NEOCH in confinement, but the present coefficient-family
mapping does not yet claim to convert CDBM's full model-resolved outputs.

### DIII-D `147634T16`

The supplied source contains `t16_pt.nml` (MMM and NEOCH enabled in the
confinement-region text), **but is missing `147634T16TR.DAT`**. The converter
cannot prove the candidate template was used or determine the true radial
bounds and equation choices. It saves the candidate configurations **with an
unverified association** and `unavailable_without_TR.DAT` for boundary
history. It does not pretend the 0..1 defaults are confirmed for T16, and
does not annotate an MMM model as definitively confinement-only.

## Multiple model updates and evolving radial boundaries

- Multiple `&pt_controls` blocks are stored separately with their
  `pt_update_time` selection. Only explicitly stated switches are asserted.
- Changes to `XIMIN_CONF`, `XIMAX_CONF`, X*BOUND and `LPREDICT_*` in `TR.DAT`
  `~update_time` sections are stored as cumulative time-dependent input states.
- This is a **piecewise configuration history**, **not** an alternative
  solution, and its change times are input events (not interpolated CDF time
  arrays). For time-dependent visualization, sample the latest applicable
  input state at any requested `TIME3` timestamp.
- PT_SOLVER can choose several anomalous contributions and a neoclassical
  component in one region; do not interpret a `model[].identifier` as one
  exclusive regional choice.

## Inspect and verify

```bash
make verify-pt-solver
make verify-all

python3 convert.py -i examples/data/iter/38530P94.CDF \
    -b HDF5 -new --dd-version 4.1.1 --ids core_transport
```

`make verify-pt-solver` creates
`build/pt-solver-regions/pt-solver-regions.json`, inventories **every** bundled
CDF's input provenance (including absent inputs), and checks time-dependent
boundary selection, native output/metadata separation and the actual writer
using IMAS-shaped objects. On an IMAS-equipped workstation, perform an actual
`core_transport` write/readback as a separate validation.

```python
import xml.etree.ElementTree as ET

# After reading core_transport from an existing IMAS entry:
regions = ET.fromstring(core_transport.code.parameters.value).find(
    'pt_solver_regions')
for state in regions.findall('./boundary_and_equation_timeline/state'):
    print(state.get('effective_from'))
    for interval in state.findall('./radial_domains/region'):
        print(interval.attrib)
```

## References

- [TRANSP PT_SOLVER documentation](https://transp.pppl.gov/modules/ptsolver.html):
  region definitions, equation-specific X*BOUND values, default radial limits,
  model selections/scaling and scheduled input changes.
- [TRANSP pedestal module documentation](https://transp.pppl.gov/modules/pedestal.html):
  predictive boundary conditions and their interaction with pedestal models.

**Known limitation:** This version records PT_SOLVER regional *input*
configuration. It does **not** recover historical region-specific coefficient
arrays if TRANSP saved only a globally assembled effective coefficient.
A future CDF with independently saved region/model contributions will be
necessary before separate regional output model records can be implemented
without inventing physical information.
