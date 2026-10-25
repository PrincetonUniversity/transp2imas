# Evidence sidecar format (v1)

The JSON file is validated by `curation.manifest.validate_manifest`. Unknown keys, nonfinite numbers,
inverted windows, duplicate species, inconsistent types and unsupported derivative conventions are
rejected. JSON is used deliberately to avoid an additional parser dependency. The file describes
analyst-provided evidence; references are not automatically opened or authenticated.

## Minimal document

```json
{
  "schema_version": "transp2imas.curation.inputs.v1"
}
```

This is valid, but most scientific provenance/uncertainty checks will be `not_assessed`.
`examples/curation/evidence-template.json` is a starting point, not a campaign-approved configuration.

## Assessment window

```json
"assessment_window": {
  "time_min": 0.2,
  "time_max": 0.6,
  "rho_min": 0.1,
  "rho_max": 0.8
}
```

Times are in seconds. Radius is native normalized X/XB, not a dimensional midplane radius. Bounds
are inclusive for samples. Declared-input states use half-open time intervals. Omitted bounds mean
all available samples; nothing is silently excluded near the axis, at startup or at the edge.
An empty selection is `not_assessed`, never a pass. Coordinates are not interpolated to enforce a window.

## Evidence by scientific role

Allowed roles are `electron_temperature`, `ion_temperature`, `electron_density`, `rotation`,
`composition`, and `heating`. Example (illustrative, not actual NSTX metadata):

```json
"evidence": {
  "rotation": {
    "origin": "assumed",
    "source": "analyst scenario manifest, revision 2",
    "uncertainty": "unknown",
    "roles": ["input"],
    "reliability": "unknown",
    "independent_validation": false,
    "notes": "Missing diagnostic; zero rotation is a sensitivity case, not a measurement."
  }
}
```

Origins: `measured`, `fitted`, `manual`, `automatic`, `assumed`, `inferred`, `predicted`, `missing`,
`unknown`. Uncertainty: `provided`, `unknown`, `not_applicable`. Reliability: `reviewed`, `unreliable`,
`unknown`. Usage roles: `input`, `calibration`, `validation`, `initial_condition`, `boundary_condition`.
A source reused as a constraint and claimed as an independent validator generates a warning when
its source strings match exactly. This limited check cannot prove independence of other records.

`uncertainty: provided` means the analyst declares an uncertainty description at the reference.
It does not mean covariance has been checked or propagated. `not_applicable` does not certify an
uncertainty-sensitive scientific use and currently remains unassessed for that screening dimension.

## Composition

```json
"composition": {
  "electron_density": "NE",
  "complete": true,
  "includes_fast_ions": true,
  "source": "explicit composition inventory and fast-ion accounting reference",
  "species": [
    {"density": "ND", "charge": 1},
    {"density": "NIMP_TEST", "charge_variable": "Z_TEST"}
  ]
}
```

The names above are schematic. List all positively charged populations once. Do not include a total
and its components, or reuse NE as an ion. `includes_fast_ions: true` asserts that fast-ion charge is
accounted for (including a documented case with no fast ions); it is not inferred by the evaluator.
Completeness is an analyst declaration, not established by a variable-name heuristic.

Each density must have recognized m^-3 or cm^-3 units and the same native `(TIME3, X)` dimensions
and shape as NE. A charge variable must have that same shape and dimensionless units. Exactly one of
constant `charge` or `charge_variable` is required. No cross-grid interpolation is performed. A
closure pass is not a validation of species identity, charge-state distribution, Zeff or dilution.

## Configuration reference expectations

```json
"input_expectations": [
  {
    "name": "EINJA(1)",
    "expected": 80000,
    "units": "eV",
    "atol": 1000,
    "rtol": 0,
    "reference": "approved device/campaign configuration identifier"
  }
]
```

Values are compared in the declared native namelist units; the evaluator does not convert a
configuration value or infer its units. The acceptance test is
`abs(actual-expected) <= atol + rtol*abs(expected)`. One reference applies to every selected declared
input interval. A differing value is a warning and remains unchanged. Missing/unparsed declarations
remain unknown. More complex campaign/time-specific reference policies are future work.

## Gradient uncertainty

```json
"gradient_uncertainties": {
  "TE": {
    "sigma_variable": "DTE_DX_SIGMA",
    "coordinate": "X",
    "method": "numpy_gradient_edge_order_1",
    "source": "derivative ensemble/covariance record for this exact discrete estimator"
  }
}
```

Only TE/TI are supported. The sigma array is the standard deviation of **dT/dX**, not pointwise T.
It must have the same dimensions/shape as the profile, and EV/KEV units per dimensionless X. The
algorithm uses NumPy gradient with edge_order=1; the supplied uncertainty must describe that estimator.
The reported interval check is `abs(dT/dX) <= k*sigma`. It is not automatically a calibrated
confidence probability. No independence assumption or derivative covariance is synthesized.

## Screening thresholds

```json
"thresholds": {
  "quasineutrality_relative": 0.02,
  "weak_gradient_relative": 0.001,
  "gradient_sigma_multiplier": 2,
  "newton_iterations_warning": 40,
  "residual_limits": {
    "RESL2TE": {
      "limit": 0.001,
      "rationale": "Analyst-selected example; replace with verified native residual normalization."
    }
  }
}
```

These are **example screening thresholds, not recommended universal TRANSP tolerances**. Only the
first three have defaults. Iteration/residual thresholds are off unless supplied. A residual limit
requires a rationale explaining its normalization. No PT_SOLVER input tolerance is silently equated
to a CDF residual definition. Even a residual-screening pass does not certify nonlinear convergence.

Optional top-level `notes` is free text. No general-purpose expressions, executable callbacks or
automatic corrective actions are accepted from this sidecar.
