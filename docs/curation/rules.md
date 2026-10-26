# Initial curation rules

Status applies to one rule, not the entire run. `fail` indicates a violated structural/physical
screening condition, `warning` requests review, `not_assessed` means evidence/support is missing,
and `pass` is limited to the predicate described below. No initial rule certifies a scientific use.

| Rule family | Implemented predicate | Important limitation |
|---|---|---|
| `GRID.TIME`, `GRID.TIME3` | Nonempty 1-D, finite, strictly increasing native coordinates | Duplicate/event times may need an explicit future segmentation policy; do not silently sort |
| `GRID.*.UNITS` | Time coordinate is explicitly in recognized seconds units | Unknown units are not guessed |
| `GRID.X` | Native `(TIME3, X)` radius is finite, increasing, in [0,1] | No midplane-to-flux-coordinate conversion |
| `PROFILE.*.FINITE` | Selected NE/TE/TI samples are not fill/missing/nonfinite | Missing measurements remain unavailable, not zeros |
| `PROFILE.*.NEGATIVE` | Selected finite density/temperature is nonnegative | No clipping or correction |
| `PROFILE.*.ZERO` | Report zero transport denominators for review | Startup/off-plasma samples are not automatically discarded |
| `PROFILE.*.GRID_MISMATCH` | Profile shape agrees with explicit X | Does not infer correspondence from equal lengths |
| `PROFILE.*.UNITS` | Known density and eV/keV temperature units | Other conventions require an explicit adapter |
| `FIT.*.WEAK_GRADIENT` | Small native dT/dX relative to profile scale | An identifiability warning, not a rejection of flat-temperature physics |
| `FIT.*.SIGMA_VALIDITY` | Supplied derivative sigma is finite and nonnegative | No covariance estimation |
| `FIT.*.GRADIENT_UNCERTAINTY` | Check whether declared derivative uncertainty includes zero | Uses the supplied exact-estimator uncertainty only |
| `EVIDENCE.*`, `UNCERTAINTY.*` | Record origin, reference, reliability, assumption and uncertainty declarations, with `coverage` of the record's time applicability on the quantity's native clock | The external evidence is not authenticated; disjoint or unplaceable coverage is `not_assessed`, partial coverage downgrades a pass to a warning |
| `INDEPENDENCE.*` | Shared lineage between an independent validator and constraining evidence: same node, same SHA-256 or source string anywhere in the ancestries; selectors separate a shared container from a shared resource; usage is scoped in time | Only declared dependencies are visible; no shared lineage is not proof of independence; a shared container with different selectors is a `pass` that says so |
| `ROTATION.AVAILABILITY` | Recognized rotation output present/absent | Presence does not imply reliable measured rotation |
| `COMPOSITION.QUASINEUTRALITY` | Charge closure of an explicitly complete, grid/unit-compatible inventory | Closure does not identify a unique composition or validate Zeff |
| `INPUT.TIMELINE` | Declared assignments in file order (f90nml-based adapter), update intervals, parser kind/version, gaps and unresolved constructs | Not the resolved TRANSP configuration: defaults, omitted bounds, preprocessing and model-execution semantics are not applied |
| `MODEL.DENSITY_CURRENT_COUPLING` | Review concurrent declared density prediction and magnetic diffusion | This is not a forbidden combination |
| `MODEL.ROTATION_CONTROLS` | Review simultaneous declared NLVWNC/NLIVPO choices | Verify compatibility/defaults against actual TRANSP version |
| `HEATING.FRACTIONS` | Explicit scalar FFULLA/FHALFA in [0,1] with sum <=1 | Current fractions, not power fractions; no general beam geometry verification |
| `HEATING.REFERENCE.*` | Native scalar declarations versus analyst reference/tolerances | Reference is not independently authenticated |
| `NUMERICS.*.FINITE` | Runtime-channel availability on its native grid | All-zero output may be inactive/unused |
| `NUMERICS.ITERATION_DOMAIN` | NPTRITER is nonnegative and integer-valued | Count zero does not prove convergence |
| `NUMERICS.ITERATION_LIMIT` | Optional analyst screening threshold | Not the solver's iteration limit unless independently established |
| `NUMERICS.RESL2*.DOMAIN` | L2 residual nonnegative | Residual definitions/normalizations are not inferred |
| `NUMERICS.RESL2*.THRESHOLD` | Optional threshold with supplied rationale | Not a complete stopping-condition check |
| `RUN.IDENTITY` | Declared device/shot equal the CDF's own `device`/`shot` attributes | Not assessed if the CDF states neither; a match does not authenticate the evidence |
| `INPUT.RESOLVED.CODE` | Declared code version equals the CDF's `TRANSP_version` | The attribute is the build's, not proof of the executable |
| `INPUT.RESOLVED.FILES` | Declared input files hashed (SHA-256) and compared | A match proves identity of the file, not that TRANSP read it |
| `INPUT.RESOLVED.CONTROL.<NAME>` | Resolved runtime control compared with the declaration over every atomic half-open piece of the assessed support; `<NAME>` may be an indexed element | Gaps, unknown initial states and undeclared pieces are `not_assessed`; a declaration-only control cannot confirm the run; no type coercion |
| `EVENTS.INVENTORY` | Declared events with origin, padding, review state and membership; declarations that matched nothing, can act on nothing, or lie outside the source time; mandatory rules an event cannot waive | Membership is not detection; exclusions are policy and never change raw findings |
| `EVENTS.DETECT.<id>` | Optional abrupt-drop candidates on an explicitly configured, unit-checked diagnostic | Scores are uncalibrated; no candidates is not evidence of an event-free plasma |
| `FIT.*.GRADIENT_COVARIANCE` | Gradient uncertainty `D C D^T` propagated from a supplied fixed-grid covariance, flagged where `|g| <= k sigma` | Conditional on the supplied reconstruction and grid; radial-coordinate uncertainty and cross-time correlation are not included |
| `RULE.<id>.EXECUTION` | A rule raised an exception (an execution error, not missing evidence) | Added by the engine; the CLI exits 2 |
| `HEATING.REFERENCE.*` (v2) | Applicability, time restriction and reviewed override applied before comparison | Unknown applicability is `not_assessed`; an override is shown, never silent |
| `NUMERICS.CONVERGENCE` | Explicitly unassessed in stage 1 | Requires equation activity, stopping conditions and sensitivity evidence |

## Local masks and weak-gradient criterion

For each valid native profile row, the first implementation computes

    g = numpy.gradient(T_eV, X, edge_order=1)
    g_small = weak_gradient_relative * max(abs(T_eV)) / (max(X) - min(X))

and flags `abs(g) <= g_small`. This is a conditioning heuristic, not a turbulence threshold. A row
with unavailable temperature or nonmonotonic radius is not differentiated. Inspect evaluated masks
and coverage before interpreting an unflagged sample. Near-axis symmetry and intentional flat profiles
can legitimately trigger this warning; select a scientifically justified analysis window.

Charge closure uses `abs(sum(Z_s*n_s)-NE)/NE`, only where all required quantities are valid and NE>0.
The default 2% tolerance is a screening convenience, not a confidence interval.

## Native runtime inventory

Explicit names: `NPTRITER`, `PT_CPTIME`, `PT_WTIME`, `RESL2TE`, `RESL2TI`, `RESL2OMG`, `UPWIND_TE`,
`UPWIND_TI`, `UPWIND_MO`. This is a curation inventory, **not an extension of the converter's typed
IMAS numerical mappings**. The old archived examples contain residual/upwind channels even where
some are inactive/all-zero. Upwind profiles can live on XB rather than X.

Do not infer transport-solver behavior from unrelated iteration channels such as ISOITER/RZITER.
No broad wildcard is used to treat arbitrary diagnostics as PT_SOLVER output.

## Safety and limitations of the declared-input parser

The parser supports scalar numbers (including D exponents), logicals, quoted strings, comments
outside quotes, indexed scalar keys and `~update_time` snapshots. Unsupported statements are recorded
as gaps; a superseded value that becomes unparseable is removed, not silently reused. Nonincreasing
or invalid update times disable temporal model/reference rules.

Arrays/repeat notation, expressions, multiple assignments on one physical line, includes, derived
model-default resolution and PT_TEMPLATE execution semantics are not resolved in this initial
parser. Their presence is evidence of an incomplete reconstruction, not of an invalid TRANSP run.
The full original inputs remain referenced by path/hash; they are never rewritten.


## Rule ids, versions and the ledger

Each rule family above is produced by a registered rule (see [architecture.md](architecture.md)) with a stable id and
version; every finding records the `rule` and `rule_version` that produced it. A finding only says what a rule
observed. Whether a rule ran at all is in the report's `execution_ledger`, and a rule that did not execute leaves its
checks unassessed.

## Policy exclusions are not finding statuses

A declared event can make flagged samples ineligible for an application (`policy_exclusions`, `summary.excused_findings`)
but it never turns a `fail` or `warning` into a `pass`. Source-integrity rules (`GRID.*`, `PROFILE.*.GRID_MISMATCH`,
`PROFILE.*.UNITS`, `*.READ`, `RUN.IDENTITY`) cannot be excluded at all.
