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
| `EVIDENCE.*`, `UNCERTAINTY.*` | Record origin, reference, reliability, assumption and uncertainty declarations | The external evidence is not authenticated |
| `INDEPENDENCE.*` | Exact-source-string reuse as both constraint and independent validator | No full dependency graph yet |
| `ROTATION.AVAILABILITY` | Recognized rotation output present/absent | Presence does not imply reliable measured rotation |
| `COMPOSITION.QUASINEUTRALITY` | Charge closure of an explicitly complete, grid/unit-compatible inventory | Closure does not identify a unique composition or validate Zeff |
| `INPUT.TIMELINE` | Explicit scalar declarations and update intervals, with parser gaps | Not a full Fortran/preprocessor/runtime-default resolver |
| `MODEL.DENSITY_CURRENT_COUPLING` | Review concurrent declared density prediction and magnetic diffusion | This is not a forbidden combination |
| `MODEL.ROTATION_CONTROLS` | Review simultaneous declared NLVWNC/NLIVPO choices | Verify compatibility/defaults against actual TRANSP version |
| `HEATING.FRACTIONS` | Explicit scalar FFULLA/FHALFA in [0,1] with sum <=1 | Current fractions, not power fractions; no general beam geometry verification |
| `HEATING.REFERENCE.*` | Native scalar declarations versus analyst reference/tolerances | Reference is not independently authenticated |
| `NUMERICS.*.FINITE` | Runtime-channel availability on its native grid | All-zero output may be inactive/unused |
| `NUMERICS.ITERATION_DOMAIN` | NPTRITER is nonnegative and integer-valued | Count zero does not prove convergence |
| `NUMERICS.ITERATION_LIMIT` | Optional analyst screening threshold | Not the solver's iteration limit unless independently established |
| `NUMERICS.RESL2*.DOMAIN` | L2 residual nonnegative | Residual definitions/normalizations are not inferred |
| `NUMERICS.RESL2*.THRESHOLD` | Optional threshold with supplied rationale | Not a complete stopping-condition check |
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
