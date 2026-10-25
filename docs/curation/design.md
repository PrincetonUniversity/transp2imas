# Scientific curation design and ideas

## Goal and boundaries

Turn a TRANSP result into a defensible **scientific usability record**: which quantities can support
which conclusions, in which time/radial regions, under which assumptions, and with which unresolved
uncertainties. A run can be internally consistent and numerically stable but experimentally weakly
constrained. Conversely, missing rotation need not invalidate every use of a discharge.

This design follows the project discussion about experimental workflow inputs, manually entered beam
and heating parameters, fitting errors, missing/unreliable NSTX rotation, uncertain impurity content,
and consistency of density prediction with magnetic diffusion. It is a development plan, not a claim
that those capabilities are already implemented or that a particular discharge is wrong.

## Separate five kinds of evidence

1. **Origin:** measured, fitted, automatic, manually entered, assumed, inferred or predicted. Automatic
   retrieval is not proof of accuracy; manual entry with an authoritative configuration can be reliable.
2. **Uncertainty:** measurement noise, shared calibration shifts, fitting/regularization choices,
   radial registration, model discrepancy and missing-input assumptions. Missing information stays unknown.
3. **Internal consistency:** units, coordinate/time alignment, particle/energy/momentum/current accounting,
   quasineutrality and compatible equation/model inputs.
4. **Numerical reliability:** actual residuals, stopping conditions, rejected steps, update cadence and
   sensitivity to resolution, time step, tolerances and stabilization.
5. **Applicability:** separate policies for heat transport, momentum transport, current evolution,
   turbulence-model validation, digital-twin use and surrogate training.

Do not merge these into a single opaque quality score. Keep pass/fail/warning/unknown per rule and
preserve local validity masks. The first implementation uses preliminary application summaries, not
calibrated acceptance policies.

## Input evidence package

Capture the original diagnostic or configuration reference, revision/calibration, coordinates, units,
valid interval, transformation/fitting recipe, original value and any override with justification.
Record whether a quantity is an input, initial/boundary condition, calibration target or independent
validation observation. Shared diagnostic evidence cannot become independent simply by appearing in a
second IDS or in a different fit.

A CDF and TR.DAT can support retrospective checks but usually cannot recover rejected diagnostic
points, fit covariance, calibration uncertainty or every user override. Extend preprocessing workflows
to emit a small evidence package alongside TRANSP inputs. The current JSON sidecar starts that contract.

Never silently modify source data. Corrections must generate a new revision with parent hashes,
changed values, reason, author/reviewer, and an audit trail. An uncertainty realization is a separate
run/member, not a replacement of the baseline.

## Heating and beam configuration

Maintain device/campaign-specific authoritative configurations for source identity, geometry,
injection direction, species, beam-energy component definitions, RF launchers and power reference
locations. Compare **effective inputs actually used**, not just an unexpanded template.

Check source-specific power/timing, energy and units, current versus particle versus power fractions,
source identities, co/counter conventions, duplicate source contributions and configuration applicability.
For example, FFULLA/FHALFA describe current fractions; no equality with power fractions should be assumed.
The initial implementation validates explicit scalar fractions and analyst-supplied native-unit
expectations; it does not yet authenticate campaign configurations or resolve all TRANSP input syntax.

Do not force injected power to equal instantaneous thermal heating. Source-accounting tests need the
appropriate fast-particle storage and loss terms. Likewise, do not tune beam inputs merely to make a
preferred transport model agree: distinguish calibration from independent model validation.

## Profile fits, gradients and identifiability

Keep measurement points, uncertainty, exclusions, fit method, smoothing/regularization, boundary
constraints, covariance or an ensemble of profiles, and diagnostics used. Test the fit against points,
against alternative defensible fitting choices, and against the derived gradients/time derivatives.

For an illustrative local convention, with temperature in energy units,

    chi_eff = -q_conductive / (n * dT/dr)
    particle_flux = -D * dn/dr + V * n

Weak/uncertain denominators can make inferred coefficients ill-conditioned even when fluxes are useful.
One particle-flux relation does not independently identify both D and V. Do not remove large or negative
coefficients merely to obtain attractive results; check sign conventions, convection subtraction,
uncertainty and the suitability of a purely diffusive interpretation.

A future profile-ensemble implementation should propagate correlated measurement, calibration,
registration and fitting uncertainty together. It should report quantiles/sign uncertainty for ratios
whose denominator can change sign rather than impose a symmetric error bar. Pointwise errors are not
an adequate surrogate for derivative covariance. The initial gradient rule accepts a supplied
uncertainty for its exact discrete derivative; it does not manufacture that uncertainty.

## Rotation assumptions

Distinguish measured, partially measured, proxy-based, model-predicted and deliberately zero-rotation
scenarios. An OMEGA array proves neither measurement availability nor diagnostic reliability. A zero
array is not evidence that missing rotation was measured to be zero; conversely, a supported measured
zero should not be rejected automatically.

Preserve species, velocity/angular-frequency definition, reference radius and sign conventions.
Construct scenario ensembles constrained by available observations/torque information. Determine which
conclusions are sensitive to rotation assumptions; do not impose a universal statement that rotation
never matters for a particular transport channel.

## Composition and charge closure

Quasineutrality and Zeff constrain different composition moments; they do not generally identify a
unique multi-impurity/charge-state distribution. A quasineutral run can retain large uncertainties in
dilution, radiation, resistivity, neoclassical contributions and species-resolved transport.

Use a joint composition representation. Avoid independent perturbations of electron, main-ion,
impurity and fast-ion densities followed by undocumented repairs. Prevent double counting of total and
component densities. Distinguish diagnostic Zeff, composition-derived Zeff and effective resistivity
adjustments. Charge-state distributions may require first and second charge moments separately.

The current charge-closure rule operates only on an analyst-declared complete inventory including
fast-ion accounting. It does not discover an authoritative inventory, prove its completeness, infer
Zeff uncertainty, or validate the composition merely because charge closes.

## Equation and model-state consistency

Build a version-aware resolved timeline of prescribed/evolved/inferred quantities, boundary conditions,
initial conditions, source-update cadence, regions, selected contributions and channel multipliers.
Keep model switches and disabled intervals; do not summarize the entire run by its final selection.

Density prediction plus magnetic diffusion is **not prohibited**. Review closure among evolving electron,
main-ion/impurity/fast-ion populations, charge/Zeff constraints, resistivity/bootstrap/noninductive
current, plasma-current/loop-voltage conditions, equilibrium updates and source-update cadence.
A measured q profile might be an initial condition or validation reference, not a simultaneous imposed
constraint. A frozen geometry/model can be an explicit approximation rather than an illegal setup.

Legal combinations and defaults must be version-verified against TRANSP. In particular, confirm whether
simultaneously enabled models replace or add to each other; do not certify this solely from a converter
parser's convention. A full model validator needs the executable version and resolved runtime settings.
The initial parser intentionally reports unsupported syntax and does not fill missing switches with false.

## Balance closure versus validation

A flux inferred from balance can close that same balance by construction. Call that internal consistency,
not independent experimental validation. Comparisons using shared inputs must account for correlated
uncertainty. Report agreement and precision separately: very large error bars do not turn uninformative
results into high-quality validation.

Useful independent tests may include line-integrated diagnostics, radiation, magnetic information,
neutrons and other synthetic diagnostics, provided they were not already used to constrain the result.
Record calibration/validation roles and shared-source dependencies explicitly.

## Numerical and local validity

Keep runtime counts and timers on native TIME, not the surrounding TIME3 coordinate. Record residual
and upwind-channel definitions, equation activity and update cadence. The older example files contain
RESL2TE/RESL2TI/RESL2OMG and UPWIND_TE/UPWIND_TI/UPWIND_MO; their presence does not establish convergence,
and an upwind adjustment is not automatically a numerical heat flux. Zero-filled inactive channels
must not be interpreted as successful equation solves.

A later stage should include time-step/tolerance/stabilization reruns, physical-versus-numerical flux
contributions where actually available, and flags for startup, termination, sawteeth, ELMs, pellets,
heating changes and diagnostic outages. Native coordinates and explicit evaluated masks are retained
in the first stage; automatic event segmentation is not yet implemented.

## Integration and representation

Keep curation independent of database creation and IDS writes. Proposed progression:

    NativeSource + input/evidence records
        -> declared/resolved state timeline
        -> deterministic rule engine
        -> optional uncertainty/scenario runner
        -> application policy + review decision
        -> JSON/Markdown/masks and eventual IMAS representation

The current `curation/` package implements the source-screening part. Future work should share stable
source/provenance interfaces with conversion services while avoiding changes to scientific meanings.
Use suitable typed IMAS validity/fitting/numerics structures only after checking DD support/semantics;
retain richer curation findings and ensemble identities in a versioned sidecar/manifest.

Retain occurrence 0 for combined transport and regional occurrences for model views. Do not reuse
regional occurrence numbers as uncertainty-ensemble identities. Track ensemble members with separate
run/dataset identities and a parent manifest.

## Other useful extensions

Campaign anomaly detection can prioritize manual review, not reject unusual physics. Repeated NUBEAM
seeds can separate sampling variation from parameter sensitivity. Radial/equilibrium registration
uncertainty deserves explicit perturbations. Root-cause ranking should be based on dependency/sensitivity
evidence rather than confident guesses. For machine-learning datasets, preserve measured/inferred/
predicted labels and group related reruns/realizations in train/test splits to avoid leakage.
