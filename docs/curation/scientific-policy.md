# Scientific policy and guardrails

What the curation and conversion code will and will not do with uncertain information. These are policy
statements backed by tests; they are not a validation of any run.

## Missing data

- A masked or fill-value sample is **missing**, never a physical number. `NativeSource` decodes `_FillValue`,
  `missing_value`, `valid_range`/`valid_min`/`valid_max` and the netCDF default fill explicitly, applies packing
  once, and judges `valid_range` on the packed value. The converter reads through netCDF4's own masking and treats
  a masked sample that reaches a writer as a named broken mapping.
- **Missing rotation is not zero.** A run without rotation output gets `ROTATION.AVAILABILITY: warning`;
  zero rotation is a scenario that must be declared as such.
- A quantity that is absent is different from one that is present and zero (for example `TQRPL` absent versus all
  zero). The NBI torque total is complete only when each required term is present; an absent term withholds the
  total instead of summing the rest.
- `not_assessed` is not `pass` and is not an uncertainty of zero.

## Physics that is not decided by a screen

- **Density prediction with magnetic diffusion** (`MODEL.DENSITY_CURRENT_COUPLING`) is a review flag: the
  combination is allowed. It never invalidates a run.
- **Closure is not validation.** Charge closure, a power-fraction sum, zero residuals or a NPTRITER of zero do not
  establish correctness; each message says so.
- **No clipping.** A negative density or temperature, or a negative diffusivity, is reported unchanged. No
  quantity is repaired by curation.
- **Regional occurrence numbers are model views**, not uncertainty-ensemble identities.
- **PT_SOLVER model semantics** (whether enabling a model in a category replaces or adds to the previous one) are
  *unresolved*: no TRANSP source or documentation was available. The exported timeline assumes replacement and says
  so (`execution_semantics="unresolved"`, `interpretation="replacement_assumed_unverified"`), and the raw
  `%active` declarations are exported beside it so a reader can apply the correct rule.

## Evidence

- Independence is lineage-based and one-sided: shared lineage is a warning, absence of shared lineage is not proof.
- Applicability that cannot be established is *unknown*, not *applies*. An override needs a value, a reason and a
  reviewer; the original value stays in the report.
- Events are declared by the analyst and are applied as exclusions with counts; they are not detected.
- A resolved control must say whether it is a resolved runtime value; a declaration cannot confirm the run.

## Conversion publication

Writes to an IMAS entry are sequential, not transactional. When a backend write fails, `ImasStore` reports what was
requested, accepted, failed and left unattempted; accepted writes stay accepted and there is no rollback. The journal
records `publication.status` (`none`, `partial`, `complete`) and a partial conversion does not print the completion
marker. Remove or rebuild the entry yourself after a partial publication.


## Added with the rule/policy/schema work

- **Raw findings are never rewritten by an exclusion.** A declared event can make samples ineligible for an application;
  the finding stays as observed and the exclusion is a separate, counted record. Source-integrity checks cannot be
  excluded at all, and an automatic candidate or rejected event excludes nothing.
- **A check that did not run is not a check that passed.** Rules disabled by policy, rules whose prerequisites did not
  execute, and rules that raised leave every application out of `screening_only_not_certified`. An exception in a rule
  is an execution error, never `not_assessed`.
- **Evidence covers only the time it applies to.** A record whose time applicability does not meet the quantity's native
  support supplies no coverage, and partial coverage is reported with its unassessed intervals.
- **Aliases are not independence.** Different strings or hashes never prove independence; a shared container with
  different selectors is reported as a distinction, not as dependence and not as independence.
- **A conflicting declared identity is not trusted, and the file's value is not silently substituted.** The item becomes
  unknown and authorizes no reference.
- **Unresolved array updates leave nothing stale.** An update the reader cannot resolve invalidates the elements it may
  have changed.
- **Uncertainty scope is labelled.** A supplied covariance gives profile/gradient uncertainty conditional on the supplied
  reconstruction and grid; it says what it does not include (cross-time correlation, radial-coordinate uncertainty, any
  source not listed). Missing uncertainty stays unknown. Detector scores are uncalibrated and an absence of detected
  events is not evidence of an event-free plasma.
- **A schema-valid sidecar is not authenticated evidence.**
