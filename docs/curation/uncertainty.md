# Profile and gradient uncertainty (`curation/uncertainty.py`)

A bounded first increment. It propagates an explicitly supplied covariance, or
an explicitly supplied sample ensemble, through the estimator curation already
uses for gradients. It does not rerun TRANSP.

## Three tiers (only the first two are provided)

1. **Profile/gradient uncertainty conditional on the supplied reconstruction.**
   Result `scope` is `profile_gradient_conditional_on_supplied_reconstruction`.
   The profile covariance (or samples) is an input; the gradient uncertainty is
   what follows from it for a fixed estimator and a fixed radial grid.
2. **Frozen-physics sensitivity of a defined quantity.** Result `scope` is
   `frozen_physics_sensitivity_of_defined_quantity` (`ratio_statistics`). It
   describes the spread of a named ratio over supplied joint members. It is not
   a probability that a transport coefficient has a given value.
3. **Self-consistent TRANSP ensembles from independent reruns. NOT provided.**
   Sources, equilibrium, and transport feedback would change between members.
   Nothing here should be described as such an ensemble.

Missing uncertainty is UNKNOWN: a missing covariance raises (`validate_covariance`)
or yields `None`/NaN with `evaluated=False` for that slice (`propagate_gradient`).
It is never replaced by zeros. Synthetic data support no calibration claim.

## Fixed-grid covariance propagation

`derivative_operator(x)` returns the exact matrix D with
`D @ y == numpy.gradient(y, x, edge_order=1)` (built by applying `numpy.gradient`
to identity columns; tested to ~1e-13 on random nonuniform grids). With
`unit_scale` s:

    mean_gradient = s * D @ mean_profile
    C_gradient    = s^2 * D @ C_profile @ D.T

`unit_scale` converts both profile and sigma (EV to KEV: 1e-3, covariance 1e-6).
The gradient stays per unit of the coordinate as given. One operator is built
per time slice because the grid may differ per slice. The edge points use the
first-order one-sided difference, as in `check_gradient`.

Each `CovarianceSpec` carries quantity, units, coordinate, covariance, ordering
(`ascending_coordinate` only), estimator (must equal
`numpy.gradient(y, x, edge_order=1)`), provenance, and `included_sources`.

### Validation (no silent repair)

- shape `(n, n)` for a profile of length `n`; `n <= MAX_DENSE_POINTS = 4096`;
- finite entries;
- symmetry: `max|C - C^T| <= 1e-10 * max|C|` (`SYMMETRY_RTOL`);
- positive semidefinite: `lambda_min >= -1e-10 * max|lambda|` (`PSD_RTOL`).
  A singular PSD matrix (zero eigenvalues) is valid.

Violations raise `ValueError` with the reason. Nothing is symmetrized or
clipped. The only adjustment is opt-in `regularization=r` with
`r <= 1e-6`: it adds `r * mean(diag C) * I` after the checks pass and is
reported in the validation record. It cannot rescue an invalid matrix.

### What the formula excludes

- **Shared additive calibration offset.** `C = sigma^2 * 1 1^T` gives zero
  gradient variance because `D @ 1 = 0`; independent per-point terms do not.
  This is exact for an additive offset only. A multiplicative gain error does
  affect the gradient and must be represented in the covariance.
- **Uncertain radial coordinate.** The formula treats x as exact. Coordinate or
  registration uncertainty needs an augmented model or a sample ensemble with
  per-member coordinates (`gradient_from_samples` accepts `x` of shape
  `(members, points)`).
- **Cross-time and cross-diagnostic correlation.** One covariance per slice
  means correlation between slices is unknown and unmodeled; independence is not
  assumed, and the result records `cross_time = 'unmodeled'`. The `cross_time`
  argument and time-stacked covariances raise `NotImplementedError`.

No dense time-by-radius matrix is ever built.

## Samples

- `gradient_from_samples(samples, x, member_ids=None)` applies the estimator to
  each member and returns mean, covariance (unbiased, `None` for fewer than two
  accepted members), and quantiles. Members with non-finite values or
  non-strictly-increasing coordinates are rejected and recorded with member id
  and reason. Coordinates are never sorted and temperatures never clipped. If
  all members are rejected it raises.
- `sample_from_covariance(mean, cov, n, seed, distribution='gaussian', ...)`:
  covariance alone does not determine a distribution, so the caller must declare
  it; only `'gaussian'` is implemented. `seed` is a required integer
  (`numpy.random.Generator(PCG64)`). Same seed gives bit-identical samples, and
  the result does not depend on `batch_size` (tested). The factorization is an
  eigen-decomposition; eigenvalues `<= 1e-10 * lambda_max` are set to zero (count
  recorded), so variance below that relative level is dropped. Limits:
  `n <= max_draws` (default 1e6) and `n * points * 8 <= max_bytes` (default
  256 MiB). The record holds algorithm, seed, numpy version, limits, validation,
  and failed members.
- Monte Carlo test tolerance: with a fixed seed and n = 40000, every entry of the
  sample gradient covariance must lie within 5 standard errors of `D C D^T`, where
  for Gaussian members `var(S_ij) ~ (C_ij^2 + C_ii C_jj)/(n-1)`; means within 5
  standard errors of the mean.
- `calibration_offset_scenarios(profile, {name: {'kind', 'value', 'source'}})`
  makes explicit additive (`y + v`) or multiplicative (`y * v`) scenarios. A
  scenario already listed in `included_sources` raises `ValueError` (no double
  counting). They are named scenarios, not a statistical ensemble.
- `ratio_statistics(num, den, num_units, den_units, definition, eps)` uses paired
  members (ids must match when given) and requires a non-empty definition and
  units. It reports 5/50/95 percentiles of finite ratios, the undefined fraction
  (non-finite numerator or denominator, or denominator exactly zero, with ids),
  the near-zero denominator fraction (`|den| <= eps`, over members with finite
  denominator), and sign ambiguity (`sign_ambiguous` and the minority-sign
  fraction). Extreme finite ratios are kept; nothing is clipped.

## Outputs

Every result has `scope`, `as_dict()` (JSON-safe, NaN becomes null), and
`save_npz(path)`, which writes arrays without pickle plus `meta_json` and
`checksums_json` (sha256 of each array's bytes).

## Synthetic worked example (synthetic, not diagnostic data)

    import numpy as np
    from curation import uncertainty as U

    x = np.array([0.0, 0.2, 0.5, 0.7, 1.0])      # synthetic grid
    y = 2.0 - 1.5 * x**2                          # synthetic profile, keV
    cov = 0.05**2 * np.exp(-abs(x[:, None] - x[None, :]) / 0.4)   # synthetic
    spec = U.CovarianceSpec('synthetic_te', 'keV', x, cov, U.ORDERING,
                            U.ESTIMATOR, 'synthetic, not diagnostic data', ())
    res = U.propagate_gradient([spec], y[None], x[None])
    print(res.mean_gradient, res.std_gradient, res.scope)

`U.synthetic_example()` returns the same thing as a dict. Tests:
`python examples/verify_curation_uncertainty.py`.

## Use inside a curation run

The engine's optional rule `uncertainty.gradient_covariance` (finding `FIT.TE.GRADIENT_COVARIANCE` or `FIT.TI...`)
propagates a covariance supplied in the sidecar (`uncertainty.covariances`, evidence format v3):

* the resource is an `.npz` with `time` (exactly the source `TIME3`, no interpolation) and `covariance` of shape
  `(nt, nr, nr)`, one matrix per time slice on the source's own `X` grid; its SHA-256 is checked when the sidecar is
  loaded and the file is protected from being overwritten by the report;
* the covariance units (`EV` or `KEV`) are converted to the profile's own units by the square of the factor;
* a malformed covariance (wrong shape, asymmetric, not positive semidefinite, nonfinite, time-stacked, wrong time axis,
  unknown units) is an `error` finding, never repaired;
* the result is flagged where `|dT/dX| <= k sigma_g` (`k` is `gradient_sigma_multiplier`); `report.uncertainty` records
  the scope `profile_gradient_conditional_on_supplied_reconstruction`, `cross_time: unmodeled`, and the propagated
  `sigma` array is stored in `curation.npz` as `uncertainty_<id>__std_gradient`;
* the finding lists what is **not** included: radial-coordinate uncertainty, cross-time correlation, and any source not
  named in `included_sources`.

When no covariance is supplied nothing is invented: `UNCERTAINTY.*` stays `not_assessed`. The existing supplied-derivative
sigma path (`gradient_uncertainties`) is unchanged and independent.
