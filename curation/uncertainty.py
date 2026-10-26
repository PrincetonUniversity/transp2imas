#!/usr/bin/python
"""Bounded profile/gradient uncertainty propagation for curation.

Scope: uncertainty of a profile gradient CONDITIONAL on a supplied reconstruction
(covariance or samples), and frozen-physics sensitivity of explicitly defined
ratios.  This is NOT a self-consistent TRANSP ensemble: no rerun feedback on
sources, equilibrium or transport is represented.  Covariance is per time slice;
cross-time and cross-diagnostic correlation is unknown/unmodeled (NotImplementedError
if supplied).  The fixed-grid formula assumes the radial coordinate is exact; an
uncertain coordinate needs an augmented model or ensemble.  Missing uncertainty is
UNKNOWN (None or error), never zero.  Inputs are never mutated or repaired.
"""

from dataclasses import dataclass
import hashlib
import json

import numpy as np

SCOPE_GRADIENT = 'profile_gradient_conditional_on_supplied_reconstruction'
SCOPE_RATIO = 'frozen_physics_sensitivity_of_defined_quantity'
ESTIMATOR = 'numpy.gradient(y, x, edge_order=1)'
ORDERING = 'ascending_coordinate'
MAX_DENSE_POINTS = 4096          # per-slice dense covariance limit
SYMMETRY_RTOL = 1e-10            # max|C-C^T| <= rtol * max|C|
PSD_RTOL = 1e-10                 # lambda_min >= -rtol * max|lambda|
MAX_REGULARIZATION = 1e-6        # relative to mean(diag C)
DEFAULT_MAX_DRAWS = 1_000_000
DEFAULT_MAX_BYTES = 256 * 1024 * 1024


@dataclass
class CovarianceSpec:
    """Explicit covariance of one profile on one radial grid (None = unknown)."""
    quantity: str
    units: str
    coordinate: np.ndarray
    covariance: object
    ordering: str
    estimator: str
    provenance: str
    included_sources: tuple = ()


def _sha256(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def _write_npz(path, arrays, meta):
    checksums = {key: _sha256(value) for key, value in arrays.items()}
    np.savez(path, **arrays, meta_json=np.array(json.dumps(_jsonable(meta))),
             checksums_json=np.array(json.dumps(checksums)))
    return checksums


def validate_covariance(covariance, n, symmetry_rtol=SYMMETRY_RTOL, psd_rtol=PSD_RTOL,
                        regularization=None):
    """Return (covariance, report); raise ValueError on any defect.

    No symmetrization or eigenvalue clipping is done.  Only if the caller passes
    `regularization` (relative, <= MAX_REGULARIZATION) is C + reg*mean(diag C)*I
    returned, after the checks, and the amount is reported.
    """
    if covariance is None:
        raise ValueError('covariance missing: uncertainty is unknown, not zero')
    c = np.asarray(covariance, dtype=float)
    if n > MAX_DENSE_POINTS:
        raise ValueError(f'{n} points exceed MAX_DENSE_POINTS={MAX_DENSE_POINTS}')
    if c.shape != (n, n):
        raise ValueError(f'covariance shape {c.shape} does not match profile length {n}')
    if not np.all(np.isfinite(c)):
        raise ValueError('covariance contains non-finite entries')
    scale = float(np.max(np.abs(c))) if n else 0.0
    asym = float(np.max(np.abs(c - c.T))) if n else 0.0
    if asym > symmetry_rtol * scale:
        raise ValueError(f'covariance not symmetric: max|C-C^T|={asym:.3e} exceeds '
                         f'{symmetry_rtol:.1e}*max|C|; not repaired')
    eig = np.linalg.eigvalsh(c)
    emax = float(np.max(np.abs(eig))) if n else 0.0
    if float(eig.min()) < -psd_rtol * emax:
        raise ValueError(f'covariance not positive semidefinite: lambda_min={eig.min():.3e}, '
                         f'tolerance {psd_rtol:.1e}*{emax:.3e}; not repaired')
    report = {'symmetry_rtol': symmetry_rtol, 'psd_rtol': psd_rtol, 'asymmetry_abs': asym,
              'eigenvalue_min': float(eig.min()), 'eigenvalue_max': float(eig.max()),
              'regularization_relative': 0.0, 'regularization_absolute': 0.0}
    if regularization:
        if not 0 < regularization <= MAX_REGULARIZATION:
            raise ValueError(f'regularization must be in (0, {MAX_REGULARIZATION}]')
        absolute = float(regularization * np.mean(np.diag(c)))
        c = c + absolute * np.eye(n)
        report['regularization_relative'] = float(regularization)
        report['regularization_absolute'] = absolute
    return c, report


def covariance_from_samples(samples):
    """Unbiased sample covariance (ddof=1) of (members, points); needs >= 2 finite members."""
    s = np.asarray(samples, dtype=float)
    if s.ndim != 2 or s.shape[0] < 2:
        raise ValueError('need a (members>=2, points) sample array')
    if s.shape[1] > MAX_DENSE_POINTS:
        raise ValueError(f'{s.shape[1]} points exceed MAX_DENSE_POINTS={MAX_DENSE_POINTS}')
    if not np.all(np.isfinite(s)):
        raise ValueError('samples contain non-finite values; use gradient_from_samples to record rejections')
    return np.cov(s, rowvar=False, ddof=1).reshape(s.shape[1], s.shape[1])


def _check_grid(x):
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or x.size < 2:
        raise ValueError('coordinate must be 1-D with at least 2 points')
    if x.size > MAX_DENSE_POINTS:
        raise ValueError(f'{x.size} points exceed MAX_DENSE_POINTS={MAX_DENSE_POINTS}')
    if not np.all(np.isfinite(x)) or not np.all(np.diff(x) > 0):
        raise ValueError('coordinate must be finite and strictly increasing (not sorted here)')
    return x


def derivative_operator(x, edge_order=1):
    """Exact matrix D with D @ y == numpy.gradient(y, x, edge_order=1) for any y."""
    if edge_order != 1:
        raise ValueError('only edge_order=1 matches the curation estimator')
    x = _check_grid(x)
    return np.gradient(np.eye(x.size), x, axis=0, edge_order=1)


@dataclass
class GradientResult:
    """Fixed-grid gradient mean/covariance per slice; unknown slices have covariance None."""
    quantity: str
    units: str
    gradient_units: str
    coordinate: np.ndarray
    mean_gradient: np.ndarray
    covariance: list
    std_gradient: np.ndarray
    evaluated: np.ndarray
    unit_scale: float
    estimator: str
    provenance: list
    included_sources: list
    validation: list
    cross_time: str = 'unmodeled'
    scope: str = SCOPE_GRADIENT

    def as_dict(self):
        return _jsonable({
            'scope': self.scope, 'quantity': self.quantity, 'units': self.units,
            'gradient_units': self.gradient_units, 'unit_scale': self.unit_scale,
            'estimator': self.estimator, 'cross_time': self.cross_time,
            'provenance': self.provenance, 'included_sources': self.included_sources,
            'unknown_slices': [i for i, c in enumerate(self.covariance) if c is None],
            'validation': self.validation, 'mean_gradient': self.mean_gradient,
            'std_gradient': self.std_gradient, 'evaluated': self.evaluated})

    def save_npz(self, path):
        """Write arrays (no pickle) with sha256 checksums; returns the checksum dict."""
        arrays = {'coordinate': self.coordinate, 'mean_gradient': self.mean_gradient,
                  'std_gradient': self.std_gradient, 'evaluated': self.evaluated,
                  'units': np.array(self.units), 'gradient_units': np.array(self.gradient_units)}
        for i, c in enumerate(self.covariance):
            if c is not None:
                arrays[f'covariance_t{i:04d}'] = c
        return _write_npz(path, arrays, self.as_dict())


def propagate_gradient(spec_per_slice, profile, x, unit_scale=1.0, cross_time=None,
                       coordinate_units='coordinate unit'):
    """C_g = s^2 D C D^T and g = s D y per slice, with D the numpy.gradient operator.

    `profile` and `x` are (time, radius); one CovarianceSpec (or None = unknown) per
    slice.  `unit_scale` converts profile and sigma units (EV->KEV: 1e-3, covariance
    1e-6); the gradient stays per unit of the coordinate as given.  Slices with a
    None covariance get NaN and evaluated=False.  Cross-time correlation is not
    modeled: any `cross_time` argument or time-stacked covariance raises
    NotImplementedError.  A shared additive offset (C = s^2 11^T) yields exactly zero
    gradient variance; an uncertain radial coordinate is NOT covered.
    """
    if cross_time is not None:
        raise NotImplementedError('cross-time correlation is not implemented; per-slice '
                                  'covariances do not assume independence, they leave it unmodeled')
    y = np.asarray(profile, dtype=float)
    x = np.asarray(x, dtype=float)
    if y.ndim != 2 or x.shape != y.shape:
        raise ValueError('profile and x must both be (time, radius)')
    nt, nr = y.shape
    if len(spec_per_slice) != nt:
        raise ValueError(f'need one covariance spec per time slice ({nt}), got {len(spec_per_slice)}')
    if not np.isfinite(unit_scale) or unit_scale <= 0:
        raise ValueError('unit_scale must be finite and positive')
    mean = np.full(y.shape, np.nan)
    std = np.full(y.shape, np.nan)
    evaluated = np.zeros(y.shape, dtype=bool)
    covs, validation, provenance, sources = [], [], [], []
    names = set()
    for t, spec in enumerate(spec_per_slice):
        if spec is None or spec.covariance is None:
            covs.append(None)
            validation.append(None)
            continue
        if nt > 1 and np.shape(spec.covariance) == (nt * nr, nt * nr):
            raise NotImplementedError('time-stacked covariance (cross-time correlation) is not implemented')
        if spec.estimator != ESTIMATOR:
            raise ValueError(f'slice {t}: estimator {spec.estimator!r} is not {ESTIMATOR!r}')
        if spec.ordering != ORDERING:
            raise ValueError(f'slice {t}: ordering {spec.ordering!r} is not {ORDERING!r}')
        if not np.array_equal(np.asarray(spec.coordinate, dtype=float), x[t]):
            raise ValueError(f'slice {t}: spec coordinate differs from x')
        if not np.all(np.isfinite(y[t])):
            raise ValueError(f'slice {t}: profile has non-finite values')
        c, report = validate_covariance(spec.covariance, nr)
        d = derivative_operator(x[t])
        mean[t] = unit_scale * (d @ y[t])
        cg = unit_scale ** 2 * (d @ c @ d.T)
        std[t] = np.sqrt(np.clip(np.diag(cg), 0.0, None))  # diag >= 0 up to rounding only
        evaluated[t] = True
        covs.append(cg)
        validation.append(report)
        if spec.provenance not in provenance:
            provenance.append(spec.provenance)
        sources.extend(s for s in spec.included_sources if s not in sources)
        names.add((spec.quantity, spec.units))
    if len(names) > 1:
        raise ValueError(f'mixed quantity/units across slices: {sorted(names)}')
    quantity, units = next(iter(names), ('unknown', 'unknown'))
    return GradientResult(quantity, units, f'{units}/{coordinate_units}', x.copy(), mean, covs, std,
                          evaluated, float(unit_scale), ESTIMATOR, provenance, sources, validation)


@dataclass
class SampleSet:
    """Gaussian samples with a full reproducibility record."""
    samples: np.ndarray
    member_ids: np.ndarray
    record: dict
    scope: str = SCOPE_GRADIENT

    def as_dict(self):
        return _jsonable({'scope': self.scope, 'shape': list(self.samples.shape), **self.record})

    def save_npz(self, path):
        arrays = {'samples': self.samples, 'member_ids': self.member_ids}
        return _write_npz(path, arrays, self.as_dict())


def sample_from_covariance(mean, cov, n, seed, distribution=None, max_draws=DEFAULT_MAX_DRAWS,
                           max_bytes=DEFAULT_MAX_BYTES, batch_size=4096, regularization=None):
    """Draw n members mean + L z with z ~ N(0, I), L from the eigen-decomposition of C.

    The distribution must be declared by the caller (covariance does not fix it);
    only 'gaussian' is implemented.  `seed` is a required non-negative integer
    (numpy PCG64).  Eigenvalues <= PSD_RTOL*lambda_max (rounding noise or
    zero; variance below ~1e-10 of the largest is dropped) are set to zero in the
    factor and counted in the record.  Normals are drawn in row-major
    batches from one stream, so results do not depend on `batch_size` (tested).
    """
    if distribution != 'gaussian':
        raise ValueError(f'distribution {distribution!r} not implemented; declare distribution="gaussian"')
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError('seed must be a non-negative integer')
    mean = np.asarray(mean, dtype=float)
    if mean.ndim != 1 or not np.all(np.isfinite(mean)):
        raise ValueError('mean must be a finite 1-D array')
    if not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError('n must be a positive integer')
    if n > max_draws:
        raise ValueError(f'n={n} exceeds max_draws={max_draws}')
    if n * mean.size * 8 > max_bytes:
        raise ValueError(f'{n}x{mean.size} float64 samples exceed max_bytes={max_bytes}')
    if batch_size < 1:
        raise ValueError('batch_size must be >= 1')
    c, report = validate_covariance(cov, mean.size, regularization=regularization)
    w, v = np.linalg.eigh(c)
    keep = w > PSD_RTOL * max(float(w.max()), 0.0)
    clipped = int(np.sum(~keep))
    factor = v * np.sqrt(np.where(keep, w, 0.0))
    rng = np.random.Generator(np.random.PCG64(int(seed)))
    out = np.empty((n, mean.size))
    for start in range(0, n, batch_size):
        stop = min(start + batch_size, n)
        out[start:stop] = mean + rng.standard_normal((stop - start, mean.size)) @ factor.T
    ids = np.arange(n)
    bad = ~np.all(np.isfinite(out), axis=1)
    failed = [{'member_id': int(i), 'reason': 'non-finite draw'} for i in ids[bad]]
    record = {'distribution': distribution, 'algorithm': 'numpy.random.Generator(PCG64)',
              'seed': int(seed), 'numpy_version': np.__version__, 'n': int(n),
              'batch_size': int(batch_size), 'max_draws': int(max_draws), 'max_bytes': int(max_bytes),
              'factorization': 'eigh', 'eigenvalues_set_to_zero': clipped, 'validation': report,
              'failed_members': failed, 'member_ids': f'0..{n - 1}'}
    return SampleSet(out, ids, record)


@dataclass
class SampleGradientResult:
    """Sample-based gradient statistics; rejected members are recorded, never dropped silently."""
    mean_gradient: np.ndarray
    covariance: object
    quantile_levels: tuple
    quantiles: np.ndarray
    gradients: np.ndarray
    member_ids_used: np.ndarray
    rejected: list
    unit_scale: float
    estimator: str = ESTIMATOR
    scope: str = SCOPE_GRADIENT

    def as_dict(self):
        return _jsonable({
            'scope': self.scope, 'estimator': self.estimator, 'unit_scale': self.unit_scale,
            'n_used': int(self.member_ids_used.size), 'rejected': self.rejected,
            'covariance_known': self.covariance is not None,
            'quantile_levels': self.quantile_levels, 'mean_gradient': self.mean_gradient,
            'quantiles': self.quantiles})

    def save_npz(self, path):
        arrays = {'mean_gradient': self.mean_gradient, 'quantiles': self.quantiles,
                  'gradients': self.gradients, 'member_ids_used': self.member_ids_used}
        if self.covariance is not None:
            arrays['covariance'] = self.covariance
        return _write_npz(path, arrays, self.as_dict())


def gradient_from_samples(samples, x, member_ids=None, quantiles=(0.05, 0.5, 0.95), unit_scale=1.0):
    """Apply numpy.gradient(edge_order=1) to each member of (members, points).

    `x` is (points,) shared or (members, points) per member (e.g. displaced
    coordinates); it is never sorted.  Members with non-finite values or
    non-strictly-increasing coordinates are rejected and recorded with id and
    reason.  Covariance is None with fewer than 2 accepted members (unknown).
    """
    s = np.asarray(samples, dtype=float)
    if s.ndim != 2:
        raise ValueError('samples must be (members, points)')
    m, n = s.shape
    if n < 2 or n > MAX_DENSE_POINTS:
        raise ValueError(f'points must be in [2, {MAX_DENSE_POINTS}]')
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        x = _check_grid(x)
        if x.size != n:
            raise ValueError('x length does not match samples')
    elif x.shape != (m, n):
        raise ValueError('x must be (points,) or (members, points)')
    ids = np.arange(m) if member_ids is None else np.asarray(member_ids)
    if ids.shape != (m,) or len(set(ids.tolist())) != m:
        raise ValueError('member_ids must be unique, one per member')
    rejected, keep = [], []
    for i in range(m):
        xi = x if x.ndim == 1 else x[i]
        if not np.all(np.isfinite(s[i])):
            reason = 'non-finite profile value'
        elif not np.all(np.isfinite(xi)):
            reason = 'non-finite coordinate'
        elif not np.all(np.diff(xi) > 0):
            reason = 'coordinate not strictly increasing'
        else:
            keep.append(i)
            continue
        rejected.append({'member_id': ids[i].item(), 'reason': reason})
    if not keep:
        raise ValueError(f'all {m} members rejected: {rejected[:3]}')
    if x.ndim == 1:
        g = np.gradient(s[keep], x, axis=1, edge_order=1)
    else:
        g = np.array([np.gradient(s[i], x[i], edge_order=1) for i in keep])
    g = unit_scale * g
    cov = np.cov(g, rowvar=False, ddof=1).reshape(n, n) if len(keep) > 1 else None
    levels = tuple(float(q) for q in quantiles)
    return SampleGradientResult(g.mean(axis=0), cov, levels, np.quantile(g, levels, axis=0), g,
                                ids[keep], rejected, float(unit_scale))


@dataclass
class ScenarioSet:
    """Named explicit profile scenarios (not a statistical ensemble)."""
    scenarios: dict
    definitions: dict
    scope: str = SCOPE_GRADIENT

    def as_dict(self):
        return _jsonable({'scope': self.scope, 'definitions': self.definitions,
                          'scenarios': self.scenarios})


def calibration_offset_scenarios(profile, offsets, included_sources=()):
    """Named scenarios {name: {'kind': 'additive'|'multiplicative', 'value': v}}.

    additive: y + v; multiplicative: y * v.  A scenario whose name or optional
    'source' is already in `included_sources` is refused (double counting).
    """
    y = np.asarray(profile, dtype=float)
    included = set(included_sources)
    scenarios, definitions = {}, {}
    for name, spec in offsets.items():
        source = spec.get('source', name)
        if name in included or source in included:
            raise ValueError(f'scenario {name!r} is already in the covariance included_sources')
        kind, value = spec.get('kind'), spec.get('value')
        if kind not in ('additive', 'multiplicative') or not np.isfinite(value):
            raise ValueError(f'scenario {name!r}: kind must be additive/multiplicative with finite value')
        scenarios[name] = y + value if kind == 'additive' else y * value
        definitions[name] = {'kind': kind, 'value': float(value), 'source': source}
    return ScenarioSet(scenarios, definitions)


@dataclass
class RatioResult:
    """Paired-sample ratio statistics; extreme finite values are kept, nothing is clipped."""
    definition: str
    units: str
    quantile_levels: tuple
    quantiles: object
    n_members: int
    n_defined: int
    undefined_fraction: float
    near_zero_denominator_fraction: float
    sign_ambiguity_fraction: float
    sign_ambiguous: bool
    eps: float
    undefined_ids: list
    scope: str = SCOPE_RATIO

    def as_dict(self):
        return _jsonable(dict(self.__dict__))


def ratio_statistics(numerator_samples, denominator_samples, numerator_units, denominator_units,
                     definition, eps, quantiles=(5, 50, 95), numerator_ids=None, denominator_ids=None):
    """Quantiles of num/den over JOINT members (same ids), with failure accounting.

    Undefined: non-finite numerator or denominator, or denominator exactly 0.
    near_zero fraction: |den| <= eps over all members with finite den.
    sign_ambiguity_fraction: minority-sign share among members with nonzero finite
    den; sign_ambiguous is True if both signs occur.  Frozen-physics sensitivity of
    the defined quantity only; not a TRANSP rerun ensemble.
    """
    if not str(definition).strip() or not numerator_units or not denominator_units:
        raise ValueError('definition and both units are required')
    num = np.asarray(numerator_samples, dtype=float)
    den = np.asarray(denominator_samples, dtype=float)
    if num.ndim != 1 or num.shape != den.shape:
        raise ValueError('numerator and denominator must be paired 1-D samples of equal length')
    if (numerator_ids is None) != (denominator_ids is None):
        raise ValueError('give both member id arrays or neither')
    if numerator_ids is not None and not np.array_equal(numerator_ids, denominator_ids):
        raise ValueError('numerator and denominator member ids differ; pairing is required')
    if not np.isfinite(eps) or eps < 0:
        raise ValueError('eps must be finite and >= 0')
    ids = np.arange(num.size) if numerator_ids is None else np.asarray(numerator_ids)
    finite_den = np.isfinite(den)
    ratio = np.full(num.shape, np.nan)
    can = np.isfinite(num) & finite_den & (den != 0)
    ratio[can] = num[can] / den[can]
    ok = np.isfinite(ratio)
    levels = tuple(float(q) for q in quantiles)
    q = np.percentile(ratio[ok], levels) if ok.any() else None
    nonzero = finite_den & (den != 0)
    pos, neg = int(np.sum(den[nonzero] > 0)), int(np.sum(den[nonzero] < 0))
    minority = min(pos, neg) / (pos + neg) if pos + neg else 0.0
    near = float(np.sum(np.abs(den[finite_den]) <= eps) / num.size) if num.size else 0.0
    n_def = int(ok.sum())
    undefined_ids = [i.item() for i in ids[~ok]]
    return RatioResult(str(definition), f'({numerator_units})/({denominator_units})', levels, q,
                       int(num.size), n_def, 1.0 - n_def / num.size if num.size else 1.0, near,
                       float(minority), bool(pos and neg), float(eps), undefined_ids)


def synthetic_example():
    """Small SYNTHETIC example (made-up numbers; not diagnostic data, no calibration claim)."""
    x = np.array([0.0, 0.2, 0.5, 0.7, 1.0])
    y = 2.0 - 1.5 * x ** 2
    sigma = 0.05 * np.ones(5)
    corr = np.exp(-np.abs(x[:, None] - x[None, :]) / 0.4)
    cov = np.outer(sigma, sigma) * corr
    spec = CovarianceSpec('synthetic_temperature', 'keV', x, cov, ORDERING, ESTIMATOR,
                          'synthetic, not diagnostic data', ())
    res = propagate_gradient([spec], y[None], x[None])
    return {'label': 'synthetic, not diagnostic data', 'x': x, 'profile': y,
            'covariance': cov, 'result': res}
