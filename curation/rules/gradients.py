#!/usr/bin/python
"""Weak-gradient identifiability screening and the supplied derivative-uncertainty path."""

import numpy as np

from curation import uncertainty
from curation.results import ROLE_APPS
from curation.source import unit_key



def check_gradient(ctx, rec, name):
    apps = ('heat_transport', 'surrogate_training')
    ident = 'FIT.' + name
    if not ctx.source.has('X'):
        rec.unavailable(ident + '.GRADIENT', 'No native radial grid.', [name, 'X'], apps)
        return
    data, x = ctx.source.array(name), ctx.source.array('X')
    if x.shape != data.shape or not ctx.normalized_radius('X'):
        rec.unavailable(ident + '.GRADIENT', 'Profile and normalized radial grid are not aligned.', [name, 'X'], apps)
        return
    scale = 1000.0 if unit_key(ctx.source.meta(name)['units']) == 'KEV' else 1.0
    y = data * scale
    gradient = np.full(y.shape, np.nan)
    limit = np.full((len(y), 1), np.nan)
    for row in range(len(y)):
        if (y.shape[1] >= 2 and np.all(np.isfinite(y[row])) and np.all(np.isfinite(x[row]))
                and np.all(np.diff(x[row]) > 0)):
            gradient[row] = np.gradient(y[row], x[row], edge_order=1)
            limit[row, 0] = (ctx.thresholds['weak_gradient_relative'] * np.max(np.abs(y[row]))
                             / (x[row, -1] - x[row, 0]))
    evaluated = rec.selection(name) & np.isfinite(gradient)
    weak = evaluated & (np.abs(gradient) <= limit)
    rec.samples(ident + '.WEAK_GRADIENT', name, evaluated, weak,
                  'Small dT/dX makes a local diffusive interpretation ill-conditioned; this is not a bad-profile verdict.',
                  applications=apps, metrics={'relative_screening_threshold': ctx.thresholds['weak_gradient_relative']},
                  recommendation='Check flux and coefficient validity separately; assess gradient uncertainty.')
    spec = ctx.manifest['gradient_uncertainties'].get(name)
    if spec is None:
        rec.unavailable(ident + '.GRADIENT_UNCERTAINTY',
                          'No supplied uncertainty for this derivative; none is inferred from pointwise errors.', [name], apps)
        return
    sigma_name = spec['sigma_variable']
    if not ctx.source.has(sigma_name):
        rec.unavailable(ident + '.GRADIENT_UNCERTAINTY', 'Declared derivative uncertainty is absent.', [sigma_name], apps)
        return
    sigma = ctx.source.array(sigma_name)
    units = unit_key(ctx.source.meta(sigma_name)['units'])
    if (sigma.shape != y.shape or ctx.source.meta(sigma_name)['dimensions'] != ctx.source.meta(name)['dimensions']
            or units not in ('EV', 'KEV')):
        rec.unavailable(ident + '.GRADIENT_UNCERTAINTY',
                          'Derivative uncertainty must share native dimensions and have eV/keV units per dimensionless X.',
                          [name, sigma_name], apps)
        return
    sigma = sigma * (1000.0 if units == 'KEV' else 1.0)
    bad_sigma = evaluated & (~np.isfinite(sigma) | (sigma < 0))
    rec.samples(ident + '.SIGMA_VALIDITY', name, evaluated, bad_sigma,
                  'Derivative standard deviations must be finite and nonnegative.', severity='fail', applications=apps)
    usable = evaluated & np.isfinite(sigma) & (sigma >= 0)
    k = ctx.thresholds['gradient_sigma_multiplier']
    finding = rec.samples(ident + '.GRADIENT_UNCERTAINTY', name, usable,
                           usable & (np.abs(gradient) <= k * sigma),
                           'The declared gradient uncertainty interval includes zero; coefficient sign/magnitude is not resolved.',
                           applications=apps, metrics={'sigma_multiplier': k})
    finding.evidence = dict(spec)


# profile units a supplied covariance may be stated in, as factors to the eV the screening works in
_TEMPERATURE_FACTOR = {'EV': 1.0, 'KEV': 1.0e3}
MAX_COVARIANCE_BYTES = 64 * 1024 * 1024


def _unit_factor(units, where):
    key = unit_key(units)
    if key not in _TEMPERATURE_FACTOR:
        raise ValueError(f'{where}: temperature units {units!r} are not recognized (EV or KEV)')
    return _TEMPERATURE_FACTOR[key]


def check_gradient_covariance(ctx, rec):
    """Propagate a supplied, fixed-grid profile covariance through the exact ``numpy.gradient`` operator.

    Scope: profile/gradient uncertainty *conditional on the supplied reconstruction and grids*.  The covariance is
    per time slice (cross-time correlation is unmodeled), the radial coordinate is taken as exact, and only the
    sources named in ``included_sources`` are claimed to be in it.  Nothing is repaired or regularized unless the
    resource allows a reported regularization; missing uncertainty stays unknown (no finding is invented).
    """
    k = ctx.thresholds['gradient_sigma_multiplier']
    for spec in ctx.manifest['uncertainty']['covariances']:
        name = spec['quantity']
        identifier = f'FIT.{name}.GRADIENT_COVARIANCE'
        apps = ROLE_APPS['electron_temperature' if name == 'TE' else 'ion_temperature']
        try:
            if not ctx.source.has(name) or not ctx.source.has('X'):
                raise ValueError(f'{name} or X is absent from the source; the covariance has nothing to propagate')
            profile, x = ctx.source.array(name), ctx.source.array('X')
            dims = (tuple(ctx.source.meta(name)['dimensions']), tuple(ctx.source.meta('X')['dimensions']))
            if dims != (('TIME3', 'X'), ('TIME3', 'X')) or profile.shape != x.shape:
                raise ValueError(f'{name} and X must share the native (TIME3, X) grid; no other grid is substituted')
            with np.load(spec.get('resolved_path', spec['path']), allow_pickle=False) as resource:
                time, covariance = resource['time'], resource['covariance']
            if covariance.nbytes > MAX_COVARIANCE_BYTES:
                raise ValueError(f'covariance resource {spec["id"]} exceeds the {MAX_COVARIANCE_BYTES} byte limit')
            if not np.array_equal(time, ctx.source.array('TIME3')):
                raise ValueError(f'covariance {spec["id"]}: its time axis is not the native TIME3 of the source '
                                 '(no interpolation or nearest-sample matching is done)')
            if covariance.shape != (profile.shape[0], profile.shape[1], profile.shape[1]):
                raise ValueError(f'covariance {spec["id"]}: shape {covariance.shape} is not (nt, nr, nr) = '
                                 f'{(profile.shape[0], profile.shape[1], profile.shape[1])}')
            native = _unit_factor(ctx.source.meta(name)['units'], name)
            stated = _unit_factor(spec['units'], f'covariance {spec["id"]}')
            factor = stated / native                                      # covariance units -> the profile's own units
            regularization = spec.get('max_regularization')
            specs = [uncertainty.CovarianceSpec(
                name, ctx.source.meta(name)['units'], x[t], covariance[t] * factor ** 2, uncertainty.ORDERING,
                uncertainty.ESTIMATOR, spec['provenance'], tuple(spec.get('included_sources', ())))
                for t in range(profile.shape[0])]
            result = uncertainty.propagate_gradient(specs, profile, x, 1.0)
        except (ValueError, NotImplementedError, OSError, KeyError) as exc:
            rec.report.add(identifier, 'error', f'Covariance could not be propagated: {exc}', variables=[name],
                           applications=list(apps))
            continue
        selected = rec.selection(name)
        slices = np.asarray(result.evaluated, dtype=bool)                  # per point: a slice without a covariance is False
        sigma = np.asarray(result.std_gradient, dtype=float)
        gradient = np.asarray(result.mean_gradient, dtype=float)
        usable = selected & slices & np.isfinite(sigma) & np.isfinite(gradient)
        finding = rec.samples(identifier, name, usable, usable & (np.abs(gradient) <= k * sigma),
                              'The gradient uncertainty propagated from the supplied covariance includes zero; the '
                              'sign and magnitude of a coefficient built on this gradient are not resolved. This '
                              'is conditional on the supplied reconstruction and grid.',
                              applications=apps, metrics={'sigma_multiplier': k, 'scope': result.scope,
                                                          'cross_time': result.cross_time})
        finding.evidence = {'resource': {'id': spec['id'], 'path': spec['path'], 'sha256': spec['sha256']},
                            'provenance': spec['provenance'], 'included_sources': list(spec.get('included_sources', [])),
                            'units_declared': spec['units'], 'scope': result.scope, 'cross_time': result.cross_time,
                            'estimator': uncertainty.ESTIMATOR, 'regularization_allowed': regularization,
                            'validation': result.validation[:12],
                            'not_included': ['radial-coordinate uncertainty', 'cross-time correlation',
                                             'any uncertainty source not listed in included_sources']}
        finite = sigma[usable]
        rec.report.uncertainty[spec['id']] = {
            'scope': result.scope, 'quantity': name, 'cross_time': result.cross_time, 'slices_evaluated': int(slices.any(axis=1).sum()),
            'slices_total': int(profile.shape[0]),
            'std_gradient_range': [float(finite.min()), float(finite.max())] if finite.size else None}
        rec.report.coordinates[f'uncertainty_{spec["id"]}__std_gradient'] = sigma
