#!/usr/bin/python
"""Charge closure of an explicitly complete, grid- and unit-compatible species inventory."""

import numpy as np

from curation.source import density_scale, unit_key



def check_composition(ctx, rec):
    ident = 'COMPOSITION.QUASINEUTRALITY'
    comp = ctx.manifest.get('composition')
    if not comp or comp.get('complete') is not True or comp.get('includes_fast_ions') is not True or not comp.get('species'):
        rec.unavailable(ident, 'An explicit complete charged-species inventory, including fast-ion accounting, is required.')
        return
    ne_name = comp.get('electron_density', 'NE')
    names = [ne_name] + [s['density'] for s in comp['species']]
    names += [s['charge_variable'] for s in comp['species'] if 'charge_variable' in s]
    missing = [name for name in names if not ctx.source.has(name)]
    if missing:
        rec.unavailable(ident, 'Declared composition variables are missing: ' + ', '.join(missing), missing)
        return
    try:
        ne = ctx.source.array(ne_name) * density_scale(ctx.source.meta(ne_name)['units'])
        dims = ctx.source.meta(ne_name)['dimensions']
        if ne.ndim != 2 or tuple(dims) != ('TIME3', 'X'):
            raise ValueError('Electron density must be on the native (TIME3, X) grid')
        charge_density = np.zeros_like(ne)
        valid = np.isfinite(ne) & (ne > 0)
        for species in comp['species']:
            name = species['density']
            density = ctx.source.array(name) * density_scale(ctx.source.meta(name)['units'])
            if density.shape != ne.shape or ctx.source.meta(name)['dimensions'] != dims:
                raise ValueError(f'{name}: native composition grids do not match; no interpolation is assumed')
            charge = species.get('charge')
            if charge is None:
                charge_name = species['charge_variable']
                charge = ctx.source.array(charge_name)
                if charge.shape != ne.shape or ctx.source.meta(charge_name)['dimensions'] != dims:
                    raise ValueError(f'{charge_name}: charge-state grid mismatch')
                if unit_key(ctx.source.meta(charge_name)['units']) not in ('', '1', '-', 'DIMENSIONLESS'):
                    raise ValueError(f'{charge_name}: charge-state units are not dimensionless')
            valid &= np.isfinite(density) & (density >= 0) & np.isfinite(charge) & (charge > 0)
            charge_density += density * charge
        selected = rec.selection(ne_name)
        rec.samples(ident + '.INPUTS', ne_name, selected, selected & ~valid,
                      'Some composition samples are unavailable/nonphysical and cannot be tested for charge closure.')
        valid &= selected
        relative = np.full(ne.shape, np.nan)
        np.divide(np.abs(charge_density - ne), ne, out=relative, where=valid)
        tolerance = ctx.thresholds['quasineutrality_relative']
        metrics = {'relative_tolerance': tolerance, 'tolerance_kind': 'screening_not_statistical'}
        if valid.any():
            metrics['max_relative_residual'] = float(np.max(relative[valid]))
        finding = rec.samples(ident, ne_name, valid, valid & (relative > tolerance),
                               'Charge closure tested for the declared inventory; closure does not validate composition or Zeff.',
                               metrics=metrics)
        finding.variables = names
        finding.evidence = {'source': comp.get('source', ''), 'inventory_is_analyst_declared': True}
    except (ValueError, TypeError, OSError, RuntimeError) as exc:
        rec.source_problem(ident, exc, names)
