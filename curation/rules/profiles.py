#!/usr/bin/python
"""Density and temperature profile screening on their native (TIME3, X) grid."""

import numpy as np

from curation.results import ROLE_APPS
from curation.rules.gradients import check_gradient
from curation.source import density_scale, unit_key

PROFILES = {'NE': 'electron_density', 'TE': 'electron_temperature', 'TI': 'ion_temperature'}

def check_profiles(ctx, rec):
    for name, role in PROFILES.items():
        apps = ROLE_APPS[role]
        if not ctx.source.has(name):
            rec.unavailable('PROFILE.' + name, 'Profile is absent; no substitute is assumed.', [name], apps)
            continue
        try:
            data = ctx.source.array(name)
            if data.ndim != 2 or tuple(ctx.source.meta(name)['dimensions']) != ('TIME3', 'X'):
                raise ValueError('Expected a native (TIME3, X) profile; no shape-based reinterpretation')
            if ctx.source.has('X') and data.shape != ctx.source.array('X').shape:
                rec.report.add('PROFILE.' + name + '.GRID_MISMATCH', 'fail',
                                'Profile shape disagrees with the explicit native X grid.',
                                variables=[name, 'X'], applications=list(apps))
                continue
            selected = rec.selection(name)
            rec.samples('PROFILE.' + name + '.FINITE', name, selected, selected & ~np.isfinite(data),
                          'Missing, masked or nonfinite profile samples require review.', applications=apps)
            finite = selected & np.isfinite(data)
            rec.samples('PROFILE.' + name + '.NEGATIVE', name, finite, finite & (data < 0),
                          'Negative density/temperature samples are not physically admissible.',
                          severity='fail', applications=apps)
            rec.samples('PROFILE.' + name + '.ZERO', name, finite, finite & (data == 0),
                          'Zero values may be startup/off-plasma samples; do not use them as transport denominators.',
                          applications=apps)
            units = ctx.source.meta(name)['units']
            if name == 'NE':
                density_scale(units)
            elif unit_key(units) not in ('EV', 'KEV'):
                raise ValueError(f'{name}: unknown temperature units {units!r}')
            rec.report.add('PROFILE.' + name + '.UNITS', 'pass', 'Native units are recognized.',
                            variables=[name], applications=list(apps), evidence={'units': units})
            if name in ('TE', 'TI'):
                check_gradient(ctx, rec, name)
        except (ValueError, TypeError, OSError, RuntimeError) as exc:
            rec.source_problem('PROFILE.' + name + '.ASSESSMENT', exc, [name], apps)
