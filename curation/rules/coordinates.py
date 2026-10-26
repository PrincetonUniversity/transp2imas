#!/usr/bin/python
"""Native time and radial coordinate screening: source integrity, mandatory."""

import numpy as np

from curation.source import unit_key



def check_coordinates(ctx, rec):
    for name in ('TIME', 'TIME3'):
        identifier = 'GRID.' + name
        if not ctx.source.has(name):
            rec.unavailable(identifier, 'Native time coordinate is absent.', [name])
            continue
        try:
            values = ctx.source.array(name)
            units = unit_key(ctx.source.meta(name)['units'])
            good = (values.ndim == 1 and values.size > 0 and np.all(np.isfinite(values))
                    and np.all(np.diff(values) > 0))
            rec.report.add(identifier, 'pass' if good else 'fail',
                            'Native time coordinate must be finite and strictly increasing.', variables=[name],
                            metrics={'sample_count': int(values.size)},
                            evidence={'units': ctx.source.meta(name)['units']})
            if units not in ('S', 'SEC', 'SECS', 'SECOND', 'SECONDS'):
                rec.unavailable(identifier + '.UNITS', 'Time units are not recognized as seconds.', [name])
        except (ValueError, TypeError, OSError, RuntimeError) as exc:
            rec.source_problem(identifier, exc, [name])
    if not ctx.source.has('X'):
        rec.unavailable('GRID.X', 'No zone-centre normalized-radius grid.', ['X'])
        return
    try:
        x = ctx.source.array('X')
        if x.ndim != 2 or tuple(ctx.source.meta('X')['dimensions']) != ('TIME3', 'X'):
            raise ValueError('X is not a native (TIME3, X) grid')
        if not ctx.normalized_radius('X'):
            raise ValueError('X units do not establish a normalized radius')
        selected = rec.selection('X')
        invalid = ~np.isfinite(x) | (x < -1e-6) | (x > 1 + 1e-6)
        invalid[:, 1:] |= np.diff(x, axis=1) <= 0
        rec.samples('GRID.X', 'X', selected, selected & invalid,
                      'Normalized X must be finite, increasing and within [0, 1].', severity='fail')
        rec.report.coordinates['coordinate_X'] = x
    except (ValueError, TypeError, OSError, RuntimeError) as exc:
        rec.source_problem('GRID.X', exc, ['X'])
