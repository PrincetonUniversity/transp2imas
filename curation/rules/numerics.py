#!/usr/bin/python
"""Runtime solver channels on their native grids; no convergence inference."""

import numpy as np

SOLVER_CHANNELS = ('NPTRITER', 'PT_CPTIME', 'PT_WTIME', 'RESL2TE', 'RESL2TI', 'RESL2OMG',
                   'UPWIND_TE', 'UPWIND_TI', 'UPWIND_MO')

def check_numerics(ctx, rec):
    present = []
    for name in SOLVER_CHANNELS:
        if not ctx.source.has(name):
            continue
        try:
            values = ctx.source.array(name)
            selected = rec.selection(name)
            finite = selected & np.isfinite(values)
            entry = ctx.source.meta(name)
            entry.update({'selected_samples': int(selected.sum()), 'finite_samples': int(finite.sum()),
                          'nonzero_samples': int(np.count_nonzero(values[finite]))})
            if finite.any():
                entry.update({'minimum': float(values[finite].min()), 'maximum': float(values[finite].max())})
            rec.report.channel_inventory[name] = entry
            present.append(name)
            rec.samples('NUMERICS.' + name + '.FINITE', name, selected, selected & ~np.isfinite(values),
                          'Runtime diagnostic contains unavailable/nonfinite samples; no convergence inference.')
            if name == 'NPTRITER':
                invalid = (values < 0) | (np.abs(values - np.round(values)) > 1e-6)
                rec.samples('NUMERICS.ITERATION_DOMAIN', name, finite, finite & invalid,
                              'Iteration counts should be nonnegative integers; zero is not proof of convergence.', severity='fail')
                limit = ctx.thresholds.get('newton_iterations_warning')
                if limit is not None:
                    rec.samples('NUMERICS.ITERATION_LIMIT', name, finite, finite & (values > limit),
                                  'Iteration count exceeds the analyst screening threshold; this is not the solver stopping criterion.',
                                  metrics={'threshold': limit})
            elif name.startswith('RESL2'):
                rec.samples('NUMERICS.' + name + '.DOMAIN', name, finite, finite & (values < 0),
                              'An L2 residual should not be negative.', severity='fail')
                spec = ctx.thresholds['residual_limits'].get(name)
                if spec:
                    finding = rec.samples('NUMERICS.' + name + '.THRESHOLD', name, finite,
                                            finite & (values > spec['limit']),
                                            'Residual exceeds the supplied native-normalization threshold; check its rationale.',
                                            metrics={'threshold': spec['limit']})
                    finding.evidence = dict(spec)
        except (ValueError, TypeError, OSError, RuntimeError) as exc:
            rec.source_problem('NUMERICS.' + name + '.READ', exc, [name])
    rec.report.add('NUMERICS.AVAILABILITY', 'pass' if present else 'not_assessed',
                    'Runtime channels inventoried on their own coordinates; upwind adjustment is not a numerical heat flux.',
                    variables=present)
    rec.unavailable('NUMERICS.CONVERGENCE',
                      'Per-equation stopping criteria, rejected steps and time-step/stabilization sensitivity are not certified by this first-stage screen.')
