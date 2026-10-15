#!/usr/bin/python
"""Source-grounded TRANSP sawtooth crash-history -> IMAS diagnostics.

``TLASTSAW(TIME)`` is the most recent crash timestamp at each scalar output
sample, not the timestamp of that output sample.  A change of TLASTSAW
therefore identifies a newly *reported* crash, whereas the first observed
value is left-censored and is not counted as a newly detected event.

``XIQ1(TIME)`` contains the normalized xi coordinate of the q=1 surface;
negative values in the supplied CDFs are absence sentinels.  No guess is made
about ``XSAWBOUND`` ("Adhoc sawtooth boundary"), crash-trigger mechanism, or
unavailable pre/post-crash profiles.
"""

from dataclasses import dataclass

import numpy as np

# IMAS FLT missing-value convention.  Unlike a physical zero, this explicitly
# marks missing information at individual samples on the homogeneous time grid.
IMAS_MISSING_FLOAT = -9.0e40


@dataclass(frozen=True)
class SawteethPlan:
    """Native TIME-aligned diagnostics and confirmed new crash updates."""

    time: np.ndarray
    previous_crash_time: np.ndarray
    event_indices: np.ndarray
    event_times: np.ndarray
    previous_period: np.ndarray | None
    rho_tor_norm_q1: np.ndarray | None
    valid_q1_samples: int
    source_boundary_samples: int

    @property
    def events(self):
        return int(self.event_indices.size)


def _time_channel(nc, name, size):
    """Read a scalar-time channel without implicitly using TIME3."""
    variable = nc.variables.get(name)
    if variable is None:
        return None
    if tuple(variable.dimensions) != ('TIME',):
        raise ValueError(f'{name}: expected native TIME dimension, got {variable.dimensions}')
    data = np.asarray(variable[:], dtype=np.float64)
    if data.shape != (size,):
        raise ValueError(f'{name}: expected {size} samples, got {data.shape}')
    return data


def build_sawteeth_plan(nc):
    """Return a plan only when TLASTSAW reports a new crash in this CDF.

    The native TIME array is used for the IMAS homogeneous time coordinate.
    A constant initial TLASTSAW (notably ITER 38530P94) does NOT demonstrate
    a crash during the recorded time interval.
    """
    last = nc.variables.get('TLASTSAW')
    if last is None:
        return None
    tvar = nc.variables.get('TIME')
    if tvar is None:
        raise ValueError('TLASTSAW is present but native TIME is missing')
    time = np.asarray(tvar[:], dtype=np.float64)
    if time.ndim != 1 or not time.size or not np.all(np.isfinite(time)):
        raise ValueError('TIME must be a nonempty finite 1D vector')
    if np.any(np.diff(time) < 0):
        raise ValueError('TIME is not monotonically increasing')
    values = _time_channel(nc, 'TLASTSAW', time.size)
    valid = np.isfinite(values) & (values >= 0)
    # float32 timestamp roundoff can be ~1e-4 near ITER T=550 s.
    tolerances = np.maximum(2.0e-6, 4.0 * np.finfo(np.float32).eps * np.maximum(1.0, np.abs(time)))
    if np.any(valid & (values > time + tolerances)):
        raise ValueError('TLASTSAW contains a future crash timestamp')

    event_indices = []
    prior = None
    for i, (value, exists) in enumerate(zip(values, valid)):
        if not exists:
            continue
        tol = float(tolerances[i])
        if prior is None:
            prior = float(value)  # initial observed crash is NOT a new event
        elif value > prior + tol:
            event_indices.append(i)
            prior = float(value)
        elif value < prior - tol:
            raise ValueError('TLASTSAW decreases; cannot infer crashes across a reset')
    if not event_indices:
        return None
    indices = np.asarray(event_indices, dtype=np.int64)
    events = values[indices].copy()

    previous_crash = np.where(valid, values, IMAS_MISSING_FLOAT)

    # The previous period is known only after TWO observed new crash updates.
    # The initial TLASTSAW is an unverified pre-interval crash/initialization;
    # never infer the first observed period from that censored baseline.
    period = None
    if events.size >= 2:
        period = np.full(time.shape, IMAS_MISSING_FLOAT, dtype=np.float64)
        for j in range(1, len(indices)):
            stop = int(indices[j + 1]) if j + 1 < len(indices) else len(time)
            period[indices[j]:stop] = float(events[j] - events[j - 1])

    # XI (TRANSP normalized toroidal flux radius) maps to rho_tor_norm,
    # the same convention used in core_profiles.grid.rho_tor_norm.
    xiq1 = _time_channel(nc, 'XIQ1', len(time))
    valid_q1 = 0
    q1 = None
    if xiq1 is not None:
        valid_mask = np.isfinite(xiq1) & (xiq1 > 0) & (xiq1 <= 1)
        valid_q1 = int(np.count_nonzero(valid_mask))
        if valid_q1:
            q1 = np.where(valid_mask, xiq1, IMAS_MISSING_FLOAT)

    boundary = _time_channel(nc, 'XSAWBOUND', len(time))
    source_boundary_samples = (int(np.count_nonzero(np.isfinite(boundary) & (boundary > 0)))
                               if boundary is not None else 0)
    # Do NOT map this ambiguous 'adhoc boundary' into inversion/mixing radii.
    return SawteethPlan(time=time, previous_crash_time=previous_crash,
                        event_indices=indices, event_times=events,
                        previous_period=period, rho_tor_norm_q1=q1,
                        valid_q1_samples=valid_q1,
                        source_boundary_samples=source_boundary_samples)


def populate_sawteeth_ids(ids, plan):
    """Fill only validated native DD3/DD4 sawteeth diagnostic fields."""
    ids.ids_properties.homogeneous_time = 1
    ids.time = plan.time
    diag = ids.diagnostics
    diag.previous_crash_time = plan.previous_crash_time
    if plan.previous_period is not None:
        diag.previous_period = plan.previous_period
    if plan.rho_tor_norm_q1 is not None:
        diag.rho_tor_norm_q1 = plan.rho_tor_norm_q1
    # Intentionally absent: crash_trigger, previous_crash_trigger, mixing and
    # inversion radii, and post-crash profiles.  No source proves those values.
    return plan.events
