#!/usr/bin/python
"""Radial profile interpolation with a cache of the interpolation operator.

``interpolate_radial_profiles(values, x_from, x_to, axis_value)`` maps time-dependent profiles from
one radial grid to another, one time slice at a time, exactly as ``np.interp`` would
(clamping outside the source grid).  Writers repeat the *same* grid mapping for many profiles (every
impurity species and charge state of the radiation IDS, for example); the per-slice ``np.interp`` loop
was the largest cost of building that IDS.  The positions and weights depend only on the two
coordinate arrays and on whether an axis point is prepended, so they are computed once and applied to
any number of profiles with array operations.

The cache key is the *content* of the coordinates (a digest of their bytes), never their shape or
identity, so two different grids of the same size cannot share an operator and a modified array
cannot return a stale one.  The cache is bounded.  Grids that are not strictly increasing in every
slice use the plain ``np.interp`` loop, so behaviour on unusual input is unchanged.

No IMAS or NetCDF dependency.
"""

from collections import OrderedDict
import hashlib

import numpy as np

_CACHE = OrderedDict()
_CACHE_SIZE = 8
STATS = {'hits': 0, 'misses': 0, 'fallbacks': 0}


def clear_cache():
    _CACHE.clear()
    for key in STATS:
        STATS[key] = 0


def _digest(*arrays):
    h = hashlib.blake2b(digest_size=16)
    for a in arrays:
        a = np.ascontiguousarray(a, dtype=np.float64)
        h.update(str(a.shape).encode())
        h.update(a.tobytes())
    return h.digest()


def _reference(values, x_from, x_to, axis_value):
    """The original per-slice np.interp loop (also the fallback)."""
    v, xf, xt = np.asarray(values), np.asarray(x_from), np.asarray(x_to)
    out = np.empty((v.shape[0], xt.shape[1]), dtype=np.float64)
    av = None if axis_value is None else np.asarray(axis_value)
    for it in range(v.shape[0]):
        xx, yy = xf[it], v[it]
        if axis_value is not None and len(xx) and xx[0] > 0.0:
            a0 = float(av[it] if av.ndim else av)
            xx = np.concatenate(([0.0], xx))
            yy = np.concatenate(([a0], yy))
        out[it] = np.interp(xt[it], xx, yy)
    return out


class _Operator:
    """Positions and weights of one grid mapping, for ``use_axis`` True or False."""

    def __init__(self, x_from, x_to, use_axis):
        xf = np.asarray(x_from, dtype=np.float64)
        xt = np.asarray(x_to, dtype=np.float64)
        nt, n = xf.shape
        self.use_axis = use_axis
        self.valid = False
        if not (np.all(np.isfinite(xf)) and np.all(np.isfinite(xt))):
            return                                   # NaN coordinates: keep np.interp's behaviour
        # rows that get a prepended axis point (0, value at the axis) when an axis value is given
        self.prepend = (xf[:, 0] > 0.0) if use_axis else np.zeros(nt, dtype=bool)
        if use_axis:
            # every row gets one extra leading node; rows without an axis point repeat their first
            # node (same position, same value), which does not change the interpolant
            lead = np.where(self.prepend, 0.0, xf[:, 0])
            xx = np.concatenate([lead[:, None], xf], axis=1)
        else:
            xx = xf
        valid = bool(np.all(np.diff(xx, axis=1) > 0.0)) if not use_axis else bool(
            np.all(np.diff(xx[:, 1:], axis=1) > 0.0) and np.all(np.diff(xx[~self.prepend], axis=1) >= 0.0))
        if not valid or xx.shape[1] < 2:
            return
        self.valid = True
        m = xx.shape[1]
        rows = np.arange(nt)[:, None]
        # number of nodes <= target, row by row (chunked so memory stays bounded)
        count = np.empty(xt.shape, dtype=np.int64)
        chunk = max(1, int(4e6 // max(1, xt.shape[1] * m)))
        for start in range(0, nt, chunk):
            sl = slice(start, start + chunk)
            count[sl] = (xt[sl, :, None] >= xx[sl, None, :]).sum(axis=2)
        lo = np.clip(count - 1, 0, m - 2)
        hi = lo + 1
        x_lo, x_hi = xx[rows, lo], xx[rows, hi]
        span = x_hi - x_lo
        # weight of the upper node; clamped to [0, 1] outside the grid like np.interp
        with np.errstate(divide='ignore', invalid='ignore'):
            w = np.where(span > 0.0, (xt - x_lo) / np.where(span > 0.0, span, 1.0), 0.0)
        self.w = np.clip(w, 0.0, 1.0)
        self.lo, self.hi, self.rows = lo, hi, rows

    def apply(self, values, axis_value):
        v = np.asarray(values, dtype=np.float64)
        if self.use_axis:
            av = np.asarray(axis_value, dtype=np.float64)
            a0 = np.broadcast_to(av, (v.shape[0],)) if av.ndim == 0 else av.reshape(-1)[:v.shape[0]]
            lead = np.where(self.prepend, a0, v[:, 0])
            v = np.concatenate([lead[:, None], v], axis=1)
        v_lo, v_hi = v[self.rows, self.lo], v[self.rows, self.hi]
        return v_lo + (v_hi - v_lo) * self.w


def interpolate_radial_profiles(values, x_from, x_to, axis_value=None):
    """Interpolate ``values`` (time x radius) from ``x_from`` to ``x_to`` slice by slice.

    Semantics of the original ``nc2imas.interp_radial_profiles``: inputs that are not 2-D with a
    common time length are returned unchanged; with ``axis_value`` a point (0, axis value) is added
    in front of every slice whose first node is positive.
    """
    v, xf, xt = np.asarray(values), np.asarray(x_from), np.asarray(x_to)
    if v.size == 0 or xf.size == 0 or xt.size == 0 or v.ndim != 2:
        return v
    if xf.ndim != 2 or xt.ndim != 2 or v.shape[0] != xf.shape[0] or v.shape[0] != xt.shape[0]:
        return v
    if v.shape != xf.shape:
        return _reference(v, xf, xt, axis_value)        # ragged request: keep the original behaviour
    use_axis = axis_value is not None
    key = (_digest(xf, xt), use_axis)
    operator = _CACHE.get(key)
    if operator is None:
        STATS['misses'] += 1
        operator = _Operator(xf, xt, use_axis)
        _CACHE[key] = operator
        while len(_CACHE) > _CACHE_SIZE:
            _CACHE.popitem(last=False)
    else:
        STATS['hits'] += 1
        _CACHE.move_to_end(key)
    if not operator.valid or not np.all(np.isfinite(v)) or (
            use_axis and not np.all(np.isfinite(np.asarray(axis_value, dtype=np.float64)))):
        STATS['fallbacks'] += 1
        return _reference(v, xf, xt, axis_value)
    return operator.apply(v, axis_value)
