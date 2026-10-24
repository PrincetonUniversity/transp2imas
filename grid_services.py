#!/usr/bin/python
"""TRANSP time bases and radial/time utilities, as functions over explicit arguments.

TRANSP stores scalar histories on ``TIME`` and profiles on ``TIME3``.  These helpers read the two time
bases, move a 1-D history between them, interpolate in time, and handle per-zone radial quantities.
``nc2imas`` keeps the time bases as state (``timesc``, ``timep``) and exposes these under their
historical method names.
"""

import numpy as np


def read_time_bases(get_ncvar, logger):
    """``(timesc, ntimes, timep, ntimep)`` from the ``TIME`` and ``TIME3`` variables; at least two of each."""
    timesc = get_ncvar('TIME')
    ntimes = len(timesc)
    timep = get_ncvar('TIME3')
    ntimep = len(timep)
    if ntimes>1 and ntimep>1:
        coords_x = get_ncvar('X')
        sh = coords_x.shape[1] if coords_x.ndim>1 else coords_x.shape
        logger.info(f' nprof {sh}: ntimes {ntimes}'
           f' {timesc[0]} .. {timesc[-1]}')
        logger.info(f' ntimep {ntimep} {timep[0]} .. {timep[-1]}')
        logger.info(f' Dtimep {timep[1]-timep[0]};'
           f' Dtimesc {timesc[1]-timesc[0]}')
    else:
        raise Exception(f'\tInvalid times prof:{ntimep} scalar {ntimes}')
    return timesc, ntimes, timep, ntimep


def interp_time(is_valid, time_new, time, dat, ind1:int=None, ind2:int=None):
    r'''
    interpolate dat from time onto time_new
    time_new: 1d times; output dat has this size.
    dat: 2d (time, 1d)
    time: 1d times (x axis); Same size of 1st dim or values or f(x) of dat
    '''
    #if dat is 2D, 1st dim should have len of len(time)
    if is_valid(time_new) and is_valid(dat, time, ind=len(time)-1):
       if isinstance(dat, np.ndarray) and dat.ndim==2:
           if ind1 is not None and ind2 is None:
              dat = dat[ind1, :]
           if ind1 is None and ind2 is not None:
              dat = dat[:, ind2]
       # TIME and TIME3 may have *equal lengths but different sample
       # locations*.  Comparing lengths is not sufficient to skip the
       # interpolation; compare the actual coordinates instead.
       if np.array_equal(np.asarray(time_new), np.asarray(time)):
           return dat
       if np.asarray(dat).ndim != 1 or len(dat) != len(time):
           return []
       return np.interp(time_new, time, dat)
    else:
       return []


def history_on_time(nc, get_ncvar, timesc, timep, var: str, target: str = 'TIME3', scale: float = 1.0):
    """Read a 1-D TRANSP time history on an explicitly selected time grid.

    The CDF variable's declared leading dimension, not the array length,
    identifies the source grid.  Interpolation is needed whenever the
    numerical TIME and TIME3 coordinates differ, even at equal lengths.
    This routine is deliberately restricted to 1-D histories; spatial
    profiles must be handled on their native TIME3 grid.
    """
    if target not in ('TIME', 'TIME3'):
        raise ValueError(f'Unsupported TRANSP target time coordinate: {target}')
    try:
        dimensions = nc.variables[var].dimensions
    except KeyError:
        return np.array([])
    if len(dimensions) != 1 or dimensions[0] not in ('TIME', 'TIME3'):
        raise ValueError(f'{var} is not a 1-D TIME/TIME3 history: {dimensions}')
    source = np.asarray(timesc if dimensions[0] == 'TIME' else timep)
    destination = np.asarray(timesc if target == 'TIME' else timep)
    values = get_ncvar(var, scale=scale)
    if values.ndim != 1 or values.shape != source.shape:
        raise ValueError(f'{var} shape {values.shape} disagrees with {dimensions[0]} '
                         f'coordinate {source.shape}')
    if not np.all(np.isfinite(source)) or not np.all(np.diff(source) > 0):
        raise ValueError(f'{var} has an invalid {dimensions[0]} time coordinate')
    if not np.all(np.isfinite(destination)) or not np.all(np.diff(destination) > 0):
        raise ValueError(f'{var} has an invalid {target} target time coordinate')
    if np.array_equal(source, destination):
        return values
    return np.interp(destination, source, values)


def cumulative_zone_quantity(values):
    """Convert TRANSP per-zone geometry to an enclosed flux-surface quantity."""
    arr = np.asarray(values)
    if arr.ndim < 1 or arr.size == 0:
        return arr
    return np.cumsum(arr, axis=-1)


def radial_gradient(values, coordinate):
    """Return d(values)/d(coordinate) independently for each time slice.

    TRANSP profiles can contain repeated/degenerate points near the axis in
    some runs.  Use a conservative one-sided/central finite difference and
    leave a slice empty (NaN) only when the coordinate is unusable.
    """
    y = np.asarray(values, dtype=np.float64)
    x = np.asarray(coordinate, dtype=np.float64)
    if y.ndim != 2 or x.ndim != 2 or y.shape != x.shape:
        return np.array([])
    out = np.full_like(y, np.nan, dtype=np.float64)
    for it in range(y.shape[0]):
        xx = x[it]
        yy = y[it]
        good = np.isfinite(xx) & np.isfinite(yy)
        if np.count_nonzero(good) < 2:
            continue
        xxg = xx[good]
        yyg = yy[good]
        # np.gradient requires a strictly monotonic coordinate.  A tiny
        # duplicate can occur through interpolation/roundoff, so retain the
        # first occurrence and interpolate the derivative back if needed.
        keep = np.concatenate(([True], np.diff(xxg) != 0.0))
        xxu = xxg[keep]
        yyu = yyg[keep]
        if len(xxu) < 2:
            continue
        edge_order = 2 if len(xxu) > 2 else 1
        grad = np.gradient(yyu, xxu, edge_order=edge_order)
        out[it] = np.interp(xx, xxu, grad, left=grad[0], right=grad[-1])
    return out
