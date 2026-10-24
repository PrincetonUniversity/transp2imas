#!/usr/bin/python
"""Validity checks, array arithmetic and assignment of source data to IDS nodes.

These are the primitives every writer uses: ``is_valid`` (is this optional source datum present and
non-empty), ``sum_arrays`` / ``operate_arrays`` / ``np_multiply_sum`` (arithmetic that reports bad
input instead of raising), and ``fill_data`` / ``fill_value`` (assign to an IDS node).  They take the
logger and the ``MappingLedger`` explicitly, so they run without a converter, CDF or database; the
``nc2imas`` methods of the same names pass the converter's own.

Policy (see ``mapping_status``): absent optional source data is a recorded skip; data that is present
but cannot be assigned is a broken mapping, raised in strict mode and recorded as a skipped field with
the reason in lenient mode.
"""

import math

import numpy as np


def is_valid(*variables, ind=0, name='', logger=None):
    """True when every argument is present and non-empty.

    ``None``, a 0-d array, an array with fewer than ``ind + 1`` rows or an empty second axis, and a
    sequence that is too short or all falsy are invalid.  Do not call this with varargs forwarded by
    another function.  With ``ind == 1`` the reason is logged.
    """
    ii = ind + 1

    def note(message):
        if ind == 1 and logger is not None:
            logger.info(message)

    for var in variables:
        if var is None:
            note(f'{name} is None')
            return False
        elif isinstance(var, np.ndarray):
            if var.ndim == 0:                       # e.g. np.array(5)
                return False
            shape = var.shape
            if shape[0] < ii or (len(shape) > 1 and shape[1] == 0):
                note(f'{name} shape {shape}')
                return False
        elif isinstance(var, (list, tuple, set)):
            if len(var) < ii or not any(var):
                note(f'{name} len: {len(var)}')
                return False
    return True


def check_1dim(dat):
    """True for a 1-D ndarray or a flat sequence; ``None`` for any other type (as before)."""
    if isinstance(dat, np.ndarray):
        return dat.ndim == 1
    elif isinstance(dat, (list, tuple, set)):
        return np.array(dat).ndim == 1


def sum_arrays(*arrays, logger=None):
    """Element-wise sum of the non-empty array-like inputs.

    Shape and type errors are treated as invalid optional source data (logged); programming and
    process-control exceptions propagate.
    """
    result, started = [], False
    try:
        for arr in arrays:
            if len(arr) == 0:
                continue
            if started:
                result = result + arr
            else:
                result = arr
                started = True
    except (TypeError, ValueError) as exc:
        if logger is not None:
            logger.info(f' sum arrays error: {exc}')
    return result


def operate_arrays(arr1, op, arr2=None, name1='', name2='', logger=None):
    """Apply ``+ - * / sqrt`` element-wise with a local NumPy error policy; empty array on bad input.

    Division by zero gives 0 (``where=arr2 != 0``), every other floating-point error is caught.
    """
    def note(message):
        if logger is not None:
            logger.info(message)

    unary = op == 'sqrt'
    valid = is_valid(arr1, logger=logger) if unary else is_valid(arr1, arr2, logger=logger)
    if not valid:
        note(f' operate_arrays error: {name1} {op} {name2}')
        return np.array([])
    try:
        with np.errstate(divide='raise', invalid='raise', over='raise'):
            if op == '-':
                return arr1 - arr2
            if op == '+':
                return arr1 + arr2
            if op == '/':
                a1, a2 = np.asarray(arr1), np.asarray(arr2)
                dtype = np.result_type(a1.dtype, a2.dtype, np.float64)
                shape = np.broadcast_shapes(a1.shape, a2.shape)
                return np.divide(a1, a2, out=np.zeros(shape, dtype=dtype), where=a2 != 0)
            if op == '*':
                return arr1 * arr2
            if op == 'sqrt':
                if isinstance(arr1, np.ndarray):
                    return np.sqrt(arr1)
                if isinstance(arr1, list):
                    return [math.sqrt(x) for x in arr1]
                if isinstance(arr1, set):
                    return {math.sqrt(x) for x in arr1}
                if isinstance(arr1, tuple):
                    return tuple(math.sqrt(x) for x in arr1)
    except (TypeError, ValueError, FloatingPointError, OverflowError, ZeroDivisionError) as exc:
        note(f' operate_arrays error:: {name1} {op} {name2}: {exc}')
        return np.array([])
    note(f' operate_arrays error:: unsupported operation {op!r}')
    return np.array([])


def np_multiply_sum(arr1, arr2, axis, name1='', name2='', scale=None, logger=None):
    """``sum(arr1 * arr2, axis)`` times ``scale``; empty array (logged) when an input is invalid."""
    if not is_valid(arr1, arr2, logger=logger):
        if logger is not None:
            s1, s2 = len(arr1), len(arr2)
            if isinstance(arr1, np.ndarray) and isinstance(arr2, np.ndarray):
                s1, s2 = arr1.shape, arr2.shape
            logger.info(f' np_multiply_sum error: {name1}:{s1} {name2}: {s2}')
        return np.array([])
    try:
        with np.errstate(invalid='raise', over='raise'):
            result = np.sum(np.multiply(arr1, arr2), axis=axis)
            return result if scale is None else scale * result
    except (TypeError, ValueError, FloatingPointError, OverflowError) as exc:
        if logger is not None:
            logger.info(f' np_multiply_sum error:: {name1} {name2}: {exc}')
        return np.array([])


def node_path(node):
    """IDS path of ``node`` for status messages (never raises)."""
    path = getattr(getattr(node, 'metadata', None), 'path', None)
    return str(path) if path is not None else type(node).__name__


def fill_value(node, dat, rsz=-1, mult=None, logger=None):
    """Assign ``dat`` to ``node.value`` (a scalar, or an array resized to ``len(dat)`` or ``rsz``)."""
    if not is_valid(dat, logger=logger):
        return
    nresize = len(dat) if (not np.isscalar(dat) and rsz == -1) else rsz
    if nresize > 0:
        node.value.resize(nresize, refcheck=False)
    if mult is not None:
        dat = mult * dat
    if np.isscalar(dat):
        node.value = dat
    else:
        node.value[:] = dat


def fill_data(node, dat, ledger, logger, mult=None, ind=-1, ind2=-1, interp=False, times=None,
              name=None, resize=True, chk=-1):
    """Assign ``dat`` to an array node; returns True if it was written.

    ``ind`` selects a time slice (``ind`` to ``ind2`` a range) of 2-D data; ``chk`` is the radial size
    the data must have; ``interp`` maps the data from the scalar to the profile time base with
    ``times() -> (timesc, timep)``.

    Absent source data (``is_valid`` false) is optional input: the field is skipped and recorded in
    ``ledger``.  Data that is present but cannot be assigned (radial size different from ``chk``, a
    scalar for an array node, a shape or schema error) is a broken mapping, not absent data: it raises
    ``MappingError`` (strict) or is recorded as skipped with the reason (lenient).  The size mismatch
    is reported for every time slice, not only slice 0.
    """
    nm = name if name is not None else 'data'
    target = node_path(node)
    if not is_valid(dat, name=nm, logger=logger):
        ledger.skipped('no source data', target=target)
        if ind == 0:
            logger.info(f'\t\t{nm} is invalid')
        return False
    try:
        if ind >= 0 and ind2 >= 0:
            dat = dat[ind:ind2]
        elif ind >= 0:
            dat = dat[ind]

        if chk > 0 and not np.isscalar(dat) and len(dat) != chk:
            raise ValueError(f'size mismatch {len(dat)} != {chk}')
        if np.isscalar(dat):
            raise ValueError('cannot fill a scalar value into an array node')
        if mult is not None:
            dat = mult * dat
        if interp is True:
            timesc, timep = times()
            if not np.array_equal(timesc, timep):
                dat = np.interp(timep, timesc, dat)
        if (not np.isscalar(dat)) and resize:
            node.resize(len(dat), refcheck=False)
        node[:] = dat
    except (AttributeError, ValueError, TypeError, IndexError, FloatingPointError) as exc:
        ledger.broken(f'{type(exc).__name__}: {exc}', source=nm, target=target)
        return False
    return True
