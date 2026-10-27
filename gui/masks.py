#!/usr/bin/python
"""Bounded, read-only views of the evaluated/flagged/excluded masks stored beside a curation report.

The ``.npz`` is opened with ``allow_pickle=False`` and only the named arrays are read.  A mask is reduced to at most
``MAX_CELLS`` per axis for drawing, taking the highest-priority category of each block, so a small flagged region is
never averaged away.  Categories (also the legend of the page):

    0 not evaluated (outside the window or without a native coordinate)
    1 evaluated, not flagged
    2 flagged
    3 excluded by a declared event, not flagged
    4 flagged and excluded
"""

import numpy as np

MAX_CELLS = 160
CATEGORIES = ('not evaluated', 'evaluated', 'flagged', 'excluded', 'flagged and excluded')
_PRIORITY = np.array([0, 1, 3, 2, 4])          # category -> drawing priority (flagged outranks excluded outranks ok)


def categories(evaluated, flagged, excluded=None):
    evaluated, flagged = np.asarray(evaluated, dtype=bool), np.asarray(flagged, dtype=bool)
    excluded = np.zeros_like(evaluated) if excluded is None else np.asarray(excluded, dtype=bool)
    code = np.zeros(evaluated.shape, dtype=np.int8)
    code[evaluated] = 1
    code[flagged] = 2
    code[excluded & evaluated & ~flagged] = 3
    code[excluded & flagged] = 4
    return code


def reduce_grid(code, limit=MAX_CELLS):
    """Block-reduce a 2-D category grid to at most ``limit`` cells per axis; returns ``(grid, (row_step, col_step))``."""
    code = np.asarray(code)
    if code.ndim == 1:
        code = code[None, :]
    if code.ndim != 2:
        raise ValueError(f'a mask must be 1-D or 2-D to be drawn, got {code.ndim}-D')
    steps = tuple(max(1, -(-n // limit)) for n in code.shape)
    rows, cols = (-(-n // s) for n, s in zip(code.shape, steps))
    out = np.zeros((rows, cols), dtype=np.int8)
    priority = _PRIORITY[code]
    for i in range(rows):
        for j in range(cols):
            block = priority[i * steps[0]:(i + 1) * steps[0], j * steps[1]:(j + 1) * steps[1]]
            best = int(block.max())
            out[i, j] = int(np.where(_PRIORITY == best)[0][0])
    return out, steps


def _axis(arrays, name):
    """A 1-D coordinate vector for drawing, or ``None`` (a 2-D radial grid is shown by its first row)."""
    values = arrays.get(name)
    if values is None:
        return None
    values = np.asarray(values, dtype=float)
    if values.ndim == 2:
        values = values[0]
    return [None if not np.isfinite(v) else float(v) for v in values]


def mask_view(npz_path, prefix, dimensions):
    """The reduced mask of finding ``prefix`` with its axes; reads only the arrays it needs."""
    with np.load(npz_path, allow_pickle=False) as store:
        names = set(store.files)
        if prefix + '__evaluated' not in names:
            raise KeyError(prefix)
        evaluated = store[prefix + '__evaluated']
        flagged = store[prefix + '__flagged'] if prefix + '__flagged' in names else np.zeros_like(evaluated)
        excluded = store[prefix + '__excluded'] if prefix + '__excluded' in names else None
        coordinates = {n: store[n] for n in names if n.startswith('coordinate_') and
                       n[len('coordinate_'):] in dimensions}
    code = categories(evaluated, flagged, excluded)
    grid, steps = reduce_grid(code)
    counts = {name: int((code == value).sum()) for value, name in enumerate(CATEGORIES)}
    axes = []
    for dim in dimensions[:2]:
        axes.append({'name': dim, 'values': _axis(coordinates, 'coordinate_' + dim)})
    return {'prefix': prefix, 'shape': list(code.shape), 'dimensions': list(dimensions), 'axes': axes,
            'grid': grid.tolist(), 'block': list(steps), 'categories': list(CATEGORIES), 'counts': counts}


def mask_prefixes(npz_path):
    """The mask prefixes available in a report's ``.npz`` (names only; nothing is loaded)."""
    with np.load(npz_path, allow_pickle=False) as store:
        return sorted(n[:-len('__evaluated')] for n in store.files if n.endswith('__evaluated'))
