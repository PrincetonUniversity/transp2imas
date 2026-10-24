#!/usr/bin/python
"""Contour lines of a 2-D field without creating a matplotlib Figure.

``Axes.contour`` computes its lines with ``contourpy`` and only then draws them.  The equilibrium writer
needs just the vertices (``contours.allsegs[0]``) of one level for every time slice; building and closing
a Figure each time was most of the writer's run time (NSTX 142202U01: 24 of 30 s in the boundary trace).
:func:`line_contours` calls the same ``contourpy`` generator with the settings matplotlib passes, so the
vertices are identical; ``verify_equilibrium_boundary_trace.py`` checks that against ``Axes.contour``.
"""

import numpy as np


DIRECT_FROM_MATPLOTLIB = (3, 10)      # the version whose ``Axes.contour`` this mirrors and is tested against


def line_contours(x, y, z, level):
    """Vertex arrays of the contour lines of ``z`` at ``level``, as ``Axes.contour(x, y, z, [level]).allsegs[0]``.

    ``x`` and ``y`` are the 1-D grid coordinates (or 2-D meshes) and ``z`` has shape ``(len(y), len(x))``.
    Like matplotlib, non-finite values of ``z`` are masked and the level is used as given: a level outside the
    range of the data has no lines, which matplotlib reports as one empty ``(0, 2)`` array.  Matplotlib older
    than 3.10 goes through ``Axes.contour`` itself, whose handling of such levels differs between versions.
    """
    try:
        import contourpy
        import matplotlib as mpl
        direct = tuple(mpl.__version_info__[:2]) >= DIRECT_FROM_MATPLOTLIB
    except ImportError:
        direct = False
    if not direct:
        return _line_contours_with_axes(x, y, z, level)
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.ndim == 1 and y.ndim == 1:
        x, y = np.meshgrid(x, y)
    masked = np.ma.masked_invalid(np.ma.asarray(z, dtype=np.float64))
    algorithm = mpl.rcParams['contour.algorithm']
    corner_mask = False if algorithm == 'mpl2005' else mpl.rcParams['contour.corner_mask']
    generator = contourpy.contour_generator(
        x, y, masked, name=algorithm, corner_mask=corner_mask,
        line_type=contourpy.LineType.SeparateCode, fill_type=contourpy.FillType.OuterCode, chunk_size=0)
    vertices, _codes = generator.lines(level)
    return list(vertices) if len(vertices) else [np.empty((0, 2))]


def _line_contours_with_axes(x, y, z, level):
    """The same through ``Axes.contour`` on a throw-away Figure (no pyplot state)."""
    from matplotlib.figure import Figure
    axes = Figure().subplots()
    return list(axes.contour(x, y, z, levels=[level]).allsegs[0])
