#!/usr/bin/python
"""Base class of IDS writers whose construction does not write.

A writer maps one IDS.  ``build()`` creates and fills the IDS objects and *returns* them; nothing is
submitted to the database by constructing the writer or by building.  ``run()`` (called by the
conversion runner through ``run_writer``) submits what ``build()`` returned.  Because the IDS
objects are values, a writer can be built and checked without a database, and what it will write is
visible to the runner before anything is submitted.

``build()`` returns a list of IDS objects, or of ``(ids, occurrence)`` pairs for an explicit
occurrence (an empty list means "nothing to write", which the journal records as a skip).  Return a
list, not a generator: ``run()`` submits only after ``build()`` finished, so a build that fails
submits nothing.

Writers not yet migrated still do all their work in their constructor; ``run_writer`` leaves them
alone, so migrated and legacy writers share one runner.
"""

import numpy as np

from conversion_context import ConversionContext


class IdsWriter:
    """Writer interface: ``build()`` returns the IDS objects, ``run()`` submits them."""

    def __init__(self, ncimas):
        self.im = ncimas
        self.logger = ncimas.logger

    # What a writer may use of the converter, by role (see conversion_context).  Derived from ``self.im`` on
    # access and cached per converter, so a writer built without ``__init__`` or given another ``im`` in a
    # test gets matching views.  ``self.im`` itself stays for what no view offers.
    @property
    def ctx(self):
        context = self.__dict__.get('_ctx')
        if context is None or context.owner is not self.im:
            context = self._ctx = ConversionContext.of(self.im)
        return context

    @property
    def cdf(self):
        return self.ctx.cdf

    @property
    def namelist(self):
        return self.ctx.namelist

    @property
    def dd(self):
        return self.ctx.dd

    @property
    def signs(self):
        return self.ctx.signs

    @property
    def grid(self):
        return self.ctx.grid

    @property
    def fill(self):
        return self.ctx.fill

    def check_profile(self, name, array, nt, nx, grid=None):
        """``array`` if it is empty or a profile that fits ``(nt, nx)`` and its declared coordinates (see
        ``nc2imas.check_profile_shape``); under the lenient policy a profile that does not fit comes back empty."""
        dimensions = getattr(self.cdf, 'variable_dimensions', None)
        return self.fill.check_profile_shape(name, array, nt, nx,
                                             dimensions=None if dimensions is None else dimensions(name), grid=grid)

    def profile(self, name, scale=1.0):
        """A ``TIME3 x radial`` CDF profile; a radial size other than the zone grid's is a named broken mapping.

        The grid is the one of ``X`` (its shape is read once per writer).  Without this, a mis-sized profile was
        skipped by a shape comparison, or failed later in numpy broadcasting without naming the variable.
        Under the lenient policy the variable is skipped with the reason and the result is empty.
        """
        values = self.cdf.get_ncvar(name, scale)
        if values.ndim != 2:
            return values
        shape = self.__dict__.get('_x_shape')
        if shape is None:
            grid = self.cdf.get_ncvar('X')
            shape = self._x_shape = tuple(grid.shape) if grid.ndim == 2 else ()
        if len(shape) == 2:
            values = np.asarray(self.check_profile(name, values, shape[0], shape[1]))
        return values

    def build(self):
        raise NotImplementedError(f'{type(self).__name__}.build() must return the IDS objects to write')

    def run(self):
        """Submit what ``build()`` returns; the number of IDS objects submitted."""
        submitted = 0
        for item in self.build():
            ids, occurrence = item if isinstance(item, tuple) else (item, None)
            if occurrence is None:
                self.im.putdb_imas(ids)
            else:
                self.im.putdb_imas_profile(ids, occurrence)
            submitted += 1
        return submitted


def run_writer(writer):
    """Run a writer object: migrated writers submit via ``run()``, legacy ones already worked in their constructor."""
    if isinstance(writer, IdsWriter):
        return writer.run()
    return None
