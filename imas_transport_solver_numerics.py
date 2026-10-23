#!/usr/bin/python
"""Populate ``transport_solver_numerics`` with TRANSP PT_SOLVER runtime diagnostics.

Only quantities that exist in the CDF are written, as solver-specific scalar control parameters
of each ``solver_1d`` slice:

* ``NPTRITER`` ("Number Newton Iterations in PT_SOLVER") -> integer ``newton_iterations``,
* ``PT_CPTIME`` / ``PT_WTIME`` ("CPU/WALL TIME USED (PT_SOLVER)", hours) ->
  real ``cumulative_cpu_time_hours`` / ``cumulative_wall_time_hours``.  The control-parameter
  node carries no unit, so the unit is part of the name and the CDF hours are stored unchanged.

The slices are on the channels' own native CDF ``TIME`` coordinate (``ids.time``), not the
``TIME3`` grid used by the profile IDSs.  ``NPTRITER`` is one global count per time, so the
per-equation ``convergence.iterations_n`` is deliberately left empty, and no residual, tolerance
or convergence status is invented.  Channels that are zero everywhere (solver not used) are not
written.  The lossless raw arrays stay available in the ``core_transport`` code parameters.

DD 3.39.0 control parameters have no ``name`` (added in 3.42.0), so nothing is written there.
"""

import numpy as np

from pt_solver_diagnostics import channel_coordinates, decode_channel, diagnostic_xml

INTEGER_PARAMETERS = {'NPTRITER': 'newton_iterations'}
REAL_PARAMETERS = {'PT_CPTIME': 'cumulative_cpu_time_hours', 'PT_WTIME': 'cumulative_wall_time_hours'}
SOLVER_INDEX = -1            # private identifier (negative index) for the PT_SOLVER code


class imas_transport_solver_numerics:
    """Writer registered as ``transport_solver_numerics`` in ``convert.py``."""

    def __init__(self, ncimas):
        self.im = ncimas
        self.logger = ncimas.logger
        self.create_ids()

    def _named_control_parameters(self):
        """True if this DD names its control parameters (DD >= 3.42.0; DD 3.39.0 has only ``value``).

        An unnamed value cannot say which parameter it is, so nothing is written without names.
        """
        scratch = self.im.new_ids('transport_solver_numerics')
        scratch.solver_1d.resize(1)
        control = scratch.solver_1d[0].control_parameters
        control.integer0d.resize(1)
        control.real0d.resize(1)
        return hasattr(control.integer0d[0], 'name') and hasattr(control.real0d[0], 'name')

    def _channels(self, diagnostics):
        """``(time, {parameter name: values}, {parameter name: values})`` of the usable channels."""
        integers, reals, time = {}, {}, None
        for channel in diagnostics.findall('channel'):
            name = channel.get('name')
            if name not in INTEGER_PARAMETERS and name not in REAL_PARAMETERS:
                continue
            values = decode_channel(channel)
            coordinate = channel_coordinates(diagnostics, channel).get('TIME')
            if channel.get('dimensions') != 'TIME' or coordinate is None or values.shape != coordinate.shape:
                self.logger.warning(f'transport_solver_numerics: {name} is not a series on a verified TIME '
                                    'coordinate; skipped')
                continue
            if time is not None and not np.array_equal(time, coordinate):
                self.logger.warning(f'transport_solver_numerics: {name} is on a different time grid; skipped')
                continue
            if name in INTEGER_PARAMETERS:
                finite = values[np.isfinite(values)]
                if finite.size != values.size or not np.allclose(finite, np.rint(finite)):
                    self.logger.warning(f'transport_solver_numerics: {name} is not an integer count; skipped')
                    continue
                integers[INTEGER_PARAMETERS[name]] = np.rint(values).astype(np.int64)
            else:
                unit = (channel.get('units') or '').strip().upper()
                if channel.get('canonical_unit') != 'hour' or unit not in ('HOURS', 'HOUR'):
                    self.logger.warning(f'transport_solver_numerics: {name} units {unit!r} are not hours; skipped')
                    continue
                reals[REAL_PARAMETERS[name]] = values
            time = coordinate if time is None else time
        return time, integers, reals

    def create_ids(self):
        self.logger.info('\ntransport_solver_numerics')
        if not self._named_control_parameters():
            self.logger.info('transport_solver_numerics: this DD version has unnamed control parameters; not '
                             'written (the raw channels remain in the core_transport code parameters)')
            self.im.mapping.skipped('this DD version has unnamed control parameters (named from DD 3.42.0)',
                                    target='transport_solver_numerics')
            return
        diagnostics = diagnostic_xml(getattr(self.im, 'nc', None), record_inactive=False)
        time, integers, reals = self._channels(diagnostics)
        if time is None or not (integers or reals):
            self.logger.info('transport_solver_numerics: no active PT_SOLVER runtime channels; skipping')
            self.im.mapping.skipped('no active PT_SOLVER runtime channels (NPTRITER/PT_CPTIME/PT_WTIME)',
                                    target='transport_solver_numerics')
            return

        ids = self.im.new_ids('transport_solver_numerics')
        ids = self.im.fill_properties(ids)
        self.im.fill_code_properties(ids)
        ids.ids_properties.homogeneous_time = 1
        ids.time = np.asarray(time, dtype=np.float64)
        ids.ids_properties.comment = (
            f'{ids.ids_properties.comment}; TRANSP PT_SOLVER runtime diagnostics on the native CDF TIME '
            'coordinate (not TIME3); timers are cumulative hours as stored in the CDF')
        ids.solver.name = 'PT_SOLVER'
        ids.solver.index = SOLVER_INDEX
        ids.solver.description = 'TRANSP PT_SOLVER predictive transport solver'
        reference = self.im.vacuum_field_reference(time)
        if reference is not None:
            ids.vacuum_toroidal_field.r0, ids.vacuum_toroidal_field.b0 = reference

        ids.solver_1d.resize(time.size)
        for it in range(time.size):
            control = ids.solver_1d[it].control_parameters
            control.integer0d.resize(len(integers))
            for slot, (name, values) in zip(control.integer0d, integers.items()):
                slot.name = name
                slot.value = int(values[it])
            control.real0d.resize(len(reals))
            for slot, (name, values) in zip(control.real0d, reals.items()):
                slot.name = name
                slot.value = float(values[it])
        self.im.putdb_imas(ids)
        self.logger.info('transport_solver_numerics: wrote '
                         f'{", ".join(list(integers) + list(reals))} at {time.size} native TIME samples')


if __name__ == '__main__':
    print('Please use: python convert.py')
