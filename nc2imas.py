#!/usr/bin/python
"""Core TRANSP NetCDF-to-IMAS conversion helpers, version metadata, and shared mapping utilities."""

#
import os
import sys
import math
import hashlib
from pathlib import Path
import xml.etree.ElementTree as ET
from typing import List
import numpy as np
import inspect

from scipy import constants
from scipy.interpolate import interp1d
from datetime import datetime, timezone
from input_provenance import InputInventory, append_file_xml, fill_model_libraries

__version__ = "0.13.0"  # x-release-please-version
VERSION = __version__
TRANSP2IMAS_VERSION = __version__
import re
import logging
from netcdf_class import netcdf_connection
from radial_interpolation import interpolate_radial_profiles
from cocos_signs import direction_override, resolve_signs
from ids_assignment import (check_1dim, fill_data, fill_value, is_valid, node_path, np_multiply_sum,
                            operate_arrays, sum_arrays)
from code_provenance import (TRANSP_METADATA_KEYS, attach_source_provenance, build_code_parameters_xml,
                             properties_comment, source_file_records)
from imas_class import imas_connection
from conversion_context import ConversionContext
import namelist_reader
from plasma_betas import boundary_minor_radius, imas_betas
from conversion_provenance import ProvenanceRecorder
import grid_services
import ids_identity
from name_patterns import make_namelist_re_pattern
from atoms_class import *
from utilities import *


class LineNumber:
    def __str__(self):
        return str(inspect.currentframe().f_back.f_lineno)

__LINE__ = LineNumber()


def getVarName(var):
    localVars = inspect.currentframe().f_back.f_locals.items()
    return str([k for k, v in localVars if v is var][0])


class nc2imas(imas_connection, netcdf_connection):
    def __init__(self, run, backend, occur, user, machine, new_db, ncfile,
         nml_file:str=None, log:str='error', debug:int=0,
         ip_direction:str='auto', bt_direction:str='auto', dd_version:str=None):
        netcdf_connection.__init__(self, ncfile, machine)
        imas_connection.__init__(self, self.device, self.shot, run, backend, occur, user, new_db, dd_version=dd_version)
        self.nml_name = nml_file
        self.logger = None
        #namelist contents dict
        self.nmlist = None
        self.set_logger(log)
        self.logger.info(f'creating nc2imas using {ncfile} for {self.device} {self.shot} {self.run}')
        self.load_namelist()
        self._init_cocos_signs(ip_direction=ip_direction, bt_direction=bt_direction)
        self.cm2m = 0.01
        self.area_inv = 1.0e4
        self.vconv = 1.e-6
        self.aconv = 1.0e-4
        self.vol_inv = 1.e6
        self.timesc = []
        self.timep = []
        self.ntimes = 0
        self.ntimep = 0
        self.pi = constants.pi  #scipy

        self.debug = False if debug<1 else True
        self.debug_max = False if debug<2 else True
        self.get_nc_times()

    @property
    def context(self):
        """The writers' view of this converter (``conversion_context``), built on first use."""
        context = self.__dict__.get('_context')
        if context is None:
            context = self._context = ConversionContext(self)
        return context

    @property
    def _provenance(self):
        """Provenance blocks and their caches (``conversion_provenance``), built on first use."""
        recorder = self.__dict__.get('_provenance_recorder')
        if recorder is None:
            recorder = self._provenance_recorder = ProvenanceRecorder(self, TRANSP2IMAS_VERSION)
        return recorder

    def close_imas_database(self):
        """Close the IMAS entry and TRANSP CDF deterministically.

        ``netcdf_connection.__del__`` remains a best-effort safety net, but
        normal converter operation must not depend on garbage collection for
        closing the input NetCDF file.
        """
        try:
            imas_connection.close_imas_database(self)
        finally:
            self.close_ncfile()

    def __enter__(self):
        """Return the active converter for use as a context manager."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Close IMAS and NetCDF resources without suppressing exceptions."""
        try:
            self.close_imas_database()
        except (AttributeError, RuntimeError, TypeError, OSError) as close_error:
            if exc_type is None:
                raise
            if self.logger is not None:
                self.logger.error(f'Error while closing converter resources: {close_error}')
        return False


    @staticmethod
    def _parse_bool_value(value):
        """Parse TRANSP namelist logical values into Python bool."""
        return namelist_reader.parse_bool_value(value)

    def get_namelist_bool(self, name:str, default=None):
        """Return a scalar TRANSP logical from the loaded TR.DAT namelist (case-insensitive key)."""
        return namelist_reader.get_bool(self.nmlist, name, default)

    # (get_ncattr_text lives in netcdf_class: reading a CDF attribute is a CDF concern)

    def set_species_name(self, node, value):
        """Set the species/state text node for the selected DD version (see ``ids_identity``)."""
        ids_identity.set_species_name(self.dd_at_least, node, value)

    # Backward-compatible alias used by older converter modules.
    def set_name_or_label(self, node, value):
        self.set_species_name(node, value)

    def set_nuclear_charge(self, node, value):
        """Set z_n using the DD-native scalar type (see ``ids_identity``)."""
        ids_identity.set_nuclear_charge(self.dd_at_least, node, value)

    def set_ec_beam_identity(self, beam, identifier, description=None):
        """Populate EC beam identity for DD3.39 or DD4 (see ``ids_identity``)."""
        ids_identity.set_ec_beam_identity(self.dd_at_least, beam, identifier, description)

    def set_nbi_unit_identity(self, unit, identifier, description=None):
        """Populate NBI unit identity for DD3.39 or DD4 (see ``ids_identity``)."""
        ids_identity.set_nbi_unit_identity(self.dd_at_least, unit, identifier, description)

    def set_ec_mode(self, beam, o_mode_fraction, times):
        """Write EC launch mode/polarization for the target DD schema (see ``ids_identity``)."""
        ids_identity.set_ec_mode(self.dd_at_least, beam, o_mode_fraction, times)

    @staticmethod
    def _direction_override(value:str):
        return direction_override(value)

    def _init_cocos_signs(self, ip_direction:str='auto', bt_direction:str='auto'):
        """Set physical direction and target COCOS signs (policy in ``cocos_signs``).

        ``NLJCCW``/``NLBCCW`` come from the namelist unless overridden; the signs are published
        as attributes (``ip_sign``, ``bt_sign``, ``psi_sign`` ...) and as ``self.signs``.
        """
        signs = resolve_signs(ip_direction, bt_direction,
                              lambda name: self.get_namelist_bool(name, default=None),
                              self.dd_major_version)
        self.signs = signs
        for name, value in signs.as_attributes().items():
            setattr(self, name, value)
        if signs.missing:
            print('*** WARNING: missing TRANSP field-orientation control(s): ' + ', '.join(signs.missing))
            print('***          assuming CCW. Use --ip-direction/--bt-direction to override.')
        print(signs.summary())

    def get_nc_times(self):
        self.timesc, self.ntimes, self.timep, self.ntimep = grid_services.read_time_bases(self.get_ncvar, self.logger)


    def check_1dim(self, dat):
        return check_1dim(dat)

    #Don't call this using varargs forwarded by another function
    #Don't call this using varargs forwarded by another function
    def is_valid(self, *variables, ind:int=0, name:str=''):
        return is_valid(*variables, ind=ind, name=name, logger=self.logger)

    def sum_arrays(self, *arrays):
        """Element-wise sum of the non-empty inputs (see ``ids_assignment.sum_arrays``)."""
        return sum_arrays(*arrays, logger=self.logger)

    def operate_arrays(self, arr1, op:str, arr2=None, name1:str='', name2:str=''):
        """Apply a supported element-wise operation (see ``ids_assignment.operate_arrays``)."""
        return operate_arrays(arr1, op, arr2, name1, name2, logger=self.logger)

    def np_multiply_sum(self, arr1, arr2, axis:int, name1:str='', name2:str='', scale=None):
        return np_multiply_sum(arr1, arr2, axis, name1, name2, scale, logger=self.logger)

    def checkErrExit(self, *variables, name:str=''):
        for var in variables:
            stat = self.is_valid(var)
            if stat is False:
               self.logger.exception(f'\t\tvar {name} is invalid')
               raise Exception(f'**ERROR**: {name} is invalid: len={len(var)}: {var}')


    def interp1d_time(self, time_new, time, dat, ind1:int=None, ind2:int=None):
        r'''interpolate ``dat`` from ``time`` onto ``time_new`` (see ``grid_services.interp_time``).'''
        return grid_services.interp_time(self.is_valid, time_new, time, dat, ind1, ind2)

    def get_ncvar_on_time(self, var: str, target: str = 'TIME3', scale: float = 1.0):
        """A 1-D TRANSP history on ``TIME`` or ``TIME3`` (see ``grid_services.history_on_time``)."""
        return grid_services.history_on_time(self.nc, self.get_ncvar, self.timesc, self.timep, var, target, scale)

    #fill ids.value
    def fill_value(self, ids, dat, rsz:int=-1, mult:float=None, name:str=None):
        return fill_value(ids, dat, rsz=rsz, mult=mult, logger=self.logger)

    def fill_value_nc(self, ids, var:str, rsz:int=-1, mult:float=None):
        dat = self.get_ncvar(var)
        if dat is None:
            pass
            #self.logger.info(f'{var} not found')
        else:
            self.fill_value(ids, dat, rsz=rsz, mult=mult)

    # if var is 2D, ind must be valid. If ind<0, use this to fill ONLY 1D data
    def fill_data_nc(self, ids, var:str, mult:float=None, ind:int=-1, resize:bool=True,
         interp:bool=False):
        if not isinstance(var, str):
            raise ValueError(f'variable {var} is not string')
        dat = self.get_ncvar(var, ind=ind)
        if not self.is_valid(dat):
            if ind==0:
                self.logger.info(f'{var} not-found')
        else:
            if mult is not None:
                dat = mult*dat
            if interp is True and not np.array_equal(self.timesc, self.timep):
                dat = np.interp(self.timep, self.timesc, dat)
            if resize:
                ids.resize(len(dat), refcheck=False)
            ids[:] = dat


    #If ind<0, use this to fill ONLY 1D data. if var is 2D, ind must be valid.
    #Don't use this to fill scalar
    @staticmethod
    def _node_path(node):
        """IDS path of ``node`` for status messages (never raises)."""
        return node_path(node)

    def fill_data(self, ids, dat, mult:float=None, ind:int=-1, ind2:int=-1,
          interp:bool=False, name:str=None, resize:bool=True, chk:int=-1) -> bool:
        """Assign ``dat`` to an array node; returns True if it was written (policy: ``ids_assignment``).

        Absent source data is a recorded skip; present data that cannot be assigned is a broken
        mapping (``MappingError`` in strict mode, a skipped field with the reason in lenient mode).
        """
        return fill_data(ids, dat, self.mapping, self.logger, mult=mult, ind=ind, ind2=ind2,
                         interp=interp, times=lambda: (self.timesc, self.timep), name=name,
                         resize=resize, chk=chk)

    def load_namelist(self):
        loaded = namelist_reader.load_namelist(self.nml_name, self.logger)
        if self.nml_name is not None:
            self.nmlist = loaded

    def reload_namelist(self, fname:str):
        self.nml_name = fname
        if self.nmlist is not None:
            self.nmlist.clear()
        self.load_namelist()
        # A changed TR.DAT must not retain XML/provenance for the prior file.
        self._provenance.reset()


    #TODO move to utilities
    def parse_number(self, ss:str, tp):
        return namelist_reader.parse_number(ss, tp)


    #name's last character is discarded. Eg. XECECH?. () added separately
    def get_namelist_nvars(self, name:str, ndig:int=2):
        return namelist_reader.count_vars(self.nmlist, self.logger, name, ndig)


    def make_namelist_re_pattern(self, name:str, nd1:int=0, nd2:int=0, nch1:int=0,
          nch2:int=0, nd:int=0, nch:int=0) -> str:
        return make_namelist_re_pattern(name, nd1, nd2, nch1, nch2, nd, nch)

    def get_namelist_keys_re(self, name:str, cut:str='', nd1:int=0, nd2:int=0, nch1:int=0,
          nch2:int=0, nd:int=0, nch:int=0) -> List:
        return namelist_reader.keys_re(self.nmlist, name, cut, nd1, nd2, nch1, nch2, nd, nch)


    #returns scalar or list
    def get_namelist_var_(self, name:str, tp=None, scale=1., ind:int=-1):
        """Scalar or list for one namelist key; ``None`` when absent."""
        return namelist_reader.get_var_(self.nmlist, self.logger, name, tp, scale, ind)

    # Some scalar can stand for arrays, e.g: FREQECH, which causes error when
    #  an array is expected. So, call for array should handle single variable also.
    # Returns scalar for single, list for array on RHS
    def get_namelist_vars(self, name:str, tp=None, scale=1.0, ind:int=-1,
          size:int=1, resize:bool=False):
        """Scalar for a single value, array for an indexed list on the right-hand side."""
        return namelist_reader.get_vars(self.nmlist, self.logger, name, tp, scale, ind, size, resize)

    def get_namelist_var(self, name:str, tp=None, scale=1.0, ind:int=-1, size:int=1):
        """Compatibility wrapper around :meth:`get_namelist_vars`."""
        return self.get_namelist_vars(name, tp, scale, ind, size)


    def extract_str_between(self, name:str, pat1:str, pat2:str, default:str, tp=None):
        return namelist_reader.extract_str_between(self.logger, name, pat1, pat2, default, tp)

    # (fluxBdryToCenter removed: it called a non-existent ``self.valid`` and nothing used it)

    def set_logger(self, log:str = None):
        self.logger = logging.getLogger(__name__)
        levels = {'debug': logging.DEBUG, 'info': logging.INFO,
                  'error': logging.ERROR}
        log = log if log is not None else 'error'
        level = levels.get(log.lower())
        self.logger.setLevel(level=level)
        logname = 'log.out'
        if os.path.exists(logname):
            os.rename(logname, f'{logname}_')
        self.logger.addHandler(logging.FileHandler(logname, mode='w')) #'a'

    #for all times
    @staticmethod
    def cumulative_zone_quantity(values):
        """Convert TRANSP per-zone geometry to an enclosed flux-surface quantity."""
        return grid_services.cumulative_zone_quantity(values)

    @staticmethod
    def interp_radial_profiles(values, x_from, x_to, axis_value=None):
        """Interpolate time-dependent radial profiles without mixing time slices.

        ``axis_value`` optionally adds an explicit point at rho=0 before
        interpolation.  This is important for TRANSP XB data because XB starts
        at the outer boundary of the first finite-volume zone, while X contains
        the zone centre inside it.  The grid mapping is cached by coordinate
        content (see ``radial_interpolation``).
        """
        return interpolate_radial_profiles(values, x_from, x_to, axis_value)

    def reference_radius(self):
        """The fixed reference major radius ``r0`` in m: the middle of the limiter ``RLIM``, or ``None``.

        The Data Dictionary asks for "a fixed position such as the middle of the vessel".  The limiter is
        the machine's own first wall, so ``(Rmax + Rmin)/2`` of it is the same for every shot of a machine,
        unlike the magnetic-axis radius, which moves with the plasma.  ``None`` when the CDF has no usable
        limiter (fewer than three finite points).
        """
        rlim = self.get_ncvar('RLIM', self.cm2m)
        if not self.is_valid(rlim):
            return None
        points = np.asarray(rlim[0] if rlim.ndim > 1 else rlim, dtype=np.float64)
        points = points[np.isfinite(points) & (points > 0.0)]
        if points.size < 3:
            return None
        return float(0.5 * (points.max() + points.min()))

    def vacuum_field_reference(self, times=None):
        """``(r0, b0(times))`` of the vacuum toroidal field, or ``None`` without BZXR.

        ``times`` defaults to the writer time grid ``timep``; ``r0`` is :meth:`reference_radius`, else the
        mean magnetic-axis radius (needs RAXIS); see :func:`utilities.vacuum_field_reference`.
        """
        bzxR = self.get_ncvar('BZXR', self.cm2m)         # BZXR(TIME) T m, vacuum field x radius
        r0 = self.reference_radius()
        raxis = None
        if r0 is None:
            raxis = self.get_ncvar('RAXIS', self.cm2m)   # RAXIS(TIME) m
            if not self.is_valid(raxis):
                return None
        if not self.is_valid(bzxR):
            return None
        return vacuum_field_reference(raxis, bzxR, self.bt_sign,
                                      self.timep if times is None else times, self.timesc, r0=r0)

    def check_profile_shape(self, name, array, nt, nx):
        """``array`` if it is a ``(nt, nx)`` time-radius profile (or empty), else a broken mapping.

        A profile with another radial size cannot be placed on the writer's grid.  Catching it here
        names the variable, instead of failing later in some array arithmetic or being dropped.
        Returns an empty array when the lenient policy skips it.
        """
        values = np.asarray(array)
        if values.ndim == 2 and values.shape[0] == nt and values.shape[1] != nx:
            self.mapping.broken(f'radial size mismatch {values.shape[1]} != {nx} on the writer grid',
                                source=name, target=f'profile ({nt} times x {nx} radial points)')
            return np.array([])
        return array

    def vacuum_b0_magnitude(self, times=None):
        """``|b0(times)|`` (T) of the IMAS ``vacuum_toroidal_field`` node, or ``None``.

        IMAS defines the parallel current densities as ``average(J.B)/B0`` and
        ``rho_tor = sqrt(phi/(pi*B0))`` with B0 the IDS's own ``vacuum_toroidal_field/b0``,
        so those normalizations must use this field and not TRANSP's ``BZ`` (the field at a
        fixed radius that also differs from ``BZXR/r0`` by up to ~40%).  The physical sign is
        applied separately by the callers.
        """
        reference = self.vacuum_field_reference(times)
        return None if reference is None else np.abs(reference[1])

    def boundary_minor_radius(self):
        """Boundary minor radius ``(Rmax - Rmin)/2`` in m on ``TIME3`` from the ``RMCB<k>``/``RMSB<k>`` moments.

        Empty without the moments.  This is the Data Dictionary's definition of the boundary minor radius and
        agrees with the traced LCFS outline to 1-2%; ``sqrt(PAREA/pi)`` is 26-37% larger for the elongated
        bundled plasmas.
        """
        cosine, sine = [], []
        k = 0
        while f'RMCB{k}' in self.nc.variables:
            cosine.append(self.get_ncvar_on_time(f'RMCB{k}', scale=self.cm2m))
            sine.append(self.get_ncvar_on_time(f'RMSB{k}', scale=self.cm2m) if k and f'RMSB{k}' in self.nc.variables
                        else np.zeros_like(cosine[0]))
            k += 1
        if not cosine or not self.is_valid(cosine[0]):
            return np.array([])
        return boundary_minor_radius(np.array(cosine), np.array(sine))

    def plasma_betas(self):
        """``{'beta_tor', 'beta_pol', 'beta_tor_norm'}`` on ``TIME3`` by the IMAS definitions, or ``{}``.

        From the total pressure ``PTOWB`` and zone volumes ``DVOL``, the exported ``b0`` and ``r0`` (see
        :meth:`vacuum_field_reference`), the plasma current and :meth:`boundary_minor_radius`.  TRANSP's own
        ``BTEQ``, ``BETAT`` and ``BPEQ`` use other normalisations (see ``plasma_betas``), so they are not copied.
        """
        reference = self.vacuum_field_reference()
        pressure = self.get_ncvar('PTOWB')
        dvol = self.get_ncvar('DVOL', self.vconv)
        current = self.get_ncvar_on_time('PCURC')
        minor = self.boundary_minor_radius()
        if reference is None or not self.is_valid(pressure, dvol, current, minor):
            return {}
        r0, b0 = reference
        return imas_betas(pressure, dvol, b0, r0, self.ip_sign * current, minor)

    def fill_vac_tor_field(self, ids):
        reference = self.vacuum_field_reference()
        if reference is not None:
            r0, b0 = reference
            self.logger.info('\tVac_tor_field')
            ids.vacuum_toroidal_field.r0 = r0
            self.fill_data(ids.vacuum_toroidal_field.b0, b0)

    #fill ids.profiles_1d[*].grid. Resize if not already done
    def fill_profiles1d_grid(self, ids):
        # Core/source profiles in TRANSP are normally defined on X (zone centres),
        # while fluxes and SURF are on XB (zone boundaries).  Keep the IMAS grid
        # internally consistent by mapping boundary quantities to X.
        timesp = self.get_ncvar('TIME3')
        ntimep = len(timesp)
        coord_x = self.get_ncvar('X')
        coord_xb = self.get_ncvar('XB')
        rho_tor_norm = coord_x

        # Toroidal flux coordinate. TRFLX is on XB and already in Wb.
        bflux_tor_b = self.get_ncvar('TRFLX')
        bflux_tor = self.interp_radial_profiles(bflux_tor_b, coord_xb, coord_x, axis_value=0.0)
        b0 = self.vacuum_b0_magnitude(self.timep)      # rho_tor = sqrt(phi/(pi*B0)), B0 = vacuum b0
        rho_tor = []
        if self.is_valid(b0, bflux_tor):
            try:
                rho_tor = np.sqrt(np.clip(bflux_tor/(self.pi*b0[:,None]), 0.0, None))
            except (ValueError, FloatingPointError):
                rho_tor = []

        # Poloidal flux: PLFLX is an enclosed increment on XB in Wb/rad.
        psi_rel_b = self.psi_sign * self.get_ncvar('PLFLX', 2.*self.pi)
        psi_axis = self.psi_sign * self.get_ncvar('PSI0_TR', 2.*self.pi)
        if not np.array_equal(self.timesc, self.timep):
            psi_axis = self.interp1d_time(self.timep, self.timesc, psi_axis)
        psi_rel = self.interp_radial_profiles(psi_rel_b, coord_xb, coord_x, axis_value=0.0)
        psi = psi_rel + psi_axis[:,None] if self.is_valid(psi_rel, psi_axis) else psi_rel
        psi_boundary = psi_rel_b[:,-1] + psi_axis if self.is_valid(psi_rel_b, psi_axis) else []
        rho_pol_norm = []
        if self.is_valid(psi, psi_axis, psi_boundary):
            try:
                denom = psi_boundary - psi_axis
                rho_pol_norm = np.sqrt(np.clip((psi-psi_axis[:,None])/denom[:,None], 0.0, None))
            except (ValueError, TypeError, IndexError, FloatingPointError) as exc:
                self.mapping.broken(exc, source='PLFLX', target='profiles_1d/grid/rho_pol_norm')
                rho_pol_norm = []

        # Finite-zone geometry -> enclosed values at XB -> interpolate to X.
        gvol = self.get_ncvar('DVOL', self.vconv)
        volume_b = self.cumulative_zone_quantity(gvol)
        volume_enclosed = self.interp_radial_profiles(volume_b, coord_xb, coord_x, axis_value=0.0)
        grid_area = self.get_ncvar('DAREA', self.aconv)
        area_b = self.cumulative_zone_quantity(grid_area)
        area_enclosed = self.interp_radial_profiles(area_b, coord_xb, coord_x, axis_value=0.0)
        surface_b = self.get_ncvar('SURF', 1.0e-4)
        surface = self.interp_radial_profiles(surface_b, coord_xb, coord_x, axis_value=0.0)

        if len(ids.profiles_1d)==0:
            ids.profiles_1d.resize(ntimep)

        for i in range(ntimep):
            gd = ids.profiles_1d[i].grid
            self.fill_data(gd.rho_tor_norm, rho_tor_norm, ind=i)
            self.fill_data(gd.rho_tor, rho_tor, ind=i)
            self.fill_data(gd.rho_pol_norm, rho_pol_norm, ind=i)
            self.fill_data(gd.psi, psi, ind=i)
            if self.is_valid(psi_axis, ind=i) and hasattr(gd, 'psi_magnetic_axis'):
                gd.psi_magnetic_axis = float(psi_axis[i])
            if self.is_valid(psi_boundary, ind=i) and hasattr(gd, 'psi_boundary'):
                gd.psi_boundary = float(psi_boundary[i])
            self.fill_data(gd.volume, volume_enclosed, ind=i)
            self.fill_data(gd.area, area_enclosed, ind=i)
            self.fill_data(gd.surface, surface, ind=i)

    def _provenance_file_records(self):
        """Portable source-file records for IDS-level provenance (see ``ProvenanceRecorder``)."""
        return self._provenance.file_records()

    def _fill_source_provenance(self, ids):
        """Attach source CDF/TR.DAT provenance to the whole IDS (see ``ProvenanceRecorder``)."""
        return self._provenance.fill_source_provenance(ids)

    def fill_properties(self, ids):
        """Fill ``ids_properties`` (comment, creation date, name, provenance); see ``ProvenanceRecorder``."""
        return self._provenance.fill_properties(ids)

    @staticmethod
    def radial_gradient(values, coordinate):
        """d(values)/d(coordinate) per time slice (see ``grid_services.radial_gradient``)."""
        return grid_services.radial_gradient(values, coordinate)

    def input_inventory(self):
        """Resolve and cache the run's TR.DAT and referenced inputs (see ``ProvenanceRecorder``)."""
        return self._provenance.input_inventory()

    def _collect_parameter_files(self):
        """Legacy API consumed by IDS-provenance summary; use new inventory."""
        return self._provenance.parameter_files()

    def fill_external_model_inputs(self, ids, family):
        """Attach model-specific XML to the appropriate IDS code libraries."""
        return self._provenance.fill_external_model_inputs(ids, family)

    def build_code_parameters_xml(self):
        """XML for IMAS code.parameters; cached (see ``ProvenanceRecorder``)."""
        return self._provenance.build_code_parameters_xml()

    def fill_code_properties(self, ids):
        """Fill common code provenance, if this IDS contains a code node."""
        return self._provenance.fill_code_properties(ids)


if __name__ == '__main__':
    print("Please use: python convert.py")
