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

__version__ = "0.12.1"  # x-release-please-version
VERSION = __version__
TRANSP2IMAS_VERSION = __version__
import re
import logging
from netcdf_class import netcdf_connection
from imas_class import imas_connection
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
        self._code_parameters_xml_cache = None
        self._provenance_file_records_cache = None
        self._input_inventory_cache = None

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
        if isinstance(value, (bool, np.bool_)):
            return bool(value)
        if isinstance(value, (int, np.integer)):
            return bool(value)
        if value is None:
            return None
        text = str(value).strip().lower().replace('.', '')
        if text in ('true', 't', 'yes', 'y', '1'):
            return True
        if text in ('false', 'f', 'no', 'n', '0'):
            return False
        return None

    def get_namelist_bool(self, name:str, default=None):
        """Return a scalar TRANSP logical from the loaded TR.DAT namelist.

        ``load_namelist`` preserves the spelling/case of keys.  Field-orientation
        controls are important enough that lookup here is deliberately
        case-insensitive.
        """
        if not self.nmlist:
            return default
        key = next((k for k in self.nmlist if k.strip().upper() == name.upper()), None)
        if key is None:
            return default
        value = self._parse_bool_value(self.nmlist[key])
        return default if value is None else value

    def get_ncattr_text(self, name:str, default:str=''):
        """Return a CDF attribute as stripped text.

        netCDF readers may expose character attributes as ``str``, ``bytes``,
        NumPy byte/string scalars, or small character arrays.  Metadata writers
        should not depend on a reader-private decoding helper, so normalize the
        common representations here.
        """
        try:
            value = self.get_ncattr(name)
        except (AttributeError, KeyError, RuntimeError):      # closed/unreadable dataset: attribute not available
            return default
        if value is None:
            return default
        if isinstance(value, (bytes, np.bytes_)):
            return bytes(value).decode('ISO-8859-1', errors='replace').strip()
        if isinstance(value, np.ndarray):
            if value.size == 0:
                return default
            if value.size == 1:
                value = value.reshape(-1)[0]
                if isinstance(value, (bytes, np.bytes_)):
                    return bytes(value).decode('ISO-8859-1', errors='replace').strip()
                return str(value).strip()
            if value.dtype.kind in ('S', 'U'):
                parts = []
                for item in value.reshape(-1):
                    if isinstance(item, (bytes, np.bytes_)):
                        parts.append(bytes(item).decode('ISO-8859-1', errors='replace'))
                    else:
                        parts.append(str(item))
                return ''.join(parts).strip()
        return str(value).strip()

    def set_species_name(self, node, value):
        """Set the species/state text node for the selected DD version.

        DD 3.39 uses ``label``.  The species/state nodes were renamed to
        ``name`` in late DD3 (3.42) and remain ``name`` in DD4.
        """
        text = str(value)
        if self.dd_at_least(3, 42, 0):
            node.name = text
        else:
            node.label = text

    # Backward-compatible alias used by older converter modules.
    def set_name_or_label(self, node, value):
        self.set_species_name(node, value)

    def set_nuclear_charge(self, node, value):
        """Set z_n using the DD-native scalar type.

        z_n changed from FLT_0D to INT_0D in DD4.
        """
        z = float(value)
        node.z_n = int(round(z)) if self.dd_at_least(4, 0, 0) else z

    def set_ec_beam_identity(self, beam, identifier, description=None):
        """Populate EC beam identity without mixing DD3/DD4 semantics.

        DD3.39: ``identifier`` is the short identifier and ``name`` is the
        descriptive beam name.  DD4: ``name`` is the short identifier and
        ``description`` carries descriptive text.
        """
        identifier = str(identifier)
        description = str(description if description is not None else identifier)
        if self.dd_at_least(4, 0, 0):
            beam.name = identifier
            beam.description = description
        else:
            beam.identifier = identifier
            beam.name = description

    def set_nbi_unit_identity(self, unit, identifier, description=None):
        """Populate NBI unit identity for DD3.39 or DD4."""
        identifier = str(identifier)
        description = str(description if description is not None else identifier)
        if self.dd_at_least(4, 0, 0):
            unit.name = identifier
            unit.description = description
        else:
            unit.identifier = identifier
            unit.name = description

    def set_ec_mode(self, beam, o_mode_fraction, times):
        """Write EC launch mode/polarization using the target DD schema.

        DD3 uses the integer ``mode`` node and, where available, the legacy
        beam-level ``o_mode_fraction`` array.  In DD4.1+, ``mode`` is removed
        and the fraction is stored under ``beam/polarization``.
        """
        frac = float(np.clip(float(o_mode_fraction), 0.0, 1.0))
        times = np.asarray(times, dtype=np.float64)
        vals = np.full(len(times), frac, dtype=np.float64)
        if self.dd_at_least(4, 1, 0):
            beam.polarization.time = times
            beam.polarization.o_mode_fraction = vals
        elif self.dd_at_least(4, 0, 0):
            # DD4.0 was transitional.  Prefer polarization when available;
            # fall back to the pre-4.1 beam-level fraction if present.
            if hasattr(beam, 'polarization') and hasattr(beam.polarization, 'o_mode_fraction'):
                if hasattr(beam.polarization, 'time'):
                    beam.polarization.time = times
                beam.polarization.o_mode_fraction = vals
            elif hasattr(beam, 'o_mode_fraction'):
                beam.o_mode_fraction = vals
        else:
            # DD3 mode is the *main* excited wave mode.  For a mixed launch,
            # infer the dominant mode from the O-mode power fraction.
            beam.mode = 1 if frac >= 0.5 else -1
            if hasattr(beam, 'o_mode_fraction'):
                beam.o_mode_fraction = vals

    @staticmethod
    def _direction_override(value:str):
        if value is None:
            return None
        text = str(value).strip().lower()
        if text == 'ccw':
            return True
        if text == 'cw':
            return False
        return None

    def _init_cocos_signs(self, ip_direction:str='auto', bt_direction:str='auto'):
        """Set physical direction and target COCOS signs.

        TRANSP stores the field-orientation controls separately from many CDF
        quantities.  ``NLJCCW`` says whether the plasma current is
        counter-clockwise from above and ``NLBCCW`` says the same for the
        toroidal magnetic field.  IMAS uses +phi counter-clockwise from above.

        IMAS DD3 uses COCOS 11; DD4 uses COCOS 17.  COCOS 11 and 17 have the
        same cylindrical/toroidal orientation and q sign convention, but the
        poloidal-flux sign is reversed.  In particular, for positive Ip,
        COCOS17 psi decreases from magnetic axis to boundary.
        """
        jccw = self._direction_override(ip_direction)
        bccw = self._direction_override(bt_direction)
        if jccw is None:
            jccw = self.get_namelist_bool('NLJCCW', default=None)
        if bccw is None:
            bccw = self.get_namelist_bool('NLBCCW', default=None)

        missing = []
        if jccw is None:
            missing.append('NLJCCW')
            jccw = True
        if bccw is None:
            missing.append('NLBCCW')
            bccw = True

        self.nljccw = bool(jccw)
        self.nlbccw = bool(bccw)
        self.ip_sign = 1.0 if self.nljccw else -1.0
        self.bt_sign = 1.0 if self.nlbccw else -1.0

        try:
            dd_major = int(str(self.dd_major_version).split('.')[0])
        except (ValueError, TypeError, AttributeError):
            dd_major = 4
        self.target_cocos = 17 if dd_major >= 4 else 11

        # Physical/vector quantities in IMAS +phi (counter-clockwise) direction.
        self.jphi_sign = self.ip_sign
        self.jparallel_sign = self.ip_sign
        self.phi_sign = self.bt_sign
        self.f_sign = self.bt_sign
        self.q_sign = self.ip_sign * self.bt_sign
        self.vloop_sign = self.ip_sign
        self.toroidal_rotation_sign = self.ip_sign
        self.poloidal_field_sign = self.ip_sign

        # COCOS 11 -> 17 reverses the sign of full poloidal flux while leaving
        # Ip, Bphi and q physical sign definitions unchanged.
        self.psi_sign = self.ip_sign if self.target_cocos == 11 else -self.ip_sign

        # TRANSP moments use the geometrical atan2(Z-Zaxis,R-Raxis) angle
        # (counter-clockwise in the R-Z plane).  DD4/COCOS17 uses
        # theta=-atan2(...), clockwise.  This multiplier maps target theta to
        # the raw TRANSP moment angle when evaluating the Fourier surfaces.
        self.theta_raw_from_target = 1.0 if self.target_cocos == 11 else -1.0

        if missing:
            print('*** WARNING: missing TRANSP field-orientation control(s): ' + ', '.join(missing))
            print('***          assuming CCW. Use --ip-direction/--bt-direction to override.')

        print('COCOS/signs:'
              f' target={self.target_cocos}'
              f' NLJCCW={self.nljccw} (Ip sign {self.ip_sign:+.0f})'
              f' NLBCCW={self.nlbccw} (Bphi sign {self.bt_sign:+.0f})'
              f' q sign {self.q_sign:+.0f}'
              f' psi sign {self.psi_sign:+.0f}')

    def get_nc_times(self):
        self.timesc = self.get_ncvar('TIME')
        self.ntimes = len(self.timesc)
        self.timep = self.get_ncvar('TIME3')
        self.ntimep = len(self.timep)
        if self.ntimes>1 and self.ntimep>1:
            coords_x = self.get_ncvar('X')
            sh = coords_x.shape[1] if coords_x.ndim>1 else coords_x.shape
            self.logger.info(f' nprof {sh}: ntimes {self.ntimes}'
               f' {self.timesc[0]} .. {self.timesc[-1]}')
            self.logger.info(f' ntimep {self.ntimep} {self.timep[0]} .. {self.timep[-1]}')
            self.logger.info(f' Dtimep {self.timep[1]-self.timep[0]};'
               f' Dtimesc {self.timesc[1]-self.timesc[0]}')
        else:
            raise Exception(f'\tInvalid times prof:{self.ntimep} scalar {self.ntimes}')


    def check_1dim(self, dat):
         if isinstance(dat, np.ndarray):
            return dat.ndim==1
         elif isinstance(dat, (list, tuple, set)):
            return np.array(dat).ndim==1


    #Don't call this using varargs forwarded by another function
    def is_valid(self, *variables, ind:int=0, name:str=''):
        ii = ind+1
        for var in variables:
          if var is None:
              if ind==1:
                  self.logger.info(f'{name} is None')
              return False
          elif isinstance(var, np.ndarray):
              #ndim=0; eg:np.array(5)
              if var.ndim==0:
                  return False
              sh = var.shape
              if sh[0]<ii or (len(sh)>1 and sh[1]==0):
                  if ind==1:
                      self.logger.info(f'{name} shape {sh}')
                  return False
          elif isinstance(var, (list, tuple, set)):
              if len(var)<ii or not any(var):
                  if ind==1:
                      self.logger.info(f'{name} len: {len(var)}')
                  return False
        return True

    def sum_arrays(self, *arrays):
        """Return the element-wise sum of non-empty array-like inputs.

        Shape/type errors are treated as invalid optional source data, while
        programming and process-control exceptions are allowed to propagate.
        """
        res = []
        init = False
        try:
            for arr in arrays:
                if len(arr) == 0:
                    continue
                if init:
                    res = res + arr
                else:
                    res = arr
                    init = True
        except (TypeError, ValueError) as exc:
            self.logger.info(f' sum arrays error: {exc}')
        return res


    def operate_arrays(self, arr1, op:str, arr2=None, name1:str='', name2:str=''):
        """Apply a supported element-wise operation with local NumPy error policy."""
        unary = op == 'sqrt'
        valid = self.is_valid(arr1) if unary else self.is_valid(arr1, arr2)
        if not valid:
            self.logger.info(f' operate_arrays error: {name1} {op} {name2}')
            return np.array([])

        try:
            with np.errstate(divide='raise', invalid='raise', over='raise'):
                if op == '-':
                    return arr1 - arr2
                if op == '+':
                    return arr1 + arr2
                if op == '/':
                    a1 = np.asarray(arr1)
                    a2 = np.asarray(arr2)
                    dtype = np.result_type(a1.dtype, a2.dtype, np.float64)
                    shape = np.broadcast_shapes(a1.shape, a2.shape)
                    return np.divide(
                        a1, a2,
                        out=np.zeros(shape, dtype=dtype),
                        where=a2 != 0,
                    )
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
            self.logger.info(f' operate_arrays error:: {name1} {op} {name2}: {exc}')
            return np.array([])

        self.logger.info(f' operate_arrays error:: unsupported operation {op!r}')
        return np.array([])

    def np_multiply_sum(self, arr1, arr2, axis:int, name1:str='', name2:str='', scale=None):
        if not self.is_valid(arr1, arr2):
            mes = ' np_multiply_sum error'
            s1, s2 = len(arr1), len(arr2)
            if isinstance(arr1, np.ndarray) and isinstance(arr2, np.ndarray):
               s1, s2 = arr1.shape, arr2.shape
            self.logger.info(f'{mes}: {name1}:{s1} {name2}: {s2}')
            return np.array([])
        try:
            with np.errstate(invalid='raise', over='raise'):
                result = np.sum(np.multiply(arr1, arr2), axis=axis)
                return result if scale is None else scale * result
        except (TypeError, ValueError, FloatingPointError, OverflowError) as exc:
            self.logger.info(f' np_multiply_sum error:: {name1} {name2}: {exc}')
            return np.array([])


    def checkErrExit(self, *variables, name:str=''):
        for var in variables:
            stat = self.is_valid(var)
            if stat is False:
               self.logger.exception(f'\t\tvar {name} is invalid')
               raise Exception(f'**ERROR**: {name} is invalid: len={len(var)}: {var}')


    def interp1d_time(self, time_new, time, dat, ind1:int=None, ind2:int=None):
        r'''
        interpolate dat from time onto time_new
        time_new: 1d times; output dat has this size.
        dat: 2d (time, 1d)
        time: 1d times (x axis); Same size of 1st dim or values or f(x) of dat
        '''
        #if dat is 2D, 1st dim should have len of len(time)
        if self.is_valid(time_new) and self.is_valid(dat, time, ind=len(time)-1):
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

    def get_ncvar_on_time(self, var: str, target: str = 'TIME3', scale: float = 1.0):
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
            dimensions = self.nc.variables[var].dimensions
        except KeyError:
            return np.array([])
        if len(dimensions) != 1 or dimensions[0] not in ('TIME', 'TIME3'):
            raise ValueError(f'{var} is not a 1-D TIME/TIME3 history: {dimensions}')
        source = np.asarray(self.timesc if dimensions[0] == 'TIME' else self.timep)
        destination = np.asarray(self.timesc if target == 'TIME' else self.timep)
        values = self.get_ncvar(var, scale=scale)
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

    #fill ids.value
    def fill_value(self, ids, dat, rsz:int=-1, mult:float=None, name:str=None):
        name_ = name if name else 'data'
        if not self.is_valid(dat):
            #self.logger.info(f'{name} is None')
            return

        nresize = len(dat) if (not np.isscalar(dat) and rsz==-1) else rsz

        if nresize > 0:
            ids.value.resize(nresize, refcheck=False)
        if mult is not None:
            dat = mult*dat
        if np.isscalar(dat):
            ids.value = dat
        else:
            ids.value[:] = dat
        #if self.src:
        #    ids.source = self.src

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
        path = getattr(getattr(node, 'metadata', None), 'path', None)
        return str(path) if path is not None else type(node).__name__

    def fill_data(self, ids, dat, mult:float=None, ind:int=-1, ind2:int=-1,
          interp:bool=False, name:str=None, resize:bool=True, chk:int=-1) -> bool:
        """Assign ``dat`` to an array node; returns True if it was written.

        Absent source data (``is_valid`` false) is optional input: the field is skipped and recorded in
        the mapping ledger.  Data that is present but cannot be assigned (radial size different from
        ``chk``, a scalar for an array node, a shape or schema error) is a broken mapping, not absent
        data: it raises ``MappingError`` (strict) or is recorded as skipped with the reason (lenient).
        Before, a size mismatch was only reported for time slice 0 and silently skipped for the others.
        """
        nm = name if name is not None else 'data'
        target = self._node_path(ids)
        if not self.is_valid(dat, name=nm):
            self.mapping.skipped('no source data', target=target)
            if ind==0:
                self.logger.info(f'\t\t{nm} is invalid')
            return False
        try:
            if (ind>= 0 and ind2>=0):
               dat = dat[ind:ind2]
            elif ind>= 0:
               dat = dat[ind]

            if chk>0 and not np.isscalar(dat) and len(dat)!=chk:
               raise ValueError(f'size mismatch {len(dat)} != {chk}')
            if np.isscalar(dat):
               raise ValueError('cannot fill a scalar value into an array node')
            if mult is not None:
                dat = mult*dat
            if interp is True and not np.array_equal(self.timesc, self.timep):
                dat = np.interp(self.timep, self.timesc, dat)
            if (not np.isscalar(dat)) and resize:
                ids.resize(len(dat), refcheck=False)
            ids[:] = dat
        except (AttributeError, ValueError, TypeError, IndexError, FloatingPointError) as exc:
            self.mapping.broken(f'{type(exc).__name__}: {exc}', source=nm, target=target)
            return False
        return True

    def load_namelist(self):
        if self.nml_name is None:
            print('**Error: namelist file not set')
            return
        self.nmlist = {}
        if not check_file_exists(self.nml_name):
            print(f'**WARNING: namelist file {self.nml_name} not found')
            return

        with open(self.nml_name, encoding = 'ISO-8859-1') as nmf:
            self.logger.info(f'loading {self.nml_name}')
            for line in nmf:
                lns = line.strip()
                if lns and (lns[0] not in ['!', '&', '/']):
                    lns = lns.split('!',1)[0]
                    kv = lns.split('=',1)
                    if len(kv) == 2:
                        val = kv[1].strip()
                        #cast to int or float
                        try:
                            val = int(val)
                        except ValueError:
                            try:
                                vv = val.upper().replace('D0','')
                                val = float(vv)
                            except ValueError:
                                pass
                        self.nmlist[kv[0].strip()] = val

    def reload_namelist(self, fname:str):
        self.nml_name = fname
        if self.nmlist is not None:
            self.nmlist.clear()
        self.load_namelist()
        # A changed TR.DAT must not retain XML/provenance for the prior file.
        self._code_parameters_xml_cache = None
        self._provenance_file_records_cache = None
        self._input_inventory_cache = None


    #TODO move to utilities
    def parse_number(self, ss:str, tp):
        if not isinstance(ss, str):
            return tp(ss)

        #replace exponential 10** if any
        ss0 = ss
        res= ss.find('10**')
        if res >=0:
            ss = ss.replace('10**', 'E')
        opls = ['*', '/', '+', '-']
        nn = []
        ops = []
        p = 0

        #splitting into separate number strings
        for i,a in enumerate(ss):
            if a in opls and i>0 and ss[i-1].upper() !='E':
                nn.append(ss[p:i])
                ops.append(a)
                p = i+1

        if p<len(ss):
            nn.append(ss[p:])
        if len(ops)!=len(nn)-1:
            raise Exception(f'Error in parsing string to number {ss0}')
        n = 0
        try:
            n = tp(nn[0]) if len(nn)>0 else None
            for i,s in enumerate(ops):
                d = tp(nn[i+1])
                if s=='*':
                    n = n * d
                elif s=='/':
                    n = n/d
                elif s=='+':
                    n = n+d
                elif s=='-':
                    n = n-d
        except (ValueError, TypeError, IndexError):
            print(f'parse error: {ss0} extracted: {nn}')
        return n


    #name's last character is discarded. Eg. XECECH?. () added separately
    def get_namelist_nvars(self, name:str, ndig:int=2):
        if name[-1] != '?':
            raise Exception(f'Error: {name} needs "?" as last char')
        name = name[:-1]
        pattern = rf'^{re.escape(name)}(?:\((\d{{1,{ndig}}})\))?$'
        if self.nmlist is None:
            self.logger.info(f'**Error: {name}: namelist not found')
            return 0
        try:
            matches = [self.nmlist[k] for k in self.nmlist 
                       if re.fullmatch(pattern, k)]
            return len(matches)
        except (AttributeError, TypeError):
            self.logger.info(f'**Error: {name}: namelist not found')
            return 0


    def make_namelist_re_pattern(self, name:str, nd1:int=0, nd2:int=0, nch1:int=0,
          nch2:int=0, nd:int=0, nch:int=0) -> str:
        #only upto ?,? handled (two occurences)
        nn = name.count('?')
        st = 1 if nn>1 else -1
        pat = self.make_re_pattern(name, nd1, nd2, nch1, nch2, nd, nch, step=st)
        ii = name.find('?')

        if nn > 1:
            pat = self.make_re_pattern(fr'{pat}{name[ii+1:]}', nd1, nd2, nch1, nch2, nd, nch)
        return pat

    def get_namelist_keys_re(self, name:str, cut:str='', nd1:int=0, nd2:int=0, nch1:int=0,
          nch2:int=0, nd:int=0, nch:int=0) -> List:
        # name: regex need to be escaped \(
        pat = self.make_namelist_re_pattern(name, nd1, nd2, nch1, nch2, nd, nch)

        if pat is None:
            return []
        matches = [k for k in self.nmlist.keys() if re.search(pat, k)]
        if cut:
            matches = [s.split(cut,1)[1] for s in matches if len(s.split(cut,1))>1]
        return matches


    #returns scalar or list
    def get_namelist_var_(self, name:str, tp=None, scale=1., ind:int=-1):
        try:
            val = self.nmlist[name]
            if isinstance(val, str):
               val = re.split(r'[\s,;]+', val)
            if ind>=0:
               val = val[ind]
        except (KeyError, TypeError, IndexError, AttributeError):
            #self.logger.info(f'** {name} NOT found in name list')
            if self.nmlist is None:
                self.logger.info(f'** {self.nmlist} is not set')
            return None

        narr = 0
        try:
            if ind==-1:
               narr = len(val)
        except TypeError:
            pass

        if tp in [int, float]:
            if narr>0:
               val =  [self.parse_number(val[i], tp) for i in range(narr)]
               val = np.array(val)
            else:
               val = self.parse_number(val, tp)
            if scale != 1.0:
               val = scale*val
        return val

    # Some scalar can stand for arrays, e.g: FREQECH, which causes error when
    #  an array is expected. So, call for array should handle single variable also.
    # Returns scalar for single, list for array on RHS
    def get_namelist_vars(self, name:str, tp=None, scale=1.0, ind:int=-1,
          size:int=1, resize:bool=False):
        #name: don't escape (); only regex needs \(
        # search for single or without () for array on rhs
        res = self.get_namelist_var_(name, tp, scale, ind=ind)
        n = 0
        # search for index
        if ind>=0 and res is None:
           namei = f'{name}({ind+1})'
           res = self.get_namelist_var_(namei, tp, scale)
           n = 1
        # search for array
        elif res is None:
           n = self.get_namelist_nvars(f'{name}?')
           res = []
           for i in range(n):
              namei = f'{name}({i+1})'
              val = self.get_namelist_var_(namei, tp, scale)
              if val is not None:
                 res.append(val)
           res = np.array(res)

        # resize only if scalar, not an array to different size
        if res is not None and resize:
           try:
               n = len(res)
           except TypeError:
               pass
           if n<=1:
               res = np.full(size, res)

        if size>1 and res is not None and len(np.array(res))!= size:
            raise Exception(f'resize failed: size of {name} != {size}')
        return res

    def get_namelist_var(self, name:str, tp=None, scale=1.0, ind:int=-1, size:int=1):
        """Compatibility wrapper around :meth:`get_namelist_vars`."""
        return self.get_namelist_vars(name, tp, scale, ind, size)


    def extract_str_between(self, name:str, pat1:str, pat2:str, default:str, tp=None):
        res = None
        try:
            li = name.find(pat1) + len(pat1) + 1
            #can have multiple occurences. eg: ','
            ris = [s.start() for s in re.finditer(pat2, name)]
            ri = [i for i in ris if i>li][0]
            res = name[li:ri] if ri>li+1 else default
            if tp is not None:
                res = tp(res)
        except (IndexError, ValueError, TypeError):
            self.logger.debug(f'extract error: {name} doesnot have {pat1} and {pat2}')
        return res

    def fluxBdryToCenter(self, data, dat0=None):
        dat = []
        if not self.valid(data):
            return dat

        dat = dat.copy()
        #average of consecutive columns, except first
        dat[:,1:] += dat[:,:-1] #add previous column
        if self.valid(dat0):
            dat[:,0] += dat0
        else:
            dat[:,0] *= 2.
        return dat/2.

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
        arr = np.asarray(values)
        if arr.ndim < 1 or arr.size == 0:
            return arr
        return np.cumsum(arr, axis=-1)

    @staticmethod
    def interp_radial_profiles(values, x_from, x_to, axis_value=None):
        """Interpolate time-dependent radial profiles without mixing time slices.

        ``axis_value`` optionally adds an explicit point at rho=0 before
        interpolation.  This is important for TRANSP XB data because XB starts
        at the outer boundary of the first finite-volume zone, while X contains
        the zone centre inside it.
        """
        v = np.asarray(values)
        xf = np.asarray(x_from)
        xt = np.asarray(x_to)
        if v.size == 0 or xf.size == 0 or xt.size == 0 or v.ndim != 2:
            return v
        if xf.ndim != 2 or xt.ndim != 2 or v.shape[0] != xf.shape[0] or v.shape[0] != xt.shape[0]:
            return v
        out = np.empty((v.shape[0], xt.shape[1]), dtype=np.float64)
        av = None if axis_value is None else np.asarray(axis_value)
        for it in range(v.shape[0]):
            xx = xf[it]
            yy = v[it]
            if axis_value is not None and len(xx) and xx[0] > 0.0:
                a0 = float(av[it] if av.ndim else av)
                xx = np.concatenate(([0.0], xx))
                yy = np.concatenate(([a0], yy))
            out[it] = np.interp(xt[it], xx, yy)
        return out

    def vacuum_field_reference(self, times=None):
        """``(r0, b0(times))`` of the vacuum toroidal field, or ``None`` without RAXIS/BZXR.

        ``times`` defaults to the writer time grid ``timep``; see
        :func:`utilities.vacuum_field_reference` for the convention.
        """
        raxis = self.get_ncvar('RAXIS', self.cm2m)       # RAXIS(TIME) m
        bzxR = self.get_ncvar('BZXR', self.cm2m)         # BZXR(TIME) T m, vacuum field x radius
        if not self.is_valid(raxis, bzxR):
            return None
        return vacuum_field_reference(raxis, bzxR, self.bt_sign,
                                      self.timep if times is None else times, self.timesc)

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

    @staticmethod
    def _iso8601_file_timestamp(path):
        """Return a UTC ISO-8601 timestamp for a local provenance file."""
        try:
            ts = Path(path).stat().st_mtime
            return datetime.fromtimestamp(ts, tz=timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
        except (OSError, ValueError, OverflowError):
            return ''

    def _provenance_file_records(self):
        """Return portable source-file records for IDS-level provenance.

        Keep the whole-IDS provenance intentionally small (the DD caps the
        reference list at 10). Full file contents and recursively referenced
        templates remain preserved in ``code.parameters``.
        """
        if self._provenance_file_records_cache is not None:
            return list(self._provenance_file_records_cache)

        records = []
        candidates = [('TRANSP CDF', self.ncfile)]
        if self.nml_name:
            candidates.append(('TRANSP TR.DAT', self.nml_name))
        # Include resolved model/template files when present, up to the DD cap.
        try:
            for role, key, filename, resolved in self._collect_parameter_files():
                if role == 'TRDAT':
                    continue
                if resolved is not None:
                    candidates.append((f'{role}:{key}', str(resolved)))
        except Exception as exc:          # best-effort provenance: never abort, but never silent either
            self.mapping.skipped(f'could not collect model/template files: {type(exc).__name__}: {exc}',
                                 target='ids_properties/provenance')
        seen = set()
        for role, filename in candidates:
            if not filename:
                continue
            path = Path(filename).expanduser()
            try:
                resolved = path.resolve()
            except (OSError, RuntimeError):
                resolved = path
            key = str(resolved)
            if key in seen:
                continue
            seen.add(key)
            digest = ''
            try:
                import hashlib
                h = hashlib.sha256()
                with open(resolved, 'rb') as stream:
                    for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                        h.update(chunk)
                digest = h.hexdigest()
            except OSError:
                pass            # unreadable provenance file: no digest (recorded as such below)
            # Avoid embedding machine-specific absolute paths in the portable
            # provenance string. The exact resolved path remains in code.parameters.
            name = f'{role}:file:{Path(filename).name}'
            if digest:
                name += f'#sha256={digest}'
            records.append((name, self._iso8601_file_timestamp(resolved)))
            if len(records) >= 10:
                break
        self._provenance_file_records_cache = list(records)
        return list(records)

    def _fill_source_provenance(self, ids):
        """Attach source CDF/TR.DAT provenance to the whole IDS.

        DD <=3.41 uses ``provenance/node/sources``. DD >=3.42 uses the
        ``reference`` AoS with name/timestamp. The full original inputs are
        additionally embedded in ``code.parameters`` wherever that node exists.
        """
        try:
            provenance = ids.ids_properties.provenance
            provenance.node.resize(1)
            node = provenance.node[0]
            if hasattr(node, 'path'):
                node.path = ''
            records = self._provenance_file_records()
            if hasattr(node, 'reference'):
                node.reference.resize(len(records))
                for i, (name, timestamp) in enumerate(records):
                    ref = node.reference[i]
                    if hasattr(ref, 'name'):
                        ref.name = name
                    if timestamp and hasattr(ref, 'timestamp'):
                        ref.timestamp = timestamp
            elif hasattr(node, 'sources'):
                node.sources = [name for name, _ in records]
        except Exception as exc:          # best-effort provenance: never abort, but never silent either
            self.logger.warning(f'Could not attach IDS provenance: {exc}')
            self.mapping.skipped(f'could not attach IDS provenance: {type(exc).__name__}: {exc}',
                                 target='ids_properties/provenance')
        return ids

    def fill_properties(self, ids):
        runid, device, year, cdf_date, build_date, tglf, glf = self.get_nc_global_atts()
        comment = f'TRANSP results for run {runid} {device} year: {year} CDF date: {cdf_date}'
        cocos = (f' DD={self.dd_version}; COCOS={self.target_cocos};'
                 f' NLJCCW={self.nljccw}; NLBCCW={self.nlbccw};'
                 f' Ip_sign={self.ip_sign:+.0f}; Bphi_sign={self.bt_sign:+.0f};'
                 f' q_sign={self.q_sign:+.0f}; psi_sign={self.psi_sign:+.0f}')
        ids.ids_properties.comment = comment + f' build: {build_date} tglf: {tglf} glf23: {glf}' + cocos
        ids.ids_properties.creation_date = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
        if hasattr(ids.ids_properties, 'name'):
            ids.ids_properties.name = f'TRANSP {runid}'
        ids = self._fill_source_provenance(ids)

        # Do not fill the legacy top-level fields
        #   ids_properties.data_dictionary
        #   ids_properties.access_layer
        # They are not present in current IMAS Data Dictionary versions (including
        # DD 4.1.1).  In modern IMAS they live below ids_properties.version_put.
        # More importantly, IMAS-Python/DBEntry.put() fills version_put from the
        # Data Dictionary and backend actually used for the write.  Leaving this
        # metadata to DBEntry therefore avoids stale or inconsistent version data
        # and remains compatible with older DD versions as well.
        return ids

    @staticmethod
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

    @staticmethod
    def _xml_safe_text(text):
        """Remove characters forbidden by XML 1.0 while preserving file text."""
        if text is None:
            return ''
        return ''.join(ch if (ch in '\t\n\r' or ord(ch) >= 0x20) else '\ufffd'
                       for ch in str(text))

    @staticmethod
    def _clean_namelist_path(value):
        if value is None:
            return ''
        text = str(value).strip().rstrip(',').strip()
        if len(text) >= 2 and text[0] in ("'", '"') and text[-1] == text[0]:
            text = text[1:-1]
        return text.strip()

    def _resolve_input_file(self, filename):
        """Resolve a TRANSP/template path relative to the TRANSP namelist."""
        if not filename:
            return None
        candidate = Path(os.path.expandvars(str(filename))).expanduser()
        if candidate.is_absolute() and candidate.exists():
            return candidate.resolve()
        roots = []
        if self.nml_name:
            try:
                roots.append(Path(self.nml_name).expanduser().resolve().parent)
            except (OSError, RuntimeError):
                roots.append(Path(self.nml_name).expanduser().parent)
        roots.append(Path.cwd())
        for root in roots:
            path = (root / candidate).resolve()
            if path.exists():
                return path
        # Return the most meaningful unresolved location for provenance.
        return (roots[0] / candidate).resolve() if roots else candidate

    def _read_parameter_file(self, path, max_bytes=10 * 1024 * 1024):
        """Read a text parameter file and return (text, sha256, status)."""
        try:
            path = Path(path)
            if not path.exists() or not path.is_file():
                return '', '', 'missing'
            size = path.stat().st_size
            if size > max_bytes:
                return '', '', f'too_large:{size}'
            raw = path.read_bytes()
            digest = hashlib.sha256(raw).hexdigest()
            try:
                text = raw.decode('utf-8')
            except UnicodeDecodeError:
                text = raw.decode('ISO-8859-1')
            return self._xml_safe_text(text), digest, 'embedded'
        except Exception as exc:
            return '', '', f'error:{type(exc).__name__}'

    @staticmethod
    def _append_namelist_xml(parent, text):
        """Add a lightweight structured Fortran-namelist representation.

        The original text is retained separately and remains authoritative.
        This parser is intentionally conservative: it records groups and
        assignment strings without trying to evaluate Fortran expressions.
        Continuation text is appended to the preceding assignment.
        """
        parsed = ET.SubElement(parent, 'parsed_namelist')
        group_name = 'GLOBAL'
        group = ET.SubElement(parsed, 'group', {'name': group_name})
        last_parameter = None

        def strip_comment(line):
            # Remove ! comments while respecting simple quoted strings.
            quote = None
            out = []
            for ch in line:
                if ch in ("'", '"'):
                    if quote is None:
                        quote = ch
                    elif quote == ch:
                        quote = None
                if ch == '!' and quote is None:
                    break
                out.append(ch)
            return ''.join(out)

        for raw_line in str(text).splitlines():
            line = strip_comment(raw_line).strip()
            if not line:
                continue
            if line.startswith('&'):
                group_name = line[1:].strip().split()[0] or 'GLOBAL'
                group = ET.SubElement(parsed, 'group', {'name': group_name})
                last_parameter = None
                continue
            if line == '/' or line.startswith('/ '):
                group_name = 'GLOBAL'
                group = ET.SubElement(parsed, 'group', {'name': group_name})
                last_parameter = None
                continue
            if '=' in line:
                name, value = line.split('=', 1)
                name = name.strip()
                if not name:
                    continue
                last_parameter = ET.SubElement(group, 'parameter', {
                    'name': name,
                    'value': value.strip(),
                })
            elif last_parameter is not None:
                prior = last_parameter.get('value', '')
                last_parameter.set('value', (prior + ' ' + line).strip())

    def input_inventory(self):
        """Resolve and cache the run's TR.DAT and recursively referenced inputs.

        Missing or oversized files remain explicit provenance records; they are
        *not* silently counted as archived inputs.
        """
        if self._input_inventory_cache is None:
            self._input_inventory_cache = InputInventory(
                self.nml_name, runid=self.get_ncattr_text('Runid', default=''))
            for item in self._input_inventory_cache.missing:
                message = (f'Source input NOT embedded ({item.status}): '
                           f'{item.key} -> {item.path}')
                self.logger.warning(message)
                print(f'*** WARNING: {message}')
            if self._input_inventory_cache.truncated:
                message = ('Source input discovery exceeded safety limits; '
                           'input inventory is incomplete')
                self.logger.warning(message)
                print(f'*** WARNING: {message}')
        return self._input_inventory_cache

    def _collect_parameter_files(self):
        """Legacy API consumed by IDS-provenance summary; use new inventory."""
        return [(item.role, item.key, item.reference, item.path)
                for item in self.input_inventory().records]

    def fill_external_model_inputs(self, ids, family):
        """Attach model-specific XML to the appropriate IDS code libraries."""
        if not hasattr(ids, 'code'):
            return []
        return fill_model_libraries(ids.code, self.input_inventory(), family)

    def build_code_parameters_xml(self):
        """Build XML provenance for IMAS code.parameters.

        IMAS defines code.parameters as code-specific parameters in XML format.
        Preserve the complete TRANSP input and referenced model/template input
        files verbatim (XML escaped), together with hashes and CDF metadata.
        """
        if self._code_parameters_xml_cache is not None:
            return self._code_parameters_xml_cache

        root = ET.Element('transp2imas_parameters', {
            'format_version': '1',
            'target_dd': str(self.dd_version),
            'target_cocos': str(self.target_cocos),
        })
        run = ET.SubElement(root, 'run', {
            'device': str(self.device),
            'shot': str(self.shot),
            'imas_run': str(self.run),
            'occurrence': str(self.occur),
        })
        run.set('source_cdf', str(self.ncfile))
        ET.SubElement(root, 'converter', {
            'name': 'transp2imas',
            'version': TRANSP2IMAS_VERSION,
        })

        metadata = ET.SubElement(root, 'transp_metadata')
        for key in ('Runid', 'CDF_date', 'BUILD_date', 'TRANSP_version',
                    'TRANSP_DOI', 'TGLF_hash', 'GLF23_hash'):
            try:
                value = self.get_ncattr(key)
            except (AttributeError, KeyError, RuntimeError):
                value = None
            if value is None:
                continue
            if isinstance(value, bytes):
                value = value.decode('ISO-8859-1', errors='replace')
            ET.SubElement(metadata, 'item', {'name': key, 'value': str(value).strip()})

        signs = ET.SubElement(root, 'coordinate_conventions')
        for key, value in (
            ('NLJCCW', self.nljccw), ('NLBCCW', self.nlbccw),
            ('ip_sign', self.ip_sign), ('bphi_sign', self.bt_sign),
            ('q_sign', self.q_sign), ('psi_sign', self.psi_sign),
            ('phi_sign', self.phi_sign)):
            ET.SubElement(signs, 'item', {'name': key, 'value': str(value)})

        inventory = self.input_inventory()
        inputs = ET.SubElement(root, 'input_files', {
            'complete': str(not (inventory.missing or inventory.truncated)).lower(),
            'discovery_truncated': str(inventory.truncated).lower(),
        })
        for record in inventory.records:
            append_file_xml(inputs, record)

        xml = ET.tostring(root, encoding='unicode')
        self._code_parameters_xml_cache = xml
        return xml

    def fill_code_properties(self, ids):
        """Fill common code provenance, if this IDS contains a code node."""
        if not hasattr(ids, 'code'):
            return ids
        ids.code.name = 'transp2imas'
        if hasattr(ids.code, 'description'):
            ids.code.description = 'TRANSP CDF to IMAS converter'
        if hasattr(ids.code, 'version'):
            ids.code.version = TRANSP2IMAS_VERSION
        if hasattr(ids.code, 'repository'):
            ids.code.repository = 'https://github.com/PrincetonUniversity/transp2imas'
        if hasattr(ids.code, 'parameters'):
            ids.code.parameters = self.build_code_parameters_xml()
        return ids


if __name__ == '__main__':
    print("Please use: python convert.py")
