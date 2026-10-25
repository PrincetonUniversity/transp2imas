#!/usr/bin/python
"""Native CDF access for curation, independent of conversion and IMAS lifecycle.

No interpolation, unit guessing, profile repair, or in-place source modification.
Small, copied arrays are cached read-only; fill/missing values remain unavailable.
"""

from collections import OrderedDict
from contextlib import contextmanager
import hashlib
from pathlib import Path

import numpy as np


def text(value):
    if isinstance(value, bytes):
        return value.decode('utf-8', errors='replace').strip()
    return str(value).strip()


def digest_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


@contextmanager
def open_cdf(path):
    """Prefer netCDF4; use SciPy only for classic CDF1/2 when netCDF4 is absent.

    Corrupt files are not silently retried with another parser. SciPy's classic
    reader does not apply scale_factor/add_offset; NativeSource does that itself.
    """
    path = Path(path)
    try:
        from netCDF4 import Dataset
    except ImportError:
        with path.open('rb') as stream:
            signature = stream.read(4)
        if signature not in (b'CDF\x01', b'CDF\x02'):
            raise RuntimeError('netCDF4 is required for NetCDF4/HDF5 or CDF5 files')
        from scipy.io import netcdf_file
        handle = netcdf_file(str(path), 'r', mmap=False)
        backend = 'scipy.io.netcdf_file'
        packed = True
    else:
        handle = Dataset(str(path), 'r')
        # Use one explicit implementation of masking/scaling for both backends.
        handle.set_auto_maskandscale(False)
        backend = 'netCDF4'
        packed = True
    try:
        yield NativeSource(handle, path=path, backend=backend, packed=packed)
    finally:
        handle.close()


class SourceDataError(ValueError):
    """The source variable exists but cannot be decoded: corrupt or unsupported data, not absent evidence."""


def _numeric_attribute(name, variable, attr):
    """The attribute as a flat float array; metadata that is not numeric is corrupt data, not a missing input."""
    try:
        return np.asarray(getattr(variable, attr), dtype=float).ravel()
    except (TypeError, ValueError):
        raise SourceDataError(f'{name}: invalid numeric {attr}')


def _scalar_attribute(name, variable, attr, default):
    if not hasattr(variable, attr):
        return default
    values = _numeric_attribute(name, variable, attr)
    if values.size != 1:
        raise SourceDataError(f'{name}: {attr} must be a single number')
    return float(values[0])


class NativeSource:
    """A read-only view over any dataset exposing a ``variables`` mapping.

    Tests can supply an in-memory dataset. The byte limits apply to decoded
    variables and the cache, not the backend's own file buffering.
    """

    def __init__(self, dataset, *, path=None, backend='memory', packed=False,
                 max_variable_bytes=64 * 1024 * 1024, cache_bytes=32 * 1024 * 1024):
        self.dataset = dataset
        self.path = Path(path) if path is not None else None
        self.backend = backend
        self.packed = packed
        self.max_variable_bytes = max_variable_bytes
        self.cache_bytes = cache_bytes
        self._cache = OrderedDict()
        self._cache_size = 0
        self.read_counts = {}

    def has(self, name):
        return name in self.dataset.variables

    def meta(self, name):
        variable = self.dataset.variables[name]
        return {'name': name, 'shape': list(variable.shape),
                'dimensions': list(variable.dimensions),
                'units': text(getattr(variable, 'units', '')),
                'long_name': text(getattr(variable, 'long_name', ''))}

    def array(self, name):
        if name in self._cache:
            self._cache.move_to_end(name)
            return self._cache[name]
        variable = self.dataset.variables[name]
        count = 1
        for length in variable.shape:
            count *= int(length)
        if count * 8 > self.max_variable_bytes:
            raise SourceDataError(f'{name}: decoded variable exceeds curation memory limit')
        raw = variable[...]  # works for scalar and array variables
        if np.asanyarray(raw).dtype.kind not in 'biuf':
            raise SourceDataError(f'{name}: curation requires real numeric data')
        values = np.array(np.ma.getdata(raw), dtype=float, copy=True)
        mask = np.ma.getmaskarray(raw).copy()
        for attr in ('_FillValue', 'missing_value'):
            if hasattr(variable, attr):
                for fill in np.asarray(getattr(variable, attr)).ravel():
                    try:
                        mask |= values == float(fill)
                    except (TypeError, ValueError):
                        raise SourceDataError(f'{name}: invalid numeric {attr}')
        if hasattr(variable, 'valid_range'):
            bounds = _numeric_attribute(name, variable, 'valid_range')
            if len(bounds) != 2 or not np.all(np.isfinite(bounds)) or bounds[0] > bounds[1]:
                raise SourceDataError(f'{name}: invalid valid_range metadata')
            mask |= (values < bounds[0]) | (values > bounds[1])
        for attr, compare in (('valid_min', np.less), ('valid_max', np.greater)):
            if hasattr(variable, attr):
                boundary = _scalar_attribute(name, variable, attr, np.nan)
                if not np.isfinite(boundary):
                    raise SourceDataError(f'{name}: nonfinite {attr}')
                mask |= compare(values, boundary)
        if self.packed:
            if text(getattr(variable, '_Unsigned', '')).lower() == 'true':
                raise SourceDataError(f'{name}: packed unsigned reinterpretation is not implemented')
            factor = _scalar_attribute(name, variable, 'scale_factor', 1.0)
            offset = _scalar_attribute(name, variable, 'add_offset', 0.0)
            if not np.isfinite(factor) or not np.isfinite(offset):
                raise SourceDataError(f'{name}: nonfinite packing metadata')
            values = values * factor + offset
        values[mask] = np.nan
        values.setflags(write=False)
        self.read_counts[name] = self.read_counts.get(name, 0) + 1
        if values.nbytes <= self.cache_bytes:
            while self._cache and self._cache_size + values.nbytes > self.cache_bytes:
                _, old = self._cache.popitem(last=False)
                self._cache_size -= old.nbytes
            self._cache[name] = values
            self._cache_size += values.nbytes
        return values

    def provenance(self):
        if self.path is None:
            return {'path': None, 'backend': self.backend, 'sha256': None,
                    'note': 'in-memory source; not a file-backed run'}
        return {'path': str(self.path.resolve()), 'backend': self.backend,
                'size_bytes': self.path.stat().st_size,
                'sha256': digest_file(self.path)}


def unit_key(units):
    return text(units).upper().replace(' ', '').replace('**', '^')


def density_scale(units):
    """Explicit native-density unit conversion; never infer units from magnitude."""
    key = unit_key(units)
    if key in ('N/CM^3', '#/CM^3', 'CM^-3', 'CM-3', '1/CM^3', 'N/CM3', '#/CM3'):
        return 1.0e6
    if key in ('N/M^3', '#/M^3', 'M^-3', 'M-3', '1/M^3', 'N/M3', '#/M3'):
        return 1.0
    raise ValueError(f'unsupported or absent density units: {units!r}')
