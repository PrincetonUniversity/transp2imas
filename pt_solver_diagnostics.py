#!/usr/bin/python
"""Lossless discovery of TRANSP PT_SOLVER runtime diagnostic CDF arrays.

This records raw native CDF numbers, not standardized IMAS transport coefficients.
Only the explicitly verified channel names in ``KNOWN_CHANNELS`` are accepted (``NPTRITER``
and the ``PT_*`` timers), so ordinary plasma/physics channels are not inadvertently labeled
as solver diagnostics.

Every channel keeps the values of its *own* native coordinate (``TIME`` for all of them).  That
is not the ``TIME3`` grid of the surrounding IDS: the two differ in length in some CDFs
(DIII-D 147634T13 has 4 ``TIME`` and 7 ``TIME3`` samples), so index ``i`` of a channel is
not index ``i`` of the IDS time base.

Discovery does not depend on the PT_SOLVER input provenance: a CDF may hold solver output
without an available model namelist.  Without a PT plan a channel is recorded only if it has
nonzero samples, because the channels exist (all zero) in every CDF, including runs that
never used the solver.
"""
import base64
import hashlib
import xml.etree.ElementTree as ET
import zlib

import numpy as np

# Verified against long_name metadata in the 11 bundled CDF examples.
# PT_TIME is an output-definition scalar (value 1), not a runtime trace.
# Other iteration counts (ISOITER, RZITER) refer to different solvers.
KNOWN_CHANNELS = {
    'NPTRITER': ('newton_iterations', 'count'),
    'PT_CPTIME': ('cumulative_cpu_time', 'hour'),
    'PT_WTIME': ('cumulative_wall_time', 'hour'),
}
# Allow explicit future additions only after checking source CDF documentation.
EXPECTED_UNAVAILABLE = ('newton_residual_norm', 'peclet_stabilization_correction',
                        'nonlinear_convergence_status')

_MAX_BYTES_PER_CHANNEL = 16 * 1024 * 1024
_MAX_TOTAL_BYTES = 64 * 1024 * 1024


def diagnostic_xml(nc, *, max_bytes=_MAX_TOTAL_BYTES, record_inactive=True):
    """Return self-contained XML (metadata + reversible numeric arrays).

    Values stored via base64(zlib(little-endian float64 C-order bytes)).
    Empty or nonfinite-only variables are omitted; mixed nonfinite retained.
    Arrays exceeding limits are reported as skipped without loading into XML.  Each channel
    lists, per dimension, the CDF coordinate variable of that name (``coordinates``); the
    coordinate values are stored once as ``<coordinate>`` elements, in the same encoding.
    All-zero channels are recorded only with ``record_inactive`` (a PT plan exists).
    """
    root = ET.Element('pt_solver_runtime_diagnostics', {
        'kind': 'native_CDF_raw_arrays',
        'schema': 'transp2imas.pt_solver_diagnostics.v1',
        'availability': 'only_matching_channels_found_in_CDF',
        'not_standardized_imas_fields': 'true',
        'selection': 'verified_explicit_CDF_names',
    })
    variables = getattr(nc, 'variables', None)
    if variables is None:
        root.set('status', 'CDF_variable_inventory_unavailable')
        return root
    used = 0
    coordinates = {}          # coordinate name -> (values, units)
    for name in sorted(variables):
        if str(name).upper() not in KNOWN_CHANNELS:
            continue
        variable = variables[name]
        try:
            shape = tuple(int(n) for n in variable.shape)
            dtype = np.dtype(getattr(variable, 'dtype', getattr(getattr(variable, 'data', None), 'dtype', None)))
            if dtype.kind not in 'ifub' or len(shape) > 3:
                raise ValueError('unsupported dtype/dimensionality')
            estimated = int(np.prod(shape, dtype=np.int64)) * 8
            if estimated > _MAX_BYTES_PER_CHANNEL or used + estimated > max_bytes:
                raise ValueError('per-channel or total size limit')
            raw = variable[:] if shape else variable.getValue()
            if np.ma.isMaskedArray(raw):
                raw = raw.filled(np.nan)
            values = np.asarray(raw, dtype='<f8')
            if not values.size or not np.isfinite(values).any():
                raise ValueError('no finite samples')
            if not record_inactive and not np.any(values[np.isfinite(values)]):
                raise ValueError('all finite samples are zero (solver inactive or channel unused)')
            dims = tuple(str(d) for d in getattr(variable, 'dimensions', ()))
            channel_coordinates = []
            for axis, dim in enumerate(dims):
                coordinate = _coordinate(variables, dim, shape[axis] if axis < len(shape) else None)
                if coordinate is not None:
                    coordinates.setdefault(dim, coordinate)
                channel_coordinates.append(dim if coordinate is not None else '')
            data = values.tobytes(order='C')
            kind, canonical_unit = KNOWN_CHANNELS[str(name).upper()]
            attrs = {'name': str(name), 'quantity': kind, 'canonical_unit': canonical_unit, 'shape': ','.join(map(str,shape)),
                     'source_dtype': str(dtype), 'stored_dtype': 'little_endian_float64',
                     'dimensions': ','.join(dims),
                     'coordinates': ','.join(channel_coordinates),
                     'encoding': 'base64+zlib', 'sha256': hashlib.sha256(data).hexdigest(),
                     'finite_count': str(int(np.isfinite(values).sum()))}
            for key in ('units', 'long_name'):
                value = getattr(variable, key, None)
                if isinstance(value, bytes):      # SciPy returns NetCDF text attributes as bytes
                    value = value.decode('utf8', 'replace')
                if value is not None:
                    attrs[key] = str(value).strip()
            ET.SubElement(root, 'channel', attrs).text = base64.b64encode(zlib.compress(data, 9)).decode('ascii')
            used += len(data)
        except (ValueError, TypeError, OverflowError, MemoryError, AttributeError) as error:
            ET.SubElement(root, 'skipped', {'name': str(name), 'reason': str(error)})
    for dim, (values, units) in sorted(coordinates.items()):
        data = values.tobytes(order='C')
        attrs = {'name': dim, 'shape': str(values.size), 'stored_dtype': 'little_endian_float64',
                 'encoding': 'base64+zlib', 'sha256': hashlib.sha256(data).hexdigest()}
        if units:
            attrs['units'] = units
        ET.SubElement(root, 'coordinate', attrs).text = base64.b64encode(zlib.compress(data, 9)).decode('ascii')
    root.set('count', str(len(root.findall('channel'))))
    for quantity in EXPECTED_UNAVAILABLE:
        ET.SubElement(root, 'unavailable', {'quantity': quantity, 'reason': 'not_mapped_from_verified_CDF_inventory'})
    if not root.findall('channel'):
        root.set('status', 'no_matching_runtime_channels')
    return root


def _coordinate(variables, name, length):
    """Values and units of the 1-D CDF variable ``name`` if it can be the coordinate of an axis.

    Returns ``None`` when there is no such variable or its length is not the axis length.
    """
    variable = variables.get(name) if hasattr(variables, 'get') else None
    if variable is None:
        return None
    try:
        values = np.asarray(variable[:], dtype='<f8')
    except (ValueError, TypeError):
        return None
    if values.ndim != 1 or values.size != length or not np.isfinite(values).all():
        return None
    units = getattr(variable, 'units', None)
    if isinstance(units, bytes):
        units = units.decode('utf8', 'replace')
    return values, ('' if units is None else str(units).strip())


def append_runtime_diagnostics_xml(parameters, nc, *, record_inactive=True):
    """Append the runtime-diagnostics element to the ``code.parameters`` XML.

    With ``record_inactive`` (a PT_SOLVER plan exists) the element is always written, with
    an explicit status when nothing matched.  Without it, nothing is added unless at least one
    channel with nonzero data was found, and ``parameters`` is returned unchanged otherwise.
    """
    if hasattr(parameters, 'value'):
        parameters = parameters.value
    if isinstance(parameters, bytes):
        parameters = parameters.decode('utf8')
    diagnostics = diagnostic_xml(nc, record_inactive=record_inactive)
    if not record_inactive and not diagnostics.findall('channel'):
        return parameters
    root = ET.fromstring(parameters) if parameters else ET.Element('transp2imas_parameters')
    for prev in root.findall('pt_solver_runtime_diagnostics'):
        root.remove(prev)
    root.append(diagnostics)
    return ET.tostring(root, encoding='unicode')


def channel_coordinates(diagnostics, channel):
    """``{dimension: values}`` of the coordinates recorded for ``channel`` (checksum verified)."""
    stored = {c.attrib['name']: c for c in diagnostics.findall('coordinate')}
    found = {}
    for dim in (d for d in channel.attrib.get('coordinates', '').split(',') if d):
        found[dim] = decode_channel(stored[dim])
    return found


def decode_channel(channel):
    """Reconstruct a channel exactly as stored, verifying the SHA-256 digest."""
    data = zlib.decompress(base64.b64decode(channel.text))
    if hashlib.sha256(data).hexdigest() != channel.attrib['sha256']:
        raise ValueError('CDF diagnostic checksum mismatch')
    shape = tuple(int(v) for v in channel.attrib['shape'].split(',') if v)
    return np.frombuffer(data, dtype='<f8').reshape(shape)
