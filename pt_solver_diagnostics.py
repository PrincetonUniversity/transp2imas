#!/usr/bin/python
"""Lossless, opt-in discovery of TRANSP PT_SOLVER runtime diagnostic CDF arrays.

This records raw native CDF numbers, not standardized IMAS transport coefficients.
Only explicitly PT_SOLVER-prefixed numeric variables are accepted, so ordinary
plasma/physics channels are not inadvertently labeled as solver diagnostics.
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


def diagnostic_xml(nc, *, max_bytes=_MAX_TOTAL_BYTES):
    """Return self-contained XML (metadata + reversible numeric arrays).

    Values stored via base64(zlib(little-endian float64 C-order bytes)).
    Empty or nonfinite-only variables are omitted; mixed nonfinite retained.
    Arrays exceeding limits are reported as skipped without loading into XML.
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
            data = values.tobytes(order='C')
            kind, canonical_unit = KNOWN_CHANNELS[str(name).upper()]
            attrs = {'name': str(name), 'quantity': kind, 'canonical_unit': canonical_unit, 'shape': ','.join(map(str,shape)),
                     'source_dtype': str(dtype), 'stored_dtype': 'little_endian_float64',
                     'dimensions': ','.join(map(str, getattr(variable, 'dimensions', ()))),
                     'encoding': 'base64+zlib', 'sha256': hashlib.sha256(data).hexdigest(),
                     'finite_count': str(int(np.isfinite(values).sum()))}
            for key in ('units', 'long_name'):
                value = getattr(variable, key, None)
                if value is not None:
                    attrs[key] = str(value)
            ET.SubElement(root, 'channel', attrs).text = base64.b64encode(zlib.compress(data, 9)).decode('ascii')
            used += len(data)
        except (ValueError, TypeError, OverflowError, MemoryError, AttributeError) as error:
            ET.SubElement(root, 'skipped', {'name': str(name), 'reason': str(error)})
    root.set('count', str(len(root.findall('channel'))))
    for quantity in EXPECTED_UNAVAILABLE:
        ET.SubElement(root, 'unavailable', {'quantity': quantity, 'reason': 'not_mapped_from_verified_CDF_inventory'})
    if not root.findall('channel'):
        root.set('status', 'no_matching_runtime_channels')
    return root


def append_runtime_diagnostics_xml(parameters, nc):
    if hasattr(parameters, 'value'):
        parameters = parameters.value
    if isinstance(parameters, bytes):
        parameters = parameters.decode('utf8')
    root = ET.fromstring(parameters) if parameters else ET.Element('transp2imas_parameters')
    for prev in root.findall('pt_solver_runtime_diagnostics'):
        root.remove(prev)
    root.append(diagnostic_xml(nc))
    return ET.tostring(root, encoding='unicode')


def decode_channel(channel):
    """Reconstruct a channel exactly as stored, verifying the SHA-256 digest."""
    data = zlib.decompress(base64.b64decode(channel.text))
    if hashlib.sha256(data).hexdigest() != channel.attrib['sha256']:
        raise ValueError('CDF diagnostic checksum mismatch')
    shape = tuple(int(v) for v in channel.attrib['shape'].split(',') if v)
    return np.frombuffer(data, dtype='<f8').reshape(shape)
