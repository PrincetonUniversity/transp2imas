#!/usr/bin/python
"""TRANSP neutral-state profiles for core_profiles, separate from transport theory.

Wall, volume, and recombination neutral populations retain their individual
origin labels. This module does not create any core_transport models.
"""

import re
import xml.etree.ElementTree as ET
import numpy as np

SPECIES_ALIASES = {'3': 'HE3', '4': 'HE4'}
ORIGIN_ORDER = ('wall', 'volume', 'recombination')


def canonical_species(token):
    token = str(token).strip().upper()
    return SPECIES_ALIASES.get(token, token)


def discover_neutral_channels(names):
    """Enumerate wall, volume and recombination neutrals by actual CDF name.

    Distinguish these physically different populations; their densities must
    not be silently added when the underlying TRANSP processes may overlap.
    """
    available = set(names)
    profiles = {}
    for name in sorted(available):
        match = re.fullmatch(r'DN0([WV])([A-Z0-9]+)', name)
        if match:
            origin = 'wall' if match.group(1) == 'W' else 'volume'
            raw = match.group(2)
            prefix = match.group(1)
            profiles[(origin, canonical_species(raw))] = {
                'origin': origin, 'species': canonical_species(raw),
                'density': name,
                'temperature': f'T0{prefix}{raw}' if f'T0{prefix}{raw}' in available else None,
                'rotation': f'OM0{prefix}{raw}' if f'OM0{prefix}{raw}' in available else None,
            }
        match = re.fullmatch(r'N0V0_([A-Z0-9]+)', name)
        if match:
            raw = match.group(1)
            profiles[('recombination', canonical_species(raw))] = {
                'origin': 'recombination', 'species': canonical_species(raw),
                'density': name,
                'temperature': f'T0V0_{raw}' if f'T0V0_{raw}' in available else None,
                'rotation': f'OM0V0_{raw}' if f'OM0V0_{raw}' in available else None,
            }
    return [profiles[key] for key in sorted(
        profiles, key=lambda k: (ORIGIN_ORDER.index(k[0]), k[1]))]


NEUTRAL_UNITS = {
    'density': ('N/CM3', '1/CM3', 'CM-3'),
    'temperature': ('EV',),
    'rotation': ('RAD/SEC', 'RAD/S'),
}

def _unit(value):
    if isinstance(value, bytes):
        value = value.decode('ISO-8859-1', 'replace')
    return str(value).strip().upper().replace(' ', '').replace('**', '').replace('^', '')


def validate_profile(nc, name, kind, required_shape=None):
    """Reject dimension/unit mismatches instead of silently mis-scaling data."""
    var = nc.variables[name]
    if tuple(var.dimensions) != ('TIME3', 'X'):
        raise ValueError(f'{name}: expected (TIME3, X), got {var.dimensions}')
    found = _unit(getattr(var, 'units', ''))
    if found not in NEUTRAL_UNITS[kind]:
        raise ValueError(f'{name}: unexpected {kind} units {found!r}')
    arr = np.asarray(var[:], dtype=np.float64)
    if arr.ndim != 2 or (required_shape is not None and tuple(arr.shape) != tuple(required_shape)):
        raise ValueError(f'{name}: unexpected shape {arr.shape}, expected {required_shape}')
    if not np.isfinite(arr).all():
        raise ValueError(f'{name}: non-finite CDF data')
    return arr


def _xml_string(parameters_xml):
    """Read an IMAS IDSString0D (or a plain string) as XML text.

    IMAS-Python exposes a scalar string node through its ``value`` property.
    Passing the node itself to ElementTree, or assigning it back to another
    IDS string node, fails even though the node's value is valid XML.
    """
    if not isinstance(parameters_xml, (str, bytes, np.str_, np.bytes_)):
        parameters_xml = getattr(parameters_xml, 'value', parameters_xml)
    if isinstance(parameters_xml, np.ndarray) and parameters_xml.ndim == 0:
        parameters_xml = parameters_xml.item()
    if parameters_xml is None:
        return ''
    if isinstance(parameters_xml, (bytes, np.bytes_)):
        return parameters_xml.decode('utf-8')
    if isinstance(parameters_xml, (str, np.str_)):
        return str(parameters_xml)
    raise TypeError(
        'code.parameters must contain an XML string, got '
        f'{type(parameters_xml).__name__}')


def attach_rotation_xml(parameters_xml, neutrals, times):
    """Preserve DD4 neutral omega without pretending it is an ion rotation.

    No typed core_profiles neutral toroidal velocity/frequency node exists in
    DD 4.1.1.  The signed omega profiles are retained as explicit, indexed
    code-specific XML.  Each row follows its matching TIME3/X profile grid.

    ``parameters_xml`` may be the IMAS ``IDSString0D`` node returned by
    ``cp.code.parameters``. Always return a native Python string, including
    when there are no neutral rotations.
    """
    xml = _xml_string(parameters_xml)
    records = [n for n in neutrals if isinstance(n.omega, np.ndarray) and n.omega.size]
    if not records:
        return xml
    root = ET.fromstring(xml) if xml.strip() else ET.Element('transp2imas_parameters')
    group = ET.SubElement(root, 'neutral_toroidal_rotation', {
        'schema': 'transp2imas-neutral-omega-1',
        'units': 'rad/s', 'coordinate_system': 'IMAS positive toroidal phi',
        'time_axis': 'TIME3', 'radial_axis': 'X',
        'radial_grid_ids_path': 'core_profiles/profiles_1d/grid/rho_tor_norm',
        'note': 'Supplementary code.parameters data; DD4 core_profiles neutral has no rotation field',
    })
    t = ET.SubElement(group, 'time', {'units': 's'})
    t.text = ' '.join(format(float(v), '.17g') for v in times)
    for n in records:
        arr = np.asarray(n.omega, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[0] != len(times):
            raise ValueError(f'{n.label}: invalid rotation/time dimensions {arr.shape}')
        item = ET.SubElement(group, 'profile', {
            'neutral_name': n.label, 'source_variable': n.rotation_source,
            'rows': str(arr.shape[0]), 'cols': str(arr.shape[1]),
            'order': 'time-major-radial-inner',
            'source_sign_factor': format(float(n.rotation_sign), '.1f'),
        })
        item.text = '\n' + '\n'.join(' '.join(format(float(v), '.10g') for v in row)
                                         for row in arr) + '\n'
    return ET.tostring(root, encoding='unicode')
