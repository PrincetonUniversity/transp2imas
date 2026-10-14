#!/usr/bin/python
"""Source-grounded MMM07 component mapping for TRANSP core_transport.

MMM07 is a reduced anomalous-transport model.  Its Weiland, DRIBM/DRBM,
ETG and MTM contributions are model *outputs*, not plasma diagnostics.
CDF coefficients reside on TIME3/XB in cm^2/s (or cm/s), whereas the
reduced-model stability/intermediate quantities reside on TIME3/X.

This module deliberately has no IMAS dependency and is used by source tests.
"""

import xml.etree.ElementTree as ET

import numpy as np


# The total model must not be summed a second time with its components.
# CDF descriptions call the Weiland family "WEILAND19", while the input
# namelist uses W20 option names; preserve the literal source labels.
MODEL_SPECS = (
    {
        'name': 'mmm07_total',
        'index': -10,
        'label': 'MMM07 combined anomalous transport',
        'description': 'Combined MMM07 reduced anomalous transport coefficients',
        'coefficients': {
            'XKDMMM07': ('electrons.particles.d', 'diffusivity'),
            'XKEMMM07': ('electrons.energy.d', 'diffusivity'),
            'XKIMMM07': ('total_ion_energy.d', 'diffusivity'),
            'XPTMMM07': ('momentum_phi.d', 'diffusivity'),
            'VPHIMMM07': ('momentum_phi.v', 'velocity'),
        },
        'auxiliary': ('XKZMMM07', 'XPPMMM07'),
    },
    {
        'name': 'mmm07_weiland',
        'index': -11,
        'label': 'MMM07 Weiland (CDF Weiland19) component',
        'description': 'Weiland component of MMM anomalous heat diffusivity',
        'coefficients': {
            'XKEW19': ('electrons.energy.d', 'diffusivity'),
            'XKIW19': ('total_ion_energy.d', 'diffusivity'),
        },
        'auxiliary': (),
    },
    {
        'name': 'mmm07_dribm',
        'index': -12,
        'label': 'MMM07 drift-resistive ballooning/interchange component (DRIBM/DRBM)',
        'description': 'DRIBM/DRBM component of MMM anomalous heat diffusivity',
        'coefficients': {
            'XKEDRBM': ('electrons.energy.d', 'diffusivity'),
            'XKIDRBM': ('total_ion_energy.d', 'diffusivity'),
        },
        'auxiliary': ('KYRSDRBM',),
    },
    {
        'name': 'mmm07_etg',
        'index': -13,
        'label': 'MMM07 ETG (electron-temperature-gradient) component',
        'description': 'ETG reduced-model contribution to MMM electron heat diffusivity',
        'coefficients': {'XKEETG': ('electrons.energy.d', 'diffusivity')},
        'auxiliary': ('GAMMAETG', 'OMEGAETG', 'KYRSETG'),
    },
    {
        'name': 'mmm07_mtm',
        'index': -14,
        'label': 'MMM07 reduced microtearing (MTM) component',
        'description': 'Reduced MTM contribution to MMM electron heat diffusivity',
        'coefficients': {'XKEMTM': ('electrons.energy.d', 'diffusivity')},
        'auxiliary': ('GAMMAMTM', 'OMEGAMTM', 'KYRSMTM', 'XDBMTM'),
    },
)

# Quantity type, native radial coordinate, expected TRANSP unit, SI multiplier.
AUXILIARY_SPECS = {
    'XKZMMM07': ('mmm_chi_z_aggregate_semantics_unresolved', 'XB', 'CM**2/SEC', 1.0e-4),
    'XPPMMM07': ('poloidal_momentum_diffusivity_aggregate', 'XB', 'CM**2/SEC', 1.0e-4),
    'KYRSDRBM': ('ky_rhos_at_max_growth', 'X', '-', 1.0),
    'GAMMAETG': ('reduced_model_growth_rate', 'X', '1/SEC', 1.0),
    'OMEGAETG': ('reduced_model_mode_frequency', 'X', 'RAD/SEC', 1.0),
    'KYRSETG': ('ky_rhos_at_max_growth', 'X', '-', 1.0),
    'GAMMAMTM': ('reduced_model_growth_rate', 'X', '1/SEC', 1.0),
    'OMEGAMTM': ('reduced_model_mode_frequency', 'X', 'RAD/SEC', 1.0),
    'KYRSMTM': ('ky_rhos_at_max_growth', 'X', '-', 1.0),
    'XDBMTM': ('reduced_model_delta_b_over_b', 'X', '-', 1.0),
}

TYPED_FIELDS = {
    name: (destination, kind)
    for model in MODEL_SPECS
    for name, (destination, kind) in model['coefficients'].items()
}
ALL_FIELDS = tuple(TYPED_FIELDS) + tuple(AUXILIARY_SPECS)
MODEL_BY_NAME = {model['name']: model for model in MODEL_SPECS}


# Only accepted units are the ones established by the bundled CDF metadata.
UNIT_SETS = {
    'diffusivity': frozenset(('CM**2/SEC', 'CM^2/SEC', 'CM2/SEC')),
    'velocity': frozenset(('CM/SEC',)),
    '1/SEC': frozenset(('1/SEC', '/SEC')),
    'RAD/SEC': frozenset(('RAD/SEC',)),
    '-': frozenset(('-', '1', '')),
}


def _text(value):
    if isinstance(value, (bytes, np.bytes_)):
        return bytes(value).decode('latin1', errors='replace').strip()
    return str(value).strip()


def _array(nc, name, expected_grid, kind):
    v = nc.variables[name]
    expected = ('TIME3', expected_grid)
    if tuple(v.dimensions) != expected:
        raise ValueError(f'MMM field {name}: expected {expected}, got {v.dimensions}')
    unit = _text(getattr(v, 'units', '')).upper().replace(' ', '')
    allowed = UNIT_SETS[kind]
    if unit not in allowed:
        raise ValueError(f'MMM field {name}: unexpected units {unit!r}, expected {sorted(allowed)}')
    data = np.asarray(v[:], dtype=np.float64)
    grid_shape = np.asarray(nc.variables[expected_grid][:]).shape
    if data.shape != grid_shape:
        raise ValueError(f'MMM field {name}: shape {data.shape} != {expected_grid} {grid_shape}')
    if not np.isfinite(data).all():
        raise ValueError(f'MMM field {name}: non-finite source data')
    return data


def check_mmm_closure(fields):
    """Check the *source* MMM component decomposition (cm^2/s).

    The closure is assessed only when the total and all component channels
    needed for the equation exist.  Normalization is by max |total|, so
    roundoff near zero does not spuriously reject a physically valid profile.
    """
    comparisons = {
        'electron_heat': ('XKEMMM07', ('XKEW19', 'XKEDRBM', 'XKEETG', 'XKEMTM')),
        'ion_heat': ('XKIMMM07', ('XKIW19', 'XKIDRBM')),
    }
    result = {}
    for label, (total, components) in comparisons.items():
        missing = [n for n in (total, *components) if n not in fields]
        if missing:
            result[label] = {'status': 'INCOMPLETE', 'missing': missing}
            continue
        actual = fields[total]
        summation = np.sum(np.stack([fields[name] for name in components]), axis=0)
        difference = actual - summation
        scale = max(float(np.max(np.abs(actual))), 1.0)
        result[label] = {
            'status': 'PASS' if np.max(np.abs(difference)) <= scale * 2.0e-6 else 'FAIL',
            'max_abs_difference_cm2_s': float(np.max(np.abs(difference))),
            'max_normalized_difference': float(np.max(np.abs(difference)) / scale),
            'sources': [total, *components],
        }
    return result


def build_mmm_plan(nc):
    """Read/validate every MMM family source present, without fabricating MMM use.

    Legacy TRANSP CDFs may declare all MMM variables as zero-filled output
    placeholders, even when another anomalous model was run.  No MMM models
    are emitted unless at least one supported MMM source contains nonzero data.
    """
    variables = nc.variables
    time = np.asarray(variables['TIME3'][:], dtype=np.float64)
    x = np.asarray(variables['X'][:], dtype=np.float64)
    xb = np.asarray(variables['XB'][:], dtype=np.float64)
    if time.ndim != 1 or not np.isfinite(time).all():
        raise ValueError('MMM TIME3 is not a finite 1D time base')
    for name, grid in (('X', x), ('XB', xb)):
        if grid.ndim != 2 or grid.shape[0] != len(time) or not np.isfinite(grid).all():
            raise ValueError(f'MMM {name} is not aligned with TIME3')

    fields = {}
    for name in ALL_FIELDS:
        if name not in variables:
            continue
        if name in TYPED_FIELDS:
            kind = TYPED_FIELDS[name][1]
            fields[name] = _array(nc, name, 'XB', kind)
        else:
            _semantic, radial, unit, scale = AUXILIARY_SPECS[name]
            fields[name] = _array(nc, name, radial, 'diffusivity' if scale != 1.0 else unit)
    active = any(np.any(data != 0) for data in fields.values())
    closure = check_mmm_closure(fields) if active else {}
    return {'active': bool(active), 'time3': time, 'x': x, 'xb': xb,
            'fields': fields, 'closure': closure}


def _numbers(values):
    return ' '.join(format(float(x), '.9g') for x in np.asarray(values).ravel())


def mmm_model_xml(plan, model):
    """Preserve non-typed reduced MMM outputs in explicitly labelled XML.

    IMAS DD4 core_transport has typed transport coefficients but no native
    fields for mode growth/frequency/ky*rhos/delta-B, nor an aggregate Z heat
    or poloidal-momentum diffusivity.  This XML is attached to *the model*
    in code.parameters; it is not an experimental measurement.
    """
    root = ET.Element('transp2imas_mmm_output', {
        'format_version': '1', 'model': model['name'],
        'classification': 'reduced_anomalous_transport_model_output',
        'parent_model': 'mmm07_total' if model['name'] != 'mmm07_total' else '',
    })
    mapping = ET.SubElement(root, 'typed_coefficients')
    for source, (destination, kind) in model['coefficients'].items():
        if source not in plan['fields']:
            continue
        ET.SubElement(mapping, 'coefficient', {
            'source': source, 'destination': destination,
            'source_units': 'cm^2/s' if kind == 'diffusivity' else 'cm/s',
            'target_units': 'm^2/s' if kind == 'diffusivity' else 'm/s',
            'source_grid': 'TIME3,XB',
        })
    remaining = [name for name in model['auxiliary'] if name in plan['fields']]
    if remaining:
        supplements = ET.SubElement(root, 'supplementary_model_outputs', {
            'note': 'These are source MMM outputs, not experimentally measured diagnostics',
        })
        times = ET.SubElement(supplements, 'coordinate', {
            'name': 'TIME3', 'units': 's', 'shape': str(len(plan['time3'])),
        })
        times.text = _numbers(plan['time3'])
        for grid in ('X', 'XB'):
            if any(AUXILIARY_SPECS[name][1] == grid for name in remaining):
                a = plan[grid.lower()]
                node = ET.SubElement(supplements, 'coordinate', {
                    'name': grid, 'units': 'rho_tor_norm',
                    'shape': ' '.join(map(str, a.shape)),
                })
                node.text = _numbers(a)
        for name in remaining:
            meaning, grid, unit, scale = AUXILIARY_SPECS[name]
            data = plan['fields'][name] * scale
            node = ET.SubElement(supplements, 'quantity', {
                'source': name,
                'meaning': meaning,
                'grid': 'TIME3,' + grid,
                'units': ('m^2/s' if scale != 1 else ('1' if unit == '-' else unit.lower())),
                'shape': ' '.join(map(str, data.shape)),
            })
            node.text = _numbers(data)
    if model['name'] == 'mmm07_total' and plan['closure']:
        root_closure = ET.SubElement(root, 'source_component_closure')
        for label, result in plan['closure'].items():
            attrs = {'quantity': label, 'status': result['status']}
            if 'max_normalized_difference' in result:
                attrs['max_normalized_difference'] = format(result['max_normalized_difference'], '.9g')
            ET.SubElement(root_closure, 'check', attrs)
    return ET.tostring(root, encoding='unicode')


def mmm_summary(plan):
    """Return JSON-safe per-source audit information."""
    models = []
    if plan['active']:
        models = [model['name'] for model in MODEL_SPECS
                  if any(name in plan['fields'] for name in
                         (*model['coefficients'], *model['auxiliary']))]
    return {'active': plan['active'], 'models': models,
            'source_channels': list(plan['fields']), 'closure': plan['closure']}
