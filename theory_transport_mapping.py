#!/usr/bin/python
"""Unified model-output contracts for TRANSP core_transport (MMM07, GLF23, NC).

MMM07 is a reduced anomalous-transport model.  Its Weiland, DRIBM/DRBM,
ETG and MTM contributions are model *outputs*, not plasma diagnostics.
CDF coefficients reside on TIME3/XB in cm^2/s (or cm/s), whereas the
reduced-model stability/intermediate quantities reside on TIME3/X.

GLF23 and neoclassical output retain their own source semantics and units.
TGLF is registered but deferred until suitable test fixtures establish mapping.
This module has no IMAS dependency and is used by source tests.
"""

import re
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
        'format_version': '2', 'model': model['name'],
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
        # Coordinate *values* are stored in native IMAS arrays, never
        # duplicated inside code.parameters.  grid_d is the XB coefficient
        # grid.  For component models with reduced-mode outputs on X, the
        # otherwise-unused grid_flux stores that distinct radial coordinate.
        ET.SubElement(supplements, 'coordinate_reference', {
            'source': 'TIME3', 'ids_path': 'core_transport/time',
            'units': 's', 'length': str(len(plan['time3'])),
        })
        for grid in ('X', 'XB'):
            if any(AUXILIARY_SPECS[name][1] == grid for name in remaining):
                a = plan[grid.lower()]
                suffix = 'grid_flux' if grid == 'X' else 'grid_d'
                ET.SubElement(supplements, 'coordinate_reference', {
                    'source': grid,
                    'ids_path': ('core_transport/model[identifier.name=' +
                                 model['name'] + ']/profiles_1d[*]/' + suffix +
                                 '/rho_tor_norm'),
                    'shape': ' '.join(map(str, a.shape)),
                    'units': '1',
                })
        for name in remaining:
            meaning, grid, unit, scale = AUXILIARY_SPECS[name]
            data = plan['fields'][name] * scale
            node = ET.SubElement(supplements, 'quantity', {
                'source': name,
                'meaning': meaning,
                'source_dimensions': 'TIME3,' + grid,
                'radial_coordinate': ('grid_flux/rho_tor_norm' if grid == 'X'
                                      else 'grid_d/rho_tor_norm'),
                'time_coordinate': 'core_transport/time',
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


# Generic theory-model registry. "TGLF_DR" is not a TGLF diffusivity:
# its CDF long_name describes a scalar instability location, not a typed
# transport coefficient. Never fabricate a TGLF model from that channel.
THEORY_FAMILIES = ('neoclassical', 'mmm07', 'glf23', 'tglf')
GLF23_MODEL_SPEC = {
    'name': 'glf23',
    'index': -20,
    'label': 'GLF23 anomalous transport',
    'description': 'GLF23 predicted electron/ion heat and toroidal momentum diffusivities',
    'coefficients': {
        'XKEGLF23': ('electrons.energy.d', 'diffusivity'),
        'XKIGLF23': ('total_ion_energy.d', 'diffusivity'),
        'ETPHIGLF': ('momentum_phi.d', 'diffusivity'),
    },
    # GLF23 ion *particle* diffusivity is aggregate, not species-specific;
    # keep it in model XML until a defensible species mapping is established.
    'auxiliary': ('DIFFIGLF',),
}
GLF23_AUXILIARY_SPECS = {
    'DIFFIGLF': ('glf23_aggregate_ion_particle_diffusivity', 'XB', 'CM**2/SEC', 1.0e-4),
}
# All coefficient-model families use the same IMAS core_transport writer.
COEFFICIENT_MODEL_SPECS = {'mmm07': MODEL_SPECS, 'glf23': (GLF23_MODEL_SPEC,)}
THEORY_AUXILIARY_SPECS = {**AUXILIARY_SPECS, **GLF23_AUXILIARY_SPECS}


def build_glf23_plan(nc):
    """Read GLF23 coefficient channels on TIME3/XB, rejecting placeholders.

    The actual bundled ITER case has nonzero, independently described GLF23
    coefficients. Other CDFs declare zero-filled GLF23 arrays as placeholders.
    """
    names = GLF23_MODEL_SPEC['coefficients'].keys() | GLF23_AUXILIARY_SPECS.keys()
    fields = {}
    for name in sorted(names & nc.variables.keys()):
        fields[name] = _array(nc, name, 'XB', 'diffusivity')
    time = np.asarray(nc.variables['TIME3'][:], dtype=np.float64)
    xb = np.asarray(nc.variables['XB'][:], dtype=np.float64)
    active = any(np.any(np.abs(v) > 0) for v in fields.values())
    return {'active': bool(active), 'fields': fields, 'time3': time,
            'xb': xb, 'closure': {}}


def theory_model_xml(plan, spec, family):
    """Serialize supplemental MODEL data only, referencing native IDS grids."""
    if family == 'mmm07':
        return mmm_model_xml(plan, spec)
    if family != 'glf23':
        raise ValueError(f'unsupported coefficient-model family: {family}')
    root = ET.Element('transp2imas_theory_output', {
        'family': 'glf23', 'model': spec['name'], 'classification': 'anomalous_transport_model_output',
        'time_coordinate': 'core_transport/time',
        'radial_coordinate': 'profiles_1d/grid_d/rho_tor_norm',
    })
    mapped = ET.SubElement(root, 'typed_coefficients')
    for source, (destination, kind) in spec['coefficients'].items():
        if source in plan['fields']:
            ET.SubElement(mapped, 'coefficient', {
                'source': source, 'destination': destination,
                'source_units': 'cm^2/s', 'target_units': 'm^2/s',
                'source_dimensions': 'TIME3,XB',
            })
    for source in spec['auxiliary']:
        if source not in plan['fields']:
            continue
        arr = plan['fields'][source] * 1.0e-4
        ET.SubElement(root, 'supplementary_quantity', {
            'source': source, 'semantics': 'aggregate_ion_particle_diffusivity',
            'units': 'm^2/s', 'source_dimensions': 'TIME3,XB',
            'shape': ' '.join(map(str, arr.shape)),
        }).text = _numbers(arr)
    return ET.tostring(root, encoding='unicode')


NEOCLASSICAL_CHANNELS = {
    'particles': 'GFLNC_',
    'heat': 'QFLNC_',
    'heat_nc_classical': 'QFLNCC_',
}
NEO_AGGREGATE_TOKENS = frozenset(('I', 'X', 'TOK'))
SPECIES_ALIASES = {'3': 'HE3', '4': 'HE4'}
NEO_UNITS = {
    'particles': ('N/CM3/SEC', 'N/CM3/S'),
    'heat': ('WATTS/CM3', 'W/CM3'),
    'heat_nc_classical': ('WATTS/CM3', 'W/CM3'),
}


def discover_neoclassical_channels(names):
    """Separate real neoclassical species from total-ion/impurity aggregates."""
    available = set(names)
    result = {}
    for kind, prefix in NEOCLASSICAL_CHANNELS.items():
        channels = {}
        for name in sorted(available):
            if name.startswith(prefix):
                token = name[len(prefix):].upper()
                if re.fullmatch(r'[A-Z0-9]+', token):
                    channels[SPECIES_ALIASES.get(token, token)] = name
        result[kind] = channels
    return result


def _neoclassical_unit(value):
    if isinstance(value, bytes):
        value = value.decode('ISO-8859-1', 'replace')
    return str(value).strip().upper().replace(' ', '').replace('**', '').replace('^', '')


def validate_profile(nc, name, kind, required_shape=None):
    """Verify NC divergence is a finite TIME3/X array with physical units."""
    if kind not in NEO_UNITS:
        raise ValueError(f'not a neoclassical divergence kind: {kind}')
    var = nc.variables[name]
    if tuple(var.dimensions) != ('TIME3', 'X'):
        raise ValueError(f'{name}: expected (TIME3, X), got {var.dimensions}')
    found = _neoclassical_unit(getattr(var, 'units', ''))
    if found not in NEO_UNITS[kind]:
        raise ValueError(f'{name}: unexpected {kind} units {found!r}')
    arr = np.asarray(var[:], dtype=np.float64)
    if arr.ndim != 2 or (required_shape is not None and arr.shape != tuple(required_shape)):
        raise ValueError(f'{name}: unexpected shape {arr.shape}, expected {required_shape}')
    if not np.isfinite(arr).all():
        raise ValueError(f'{name}: non-finite CDF data')
    return arr


def divergence_to_surface_flux(div, dvol_cm3, surf_cm2, kind):
    """Integrate div(F) * zone-volume / boundary-area; convert cm^-2 to m^-2."""
    if kind not in NEOCLASSICAL_CHANNELS:
        raise ValueError(f'unknown divergence quantity: {kind}')
    div, dvol, surf = [np.asarray(a, dtype=np.float64)
                       for a in (div, dvol_cm3, surf_cm2)]
    if div.ndim != 2 or div.shape != dvol.shape or div.shape != surf.shape:
        raise ValueError(f'incompatible flux integration shapes: {div.shape}, {dvol.shape}, {surf.shape}')
    integrated = np.cumsum(div * dvol, axis=1)
    return np.divide(integrated * 1.0e4, surf, out=np.zeros_like(integrated),
                     where=np.abs(surf) > 0)


def build_theory_plans(nc):
    """Unified, testable family inventory for PT_SOLVER model outputs.

    Neoclassical source quantities are available even in interpretive runs:
    classification as *PT_SOLVER used* requires explicit input provenance
    and must not be inferred merely from CDF channel existence.
    """
    return {
        'neoclassical': discover_neoclassical_channels(nc.variables),
        'mmm07': build_mmm_plan(nc),
        'glf23': build_glf23_plan(nc),
        'tglf': {'active': False, 'status': 'WAITING_FOR_VALIDATED_TGLF_CDF',
                 'note': 'Do not infer TGLF transport from TGLF_DR or zero placeholders'},
    }
