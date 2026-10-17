#!/usr/bin/python
"""PT_SOLVER radial-region input provenance for an assembled IMAS transport solution.

A TRANSP PT_SOLVER run has ONE physical solution.  Its axial, confinement and
edge settings select theory models on portions of the normalized toroidal-flux
radius; they are neither separate core_transport occurrences nor alternative
solutions.  Distinct per-region *output profiles* are not available in the
bundled CDF cases, so this module only represents documented input settings.

The source of truth is the original TR.DAT plus a specifically referenced
PT_TEMPLATE.  Adjacent unreferenced *.nml are retained as CANDIDATES only.
"""

import math
import re
import xml.etree.ElementTree as ET

from input_provenance import assignments, strip_fortran_comment

REGIONS = ('axial', 'confinement', 'edge')
BOUNDARIES = {
    'temperature': 'XIBOUND',
    'density': 'XNBOUND',
    'toroidal_momentum': 'XPHIBOUND',
}
DEFAULTS = {
    'XIMIN_CONF': 0.0,
    'XIMAX_CONF': 1.0,
    'XIBOUND': 1.0,
    'XNBOUND': 1.0,
    'XPHIBOUND': 1.0,
}
PREDICT_FLAGS = (
    'LPREDICTIVE_MODE', 'LPREDICT_TE', 'LPREDICT_TI',
    'LPREDICT_PPHI', 'LPREDICT_NE', 'LPREDICT_NMAIN', 'LPREDICT_NIMP',
)
XANOM_CHANNELS = ('Te', 'Ti', 'Pphi', 'Ni', 'Nz', 'Ne')
# The PT_SOLVER module documentation explicitly maps these numbers.
EXB_VERSIONS = {1: 'NCEXB', 2: 'DMEXB', 3: 'TGYROEXB', 4: 'TRANSPEXB'}


def _number(raw):
    """Strict scalar Fortran number; return None for absent/expression values."""
    value = str(raw).strip().rstrip(',').strip()
    if not re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?', value):
        return None
    result = float(value.replace('D', 'E').replace('d', 'e'))
    return result if math.isfinite(result) else None


def _boolean(raw):
    value = str(raw).strip().rstrip(',').strip().upper()
    if value in ('.TRUE.', '.T.', 'T', 'TRUE'):
        return True
    if value in ('.FALSE.', '.F.', 'F', 'FALSE'):
        return False
    return None


def _xanom(raw):
    """Parse Fortran 6-entry multiplier list, including n*value notation."""
    out = []
    for token in str(raw).split(','):
        token = token.strip()
        match = re.fullmatch(r'(\d+)\s*\*\s*(.+)', token)
        count = int(match.group(1)) if match else 1
        value = _number(match.group(2) if match else token)
        if value is None or count > 6:
            return None
        out.extend([value] * count)
    return out if len(out) == 6 else None


def _trdat_states(text):
    """Snapshot updatable TR.DAT boundaries and equation flags at each update.

    TRANSP uses ~update_time to start a new input state; do not collapse all
    assignments into their final values. Defaults are from PT_SOLVER docs and
    are used ONLY when the actual TR.DAT was supplied.
    """
    controls = dict(DEFAULTS)
    origins = {name: 'TRANSP_documented_default' for name in DEFAULTS}
    snapshots = []
    at = None
    events = []
    for line in str(text).splitlines():
        source = strip_fortran_comment(line).strip()
        update = re.match(r'^~\s*update_time\s*=\s*(.+)$', source, re.I)
        if update:
            time = _number(update.group(1))
            if time is not None:
                events.append((at, dict(controls), dict(origins)))
                at = time
            continue
        for _line, key, raw in assignments(source):
            name = key.upper()
            if name in DEFAULTS:
                number = _number(raw)
                if number is not None:
                    controls[name] = number
                    origins[name] = 'TR.DAT'
            elif name in PREDICT_FLAGS:
                value = _number(raw)
                if value is not None and value.is_integer():
                    controls[name] = int(value)
                    origins[name] = 'TR.DAT'
    events.append((at, dict(controls), dict(origins)))
    # An update may contain unrelated RF/heating changes. Every phase is
    # preserved, even when no PT_SOLVER setting changes in that phase.
    for time, vals, src in events:
        snapshots.append({'start_time': time, 'values': vals, 'origins': src})
    return snapshots


def _pt_controls(text):
    """Read independent &pt_controls blocks, preserving update-time semantics.

    A block lists *explicit assignments*, not a complete model-state snapshot.
    Never infer the state of omitted models from prior blocks without checking
    the PT_SOLVER update semantics for that TRANSP version.
    """
    blocks = []
    current = None
    last_slot = None
    for line in str(text).splitlines():
        source = strip_fortran_comment(line).strip()
        if re.match(r'^&pt_controls\b', source, re.I):
            if current is not None:
                blocks.append(current)
            current = {'start_time': None, 'options': {}, 'regions': {}}
            last_slot = None
            continue
        if current is None:
            continue
        if source == '/' or re.match(r'^&end\b', source, re.I):
            blocks.append(current)
            current = None
            last_slot = None
            continue
        found_assignment = False
        for _line, key, value in assignments(source):
            found_assignment = True
            key = key.lower()
            if key == 'pt_update_time':
                current['start_time'] = _number(value)
                last_slot = None
                continue
            m = re.fullmatch(r'pt_(axial|confinement|edge)%([a-z][a-z0-9_]*)%([a-z][a-z0-9_]*)', key)
            if m:
                region, model, option = m.groups()
                settings = current['regions'].setdefault(region, {}).setdefault(model, {})
                settings[option] = value
                last_slot = (settings, option)
            else:
                current['options'][key] = value
                last_slot = (current['options'], key)
        if not found_assignment and last_slot and source:
            # Fortran namelist continuations can split a six-channel xanom
            # list across physical lines. Preserve/parse its complete value.
            owner, name = last_slot
            owner[name] += ' ' + source.lstrip('&').strip()
    if current is not None:
        blocks.append(current)
    return blocks


def build_pt_region_plan(inventory):
    """Build documented PT_SOLVER regional input plan from archived run inputs.

    Return None if the input contains no evidence of PT_SOLVER configuration.
    A standalone model CDF channel is not evidence of a particular *region*.
    """
    records = inventory.records
    primary = next((r for r in records if r.role == 'TRANSP_INPUT'), None)
    main_ok = primary is not None and primary.status == 'embedded'
    trdat_states = _trdat_states(primary.text) if main_ok else []
    referenced = [r for r in records if r.model == 'PT_SOLVER'
                  and r.key.upper() == 'PT_TEMPLATE' and r.role == 'TEMPLATE'
                  and r.status == 'embedded']
    # Orphan inputs may be preserved from a run directory, but must never be
    # presented as the confirmed model used for the CDF simulation.
    candidates = [r for r in records if r.model == 'PT_SOLVER'
                  and r.role == 'ADJACENT_UNREFERENCED' and r.status == 'embedded']
    mode = (trdat_states[-1]['values'].get('LPREDICTIVE_MODE')
            if trdat_states else None)
    has_declared_template = any(r.key.upper() == 'PT_TEMPLATE'
                                and r.role == 'TEMPLATE' for r in records)
    if mode != 3 and not has_declared_template and not candidates:
        return None
    # Only a PT_TEMPLATE *actually linked from TR.DAT* is authoritative.
    templates = referenced if referenced else candidates
    return {
        'trdat_present': main_ok,
        'trdat_file': str(primary.path) if primary else '',
        'predictive_mode': mode,
        'boundary_states': trdat_states,
        'model_source': ('explicit_PT_TEMPLATE' if referenced else
                         'adjacent_unreferenced_candidate' if candidates else
                         'referenced_PT_TEMPLATE_missing' if has_declared_template else
                         'no_PT_TEMPLATE'),
        'trusted_model_assignment': bool(referenced),
        'templates': [
            {'path': str(r.path), 'sha256': r.sha256,
             'source_role': r.role, 'blocks': _pt_controls(r.text)}
            for r in templates
        ],
    }


def _fmt(value):
    return format(float(value), '.12g')


def _range_limits(region, values, quantity):
    outer = values[BOUNDARIES[quantity]]
    low_conf, high_conf = values['XIMIN_CONF'], values['XIMAX_CONF']
    if region == 'axial':
        low, high = 0., min(low_conf, outer)
    elif region == 'confinement':
        low, high = low_conf, min(high_conf, outer)
    else:
        low, high = high_conf, outer
    return low, high


def region_plan_xml(plan):
    """Serialize *input settings* into core_transport.code.parameters.

    CDF TIME3/XB are already in typed grid/time fields; no radial array is
    duplicated in this document. Boundary scalars are not grid coordinates:
    they describe regions of applicability for separate predictive equations.
    """
    root = ET.Element('pt_solver_regions', {
        'format_version': '1',
        'classification': 'TRANSP_PT_SOLVER_input_configuration',
        'solution_policy': 'one_core_transport_occurrence_one_assembled_solution',
        'region_representation': 'input_metadata_not_alternative_solution',
        'radial_coordinate': 'xi=sqrt(Phi_tor/Phi_tor_boundary)',
        'native_grid_reference': 'core_transport/model[*]/profiles_1d[*]/grid_d/rho_tor_norm',
        'model_source': plan['model_source'],
        'model_assignment_verified_by_reference': str(plan['trusted_model_assignment']).lower(),
    })
    if plan['trdat_file']:
        ET.SubElement(root, 'trdat', {
            'path': plan['trdat_file'],
            'status': 'embedded' if plan['trdat_present'] else 'missing',
        })
    if plan['predictive_mode'] is not None:
        ET.SubElement(root, 'predictive_mode').text = str(plan['predictive_mode'])
    timeline = ET.SubElement(root, 'boundary_and_equation_timeline')
    if not plan['boundary_states']:
        timeline.set('status', 'unavailable_without_TR.DAT')
    for phase in plan['boundary_states']:
        attrs = {'effective_from': ('initial' if phase['start_time'] is None
                                   else _fmt(phase['start_time']))}
        entry = ET.SubElement(timeline, 'state', attrs)
        vals, origins = phase['values'], phase['origins']
        for name in DEFAULTS:
            ET.SubElement(entry, 'boundary', {
                'name': name, 'xi': _fmt(vals[name]), 'source': origins[name],
            })
        for name in PREDICT_FLAGS:
            if name in vals:
                ET.SubElement(entry, 'equation_option', {
                    'name': name, 'value': str(vals[name]), 'source': origins[name],
                })
        domains = ET.SubElement(entry, 'radial_domains')
        for quantity in BOUNDARIES:
            for region in REGIONS:
                low, high = _range_limits(region, vals, quantity)
                ET.SubElement(domains, 'region', {
                    'name': 'pt_' + region, 'quantity': quantity,
                    'xi_min': _fmt(low), 'xi_max': _fmt(high),
                    'has_extent': str(high > low).lower(),
                    'outer_boundary_control': BOUNDARIES[quantity],
                })
    configs = ET.SubElement(root, 'model_configuration_timeline')
    if not plan['templates']:
        configs.set('status', 'no_readable_PT_TEMPLATE')
    for source in plan['templates']:
        inp = ET.SubElement(configs, 'input_file', {
            'path': source['path'], 'sha256': source['sha256'],
            'association': ('explicit_TRDAT_PT_TEMPLATE' if plan['trusted_model_assignment']
                            else 'adjacent_unreferenced_candidate'),
        })
        for i, block in enumerate(source['blocks']):
            at = block['start_time']
            phase = ET.SubElement(inp, 'configuration', {
                'index': str(i),
                'effective_from': ('initial' if at is None or at < 0 else _fmt(at)),
                'semantics': 'explicit_settings_not_inferred_model_state',
            })
            for region in REGIONS:
                for model, options in sorted(block['regions'].get(region, {}).items()):
                    node = ET.SubElement(phase, 'model', {
                        'region': 'pt_' + region, 'name': model.upper(),
                    })
                    for key, raw in sorted(options.items()):
                        attrib = {'key': key, 'raw': raw}
                        if key == 'active':
                            state = _boolean(raw)
                            if state is not None:
                                attrib['value'] = str(state).lower()
                        elif key == 'xanom':
                            vals = _xanom(raw)
                            if vals is not None:
                                attrib['channels'] = ','.join(XANOM_CHANNELS)
                                attrib['values'] = ','.join(_fmt(v) for v in vals)
                        elif key == 'version' and model == 'exb':
                            value = _number(raw)
                            if value is not None and value.is_integer():
                                attrib['value'] = str(int(value))
                                if int(value) in EXB_VERSIONS:
                                    attrib['scheme'] = EXB_VERSIONS[int(value)]
                        ET.SubElement(node, 'option', attrib)
            for name, raw in sorted(block['options'].items()):
                ET.SubElement(phase, 'solver_option', {'name': name, 'raw': raw})
    return root


def _text(value):
    """Support plain strings and real imas IDSString0D wrappers."""
    if hasattr(value, 'value'):
        value = value.value
    if value is None:
        return ''
    if isinstance(value, bytes):
        return value.decode('utf-8')
    if not isinstance(value, str):
        raise TypeError(f'expected IMAS XML string, got {type(value).__name__}')
    return value


def _xml_root(parameters):
    text = _text(parameters)
    return ET.fromstring(text) if text.strip() else ET.Element('transp2imas_parameters')


def append_pt_regions_xml(parameters, plan):
    """Attach PT_SOLVER region policy without overwriting original input XML."""
    root = _xml_root(parameters)
    if plan is None:
        return ET.tostring(root, encoding='unicode')
    for prior in root.findall('pt_solver_regions'):
        root.remove(prior)
    root.append(region_plan_xml(plan))
    return ET.tostring(root, encoding='unicode')


def annotate_theory_model_xml(parameters, plan, model_name):
    """Point an actual model's code.parameters to configured PT regions.

    No copy of a model coefficient array is split by region.  Input evidence
    is added ONLY for explicitly referenced PT_TEMPLATE blocks.
    """
    root = _xml_root(parameters)
    if plan is None or not plan['trusted_model_assignment']:
        return ET.tostring(root, encoding='unicode')
    name = model_name.lower()
    possible = ({'mmm'} if name.startswith('mmm07') else
                {'glf23'} if name == 'glf23' else
                {'neoch', 'neogk', 'nclass', 'gtcneo'}
                if name == 'neoclassical' else set())
    if not possible:
        return ET.tostring(root, encoding='unicode')
    applications = []
    for src in plan['templates']:
        for i, block in enumerate(src['blocks']):
            for region, models in block['regions'].items():
                for selected_model in sorted(possible & models.keys()):
                    active = _boolean(models[selected_model].get('active', ''))
                    if active is not None:
                        applications.append((region, block['start_time'], i,
                                             active, selected_model))
    if not applications:
        return ET.tostring(root, encoding='unicode')
    application = ET.SubElement(root, 'pt_solver_region_application', {
        'source': 'explicit_TRDAT_PT_TEMPLATE',
        'region_definition': 'core_transport/code/parameters/pt_solver_regions',
        'note': 'model output on native radial grid; no per-region output splitting implied',
    })
    for region, start, index, active, selected_model in applications:
        ET.SubElement(application, 'selection', {
            'region': 'pt_' + region,
            'theory_model': selected_model.upper(),
            'effective_from': ('initial' if start is None or start < 0 else _fmt(start)),
            'configuration_index': str(index),
            'active': str(active).lower(),
        })
    return ET.tostring(root, encoding='unicode')


def append_pt_regions_reference(parameters, plan):
    """Link assembled transport_solver fluxes to run-level region inputs."""
    root = _xml_root(parameters)
    if plan is not None:
        ET.SubElement(root, 'pt_solver_region_configuration_reference', {
            'path': 'core_transport/code/parameters/pt_solver_regions',
            'source': plan['model_source'],
            'model_assignment_verified_by_reference': str(plan['trusted_model_assignment']).lower(),
            'note': 'assembled predictive fluxes; no per-region coefficient arrays available',
        })
    return ET.tostring(root, encoding='unicode')
