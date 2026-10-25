#!/usr/bin/python
"""PT_SOLVER radial-region input provenance, and the region occurrences built from it.

A TRANSP PT_SOLVER run has ONE physical solution: ``core_transport`` writes it as occurrence 0, the
assembled result.  The axial, confinement and edge settings select theory models on portions of the
normalized toroidal-flux radius over time; they are not alternative solutions.  ``core_transport``
adds one further *view* occurrence per region that selected a model (``occurrence_region_models``),
carrying the selected models with their time intervals and radial validity as metadata and filling
model slots only from native CDF output channels that exist.  Distinct per-region *output profiles*
are not available in the bundled CDF cases, so nothing is relabelled or invented: a region occurrence
holds no profile that the CDF did not produce for it.

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
             'source_role': r.role, 'blocks': _pt_controls(r.text), 'text': r.text}
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


def numerical_settings_xml(plan):
    """Preserve numeric PT_SOLVER controls as structured, timed input metadata.

    All &pt_controls non-model assignments are kept (not a fragile whitelist),
    including NEWTON_ITERATES, THETA_IMPLICIT and pt_residual%*.
    These are INPUT SETTINGS, not measured Newton convergence histories.
    """
    root = ET.Element('pt_solver_numerics', {
        'scope': 'combined_solver_shared_by_all_regions',
        'kind': 'input_settings',
        'runtime_diagnostics': 'not_available_from_PT_TEMPLATE',
    })
    for source in plan['templates']:
        inp = ET.SubElement(root, 'input_file', {
            'path': source['path'], 'sha256': source['sha256'],
            'trusted': str(plan['trusted_model_assignment']).lower(),
        })
        for index, block in enumerate(source['blocks']):
            at = block['start_time']
            phase = ET.SubElement(inp, 'settings', {
                'index': str(index),
                'effective_from': 'initial' if at is None or at < 0 else _fmt(at),
                'semantics': 'explicit_namelist_assignments',
            })
            for name, raw in sorted(block['options'].items()):
                attrs = {'name': name, 'raw': raw}
                parsed = _number(raw)
                flag = _boolean(raw)
                if parsed is not None:
                    attrs['parsed_scalar'] = _fmt(parsed)
                elif flag is not None:
                    attrs['parsed_boolean'] = str(flag).lower()
                ET.SubElement(phase, 'control', attrs)
            # Model-specific options, including E x B shear versions, are
            # numerical controls for a particular region and model.
            for region in REGIONS:
                for model, options in sorted(block['regions'].get(region, {}).items()):
                    for name, raw in sorted(options.items()):
                        ET.SubElement(phase, 'regional_model_control', {
                            'region': region, 'model': model.upper(),
                            'name': name, 'raw': raw,
                        })
    if not plan['templates']:
        root.set('status', 'no_PT_TEMPLATE')
    return root


def region_validity_xml(plan, region):
    """Exact time-varying radial validity of each un-clipped regional profile.

    Each predicted channel can have a different outer boundary.  Do not
    simplify the model's validity to one fixed [xi_min,xi_max] interval.
    """
    root = ET.Element('regional_validity', {
        'region': 'pt_' + region,
        'radial_coordinate': 'xi=sqrt(Phi_tor/Phi_tor_boundary)',
        'full_profile_storage': 'true',
        'outside_validity': 'stored_values_not_selected_by_this_region',
        'time_semantics': 'piecewise_constant_from_effective_from',
    })
    states = plan['boundary_states']
    for index, state in enumerate(states):
        at = state['start_time']
        following = states[index + 1]['start_time'] if index + 1 < len(states) else None
        phase = ET.SubElement(root, 'state', {
            'effective_from': 'initial' if at is None else _fmt(at),
            'effective_until': 'end' if following is None else _fmt(following),
        })
        for quantity in BOUNDARIES:
            low, high = _range_limits(region, state['values'], quantity)
            ET.SubElement(phase, 'domain', {
                'quantity': quantity, 'xi_min': _fmt(low),
                'xi_max': _fmt(high), 'has_extent': str(high > low).lower(),
                'outer_boundary_control': BOUNDARIES[quantity],
            })
    if not plan['boundary_states']:
        root.set('status', 'unknown_without_TR_DAT')
    return root


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
    root.append(numerical_settings_xml(plan))
    configs = ET.SubElement(root, 'model_configuration_timeline')
    if not plan['templates']:
        configs.set('status', 'no_readable_PT_TEMPLATE')
    for source in plan['templates']:
        inp = ET.SubElement(configs, 'input_file', {
            'path': source['path'], 'sha256': source['sha256'],
            'association': ('explicit_TRDAT_PT_TEMPLATE' if plan['trusted_model_assignment']
                            else 'adjacent_unreferenced_candidate'),
        })
        ET.SubElement(inp, 'namelist', {'format': 'Fortran', 'encoding': 'UTF-8'}).text = source.get('text', '')
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
    # Whole-run occurrence has no single region; timed regional bounds and
    # numerical settings are included by region_plan_xml(plan).
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


ANOMALOUS_MODELS = frozenset(('glf23', 'tglf', 'mmm', 'rlw', 'coppi',
                              'cdbm', 'user', 'paleo'))
NEOCLASSICAL_MODELS = frozenset(('neoch', 'neogk', 'nclass', 'gtcneo'))


CATEGORIES = ('anomalous', 'neoclassical')


def _start_key(block):
    """Chronological sort key of a template block: no/negative update time means initial."""
    time = block['start_time']
    return -math.inf if time is None or time < 0 else float(time)


def _model_category(name):
    return ('anomalous' if name in ANOMALOUS_MODELS else
            'neoclassical' if name in NEOCLASSICAL_MODELS else None)


def _activation_events(plan):
    """Readable ``%active`` assignments of every template block, in chronological order.

    Blocks are explicit assignments, not state snapshots, so each event is kept as written
    (an unreadable ``active`` value is ignored, never guessed).
    """
    events = []
    for source in plan['templates']:
        for block in source['blocks']:
            for region, models in block['regions'].items():
                for name, options in models.items():
                    category = _model_category(name)
                    if category is None or 'active' not in options:
                        continue
                    active = _boolean(options['active'])
                    if active is None:
                        continue
                    events.append((_start_key(block), len(events), region, category, name, active))
    events.sort(key=lambda event: (event[0], event[1]))
    return events


def region_selection_timeline(plan):
    """Time-resolved model selection ``{region: {category: [{'model', 'intervals'}]}}``.

    Active flags accumulate across updates and a model enabled in a category replaces the
    previous one in that region.  Each model lists every interval ``(start, end)`` in which it
    was selected, ``start`` None meaning "from the initial state" and ``end`` None meaning "until
    the end of the run"; a model that is disabled later keeps its earlier interval (it is not
    dropped), and a model re-enabled later gets a further interval.  Models are ordered by their
    first activation.  Nothing is inferred for regions/categories the template never mentions.
    """
    timeline = {region: {category: [] for category in CATEGORIES} for region in REGIONS}
    current = {}                                # (region, category) -> [model, start]; start -inf = initial state

    def close(key, end):
        """Close the selection held at ``key`` at ``end`` (``-inf`` the initial instant, ``+inf`` the end of the run)."""
        model, start = current.pop(key)
        if end <= start:
            return                              # replaced at the instant it was enabled, including the initial instant
        region, category = key
        entry = next((e for e in timeline[region][category] if e['model'] == model), None)
        if entry is None:
            entry = {'model': model, 'intervals': []}
            timeline[region][category].append(entry)
        entry['intervals'].append((None if start == -math.inf else start, None if end == math.inf else end))

    for time, _order, region, category, name, active in _activation_events(plan):
        key = (region, category)
        held = current.get(key)
        if active:
            if held is not None and held[0] == name:
                continue
            if held is not None:
                close(key, time)
            current[key] = [name, time]
            # keep first-activation order even if the interval closes at once
            if not any(e['model'] == name for e in timeline[region][category]):
                timeline[region][category].append({'model': name, 'intervals': []})
        elif held is not None and held[0] == name:
            close(key, time)
    for key in list(current):
        close(key, math.inf)
    for region in REGIONS:
        for category in CATEGORIES:
            timeline[region][category] = [e for e in timeline[region][category] if e['intervals']]
    return timeline


def occurrence_region_models(plan):
    """One occurrence per region that selects a model at some time, ordered axial/confinement/edge.

    Each item is ``{'region', 'selections', 'final_models'}``: ``selections`` lists every
    ``{'category', 'model', 'intervals'}`` selected during the run (see
    :func:`region_selection_timeline`), so a GLF23 -> TGLF switch keeps both models and a model
    disabled before the end is not lost; ``final_models`` is only the model still active at the
    end of the run per category.  Only an explicit, readable PT_TEMPLATE can establish the active
    model, and a missing TR.DAT cannot establish radial extent.  Omitted configurations are never
    invented.
    """
    if not plan or not plan['trusted_model_assignment'] or not plan['boundary_states']:
        return []
    timeline = region_selection_timeline(plan)
    result = []
    for region in REGIONS:
        extent = any(_range_limits(region, state['values'], quantity)[1] >
                     _range_limits(region, state['values'], quantity)[0]
                     for state in plan['boundary_states'] for quantity in BOUNDARIES)
        if not extent:
            continue
        selections = [{'category': category, 'model': entry['model'], 'intervals': list(entry['intervals'])}
                      for category in CATEGORIES for entry in timeline[region][category]]
        if not selections:
            continue
        final = {category: next((e['model'] for e in timeline[region][category]
                                 if e['intervals'][-1][1] is None), None)
                 for category in CATEGORIES}
        result.append({'region': region, 'selections': selections, 'final_models': final})
    return result


def regional_model_slots(item, neoclassical_channels_available):
    """Model slots of one regional occurrence: one per *selected* model, no placeholders.

    Anomalous models get a slot each (every model selected at some time, so a GLF23 -> TGLF
    switch keeps both).  TRANSP publishes ONE assembled neoclassical result, not one per
    neoclassical model, so when those channels exist a single assembled slot lists every selected
    neoclassical model with its intervals and the output is not relabelled with any of their names.
    Without the channels each selected neoclassical model gets a settings-only slot.
    """
    slots = [{'category': 'anomalous', 'selected': [selection]}
             for selection in item['selections'] if selection['category'] == 'anomalous']
    neoclassical = [s for s in item['selections'] if s['category'] == 'neoclassical']
    if neoclassical:
        if neoclassical_channels_available:
            slots.append({'category': 'neoclassical', 'selected': neoclassical, 'assembled': True})
        else:
            slots.extend({'category': 'neoclassical', 'selected': [s]} for s in neoclassical)
    return slots


def interval_attributes(start, end):
    """``effective_from``/``effective_until`` attributes of a selection interval."""
    return {'effective_from': 'initial' if start is None else _fmt(start),
            'effective_until': 'end' if end is None else _fmt(end)}


# What is and is not established about how TRANSP combines models that are switched on in one region and category.
# The timeline assumes replacement; no TRANSP release source or documentation was available to verify it, and nothing
# else in this repository establishes it, so the status stays unresolved and the raw declarations are preserved.
MODEL_EXECUTION_SEMANTICS = {
    'status': 'unresolved',
    'assumed_rule': 'a model enabled in a category replaces the previous one in that region',
    'evidence': 'none supplied',
    'transp_version': 'not established',
}


def region_declarations(plan, region):
    """The readable ``%active`` declarations of ``region`` as written, chronologically, with no interpretation.

    ``[{'time': None | seconds, 'category', 'model', 'active'}]`` where ``time`` None is the initial state.  Unlike
    :func:`region_selection_timeline` this applies no replacement rule, so it stays valid whatever the execution
    semantics turn out to be.
    """
    return [{'time': None if key_time == -math.inf else key_time, 'category': category, 'model': name, 'active': active}
            for key_time, _order, event_region, category, name, active in _activation_events(plan)
            if event_region == region]


def region_occurrence_xml(plan, item, occurrence):
    """Self-describing configuration for one regional IDS occurrence."""
    root = ET.Element('pt_solver_region_occurrence', {
        'occurrence': str(occurrence), 'region': item['region'],
        'radial_coordinate': 'rho_tor_norm',
        'outputs_policy': 'write_only_verified_CDF_model_channels',
        'execution_semantics': MODEL_EXECUTION_SEMANTICS['status'],
        'execution_semantics_assumed_rule': MODEL_EXECUTION_SEMANTICS['assumed_rule'],
        'execution_semantics_evidence': MODEL_EXECUTION_SEMANTICS['evidence'],
        'execution_semantics_transp_version': MODEL_EXECUTION_SEMANTICS['transp_version'],
    })
    declarations = ET.SubElement(root, 'declarations', {
        'note': 'the %active assignments exactly as written, with no replacement or addition rule applied'})
    for declaration in region_declarations(plan, item['region']):
        ET.SubElement(declarations, 'declaration', {
            'time': 'initial' if declaration['time'] is None else _fmt(declaration['time']),
            'category': declaration['category'], 'model': declaration['model'].upper(),
            'active': str(declaration['active']).lower()})
    for selection in item['selections']:
        node = ET.SubElement(root, 'model_selection', {
            'category': selection['category'], 'name': selection['model'].upper(),
            'interpretation': 'replacement_assumed_unverified'})
        for start, end in selection['intervals']:
            ET.SubElement(node, 'interval', interval_attributes(start, end))
    # This deliberately retains ALL regional boundary changes and the exact
    # executed template, not just the final model names.
    root.append(region_validity_xml(plan, item['region']))
    root.append(region_plan_xml(plan))
    return ET.tostring(root, encoding='unicode')
