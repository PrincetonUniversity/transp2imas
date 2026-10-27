#!/usr/bin/python
"""Read converted IMAS entries for plotting.  Runs in the IMAS environment, as a subprocess of the viewer.

``python -m gui.imas_worker '<json request>'`` reads one entry (HDF5 only) and prints one line ``GUI_RESULT:<json>``.
It never writes: entries are opened ``'r'``.  Requests: ``describe`` (IDSs in an entry, their time ranges and the recipes
that apply), ``nodes`` (the plottable numeric nodes of one IDS) and ``data`` (a chart payload for node paths or a recipe).

Node paths use ``/`` between names, ``[:]`` for the time index of a time-dependent array of structures (detected from the
Data Dictionary metadata, not from a name list), ``[i]`` for a fixed element and ``[*]`` for every existing element.
Chart payloads are bounded (series stride-reduced to ``MAX_POINTS``, maps to ``MAX_GRID`` per axis) and never contain
NaN: a nonfinite sample is ``null``.  Missing nodes are reported as missing, not drawn as zeros.
"""

import json
import math
from pathlib import Path
import re
import sys

import numpy as np

PREFIX = 'GUI_RESULT:'
MAX_POINTS = 1500
MAX_GRID = 200
MAX_ELEMENTS = 12                         # elements of a non-time array of structures expanded by [*] or discovery
MAX_NODES = 600
MAX_OVERLAY_TIMES = 8
MAX_FRAMES = 300                          # time slices sent at once for local scrubbing of a profile plot
DD_BY_MAJOR = {'3': '3.39.0', '4': '4.1.1'}
SEGMENT = re.compile(r'^([A-Za-z_][A-Za-z0-9_]*)(?:\[(\d+|:|\*)\])?$')


class WorkerError(ValueError):
    pass


# ---------------------------------------------------------------- pure helpers (testable without IMAS)
def parse_path(path):
    """``'profiles_1d[:]/ion[*]/temperature'`` -> ``[('profiles_1d', ':'), ('ion', '*'), ('temperature', None)]``."""
    if not isinstance(path, str) or not path or len(path) > 300:
        raise WorkerError('a node path is required')
    segments = []
    for part in path.split('/'):
        match = SEGMENT.match(part)
        if not match:
            raise WorkerError(f'bad path segment {part!r}')
        name, index = match.groups()
        segments.append((name, index if index in (None, ':', '*') else int(index)))
    return segments


def clean(values):
    """A JSON-safe nested list: nonfinite numbers become ``None`` (the page draws a gap), bools and ints stay."""
    array = np.asarray(values)
    if array.dtype.kind in 'fc':
        return np.where(np.isfinite(array.real), array.real, np.nan).astype(float).tolist() if array.ndim else _scalar(array)
    return array.tolist()


def _scalar(value):
    value = float(np.asarray(value).real)
    return value if math.isfinite(value) else None


def _none_for_nan(obj):
    if isinstance(obj, list):
        return [_none_for_nan(x) for x in obj]
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


def cleaned(values):
    return _none_for_nan(clean(values))


def stride_for(length, limit=MAX_POINTS):
    return max(1, -(-int(length) // limit))


def entry_dd_version(entry_dir, requested=None):
    """The DD version to open an entry with: the requested one, else from the ``<machine>/<major>/<shot>/<run>`` path."""
    if requested:
        return requested
    parts = Path(entry_dir).parts
    return DD_BY_MAJOR.get(parts[-3]) if len(parts) >= 3 else None


def entry_files(entry_dir):
    """``[(ids, occurrence, bytes)]`` from the ``<ids>[_<occurrence>].h5`` files of an entry directory."""
    found = []
    for path in sorted(Path(entry_dir).glob('*.h5')):
        if path.name == 'master.h5':
            continue
        match = re.fullmatch(r'(.+?)(?:_(\d+))?\.h5', path.name)
        found.append((match.group(1), int(match.group(2) or 0), path.stat().st_size))
    return found


def downsample_series(x, y, limit=MAX_POINTS):
    step = stride_for(len(y), limit)
    return cleaned(np.asarray(x, dtype=float)[::step]), cleaned(np.asarray(y, dtype=float)[::step]), step


# ---------------------------------------------------------------- IMAS access
def open_entry(entry_dir, dd_version=None):
    import imas
    path = Path(entry_dir)
    if not (path / 'master.h5').is_file():
        raise WorkerError(f'{entry_dir} is not an HDF5 IMAS entry (no master.h5); only HDF5 entries can be plotted')
    version = entry_dd_version(path, dd_version)
    kwargs = {'dd_version': version} if version else {}
    return imas.DBEntry(f'imas:hdf5?path={path}', 'r', **kwargs)


def _times(ids):
    try:
        values = np.asarray(ids.time, dtype=float)
    except (AttributeError, TypeError, ValueError):
        return np.array([])
    return values


def _is_time_aos(node):
    try:
        return any('(itime)' in str(c) for c in node.metadata.coordinates[:1])
    except (AttributeError, IndexError, TypeError):
        return False


def _numeric(node):
    """A real or integer node of at most two dimensions (strings, complex numbers and higher ranks are not plotted)."""
    try:
        data_type = node.metadata.data_type
        return getattr(data_type, 'name', str(data_type)) in ('FLT', 'INT') and node.metadata.ndim <= 2
    except AttributeError:
        return False


def _element_label(element, index):
    """A useful AoS label, preferring physical identity over generic DD ``label`` fields."""
    def text(value):
        if value is None:
            return ''
        value = getattr(value, 'value', value)
        try:
            value = str(value).strip()
        except Exception:
            return ''
        return '' if value.lower() in ('', 'value', 'none', 'unknown') else value

    # Process entries (for example radiation.process) identify their physics here.
    identifier = getattr(element, 'identifier', None)
    result = text(getattr(identifier, 'name', None))
    if not result:
        result = text(getattr(identifier, 'description', None))

    # Distribution entries often have no label or identifier. Their species is the
    # information needed to tell density_fast[0] from density_fast[1].
    if not result:
        species = getattr(element, 'species', None)
        for part in ('ion', 'neutral', 'element'):
            result = text(getattr(getattr(species, part, None), 'label', None)) or \
                     text(getattr(getattr(species, part, None), 'name', None))
            if result:
                break
        if not result:
            result = text(getattr(getattr(species, 'type', None), 'name', None))

    if not result:
        result = text(getattr(getattr(element, 'type', None), 'name', None))
    if not result:
        for attribute in ('name', 'label'):
            result = text(getattr(element, attribute, None))
            if result:
                break
    # NUBEAM may store separate distributions for the same ion species (for
    # example guiding-centre and particle-position distributions). Include
    # those descriptors so that their profile legends remain distinguishable.
    qualifiers = []
    gyro = getattr(element, 'gyro_type', None)
    if gyro is not None:
        try:
            qualifiers.append({0: 'unspecified', 1: 'particle position', 2: 'guiding centre'}.get(
                int(getattr(gyro, 'value', gyro)), f'gyro type {int(getattr(gyro, "value", gyro))}'))
        except (TypeError, ValueError):
            pass
    processes = getattr(element, 'process', None)
    try:
        process = processes[0]
        process_name = text(getattr(getattr(process, 'type', None), 'name', None))
        if process_name:
            qualifiers.append(process_name)
    except (AttributeError, IndexError, TypeError):
        pass
    details = ', '.join(([result] if result else []) + qualifiers)
    if details:
        return f'{index}: {details}'
    return str(index)


def discover(ids):
    """Plottable numeric nodes of ``ids`` (a single-time-slice IDS): ``[{path, kind, shape, units, label}]``."""
    nodes = []

    def walk(structure, trail, shown, has_time):
        for child in structure.iter_nonempty_():
            if len(nodes) >= MAX_NODES:
                return
            name = child.metadata.name
            if name in ('ids_properties', 'code', 'time') and not trail:
                continue
            kind_name = type(child).__name__
            if kind_name == 'IDSStructArray':
                if len(child) == 0:
                    continue
                if _is_time_aos(child):
                    walk(child[0], trail + [f'{name}[:]'], shown + [name], True)
                else:
                    for index in range(min(len(child), MAX_ELEMENTS)):
                        walk(child[index], trail + [f'{name}[{index}]'], shown + [f'{name}[{_element_label(child[index], index)}]'],
                             has_time)
            elif kind_name == 'IDSStructure':
                walk(child, trail + [name], shown, has_time)
            elif _numeric(child):
                ndim = child.metadata.ndim
                if ndim == 0 and not has_time:
                    continue
                kind = ('timetrace' if ndim == 0 else 'profile' if ndim == 1 else 'map') if has_time else \
                    ('timetrace' if ndim == 1 and str(child.metadata.coordinates[0]).endswith('time') else
                     'profile' if ndim == 1 else 'map')
                value = np.asarray(child.value)
                nodes.append({'path': '/'.join(trail + [name]), 'kind': kind, 'shape': list(value.shape) if value.ndim else [],
                              'units': child.metadata.units, 'label': ' / '.join(shown + [name])})

    walk(ids, [], [], False)
    return nodes


def _resolve(root, segments, time_index, element=None):
    """Follow ``segments`` from ``root``; ``[:]`` takes ``time_index`` and ``[*]`` the given ``element``."""
    node = root
    for name, index in segments:
        node = getattr(node, name)
        if index is None:
            continue
        if index == ':':
            node = node[time_index]
        elif index == '*':
            node = node[element]
        else:
            node = node[index]
    return node


def _expand(root, segments):
    """Concrete paths for ``[*]`` segments: ``[(label suffix, [segments with a number])]``."""
    results = [([], [])]
    node_options = [(root, [], [])]
    for position, (name, index) in enumerate(segments):
        expanded = []
        for node, done, labels in node_options:
            try:
                child = getattr(node, name)
            except AttributeError:
                continue
            if index == '*':
                for i in range(min(len(child), MAX_ELEMENTS)):
                    expanded.append((child[i], done + [(name, i)], labels + [_element_label(child[i], i)]))
            elif index == ':':
                expanded.append((child[0] if len(child) else None, done + [(name, ':')], labels))
            elif isinstance(index, int):
                expanded.append((child[index] if len(child) > index else None, done + [(name, index)], labels))
            else:
                expanded.append((child, done + [(name, None)], labels))
        node_options = [(n, d, l) for n, d, l in expanded if n is not None or position == len(segments) - 1]
    return [(' '.join(labels), done) for _, done, labels in node_options]


def _trace(ids, segments, times, limit):
    values, xs = [], []
    n = len(times) if len(times) else 0
    step = stride_for(n, limit)
    for i in range(0, n, step):
        try:
            values.append(float(np.asarray(_resolve(ids, segments, i).value).reshape(-1)[0]))
            xs.append(float(times[i]))
        except (IndexError, AttributeError, ValueError):
            values.append(float('nan'))
            xs.append(float(times[i]))
    return xs, values, step


def _units(value):
    """Units for an axis label; the Data Dictionary's 'mixed' (it depends on the grid type) is not a unit."""
    return '' if not value or value == 'mixed' else value


def _coordinate(node, axis):
    """Values, name and units of the node's own coordinate ``axis`` (index values if it has none)."""
    try:
        coordinate = node.coordinates[axis]
        values = np.asarray(coordinate.value, dtype=float)
        return values, coordinate.metadata.name, coordinate.metadata.units
    except (AttributeError, IndexError, ValueError, TypeError, RuntimeError):
        return None, None, None


def extract(entry, ids_name, occurrence, paths, time_index=0, all_times=False, kind=None, frames=False):
    """A chart payload for the node ``paths`` (same kind) of one IDS."""
    ids = entry.get(ids_name, occurrence, lazy=True)
    times = _times(ids)
    series, notes, missing = [], [], []
    payload = {'ids': ids_name, 'occurrence': occurrence, 'notes': notes, 'missing': missing}
    first_kind = None
    for path in paths:
        segments = parse_path(path)
        for label, concrete in _expand(ids, segments):
            time_aos = any(index == ':' for _, index in concrete)
            try:
                node = _resolve(ids, concrete, min(time_index, max(0, len(times) - 1)) if time_aos else 0)
                value = np.asarray(node.value)
            except (AttributeError, IndexError, KeyError, ValueError, TypeError, RuntimeError):
                missing.append(f'{path} [{label}]' if label else path)
                continue
            if value.size == 0 and value.ndim:
                missing.append(f'{path} [{label}]' if label else path)
                continue
            ndim = value.ndim
            this_kind = ('timetrace' if ndim == 0 else 'profile' if ndim == 1 else 'map') if time_aos else \
                ('timetrace' if ndim == 1 and str(node.metadata.coordinates[0]).endswith('time') else
                 'profile' if ndim == 1 else 'map')
            first_kind = first_kind or this_kind
            if this_kind != first_kind:
                raise WorkerError('the selected nodes are not the same kind of plot (profile, time trace, map)')
            name = (path.split('/')[-1] if not label else f'{path.split("/")[-1]} [{label}]')
            units = node.metadata.units
            if this_kind == 'timetrace':
                if time_aos:
                    x, y, step = _trace(ids, concrete, times, MAX_POINTS)
                    x_name, x_units = 'time', 's'
                else:
                    xv, x_name, x_units = _coordinate(node, 0)
                    x, y, step = downsample_series(np.arange(len(value)) if xv is None else xv, value)
                    x_name = x_name or 'index'
                series.append({'name': name, 'path': path, 'units': units, 'x': x, 'y': _none_for_nan(y), 'stride': step})
                payload.update(x_label=x_name, x_units=x_units)
            elif this_kind == 'profile':
                xv, x_name, x_units = _coordinate(node, 0)
                chosen = [min(time_index, max(0, len(times) - 1))]
                if all_times and time_aos and len(times) > 1:
                    chosen = sorted({int(i) for i in np.linspace(0, len(times) - 1, min(MAX_OVERLAY_TIMES, len(times)))})
                for ti in chosen:
                    current = _resolve(ids, concrete, ti) if time_aos else node
                    yv = np.asarray(current.value, dtype=float)
                    cx, _, _ = _coordinate(current, 0)
                    if cx is None or len(cx) != len(yv):
                        cx = np.arange(len(yv), dtype=float)
                    x, y, step = downsample_series(cx, yv)
                    tag = f' @ t={times[ti]:.4g} s' if time_aos and len(times) and len(chosen) > 1 else ''
                    series.append({'name': name + tag, 'path': path, 'units': units, 'x': x, 'y': y, 'stride': step})
                payload.update(x_label=x_name or 'index', x_units=x_units or '')
                if frames and time_aos and len(times):
                    indices = sorted({int(i) for i in np.linspace(0, len(times) - 1, min(MAX_FRAMES, len(times)))})
                    collected = payload.setdefault('frames', [{'index': i, 'time': _scalar(times[i]), 'series': []} for i in indices])
                    for frame in collected:
                        try:
                            current = _resolve(ids, concrete, frame['index'])
                            fy = np.asarray(current.value, dtype=float)
                            fx, _, _ = _coordinate(current, 0)
                            if fx is None or len(fx) != len(fy):
                                fx = np.arange(len(fy), dtype=float)
                            fxs, fys, _ = downsample_series(fx, fy)
                            frame['series'].append({'name': name, 'x': fxs, 'y': fys})
                        except (AttributeError, IndexError, ValueError, TypeError, RuntimeError):
                            frame['series'].append({'name': name, 'x': [], 'y': []})
            else:
                xv, x_name, x_units = _coordinate(node, 0)
                yv, y_name, y_units = _coordinate(node, 1)
                z = np.asarray(value, dtype=float)
                if xv is None or len(xv) != z.shape[0]:
                    xv = np.arange(z.shape[0], dtype=float)
                if yv is None or len(yv) != z.shape[1]:
                    yv = np.arange(z.shape[1], dtype=float)
                sx, sy = stride_for(z.shape[0], MAX_GRID), stride_for(z.shape[1], MAX_GRID)
                payload.update(kind='map', x=cleaned(xv[::sx]), y=cleaned(yv[::sy]), z=cleaned(z[::sx, ::sy].T),
                               x_label=x_name or 'dim1', x_units=_units(x_units), y_label=y_name or 'dim2', y_units=_units(y_units),
                               z_label=name, z_units=units, stride=[sx, sy], overlays=[], points=[])
                if time_aos and path.startswith('time_slice'):
                    _equilibrium_overlays(ids, min(time_index, max(0, len(times) - 1)), payload, entry)
                series.append({'name': name, 'path': path, 'units': units})
    if not series:
        raise WorkerError('none of the requested nodes exists in this entry: ' + ', '.join(missing or paths))
    payload['kind'] = payload.get('kind') or first_kind
    payload['series'] = [] if payload['kind'] == 'map' else series
    payload['time'] = {'index': int(min(time_index, max(0, len(times) - 1))), 'count': int(len(times)),
                       'values': cleaned(times[::stride_for(len(times), 2000)]), 'stride': stride_for(len(times), 2000),
                       'value': _scalar(times[min(time_index, len(times) - 1)]) if len(times) else None} if len(times) else None
    payload['time_dependent'] = bool(any(parse_path(p)[i][1] == ':' for p in paths for i in range(len(parse_path(p)))))
    if len(times) == 0 or not payload['time_dependent']:
        payload['time'] = payload['time'] if payload['time_dependent'] else None
    return payload


WALL_PAIRS = [['description_2d[0]/limiter/unit[*]/outline/r', 'description_2d[0]/limiter/unit[*]/outline/z']]
UNSET = 1.0e30                            # the Data Dictionary's empty-float sentinel is about -9e40


def rz_curves(entry, ids_name, occurrence, pairs):
    """Outlines in the R,Z plane for (r path, z path) pairs; ``[*]`` expands the same way for both paths."""
    ids = entry.get(ids_name, occurrence, lazy=True)
    curves, missing = [], []
    for r_path, z_path in pairs:
        r_items = _expand(ids, parse_path(r_path))
        z_items = _expand(ids, parse_path(z_path))
        for (label, r_concrete), (_, z_concrete) in zip(r_items, z_items):
            try:
                r = np.asarray(_resolve(ids, r_concrete, 0).value, dtype=float)
                z = np.asarray(_resolve(ids, z_concrete, 0).value, dtype=float)
            except (AttributeError, IndexError, KeyError, ValueError, TypeError, RuntimeError):
                r = z = np.array([])
            if r.size == 0 or r.size != z.size or r.ndim != 1:
                missing.append(f'{r_path} [{label}]' if label else r_path)
                continue
            curves.append({'name': label or ids_name, 'x': cleaned(r), 'y': cleaned(z), 'points': int(r.size),
                           'closed': bool(np.allclose([r[0], z[0]], [r[-1], z[-1]]))})
    if not curves:
        raise WorkerError('no outline in this entry: ' + ', '.join(missing or [p[0] for p in pairs]))
    return {'kind': 'rz', 'ids': ids_name, 'occurrence': occurrence, 'curves': curves, 'missing': missing, 'notes': [],
            'x_label': 'R', 'x_units': 'm', 'y_label': 'Z', 'y_units': 'm', 'series': [], 'time': None, 'time_dependent': False}


def wall_overlays(entry):
    """The wall outline(s) as overlay curves, or none when the entry has no wall (never an error: it is only context)."""
    try:
        curves = rz_curves(entry, 'wall', 0, WALL_PAIRS)['curves']
    except (WorkerError, RuntimeError, ValueError, OSError, KeyError):
        return []
    return [{'name': 'wall' if len(curves) == 1 else f'wall ({c["name"]})', 'x': c['x'], 'y': c['y']} for c in curves]


def _equilibrium_overlays(ids, index, payload, entry=None):
    """Wall, boundary outline, magnetic axis and the flux range of one time slice, when the entry has them."""
    if entry is not None:
        payload['overlays'].extend(wall_overlays(entry))
    try:
        quantities = ids.time_slice[index].global_quantities
        axis_flux, boundary_flux = float(quantities.psi_axis), float(quantities.psi_boundary)
        if all(math.isfinite(v) and abs(v) < UNSET for v in (axis_flux, boundary_flux)) and axis_flux != boundary_flux:
            payload['flux'] = {'axis': axis_flux, 'boundary': boundary_flux}
    except (AttributeError, IndexError, ValueError, TypeError):
        pass
    try:
        outline = ids.time_slice[index].boundary.outline
        r, z = np.asarray(outline.r.value, dtype=float), np.asarray(outline.z.value, dtype=float)
        if r.size and r.size == z.size:
            payload['overlays'].append({'name': 'boundary', 'x': cleaned(r), 'y': cleaned(z)})
    except (AttributeError, IndexError, ValueError):
        pass
    try:
        axis = ids.time_slice[index].global_quantities.magnetic_axis
        r, z = float(axis.r), float(axis.z)
        if math.isfinite(r) and math.isfinite(z):
            payload['points'].append({'name': 'magnetic axis', 'x': r, 'y': z})
    except (AttributeError, IndexError, ValueError, TypeError):
        pass


def describe(entry_dir, dd_version=None):
    """The IDSs of an entry with their time ranges, and the recipes whose IDS is present."""
    from gui import recipes
    entry = open_entry(entry_dir, dd_version)
    try:
        present = []
        for name, occurrence, size in entry_files(entry_dir):
            item = {'ids': name, 'occurrence': occurrence, 'bytes': size, 'n_times': None, 't_first': None, 't_last': None}
            try:
                ids = entry.get(name, occurrence, lazy=True)
                times = _times(ids)
                if len(times):
                    item.update(n_times=int(len(times)), t_first=_scalar(times[0]), t_last=_scalar(times[-1]))
            except (RuntimeError, ValueError, OSError) as exc:
                item['error'] = str(exc)
            present.append(item)
        version = entry.factory.version
    finally:
        entry.close()
    have = {p['ids'] for p in present if 'error' not in p}
    return {'entry': str(entry_dir), 'dd_version': version, 'ids': present,
            'recipes': [{k: r[k] for k in ('id', 'title', 'ids', 'kind', 'group')} for r in recipes.RECIPES if r['ids'] in have]}


def nodes(entry_dir, ids_name, occurrence=0, dd_version=None):
    """Plottable nodes of one IDS, read from its first time slice (cheap on large IDSs)."""
    from imas import ids_defs
    entry = open_entry(entry_dir, dd_version)
    try:
        lazy = entry.get(ids_name, occurrence, lazy=True)
        times = _times(lazy)
        ids = entry.get_slice(ids_name, float(times[0]), ids_defs.CLOSEST_INTERP, occurrence) if len(times) else \
            entry.get(ids_name, occurrence)
        return {'ids': ids_name, 'occurrence': occurrence, 'nodes': discover(ids)}
    finally:
        entry.close()


def data(entry_dir, request, dd_version=None):
    from gui import recipes
    entry = open_entry(entry_dir, dd_version)
    try:
        if request.get('recipe'):
            recipe = recipes.by_id(request['recipe'])
            paths, ids_name, kind = recipe['paths'], recipe['ids'], recipe['kind']
            title = recipe['title']
            if kind == 'rz':
                payload = rz_curves(entry, ids_name, int(request.get('occurrence') or 0), recipe['pairs'])
                payload['title'] = title
                return payload
        else:
            paths, ids_name, kind, title = request.get('paths') or [], request.get('ids'), None, None
        if not ids_name or not paths or not isinstance(paths, list) or len(paths) > 16:
            raise WorkerError('give a recipe, or an ids and one to sixteen node paths')
        payload = extract(entry, ids_name, int(request.get('occurrence') or 0), paths, int(request.get('time_index') or 0),
                          bool(request.get('all_times')), kind, bool(request.get('frames')))
    finally:
        entry.close()
    payload['title'] = title or ', '.join(p.split('/')[-1] for p in paths)
    return payload


def df(entry_dir, request, dd_version=None):
    from gui import fastion
    entry = open_entry(entry_dir, dd_version)
    try:
        return fastion.compute(entry, request)
    finally:
        entry.close()


def run(request):
    action = request.get('action')
    entry_dir, dd_version = request.get('entry'), request.get('dd_version')
    if not entry_dir:
        raise WorkerError('no entry given')
    if action == 'describe':
        return describe(entry_dir, dd_version)
    if action == 'nodes':
        return nodes(entry_dir, request.get('ids'), int(request.get('occurrence') or 0), dd_version)
    if action == 'data':
        return data(entry_dir, request, dd_version)
    if action == 'df':
        return df(entry_dir, request, dd_version)
    raise WorkerError(f'unknown action {action!r}')


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    try:
        result = run(json.loads(argv[0]))
        line = json.dumps(result, allow_nan=False)
        print(PREFIX + line, flush=True)
        return 0
    except (WorkerError, KeyError, RuntimeError, ValueError, OSError, ImportError, IndexError, json.JSONDecodeError) as exc:
        print(PREFIX + json.dumps({'error': f'{type(exc).__name__}: {exc}'}), flush=True)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
