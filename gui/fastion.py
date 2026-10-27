#!/usr/bin/python
"""Fast-ion distribution functions from the ``distributions`` IDS, for plotting.

The converter stores a NUBEAM distribution on a general grid description (GGD) with three spaces: space 0 is the R,Z
poloidal plane of ``zones`` (2-D faces whose ``boundary`` lists edge indices, edges list node indices, all 1-based; each
face's ``geometry`` is [row, zone in row, volume in m^3]), space 1 the energy bin edges (eV) and space 2 the pitch bin
edges.  ``expansion_fd3v[0].grid_subset[0].values`` is the number of fast ions per bin and per m^3 (``fd3v``), flattened
with the zone index varying fastest, then energy, then pitch: ``reshape((zones, energy, pitch), order='F')``.  That layout
is not guessed: the volume-weighted total of the array equals ``global_quantities.particles_n``, and the volume-weighted
row averages equal ``profiles_1d.density_fast`` (see ``examples/verify_gui_plots.py``).  A file whose size does not
equal zones x energy x pitch is refused rather than reshaped.

The functions up to ``compute`` are pure (numpy only) so they are unit tested without an IMAS installation.
"""

import numpy as np

from gui import imas_worker as worker

GYRO_TYPES = {0: 'unspecified', 1: 'particle position', 2: 'guiding centre'}


def fortran_cube(values, n_zones, n_energy, n_pitch):
    """``values`` as an array indexed [zone, energy, pitch]; the size must match exactly."""
    values = np.asarray(values, dtype=float)
    if values.size != n_zones * n_energy * n_pitch:
        raise worker.WorkerError(f'the distribution has {values.size} values but the grid has {n_zones} zones x '
                                 f'{n_energy} energy x {n_pitch} pitch bins = {n_zones * n_energy * n_pitch}')
    return values.reshape((n_zones, n_energy, n_pitch), order='F')


def centers(edges):
    edges = np.asarray(edges, dtype=float)
    return 0.5 * (edges[:-1] + edges[1:])


def window(cube, energy_edges, pitch_edges, e_min=None, e_max=None, p_min=None, p_max=None):
    """``cube`` with the bins whose centre lies outside the energy / pitch window set to zero (a copy)."""
    keep_e = np.ones(len(energy_edges) - 1, dtype=bool)
    keep_p = np.ones(len(pitch_edges) - 1, dtype=bool)
    ec, pc = centers(energy_edges), centers(pitch_edges)
    if e_min is not None:
        keep_e &= ec >= e_min
    if e_max is not None:
        keep_e &= ec <= e_max
    if p_min is not None:
        keep_p &= pc >= p_min
    if p_max is not None:
        keep_p &= pc <= p_max
    return cube * keep_e[None, :, None] * keep_p[None, None, :]


def reduce(cube, volumes, energy_edges, pitch_edges, what, zone=None):
    """The requested reduction of ``cube`` (fast ions per bin per m^3, indexed [zone, energy, pitch]).

    All zones: volume-weighted, so the sum of everything is the number of fast ions.  One zone: the stored per-m^3 values.
    ``ep`` is [energy, pitch]; ``energy`` and ``pitch`` are spectra per unit energy (eV) and per unit pitch; ``zones`` the
    density per zone in m^-3 (all energies and pitches in the window).
    """
    de = np.diff(np.asarray(energy_edges, dtype=float))
    dp = np.diff(np.asarray(pitch_edges, dtype=float))
    weighted = cube[zone][None] if zone is not None else cube * np.asarray(volumes, dtype=float)[:, None, None]
    if what == 'ep':
        return weighted.sum(axis=0)
    if what == 'energy':
        return weighted.sum(axis=(0, 2)) / de
    if what == 'pitch':
        return weighted.sum(axis=(0, 1)) / dp
    if what == 'zones':
        return cube.sum(axis=(1, 2))
    raise worker.WorkerError(f'unknown reduction {what!r}')


def chain_polygon(edges):
    """Vertex ids of the closed loop formed by ``edges`` (pairs of node ids), in order; ``None`` if they do not close.

    Edge order and direction in the input do not matter.  A node with other than two edges (a pinched or open outline)
    means the zone is not a simple polygon and is not drawn.
    """
    adjacency = {}
    for a, b in edges:
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)
    if len(edges) < 3 or any(len(v) != 2 for v in adjacency.values()):
        return None
    start = edges[0][0]
    loop, previous, current = [start], None, start
    while True:
        first, second = adjacency[current]
        following = first if first != previous else second
        if previous is None:
            following = first
        if following == start:
            break
        loop.append(following)
        previous, current = current, following
        if len(loop) > len(edges):
            return None
    return loop if len(loop) == len(adjacency) else None


def polygon_area(xy):
    x, y = np.asarray(xy, dtype=float).T
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def _grid(ggd):
    spaces = ggd.grid.space
    if len(spaces) < 3:
        raise worker.WorkerError('the distribution grid does not have the expected three spaces (zones, energy, pitch)')
    zones = spaces[0].objects_per_dimension[2].object
    volumes = np.array([float(zones[i].geometry[2]) for i in range(len(zones))])
    energy = np.array([float(n.geometry[0]) for n in (spaces[1].objects_per_dimension[0].object[i]
                                                      for i in range(len(spaces[1].objects_per_dimension[0].object)))])
    pitch = np.array([float(n.geometry[0]) for n in (spaces[2].objects_per_dimension[0].object[i]
                                                     for i in range(len(spaces[2].objects_per_dimension[0].object)))])
    return spaces, volumes, energy, pitch


def zone_polygons(spaces):
    """[(vertex list of [R, Z], area in m^2 as stored)] per zone; ``None`` for a zone that is not a simple polygon."""
    space = spaces[0].objects_per_dimension
    faces, edges, nodes = space[2].object, space[1].object, space[0].object
    cache = {}

    def point(index):
        if index not in cache:
            geometry = nodes[index - 1].geometry
            cache[index] = [float(geometry[0]), float(geometry[1])]
        return cache[index]

    polygons = []
    for i in range(len(faces)):
        face = faces[i]
        pairs = []
        for k in range(len(face.boundary)):
            nodes_of_edge = edges[int(face.boundary[k].index) - 1].nodes
            pairs.append((int(nodes_of_edge[0]), int(nodes_of_edge[1])))
        loop = chain_polygon(pairs)
        polygons.append((None if loop is None else [point(n) for n in loop], float(face.measure)))
    return polygons


def smooth_zone_polygons(spaces, axis=None, samples=14):
    """Rebuild ordered MC-zone outlines and smooth their constant-flux surfaces.

    The stored GGD face boundary is an unordered edge list. Chaining those
    edges can produce a self-crossing vertex order; the face's ``nodes`` list
    is ordered and its geometry records the radial row. This follows the
    standalone plotting routine's periodic polar spline reconstruction.
    """
    from scipy.interpolate import CubicSpline

    objects = spaces[0].objects_per_dimension
    nodes = np.asarray([node.geometry for node in objects[0].object], dtype=float)[:, :2]
    faces = objects[2].object
    face_nodes = []
    rows = []
    legacy = None
    for index, face in enumerate(faces):
        ordered = [int(n) - 1 for n in getattr(face, 'nodes', [])]
        if len(ordered) < 3 or any(n < 0 or n >= len(nodes) for n in ordered):
            # Older/synthetic entries may lack the ordered node list. Retain
            # the topological polygon for those entries.
            if legacy is None:
                legacy = zone_polygons(spaces)
            fallback, _area = legacy[index]
            face_nodes.append(None if fallback is None else np.asarray(fallback, dtype=float))
        else:
            face_nodes.append(ordered)
        rows.append(max(0, int(float(face.geometry[0])) - 1))

    # If any face required the legacy fallback, return its original polygons;
    # mixing coordinate arrays and node-index lists would corrupt shared edges.
    if any(face is None or isinstance(face, np.ndarray) for face in face_nodes):
        return [(face.tolist() if isinstance(face, np.ndarray) else
                 [nodes[n].tolist() for n in face] if face is not None else None)
                for face in face_nodes]

    if axis is None:
        axis = nodes[face_nodes[0][0]] if face_nodes else np.mean(nodes, axis=0)
    axis = np.asarray(axis, dtype=float)

    level = {}
    for face, row in zip(face_nodes, rows):
        for node in face[1:3]:
            level[node] = row + 1
        if len(face) == 4:
            for node in (face[0], face[3]):
                level[node] = row
    surfaces = {}
    for surface in set(level.values()):
        surface_nodes = [node for node, row in level.items() if row == surface]
        rel = nodes[surface_nodes] - axis
        phi = np.arctan2(rel[:, 1], rel[:, 0])
        order = np.argsort(phi)
        phi, radius = phi[order], np.hypot(rel[order, 0], rel[order, 1])
        if len(phi) >= 4 and np.all(np.diff(phi) > 1e-9):
            surfaces[surface] = CubicSpline(np.r_[phi, phi[0] + 2 * np.pi],
                                            np.r_[radius, radius[0]], bc_type='periodic')

    def arc(first, second, surface):
        delta_a, delta_b = nodes[first] - axis, nodes[second] - axis
        a, b = np.arctan2(delta_a[1], delta_a[0]), np.arctan2(delta_b[1], delta_b[0])
        if surface not in surfaces:
            return np.asarray([nodes[first], nodes[second]])
        sweep = (b - a + np.pi) % (2 * np.pi) - np.pi
        angle = a + np.linspace(0, sweep, samples)
        start = surfaces[surface].x[0]
        wrapped = start + (angle - start) % (2 * np.pi)
        radius = surfaces[surface](wrapped)
        return np.column_stack((axis[0] + radius * np.cos(angle), axis[1] + radius * np.sin(angle)))

    polygons = []
    for face, row in zip(face_nodes, rows):
        if len(face) == 3:
            points = [nodes[face[0]][None, :], arc(face[1], face[2], row + 1)]
        else:
            points = [nodes[face[0]][None, :], arc(face[1], face[2], row + 1),
                      arc(face[3], face[0], row)]
        polygons.append(np.vstack(points).tolist())
    return polygons


def _magnetic_axis(ids, time_index):
    """Per-time magnetic axis when available, otherwise let the grid supply it."""
    try:
        r, z = np.asarray(ids.magnetic_axis.r, dtype=float), np.asarray(ids.magnetic_axis.z, dtype=float)
        if time_index < len(r) and time_index < len(z):
            return float(r[time_index]), float(z[time_index])
    except (AttributeError, IndexError, TypeError, ValueError):
        pass
    return None


def label(distribution, index):
    parts = [f'distribution {index}']
    try:
        parts.append(GYRO_TYPES.get(int(distribution.gyro_type), f'gyro_type {int(distribution.gyro_type)}'))
    except (AttributeError, ValueError, TypeError):
        pass
    try:
        parts.append(str(distribution.species.ion.label).strip() or 'ion')
    except (AttributeError, ValueError, TypeError):
        pass
    try:
        parts.append(str(distribution.process[0].type.name))
    except (AttributeError, IndexError, ValueError, TypeError):
        pass
    return ', '.join(p for p in parts if p)


def compute(entry, request):
    """The payload for ``{dist, time_index, what: info|ep|energy|pitch|zones, zone, e_min, e_max, p_min, p_max}``."""
    ids = entry.get('distributions', int(request.get('occurrence') or 0), lazy=True)
    n_dist = len(ids.distribution)
    if n_dist == 0:
        raise worker.WorkerError('this entry has no fast-ion distribution')
    times = worker._times(ids)
    what = request.get('what') or 'info'
    di = int(request.get('dist') or 0)
    if not 0 <= di < n_dist:
        raise worker.WorkerError(f'distribution {di} does not exist (0 to {n_dist - 1})')
    distribution = ids.distribution[di]
    n_times = len(distribution.ggd)
    ti = min(max(int(request.get('time_index') or 0), 0), max(n_times - 1, 0))
    if n_times == 0:
        raise worker.WorkerError('this distribution has no grid-based (ggd) data')
    ggd = distribution.ggd[ti]
    spaces, volumes, energy, pitch = _grid(ggd)
    n_zones, n_e, n_p = len(volumes), len(energy) - 1, len(pitch) - 1
    time_info = {'index': ti, 'count': int(len(times)), 'values': worker.cleaned(times), 'stride': 1,
                 'value': worker._scalar(times[ti]) if ti < len(times) else None}
    if what == 'info':
        return {'kind': 'df-info', 'distributions': [{'index': i, 'label': label(ids.distribution[i], i)} for i in range(n_dist)],
                'zones': n_zones, 'energy_bins': n_e, 'pitch_bins': n_p, 'energy_range': [float(energy[0]), float(energy[-1])],
                'pitch_range': [float(pitch[0]), float(pitch[-1])], 'time': time_info}
    cube = fortran_cube(ggd.expansion_fd3v[0].grid_subset[0].values, n_zones, n_e, n_p)
    cube = window(cube, energy, pitch, *(request.get(k) for k in ('e_min', 'e_max', 'p_min', 'p_max')))
    zone = request.get('zone')
    if zone is not None:
        zone = int(zone)
        if not 0 <= zone < n_zones:
            raise worker.WorkerError(f'zone {zone} does not exist (0 to {n_zones - 1})')
    result = reduce(cube, volumes, energy, pitch, what, zone)
    name = label(distribution, di) + (f', zone {zone}' if zone is not None else ', all zones')
    payload = {'ids': 'distributions', 'time': time_info, 'time_dependent': True, 'missing': [], 'notes': [], 'series': [],
               'title': name, 'zone': zone, 'dist': di}
    per_zone = zone is not None
    if what == 'ep':
        payload.update(kind='map', x=worker.cleaned(centers(pitch)), y=worker.cleaned(centers(energy)), z=worker.cleaned(result),
                       x_label='pitch', x_units='', y_label='energy', y_units='eV', z_label='fd3v' if per_zone else 'particles per bin',
                       z_units='m^-3' if per_zone else '', stride=[1, 1], overlays=[], points=[])
    elif what in ('energy', 'pitch'):
        axis = energy if what == 'energy' else pitch
        unit = 'eV' if what == 'energy' else ''
        payload.update(kind='profile', x_label='energy' if what == 'energy' else 'pitch', x_units=unit)
        payload['series'] = [{'name': 'spectrum' if not per_zone else f'spectrum, zone {zone}',
                              'path': 'distributions', 'units': ('m^-3' if per_zone else 'particles') + (' / eV' if what == 'energy' else ''),
                              'x': worker.cleaned(centers(axis)), 'y': worker.cleaned(result), 'stride': 1}]
    else:
        polygons = smooth_zone_polygons(spaces, axis=_magnetic_axis(ids, ti))
        drawn = [{'xy': xy, 'value': float(result[i]), 'zone': i} for i, xy in enumerate(polygons) if xy is not None]
        if len(drawn) < len(polygons):
            payload['notes'].append(f'{len(polygons) - len(drawn)} zones have incomplete geometry and are not drawn')
        payload.update(kind='polys', polys=drawn, x_label='R', x_units='m', y_label='Z', y_units='m',
                       z_label='fast-ion density', z_units='m^-3', overlays=worker.wall_overlays(entry), points=[])
    return payload
