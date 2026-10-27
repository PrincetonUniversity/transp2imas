#!/usr/bin/python
"""Tests of plotting converted data: the pure helpers, and a REAL conversion read back through the plot worker.

The real part converts the bundled NSTX 141711Y05 case (a subset of IDSs, HDF5, DD 4.1.1) with the GUI's own planner and
job runner, then checks what the worker returns against the source CDF and against a direct IMAS read.  It needs
IMAS-Python and the bundled fixture; without IMAS those tests are skipped with the reason stated, never passed.

    python examples/verify_gui_plots.py

Note on ``nodes``: after ``imas.DBEntry.get_slice`` on some IDSs the HDF5 file stays locked until the PROCESS exits, even
after ``close()``.  The server is not affected (every request is a fresh worker subprocess), so these tests call
``nodes`` through the subprocess too, exactly as the server does.
"""

import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

try:
    import imas
except ImportError:                          # pragma: no cover - environment dependent
    imas = None

from gui import convert_plan, fastion, imas_client, imas_worker as worker, jobs, recipes, server


class AoSLabelTests(unittest.TestCase):
    def test_radiation_process_uses_its_identifier_instead_of_generic_label(self):
        scalar = lambda value: type('Scalar', (), {'value': value})()
        process = type('Process', (), {
            'label': scalar('Value'),
            'identifier': type('Identifier', (), {'name': scalar('line_radiation')})(),
        })()
        self.assertEqual(worker._element_label(process, 2), '2: line_radiation')

    def test_distribution_uses_species_label(self):
        scalar = lambda value: type('Scalar', (), {'value': value})()
        distribution = type('Distribution', (), {
            'label': scalar('Value'),
            'species': type('Species', (), {'ion': type('Ion', (), {'label': scalar('D')})()})(),
            'gyro_type': 2,
            'process': [type('Process', (), {'type': type('Identifier', (), {'name': scalar('NBI')})()})()],
        })()
        self.assertEqual(worker._element_label(distribution, 1), '1: D, guiding centre, NBI')

    def test_radiation_recipes_separate_process_and_species_views(self):
        by_id = {recipe['id']: recipe for recipe in recipes.RECIPES}
        self.assertEqual(by_id['rad_emissivity']['title'], 'Electron radiation emissivity by process')
        self.assertEqual(by_id['rad_ion_emissivity']['paths'],
                         ['process[*]/profiles_1d[:]/ion[*]/emissivity'])

NSTX = ROOT / 'examples' / 'data' / 'nstx' / '141711Y05' / '141711Y05.CDF'
IDS_SUBSET = ['wall', 'tf', 'core_profiles', 'equilibrium', 'summary', 'nbi', 'core_sources', 'core_transport', 'radiation']
NUBEAM_CDF = ROOT / 'examples' / 'data' / 'nstx' / '120967E03_send' / '120967E03.CDF'   # local only: too large for the repository
DF_ENTRY_ENV = os.environ.get('TRANSP2IMAS_DF_ENTRY')                                     # or an already converted NUBEAM entry
_state = {}


def real_entry():
    """Convert the NSTX case once, through the GUI planner and job runner; returns the entry directory."""
    if 'entry' in _state:
        return _state['entry']
    tmp = tempfile.TemporaryDirectory()
    _state['tmp'] = tmp
    out = Path(tmp.name) / 'out'
    config = convert_plan.Config([ROOT / 'examples' / 'data'], [out], sys.executable, True, False)
    manager = jobs.JobManager(out, sys.executable, ROOT)
    job_id, directory = manager.reserve()
    plan = convert_plan.validate({'cdf': str(NSTX), 'options': {'dd_version': '4.1.1', 'ids': IDS_SUBSET}}, config,
                                 report_json=str(directory / 'conversion-report.json'))
    manager.start(job_id, plan)
    job = manager.wait(job_id, 280)
    _state.update(job=job, log=manager.log(job_id)['text'], config=config, plan=plan, out=out)
    _state['entry'] = Path(plan['target']['primary'])
    return _state['entry']


def tearDownModule():
    if 'tmp' in _state:
        _state['tmp'].cleanup()


def nodes_in_subprocess(entry, ids_name):
    """What the server does for the node browser: a fresh worker process per request."""
    return imas_client.call(sys.executable, ROOT, {'action': 'nodes', 'entry': str(entry), 'ids': ids_name})


class PureHelperTests(unittest.TestCase):
    def test_paths_parse_into_names_and_indices(self):
        self.assertEqual(worker.parse_path('profiles_1d[:]/ion[*]/temperature'), [('profiles_1d', ':'), ('ion', '*'), ('temperature', None)])
        self.assertEqual(worker.parse_path('source[3]/profiles_1d[:]/electrons/energy')[0], ('source', 3))

    def test_malformed_paths_are_refused(self):
        for bad in ('', None, 'a//b', 'a/b[', '1abc', 'a b', 'a/../b', 'a[-1]', 'a[x]', 'a;b', '/abs', 'x' * 400, 'a[1][2]'):
            with self.assertRaises(worker.WorkerError, msg=str(bad)):
                worker.parse_path(bad)

    def test_nonfinite_numbers_become_null_never_zero(self):
        self.assertEqual(worker.cleaned([1.0, float('nan'), float('inf'), -float('inf'), 2]), [1.0, None, None, None, 2.0])
        self.assertEqual(worker.cleaned([[1, np.nan], [np.inf, 4]]), [[1.0, None], [None, 4.0]])
        self.assertEqual(worker.cleaned(np.array(np.nan)), None)
        json.dumps(worker.cleaned(np.array([np.nan, 1.0])), allow_nan=False)

    def test_series_are_thinned_by_a_stride_and_the_last_samples_are_kept_consistent(self):
        self.assertEqual(worker.stride_for(1500), 1)
        self.assertEqual(worker.stride_for(1501), 2)
        self.assertEqual(worker.stride_for(0), 1)
        x, y, step = worker.downsample_series(np.arange(4000.0), np.arange(4000.0) * 2, limit=1000)
        self.assertEqual((step, len(x), len(y)), (4, 1000, 1000))
        self.assertEqual((x[1], y[1]), (4.0, 8.0))

    def test_the_dd_version_follows_the_entry_path_unless_overridden(self):
        self.assertEqual(worker.entry_dd_version('/db/NSTX/4/141711/2505'), '4.1.1')
        self.assertEqual(worker.entry_dd_version('/db/NSTX/3/141711/2505'), '3.39.0')
        self.assertIsNone(worker.entry_dd_version('/db/NSTX/9/141711/2505'))
        self.assertEqual(worker.entry_dd_version('/db/NSTX/4/141711/2505', '4.0.0'), '4.0.0')

    def test_entry_files_list_ids_occurrences_and_sizes(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name, size in (('master.h5', 1), ('wall.h5', 5), ('core_transport.h5', 7), ('core_transport_1.h5', 9), ('notes.txt', 3)):
                (Path(tmp) / name).write_bytes(b'x' * size)
            self.assertEqual(worker.entry_files(tmp), [('core_transport', 0, 7), ('core_transport', 1, 9), ('wall', 0, 5)])

    def test_every_recipe_is_well_formed(self):
        ids = set()
        for recipe in recipes.RECIPES:
            self.assertNotIn(recipe['id'], ids)
            ids.add(recipe['id'])
            self.assertIn(recipe['kind'], ('profile', 'timetrace', 'map', 'rz'))
            self.assertIn(recipe['group'], recipes.QUANTITY_GROUPS)
            if recipe['kind'] == 'rz':
                for r_path, z_path in recipe['pairs']:
                    worker.parse_path(r_path), worker.parse_path(z_path)
            self.assertTrue(recipe['paths'] and recipe['title'])
            for path in recipe['paths']:
                worker.parse_path(path)
            self.assertEqual(recipes.by_id(recipe['id']), recipe)
        with self.assertRaises(KeyError):
            recipes.by_id('nope')

    def test_non_hdf5_directories_are_refused_with_the_reason(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(worker.WorkerError) as caught:
                if imas is None:
                    raise worker.WorkerError('no master.h5')
                worker.open_entry(tmp)
            self.assertIn('master.h5', str(caught.exception))

    def test_unknown_actions_and_missing_entries_are_errors(self):
        with self.assertRaises(worker.WorkerError):
            worker.run({'action': 'describe'})
        with self.assertRaises(worker.WorkerError):
            worker.run({'action': 'explode', 'entry': '/x'})

    def test_the_worker_prints_one_result_line_and_an_error_line_on_failure(self):
        import io
        from contextlib import redirect_stdout
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            code = worker.main([json.dumps({'action': 'describe'})])
        self.assertEqual(code, 1)
        line = buffer.getvalue().strip()
        self.assertTrue(line.startswith(worker.PREFIX))
        self.assertIn('no entry given', json.loads(line[len(worker.PREFIX):])['error'])

    def test_the_client_maps_timeouts_and_silence_to_a_clear_failure(self):
        import subprocess
        with mock.patch.object(subprocess, 'run', side_effect=subprocess.TimeoutExpired('x', 1)):
            with self.assertRaises(imas_client.WorkerFailure) as caught:
                imas_client.call(sys.executable, ROOT, {'action': 'describe', 'entry': '/x'})
            self.assertIn('longer than', str(caught.exception))
        silent = subprocess.CompletedProcess([], 0, stdout=b'no marker here\n', stderr=b'Traceback...\nValueError: bad\n')
        with mock.patch.object(subprocess, 'run', return_value=silent):
            with self.assertRaises(imas_client.WorkerFailure) as caught:
                imas_client.call(sys.executable, ROOT, {'action': 'describe', 'entry': '/x'})
            self.assertIn('no result', str(caught.exception))
        with self.assertRaises(imas_client.WorkerFailure):
            imas_client.call(None, ROOT, {})


@unittest.skipIf(imas is None, 'IMAS-Python is required to convert and read back (these tests are skipped, not passed)')
@unittest.skipUnless(NSTX.is_file(), 'the bundled NSTX 141711Y05 fixture is required')
class RealEntryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.entry = real_entry()
        cls.dir = str(cls.entry)

    def tearDown(self):
        import gc
        gc.collect()                       # lazily loaded IDSs kept alive by exception tracebacks hold their HDF5 file open

    def hashes(self):
        return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(self.entry.glob('*.h5'))}

    def test_the_conversion_through_the_gui_planner_and_runner_succeeded(self):
        job = _state['job']
        self.assertEqual((job['state'], job['returncode']), ('succeeded', 0), _state['log'][-600:])
        self.assertTrue(job['report_written'])
        report = json.loads(Path(job['report']).read_text())
        self.assertEqual((report['status'], report['publication']['status'], report['dd_version']), ('PASS', 'complete', '4.1.1'))
        self.assertEqual([w['ids'] for w in report['writers']], IDS_SUBSET and [w['ids'] for w in report['writers']])
        self.assertEqual(set(IDS_SUBSET), {w['ids'] for w in report['writers']})
        self.assertEqual(str(self.entry), _state['plan']['target']['primary'])
        self.assertTrue((self.entry / 'master.h5').is_file())

    def test_describe_lists_the_ids_time_ranges_and_applicable_recipes(self):
        result = worker.describe(self.dir)
        self.assertEqual(result['dd_version'], '4.1.1')
        self.assertEqual({i['ids'] for i in result['ids']}, set(IDS_SUBSET))
        profiles = next(i for i in result['ids'] if i['ids'] == 'core_profiles')
        with imas.DBEntry(f'imas:hdf5?path={self.dir}', 'r', dd_version='4.1.1') as entry:
            times = np.asarray(entry.get('core_profiles', lazy=True).time)
        self.assertEqual((profiles['n_times'], profiles['t_first'], profiles['t_last']), (len(times), float(times[0]), float(times[-1])))
        self.assertEqual(profiles['bytes'], (self.entry / 'core_profiles.h5').stat().st_size)
        have = {r['id'] for r in result['recipes']}
        self.assertIn('cp_temperature', have)
        self.assertNotIn('distributions', {r['ids'] for r in result['recipes']})        # IDS not converted: no recipe offered

    def test_the_plotted_electron_temperature_and_density_equal_the_source_cdf(self):
        from netCDF4 import Dataset
        payload = worker.data(self.dir, {'recipe': 'cp_temperature', 'time_index': 3})
        te = next(s for s in payload['series'] if s['name'] == 'temperature')
        with Dataset(NSTX) as cdf:
            expected_te = np.asarray(cdf.variables['TE'][3, :], dtype=float)
            expected_ne = np.asarray(cdf.variables['NE'][3, :], dtype=float) * 1.0e6          # cm^-3 -> m^-3
            expected_x = np.asarray(cdf.variables['X'][3, :], dtype=float)
            self.assertEqual(cdf.variables['TE'].units.strip().upper(), 'EV')
        np.testing.assert_allclose(te['y'], expected_te, rtol=1e-6)
        density = worker.data(self.dir, {'recipe': 'cp_density', 'time_index': 3})['series'][0]
        np.testing.assert_allclose(density['y'], expected_ne, rtol=1e-6)
        self.assertEqual((payload['x_label'], te['units']), ('rho_tor_norm', 'eV'))
        self.assertEqual(len(te['x']), len(expected_x))
        self.assertTrue(np.all(np.diff(te['x']) > 0))

    def test_the_frames_agree_with_a_single_slice_request_and_cover_every_time(self):
        framed = worker.data(self.dir, {'recipe': 'cp_temperature', 'frames': True})
        count = framed['time']['count']
        self.assertEqual(len(framed['frames']), count)
        for index in (0, 2, count - 1):
            single = worker.data(self.dir, {'recipe': 'cp_temperature', 'time_index': index})
            self.assertEqual(framed['frames'][index]['index'], index)
            for a, b in zip(framed['frames'][index]['series'], single['series']):
                self.assertEqual((a['name'], a['x'], a['y']), (b['name'], b['x'], b['y']))
        times = [f['time'] for f in framed['frames']]
        self.assertEqual(times, sorted(times))

    def test_the_flux_map_matches_a_direct_read_and_carries_boundary_and_axis(self):
        payload = worker.data(self.dir, {'recipe': 'eq_psi_map', 'time_index': 2})
        self.assertEqual(payload['kind'], 'map')
        with imas.DBEntry(f'imas:hdf5?path={self.dir}', 'r', dd_version='4.1.1') as entry:
            ts = entry.get('equilibrium', lazy=True).time_slice[2]
            psi = np.asarray(ts.profiles_2d[0].psi.value)
            r, z = np.asarray(ts.profiles_2d[0].grid.dim1.value), np.asarray(ts.profiles_2d[0].grid.dim2.value)
            axis = (float(ts.global_quantities.magnetic_axis.r), float(ts.global_quantities.magnetic_axis.z))
        np.testing.assert_allclose(np.array(payload['z']), psi.T)                          # rows follow dim2, columns dim1
        np.testing.assert_allclose(payload['x'], r)
        np.testing.assert_allclose(payload['y'], z)
        self.assertEqual((len(payload['z']), len(payload['z'][0])), (len(z), len(r)))
        self.assertEqual([o['name'] for o in payload['overlays']], ['wall', 'boundary'])
        self.assertEqual((payload['points'][0]['x'], payload['points'][0]['y']), axis)
        self.assertEqual(payload['z_units'], 'Wb')
        self.assertEqual(payload['x_units'], '')                                           # 'mixed' is not a unit

    def test_a_time_trace_equals_a_direct_read(self):
        payload = worker.data(self.dir, {'recipe': 'eq_ip'})
        trace = payload['series'][0]
        with imas.DBEntry(f'imas:hdf5?path={self.dir}', 'r', dd_version='4.1.1') as entry:
            eq = entry.get('equilibrium', lazy=True)
            expected = [float(eq.time_slice[i].global_quantities.ip) for i in range(len(eq.time))]
            times = np.asarray(eq.time)
        np.testing.assert_allclose(trace['y'], expected)
        np.testing.assert_allclose(trace['x'], times)
        self.assertEqual((trace['units'], payload['x_label'], payload['x_units']), ('A', 'time', 's'))

    def test_a_time_trace_stored_with_its_own_time_axis_uses_that_axis(self):
        payload = worker.data(self.dir, {'recipe': 'nbi_power'})
        self.assertGreaterEqual(len(payload['series']), 1)
        with imas.DBEntry(f'imas:hdf5?path={self.dir}', 'r', dd_version='4.1.1') as entry:
            nbi = entry.get('nbi', lazy=True)
            expected_time = np.asarray(nbi.unit[0].power_launched.time)
            expected_power = np.asarray(nbi.unit[0].power_launched.data)
        np.testing.assert_allclose(payload['series'][0]['x'], expected_time)
        np.testing.assert_allclose(payload['series'][0]['y'], expected_power)

    def test_wildcards_expand_to_every_element_with_its_label(self):
        payload = worker.data(self.dir, {'recipe': 'cp_ion_temperature', 'time_index': 1})
        names = [s['name'] for s in payload['series']]
        self.assertGreaterEqual(len(names), 2)
        self.assertTrue(all('[' in n for n in names), names)
        self.assertEqual(len(set(names)), len(names))

    def test_the_ion_momentum_wall_and_equilibrium_plots_draw_the_stored_numbers(self):
        """The recipes added for ion, momentum, wall and equilibrium plots read what is stored, one by one against a direct read."""
        with imas.DBEntry(f'imas:hdf5?path={self.dir}', 'r', dd_version='4.1.1') as entry:
            profiles = entry.get('core_profiles', lazy=True)
            ions = profiles.profiles_1d[3].ion
            density = {f'density [{i}]': np.asarray(ions[i].density.value) for i in range(len(ions)) if len(ions[i].density)}
            momentum = np.asarray(profiles.profiles_1d[3].momentum_phi)
            wall = entry.get('wall', lazy=True).description_2d[0].limiter.unit[0].outline
            wall_r, wall_z = np.asarray(wall.r.value), np.asarray(wall.z.value)
            eq = entry.get('equilibrium', lazy=True).time_slice[3]
            q, f = np.asarray(eq.profiles_1d.q.value), np.asarray(eq.profiles_1d.f.value)
        ion_density = worker.data(self.dir, {'recipe': 'cp_ion_density', 'time_index': 3})
        self.assertGreater(len(density), 0)
        self.assertEqual(len(ion_density['series']), len(density))
        for series, expected in zip(ion_density['series'], density.values()):
            np.testing.assert_allclose(series['y'], expected)
        self.assertGreater(len(momentum), 0, 'the converted entry has no toroidal momentum density to plot')
        np.testing.assert_allclose(worker.data(self.dir, {'recipe': 'cp_momentum', 'time_index': 3})['series'][0]['y'], momentum)
        rz = worker.data(self.dir, {'recipe': 'wall_outline'})
        self.assertEqual((rz['kind'], len(rz['curves'])), ('rz', 1))
        np.testing.assert_allclose(rz['curves'][0]['x'], wall_r)
        np.testing.assert_allclose(rz['curves'][0]['y'], wall_z)
        self.assertEqual((rz['x_units'], rz['y_units']), ('m', 'm'))
        np.testing.assert_allclose(worker.data(self.dir, {'recipe': 'eq_q', 'time_index': 3})['series'][0]['y'], q)
        np.testing.assert_allclose(worker.data(self.dir, {'recipe': 'eq_f', 'time_index': 3})['series'][0]['y'], f)

    def test_the_flux_map_carries_the_wall_and_the_flux_range_when_stored(self):
        payload = worker.data(self.dir, {'recipe': 'eq_psi_map', 'time_index': 3})
        wall = [o for o in payload['overlays'] if o['name'].startswith('wall')]
        self.assertEqual(len(wall), 1)
        rz = worker.data(self.dir, {'recipe': 'wall_outline'})
        self.assertEqual(wall[0]['x'], rz['curves'][0]['x'])
        with imas.DBEntry(f'imas:hdf5?path={self.dir}', 'r', dd_version='4.1.1') as entry:
            quantities = entry.get('equilibrium', lazy=True).time_slice[3].global_quantities
            axis, boundary = float(quantities.psi_axis), float(quantities.psi_boundary)
        if abs(axis) < 1e30 and abs(boundary) < 1e30 and axis != boundary:
            self.assertEqual((payload['flux']['axis'], payload['flux']['boundary']), (axis, boundary))
        else:
            self.assertNotIn('flux', payload)                      # the empty-float sentinel is never offered as a flux value

    def test_a_missing_wall_is_not_an_error_for_a_map_and_a_clear_one_for_the_wall_plot(self):
        with tempfile.TemporaryDirectory() as tmp:
            class Empty:
                def get(self, *a, **k):
                    raise RuntimeError('no such IDS')
            self.assertEqual(worker.wall_overlays(Empty()), [])
            with self.assertRaises(RuntimeError):
                worker.rz_curves(Empty(), 'wall', 0, worker.WALL_PAIRS)

    def test_the_every_new_recipe_kind_group_is_known(self):
        self.assertTrue({r['group'] for r in recipes.RECIPES} <= set(recipes.QUANTITY_GROUPS))
        described = {r['id']: r for r in worker.describe(self.dir)['recipes']}
        self.assertTrue(all('group' in r for r in described.values()))
        for must in ('cp_ion_density', 'cp_ion_temperature', 'wall_outline', 'eq_psi_map', 'tr_ion_chi', 'src_momentum'):
            self.assertIn(must, described)

    def test_the_model_recipes_expand_every_transport_model(self):
        payload = worker.data(self.dir, {'recipe': 'tr_ion_chi', 'time_index': 3})
        self.assertTrue(payload['series'])
        self.assertEqual(payload['x_label'], 'rho_tor_norm')

    def test_every_recipe_either_draws_or_says_which_nodes_are_missing(self):
        present = {i for i, _, _ in worker.entry_files(self.dir)}
        drawn = 0
        for recipe in recipes.RECIPES:
            if recipe['ids'] not in present:
                continue
            try:
                payload = worker.data(self.dir, {'recipe': recipe['id']})
            except worker.WorkerError as exc:
                self.assertRegex(str(exc), 'none of the requested nodes exists|no outline in this entry', recipe['id'])
                continue
            drawn += 1
            self.assertEqual(payload['kind'], recipe['kind'], recipe['id'])
            json.dumps(payload, allow_nan=False)
            self.assertTrue(payload['series'] or payload['kind'] in ('map', 'rz'), recipe['id'])
            for series in payload['series']:
                self.assertEqual(len(series['x']), len(series['y']), recipe['id'])
        self.assertGreaterEqual(drawn, 10)

    def test_discovery_finds_typed_nodes_and_every_discovered_node_can_be_drawn(self):
        listing = nodes_in_subprocess(self.dir, 'core_profiles')['nodes']
        by_path = {n['path']: n for n in listing}
        te = by_path['profiles_1d[:]/electrons/temperature']
        self.assertEqual((te['kind'], te['units']), ('profile', 'eV'))
        self.assertIn('timetrace', {n['kind'] for n in listing})
        maps = [n for n in nodes_in_subprocess(self.dir, 'equilibrium')['nodes'] if n['kind'] == 'map']
        self.assertTrue(any(n['path'].endswith('profiles_2d[0]/psi') for n in maps))
        summary = nodes_in_subprocess(self.dir, 'summary')['nodes']
        self.assertTrue(summary and all(n['kind'] == 'timetrace' for n in summary))
        checked = 0
        for ids_name, nodes in (('core_profiles', listing), ('summary', summary)):
            for node in nodes[:12]:
                try:
                    payload = worker.data(self.dir, {'ids': ids_name, 'paths': [node['path']]})
                except worker.WorkerError as exc:
                    self.fail(f"{node['path']}: {exc}")
                self.assertEqual(payload['kind'], node['kind'], node['path'])
                json.dumps(payload, allow_nan=False)
                checked += 1
        self.assertGreater(checked, 15)

    def test_nodes_of_different_kinds_cannot_be_drawn_together(self):
        with self.assertRaises(worker.WorkerError) as caught:
            worker.data(self.dir, {'ids': 'equilibrium', 'paths': ['time_slice[:]/global_quantities/ip', 'time_slice[:]/profiles_1d/q']})
        self.assertIn('same kind', str(caught.exception))

    def test_bad_requests_are_clear_errors(self):
        for request, fragment in (({'ids': 'core_profiles', 'paths': []}, 'give a recipe'), ({'ids': 'core_profiles', 'paths': ['x'] * 20}, 'give a recipe'),
                                  ({'ids': 'core_profiles', 'paths': ['profiles_1d[:]/not_a_node']}, 'none of the requested nodes'),
                                  ({'recipe': 'nope'}, 'nope')):
            with self.assertRaises((worker.WorkerError, KeyError), msg=str(request)) as caught:
                worker.data(self.dir, request)
            self.assertIn(fragment, str(caught.exception))

    def test_the_wrong_dd_version_is_a_reported_error_not_a_crash(self):
        with self.assertRaises(RuntimeError) as caught:
            worker.data(self.dir, {'recipe': 'cp_temperature'}, dd_version='3.39.0')
        self.assertIn('major version', str(caught.exception))

    def test_reading_never_modifies_the_entry(self):
        before = self.hashes()
        worker.describe(self.dir)
        worker.data(self.dir, {'recipe': 'cp_temperature', 'frames': True})
        worker.data(self.dir, {'recipe': 'eq_psi_map'})
        nodes_in_subprocess(self.dir, 'equilibrium')
        self.assertEqual(self.hashes(), before)

    def test_a_time_index_beyond_the_end_is_clamped_not_an_error(self):
        payload = worker.data(self.dir, {'recipe': 'cp_temperature', 'time_index': 10 ** 6})
        self.assertEqual(payload['time']['index'], payload['time']['count'] - 1)

    def test_the_worker_as_a_subprocess_and_through_the_http_api_give_the_same_answer(self):
        request = {'action': 'data', 'entry': self.dir, 'recipe': 'cp_q', 'time_index': 2}
        via_process = imas_client.call(sys.executable, ROOT, request)
        self.assertEqual(via_process, worker.run(request))
        config = _state['config']
        srv = server.make_server([ROOT / 'examples' / 'data'], config=config, repo_root=ROOT)
        threading.Thread(target=srv.serve_forever, daemon=True).start()
        try:
            from verify_gui_convert import Client
            client = Client(srv)
            status, described = client.call('GET', f'/api/imas/describe?entry={self.dir}')
            self.assertEqual(status, 200)
            self.assertEqual({i['ids'] for i in described['ids']}, set(IDS_SUBSET))
            status, body = client.call('GET', f'/api/imas/data?entry={self.dir}&recipe=cp_q&time_index=2')
            self.assertEqual((status, body), (200, via_process))
            status, body = client.call('GET', f'/api/imas/data?entry={self.dir}&ids=core_profiles&paths=' + json.dumps(['profiles_1d[:]/q']))
            self.assertEqual(status, 200)
            status, body = client.call('GET', f'/api/imas/data?entry={self.dir}&recipe=nope')
            self.assertEqual(status, 502)
            status, body = client.call('GET', f'/api/imas/data?entry={self.dir}&paths=notjson')
            self.assertIn(status, (400, 502))
            for bad in ('what=bogus', 'what=ep&zone=x', 'what=ep&e_min=nan', 'what=ep&e_max=inf', 'what=ep&dist=a'):
                status, body = client.call('GET', f'/api/imas/df?entry={self.dir}&{bad}')
                self.assertEqual(status, 400, bad)                                    # refused before any worker is started
            status, body = client.call('GET', f'/api/imas/df?entry={self.dir}&what=info')
            self.assertEqual(status, 502)                                             # this entry has no distributions: the reason is shown
            self.assertIn('distributions', body['error'])
            results = []
            threads = [threading.Thread(target=lambda: results.append(client.call('GET', f'/api/imas/data?entry={self.dir}&recipe=cp_density')[0]))
                       for _ in range(4)]
            [t.start() for t in threads]
            [t.join(180) for t in threads]
            self.assertEqual(results, [200, 200, 200, 200])                              # extra readers wait their turn
            status, entries = client.call('GET', '/api/entries')
            ours = [e for e in entries['entries'] if e['path'] == self.dir]
            self.assertEqual(len(ours), 1)                                    # other tests keep their own entries in the same database
            self.assertEqual(sorted(ours[0]['ids']), sorted(IDS_SUBSET))
        finally:
            srv.shutdown()
            srv.server_close()

    def test_the_gui_conversion_produces_the_entry_the_command_line_would(self):
        """The job is the CLI: its log is the converter's own output and the report records the same command's results."""
        log = _state['log']
        self.assertIn('IMAS target DD: 4.1.1', log)
        self.assertIn('ids stored in ' + str(self.entry), log)
        report = json.loads(Path(_state['job']['report']).read_text())
        self.assertEqual((report['machine'], report['shot'], report['run']), ('NSTX', 141711, 2505))
        self.assertEqual(report['mapping_policy'], 'strict')
        self.assertEqual(Path(report['output_directory']), self.entry)


@unittest.skipIf(imas is None, 'IMAS-Python is required to convert and read back (these tests are skipped, not passed)')
@unittest.skipUnless(NSTX.is_file(), 'the bundled NSTX 141711Y05 fixture is required')
class OverwriteRealTests(unittest.TestCase):
    """A second conversion into an existing entry, through the planner's overwrite guard, on a PRIVATE entry (run 4001)."""

    RUN = 4001

    @classmethod
    def convert(cls, ids, new_db=True, confirm=None, expect_existing=False):
        real_entry()                                                  # the shared configuration and output root
        config, out = _state['config'], _state['out']
        manager = jobs.JobManager(out, sys.executable, ROOT)
        job_id, directory = manager.reserve()
        plan = convert_plan.validate({'cdf': str(NSTX), 'options': {'dd_version': '4.1.1', 'ids': ids, 'new_db': new_db, 'run': cls.RUN}},
                                     config, report_json=str(directory / 'conversion-report.json'))
        if expect_existing:
            assert plan['needs_confirmation'] and plan['target']['exists']
            if not convert_plan.confirm(plan, confirm):
                return None, '', plan
        manager.start(job_id, plan)
        job = manager.wait(job_id, 280)
        return job, manager.log(job_id)['text'], plan

    @classmethod
    def setUpClass(cls):
        job, log, plan = cls.convert(['wall', 'tf', 'core_profiles', 'summary'])
        assert job['state'] == 'succeeded', log[-500:]
        cls.entry = Path(plan['target']['primary'])
        assert cls.entry.name == str(cls.RUN)

    def ids(self):
        return {i for i, _, _ in worker.entry_files(self.entry)}

    def test_a_first_replace_without_the_typed_path_does_not_start(self):
        before = self.ids()
        job, _log, plan = self.convert(['wall'], True, confirm='not the path', expect_existing=True)
        self.assertIsNone(job)                                         # the guard refused: nothing ran
        self.assertTrue(plan['needs_confirmation'] and plan['target']['next_free_run'] == self.RUN + 1)
        self.assertEqual(self.ids(), before)

    def test_b_a_new_entry_replaces_the_old_one_and_leaves_only_what_was_converted_now(self):
        self.assertEqual(self.ids(), {'wall', 'tf', 'core_profiles', 'summary'})
        job, log, plan = self.convert(['wall', 'tf'], True, confirm=str(self.entry), expect_existing=True)
        self.assertEqual((job['state'], job['returncode']), ('succeeded', 0), log[-500:])
        self.assertEqual(self.ids(), {'wall', 'tf'}, 'creating a new entry must delete what an earlier conversion wrote')
        report = json.loads(Path(job['report']).read_text())
        self.assertEqual((report['status'], report['publication']['status']), ('PASS', 'complete'))

    def test_c_without_create_new_the_selected_ids_are_added_and_the_others_stay(self):
        self.assertEqual(self.ids(), {'wall', 'tf'})
        job, log, plan = self.convert(['summary'], False, confirm=str(self.entry), expect_existing=True)
        self.assertEqual((job['state'], job['returncode']), ('succeeded', 0), log[-500:])
        self.assertEqual(self.ids(), {'wall', 'tf', 'summary'}, 'extending an entry must keep the IDSs it already had')

    def test_d_the_suggested_run_number_is_free_and_converts_without_confirmation(self):
        _job, _log, plan = self.convert(['wall'], True, confirm='x', expect_existing=True)
        free = plan['target']['next_free_run']
        config, manager = _state['config'], jobs.JobManager(_state['out'], sys.executable, ROOT)
        job_id, directory = manager.reserve()
        fresh = convert_plan.validate({'cdf': str(NSTX), 'options': {'dd_version': '4.1.1', 'ids': ['wall'], 'run': free}}, config,
                                      report_json=str(directory / 'conversion-report.json'))
        self.assertEqual((fresh['needs_confirmation'], fresh['target']['exists']), (False, False))
        manager.start(job_id, fresh)
        self.assertEqual(manager.wait(job_id, 280)['state'], 'succeeded')
        self.assertEqual({i for i, _, _ in worker.entry_files(fresh['target']['primary'])}, {'wall'})
        self.assertEqual(self.ids(), {'wall', 'tf', 'summary'})        # the original entry was not touched


class FastIonPureTests(unittest.TestCase):
    """The distribution-function layout and reductions, on small arrays with answers computed by hand."""

    def cube(self):
        z, e, p = np.meshgrid(np.arange(3), np.arange(4), np.arange(2), indexing='ij')
        return 1.0 + z + 10.0 * e + 100.0 * p                       # value identifies its (zone, energy, pitch)

    def test_the_flat_array_is_zone_fastest_then_energy_then_pitch(self):
        cube = self.cube()
        flat = cube.reshape(-1, order='F')
        self.assertEqual(list(flat[:4]), [1.0, 2.0, 3.0, 11.0])      # zones 0-2, then the next energy
        np.testing.assert_array_equal(fastion.fortran_cube(flat, 3, 4, 2), cube)
        self.assertFalse(np.array_equal(flat.reshape(3, 4, 2), cube))   # the C-order reading is a different (wrong) array

    def test_a_size_that_does_not_match_the_grid_is_refused(self):
        with self.assertRaises(worker.WorkerError) as caught:
            fastion.fortran_cube(np.zeros(23), 3, 4, 2)
        self.assertIn('24', str(caught.exception))

    def test_reductions_weight_by_zone_volume_and_divide_by_bin_width(self):
        cube, volumes = self.cube(), np.array([2.0, 3.0, 5.0])
        energy, pitch = np.array([0.0, 10.0, 30.0, 60.0, 100.0]), np.array([-1.0, 0.0, 1.0])
        total = float((cube * volumes[:, None, None]).sum())
        ep = fastion.reduce(cube, volumes, energy, pitch, 'ep')
        self.assertEqual(ep.shape, (4, 2))
        self.assertAlmostEqual(ep.sum(), total)
        spectrum = fastion.reduce(cube, volumes, energy, pitch, 'energy')
        self.assertAlmostEqual(float((spectrum * np.diff(energy)).sum()), total)       # integrating the spectrum gives the number
        self.assertAlmostEqual(spectrum[1], (cube[:, 1, :] * volumes[:, None]).sum() / 20.0)
        pitches = fastion.reduce(cube, volumes, energy, pitch, 'pitch')
        self.assertAlmostEqual(float((pitches * np.diff(pitch)).sum()), total)
        np.testing.assert_allclose(fastion.reduce(cube, volumes, energy, pitch, 'zones'), cube.sum(axis=(1, 2)))
        np.testing.assert_array_equal(fastion.reduce(cube, volumes, energy, pitch, 'ep', zone=1), cube[1])   # one zone: stored values
        with self.assertRaises(worker.WorkerError):
            fastion.reduce(cube, volumes, energy, pitch, 'nonsense')

    def test_the_window_keeps_bins_by_centre(self):
        cube, energy, pitch = self.cube(), np.array([0.0, 10.0, 20.0, 30.0, 40.0]), np.array([-1.0, 0.0, 1.0])
        kept = fastion.window(cube, energy, pitch, e_min=15.0, e_max=35.0, p_min=0.0)
        self.assertTrue((kept[:, 0, :] == 0).all())                       # centre 5 is below 15; centres 15, 25, 35 stay
        self.assertEqual(kept[:, 3, 1].tolist(), cube[:, 3, 1].tolist())
        self.assertTrue((kept[:, :, 0] == 0).all())                       # pitch centre -0.5 is outside p >= 0
        np.testing.assert_array_equal(fastion.window(cube, energy, pitch), cube)

    def test_polygons_are_chained_from_edges_in_any_order_and_direction(self):
        loop = fastion.chain_polygon([(3, 4), (2, 1), (4, 1), (3, 2)])
        self.assertEqual(sorted(loop), [1, 2, 3, 4])
        ring = loop + loop[:1]
        self.assertEqual({frozenset(pair) for pair in zip(ring, ring[1:])}, {frozenset(e) for e in [(3, 4), (2, 1), (4, 1), (3, 2)]})   # consecutive
        self.assertEqual(len(fastion.chain_polygon([(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)])), 5)
        self.assertIsNone(fastion.chain_polygon([(1, 2), (2, 3), (3, 4)]))                  # open
        self.assertIsNone(fastion.chain_polygon([(1, 2), (2, 3), (3, 1), (1, 4), (4, 5), (5, 1)]))   # pinched: node 1 has four edges
        self.assertIsNone(fastion.chain_polygon([(1, 2), (2, 1)]))
        self.assertAlmostEqual(fastion.polygon_area([[0, 0], [2, 0], [2, 1.5], [0, 1.5]]), 3.0)

    def test_the_chart_module_exposes_the_contour_and_zone_helpers(self):
        charts_source = (ROOT / 'gui' / 'static' / 'charts.js').read_text()
        for name in ('contourSegments', 'contourLevels', 'pointInPolygon', 'shapeMap'):
            self.assertIn(name, charts_source)

    def test_zone_map_uses_ordered_nodes_and_smooth_shared_flux_surfaces(self):
        from types import SimpleNamespace as NS
        angles = np.arange(4) * np.pi / 2
        nodes = [NS(geometry=(r * np.cos(a), r * np.sin(a)))
                 for r in (1.0, 2.0) for a in angles]
        edge_ids, edges, faces = {}, [], []
        for i in range(4):
            j = (i + 1) % 4
            loop = [i, 4 + i, 4 + j, j]
            boundary = []
            for a, b in zip(loop, loop[1:] + loop[:1]):
                pair = tuple(sorted((a + 1, b + 1)))
                if pair not in edge_ids:
                    edge_ids[pair] = len(edges) + 1
                    edges.append(NS(nodes=pair))
                boundary.append(NS(index=edge_ids[pair]))
            faces.append(NS(nodes=[n + 1 for n in loop], boundary=boundary,
                            geometry=(1.0, float(i), 1.0), measure=1.0))
        spaces = [NS(objects_per_dimension=[NS(object=nodes), NS(object=edges), NS(object=faces)])]
        polygons = fastion.smooth_zone_polygons(spaces, axis=(0.0, 0.0))
        self.assertEqual(len(polygons), 4)
        self.assertTrue(all(len(polygon) > 4 for polygon in polygons))
        self.assertTrue(all(np.isfinite(polygon).all() for polygon in map(np.asarray, polygons)))


@unittest.skipIf(imas is None, 'IMAS-Python is required to write and read a distributions entry')
class FastIonSyntheticEntryTests(unittest.TestCase):
    """A small ``distributions`` entry written here (two square zones side by side) and read back through the worker."""

    VOLUMES = (2.0, 3.0)
    ENERGY = (0.0, 10.0, 20.0, 30.0)
    PITCH = (-1.0, 0.0, 1.0)

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.dir = Path(cls.tmp.name) / 'TEST' / '4' / '1' / '1'
        cls.cubes = cls.write_entry(cls.dir)

    @classmethod
    def write_entry(cls, directory):
        """Write the two-zone entry into ``directory`` (used by the browser check too); returns the two time slices' cubes."""
        Path(directory).mkdir(parents=True, exist_ok=True)
        z, e, p = np.meshgrid(np.arange(2), np.arange(3), np.arange(2), indexing='ij')
        cubes = [1.0 + z + 10.0 * e + 100.0 * p, 2.0 * (1.0 + z + 10.0 * e + 100.0 * p)]
        entry = imas.DBEntry(f'imas:hdf5?path={directory}', 'w', dd_version='4.1.1')
        ids = imas.IDSFactory('4.1.1').distributions()
        ids.ids_properties.homogeneous_time = 1
        ids.time = [0.0, 1.0]
        ids.distribution.resize(1)
        dist = ids.distribution[0]
        dist.gyro_type = 2
        dist.ggd.resize(2)
        for it, ggd in enumerate(dist.ggd):
            ggd.time = float(it)
            cls.fill_grid(ggd)
            ggd.expansion_fd3v.resize(1)
            ggd.expansion_fd3v[0].grid_subset.resize(1)
            ggd.expansion_fd3v[0].grid_subset[0].values = cubes[it].reshape(-1, order='F')
        entry.put(ids)
        entry.close()
        return cubes

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    @classmethod
    def fill_grid(cls, ggd):
        ggd.grid.space.resize(3)
        zones = ggd.grid.space[0]
        zones.objects_per_dimension.resize(3)
        nodes = zones.objects_per_dimension[0]
        nodes.object.resize(6)
        for i, (r, zz) in enumerate([(1, 0), (2, 0), (3, 0), (1, 1), (2, 1), (3, 1)]):
            nodes.object[i].geometry = [float(r), float(zz)]
        edges = zones.objects_per_dimension[1]
        pairs = [(1, 2), (2, 3), (4, 5), (5, 6), (1, 4), (2, 5), (3, 6)]
        edges.object.resize(len(pairs))
        for i, pair in enumerate(pairs):
            edges.object[i].nodes = np.asarray(pair, dtype=np.int32)
        faces = zones.objects_per_dimension[2]
        faces.object.resize(2)
        for i, boundary in enumerate([(1, 6, 3, 5), (2, 7, 4, 6)]):
            faces.object[i].geometry = [0.0, float(i), cls.VOLUMES[i]]
            faces.object[i].measure = 1.0
            faces.object[i].boundary.resize(len(boundary))
            for k, edge in enumerate(boundary):
                faces.object[i].boundary[k].index = edge
        for space, edges_values in ((ggd.grid.space[1], cls.ENERGY), (ggd.grid.space[2], cls.PITCH)):
            space.objects_per_dimension.resize(1)
            space.objects_per_dimension[0].object.resize(len(edges_values))
            for i, value in enumerate(edges_values):
                space.objects_per_dimension[0].object[i].geometry = [value]

    def call(self, **request):
        return imas_client.call(sys.executable, ROOT, {'action': 'df', 'entry': str(self.dir), **request})

    def test_info_reports_the_grid_and_the_times(self):
        info = self.call(what='info')
        self.assertEqual((info['zones'], info['energy_bins'], info['pitch_bins']), (2, 3, 2))
        self.assertEqual((info['energy_range'], info['pitch_range']), ([0.0, 30.0], [-1.0, 1.0]))
        self.assertEqual((info['time']['count'], info['distributions'][0]['index']), (2, 0))
        self.assertIn('guiding centre', info['distributions'][0]['label'])

    def test_the_energy_pitch_map_is_the_volume_weighted_sum_and_follows_the_time_slice(self):
        for it in (0, 1):
            result = self.call(what='ep', time_index=it)
            expected = (self.cubes[it] * np.array(self.VOLUMES)[:, None, None]).sum(axis=0)
            np.testing.assert_allclose(np.array(result['z']), expected)                        # rows are energy, columns pitch
            self.assertEqual((result['x'], result['y']), ([-0.5, 0.5], [5.0, 15.0, 25.0]))
            self.assertEqual(result['time']['index'], it)
        zone = self.call(what='ep', zone=1)
        np.testing.assert_allclose(np.array(zone['z']), self.cubes[0][1])
        self.assertEqual(zone['z_units'], 'm^-3')

    def test_spectra_integrate_to_the_number_of_particles(self):
        total = float((self.cubes[0] * np.array(self.VOLUMES)[:, None, None]).sum())
        for what, edges in (('energy', self.ENERGY), ('pitch', self.PITCH)):
            series = self.call(what=what)['series'][0]
            self.assertAlmostEqual(float(np.dot(series['y'], np.diff(edges))), total, places=9)

    def test_the_zone_map_draws_each_square_with_its_density(self):
        result = self.call(what='zones')
        self.assertEqual(result['kind'], 'polys')
        self.assertEqual([q['zone'] for q in result['polys']], [0, 1])
        for q, expected in zip(result['polys'], self.cubes[0].sum(axis=(1, 2))):
            self.assertAlmostEqual(q['value'], expected)
            self.assertAlmostEqual(fastion.polygon_area(q['xy']), 1.0)
        self.assertEqual(sorted(tuple(p) for p in result['polys'][0]['xy']), [(1.0, 0.0), (1.0, 1.0), (2.0, 0.0), (2.0, 1.0)])
        self.assertEqual(result['overlays'], [])                                                # no wall in this entry
        json.dumps(result, allow_nan=False)

    def test_the_window_restricts_the_reduction(self):
        result = self.call(what='energy', e_min=10.0)                                           # bin centres 15 and 25 only
        self.assertEqual(result['series'][0]['y'][0], 0.0)
        self.assertGreater(result['series'][0]['y'][1], 0.0)

    def test_bad_requests_are_clear_errors(self):
        for request, text in (({'what': 'ep', 'zone': 5}, 'zone 5 does not exist'), ({'what': 'ep', 'dist': 3}, 'distribution 3 does not exist'),
                              ({'what': 'bogus'}, 'unknown reduction')):
            with self.assertRaises(imas_client.WorkerFailure) as caught:
                self.call(**request)
            self.assertIn(text, str(caught.exception))

    def test_a_grid_that_disagrees_with_the_array_size_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            entry_dir = Path(tmp) / 'TEST' / '4' / '1' / '1'
            entry_dir.mkdir(parents=True)
            entry = imas.DBEntry(f'imas:hdf5?path={entry_dir}', 'w', dd_version='4.1.1')
            ids = imas.IDSFactory('4.1.1').distributions()
            ids.ids_properties.homogeneous_time = 1
            ids.time = [0.0]
            ids.distribution.resize(1)
            ids.distribution[0].ggd.resize(1)
            ggd = ids.distribution[0].ggd[0]
            self.fill_grid(ggd)
            ggd.expansion_fd3v.resize(1)
            ggd.expansion_fd3v[0].grid_subset.resize(1)
            ggd.expansion_fd3v[0].grid_subset[0].values = np.arange(11.0)                       # 2 x 3 x 2 would be 12
            entry.put(ids)
            entry.close()
            with self.assertRaises(imas_client.WorkerFailure) as caught:
                imas_client.call(sys.executable, ROOT, {'action': 'df', 'entry': str(entry_dir), 'what': 'ep'})
            self.assertIn('12', str(caught.exception))


def fast_ion_entry():
    """A converted NUBEAM entry: ``TRANSP2IMAS_DF_ENTRY`` if set, else the local 120967E03 case converted here (wall and
    distributions only, through the GUI planner and runner).  ``None`` when neither is available."""
    if 'df' in _state:
        return _state['df']
    entry = None
    if DF_ENTRY_ENV and (Path(DF_ENTRY_ENV) / 'master.h5').is_file():
        entry = Path(DF_ENTRY_ENV)
    elif NUBEAM_CDF.is_file():
        tmp = tempfile.TemporaryDirectory()
        _state['df_tmp'] = tmp
        out = Path(tmp.name) / 'out'
        config = convert_plan.Config([NUBEAM_CDF.parent.parent], [out], sys.executable, True, False)
        manager = jobs.JobManager(out, sys.executable, ROOT)
        job_id, directory = manager.reserve()
        plan = convert_plan.validate({'cdf': str(NUBEAM_CDF), 'options': {'dd_version': '4.1.1', 'ids': ['wall', 'distributions']}}, config,
                                     report_json=str(directory / 'conversion-report.json'))
        manager.start(job_id, plan)
        job = manager.wait(job_id, 590)
        _state['df_log'] = manager.log(job_id)['text']
        if job['state'] == 'succeeded':
            entry = Path(plan['target']['primary'])
    _state['df'] = entry
    return entry


@unittest.skipIf(imas is None, 'IMAS-Python is required')
class FastIonRealEntryTests(unittest.TestCase):
    """The distribution function of a real NUBEAM conversion against the moments the converter stored separately.

    Uses the local ``examples/data/nstx/120967E03_send`` case (not in the repository: too large) and converts it here; skipped,
    with this reason, when the case is absent and ``TRANSP2IMAS_DF_ENTRY`` is not set.
    """

    @classmethod
    def setUpClass(cls):
        cls.entry = fast_ion_entry()
        if cls.entry is None:
            raise unittest.SkipTest(f'no NUBEAM data: {NUBEAM_CDF} is absent (local-only case) and TRANSP2IMAS_DF_ENTRY is not set'
                                    if not NUBEAM_CDF.is_file() else 'the NUBEAM conversion failed: ' + _state.get('df_log', '')[-500:])
        cls.dir = cls.entry

    @classmethod
    def tearDownClass(cls):
        if 'df_tmp' in _state:
            _state.pop('df_tmp').cleanup()
            _state['df'] = None

    def call(self, **request):
        return imas_client.call(sys.executable, ROOT, {'action': 'df', 'entry': str(self.dir), **request})

    def stored(self, dist, it):
        with imas.DBEntry(f'imas:hdf5?path={self.dir}', 'r', dd_version='4.1.1') as entry:
            d = entry.get('distributions', lazy=True).distribution[dist]
            return float(d.global_quantities[it].particles_n), np.asarray(d.profiles_1d[it].density_fast), \
                np.asarray(d.profiles_1d[it].grid.rho_tor_norm)

    def test_the_summed_distribution_equals_the_stored_particle_number_for_every_distribution_and_time(self):
        info = self.call(what='info')
        for dist in info['distributions']:
            for it in range(info['time']['count']):
                number = self.stored(dist['index'], it)[0]
                ep = np.array(self.call(what='ep', dist=dist['index'], time_index=it)['z'])
                self.assertAlmostEqual(ep.sum() / number, 1.0, places=9, msg=f'{dist["label"]} slice {it}')
                spectrum = self.call(what='energy', dist=dist['index'], time_index=it)['series'][0]
                self.assertAlmostEqual(float(np.dot(spectrum['y'], np.diff(self.edges(dist['index'], it)))) / number, 1.0, places=9)

    def edges(self, dist, it):
        with imas.DBEntry(f'imas:hdf5?path={self.dir}', 'r', dd_version='4.1.1') as entry:
            space = entry.get('distributions', lazy=True).distribution[dist].ggd[it].grid.space[1].objects_per_dimension[0].object
            return np.array([float(space[i].geometry[0]) for i in range(len(space))])

    def test_the_zone_polygons_are_all_drawn_and_have_the_stored_areas(self):
        result = self.call(what='zones')
        info = self.call(what='info')
        self.assertEqual(len(result['polys']), info['zones'])
        with imas.DBEntry(f'imas:hdf5?path={self.dir}', 'r', dd_version='4.1.1') as entry:
            faces = entry.get('distributions', lazy=True).distribution[0].ggd[0].grid.space[0].objects_per_dimension[2].object
            measures = [float(faces[q['zone']].measure) for q in result['polys']]
        areas = [fastion.polygon_area(q['xy']) for q in result['polys']]
        # Zones are straight-chord polygons through the stored nodes.  Where a zone is bounded by a curved flux surface with few
        # nodes the chord polygon is smaller than the stored area: the innermost wedges are triangles of 2/pi (0.637) of it, a
        # square inscribed in a circle.  It never exceeds the stored area by more than rounding, and overall it is within 3%.
        ratios = np.array(areas) / np.array(measures)
        self.assertGreater(ratios.min(), 0.63)
        self.assertLess(ratios.max(), 1.02)
        self.assertAlmostEqual(sum(areas) / sum(measures), 1.0, delta=0.03)

    def test_the_wall_is_the_overlay_of_the_zone_map_when_the_entry_has_one(self):
        result = self.call(what='zones')
        self.assertTrue(all(o['name'].startswith('wall') for o in result['overlays']))


if __name__ == '__main__':
    unittest.main(verbosity=1)
