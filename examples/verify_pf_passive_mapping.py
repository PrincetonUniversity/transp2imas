#!/usr/bin/python
"""Portable pf_passive mapping regression (no IMAS or netCDF4 needed).

Creates a clearly synthetic CDF and checks PS_n discovery, time-coordinate
integrity, units, current-sign translation, all supported geometry types,
DD3/DD4 polygon outlines, missing-input cases and refusal of unsafe mappings.
The early non-ISolver CDF fixtures have no PS_n variables. The ITER/NSTX
ISolver examples contain PS_n and are excluded from the no-current check.
"""

from __future__ import annotations

import importlib
import json
import logging
from pathlib import Path
import sys
import tempfile
import types
from unittest import mock
import unittest
from types import SimpleNamespace

import numpy as np
from scipy.io import netcdf_file

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pf_passive_mapping import (  # noqa: E402
    build_passive_plan, load_geometry_manifest, populate_passive_ids,
    read_passive_channels, _source_scale,
)


EXAMPLE_GEOMETRY = Path(__file__).with_name('pf_passive_geometry.synthetic.json')


class Node(SimpleNamespace):
    """Simple object whose members mimic IMAS geometry structures."""


class StructureArray(list):
    def __init__(self, factory):
        super().__init__()
        self.factory = factory

    def resize(self, count):
        self.clear()
        self.extend(self.factory() for _ in range(count))


def element_node():
    return Node(
        name='', description='', identifier='', turns_with_sign=None, area=None,
        geometry=Node(
            geometry_type=None,
            rectangle=Node(), outline=Node(), oblique=Node(),
            annulus=Node(), arcs_of_circle=Node(),
            thick_line=Node(first_point=Node(), second_point=Node()),
        ),
    )


def loop_node():
    return Node(name='', description='', element=StructureArray(element_node),
                resistance=None, resistivity=None, time=None, current=None)


def ids_node():
    return Node(loop=StructureArray(loop_node))


class TestPassiveMapping(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.cdf = Path(self.tempdir.name) / 'synthetic.CDF'
        with netcdf_file(self.cdf, 'w') as nc:
            nc.createDimension('TIME', 3)
            nc.createDimension('TIME3', 2)
            t = nc.createVariable('TIME', 'd', ('TIME',))
            t[:] = [0.1, 0.2, 0.3]
            t3 = nc.createVariable('TIME3', 'd', ('TIME3',))
            t3[:] = [0.1, 0.3]
            ps1 = nc.createVariable('PS_1', 'd', ('TIME',))
            ps1[:] = [100.0, 200.0, 300.0]
            ps1.units = 'A'
            ps1.long_name = 'Current in passive structure 1 <BVFB>'
            ps2 = nc.createVariable('PS_2', 'd', ('TIME3',))
            ps2[:] = [1.0, 2.0]
            ps2.units = 'kA'
            ps2.long_name = 'Current in passive structure 2 <LOWER_RING>'
            psm = nc.createVariable('PS_MISC', 'd', ('TIME',))
            psm[:] = [0.0, 0.0, 0.0]

    def _with_cdf(self, func):
        with netcdf_file(self.cdf, 'r', mmap=False) as nc:
            return func(nc)

    def _manifest(self, fn):
        data = json.loads(EXAMPLE_GEOMETRY.read_text())
        fn(data)
        path = Path(self.tempdir.name) / 'geometry.json'
        path.write_text(json.dumps(data), encoding='utf-8')
        return path

    def test_source_discovery_time_and_units(self):
        channels = self._with_cdf(read_passive_channels)
        self.assertEqual(list(channels), ['PS_1', 'PS_2'])
        self.assertEqual(channels['PS_1'].name, 'BVFB')
        np.testing.assert_array_equal(channels['PS_1'].current_ampere, [100, 200, 300])
        np.testing.assert_array_equal(channels['PS_1'].time, [0.1, 0.2, 0.3])
        np.testing.assert_array_equal(channels['PS_2'].current_ampere, [1000, 2000])
        np.testing.assert_array_equal(channels['PS_2'].time, [0.1, 0.3])

    def test_full_geometry_and_loop_time(self):
        def run(nc):
            plan = build_passive_plan(nc, machine='SYNTHETIC', geometry_path=EXAMPLE_GEOMETRY)
            self.assertEqual(len(plan.loops), 3)
            self.assertEqual(plan.geometry_source.startswith('SYNTHETIC'), True)
            ids = ids_node()
            self.assertEqual(populate_passive_ids(ids, plan, 4), 2)
            np.testing.assert_array_equal(ids.loop[0].current, [100, 200, 300])
            np.testing.assert_array_equal(ids.loop[1].current, [-1000, -2000])
            np.testing.assert_array_equal(ids.loop[1].time, [0.1, 0.3])
            self.assertIsNone(ids.loop[2].current)
            self.assertEqual(ids.loop[0].element[0].geometry.geometry_type, 2)
            self.assertEqual(ids.loop[2].element[0].geometry.geometry_type, 5)
            outline = ids.loop[1].element[0].geometry.outline
            self.assertEqual(len(outline.r), 5)  # DD4 closed
            self.assertEqual(outline.r[0], outline.r[-1])
            self.assertEqual(outline.z[0], outline.z[-1])
            ids3 = ids_node()
            populate_passive_ids(ids3, plan, 3)
            self.assertEqual(len(ids3.loop[1].element[0].geometry.outline.r), 4)  # DD3 open
        self._with_cdf(run)

    def test_current_only_fallback(self):
        def run(nc):
            plan = build_passive_plan(nc)
            self.assertEqual(len(plan.loops), 2)
            self.assertFalse(any(loop.elements for loop in plan.loops))
            ids = ids_node()
            self.assertEqual(populate_passive_ids(ids, plan, 4), 2)
            self.assertEqual(ids.loop[0].name, 'BVFB')
            self.assertEqual(len(ids.loop[0].element), 0)
        self._with_cdf(run)

    def test_full_writer_control_flow_without_installed_imas(self):
        """Exercise the actual writer with a fake IDS factory, not just helpers."""
        def run(nc):
            class MockImas:
                device = 'SYNTHETIC'
                dd_major_version = '4'
                logger = logging.getLogger('pfpassive-smoke')

                def __init__(self):
                    self.nc = nc
                    self.put_result = None

                def new_ids(self, name):
                    assert name == 'pf_passive'
                    return Node(
                        loop=StructureArray(loop_node),
                        ids_properties=Node(comment='CDF source info', homogeneous_time=None),
                        code=Node(parameters='<source_parameters />'),
                        time=None,
                    )

                def fill_properties(self, ids):
                    return ids

                def fill_code_properties(self, ids):
                    return ids

                def get_ncvar(self, name):
                    if name not in self.nc.variables:
                        return np.array([])
                    return np.asarray(self.nc.variables[name][:], dtype=np.float64)

                def putdb_imas(self, ids):
                    self.put_result = ids

            fake_nc2imas = types.ModuleType('nc2imas')
            fake_nc2imas.nc2imas = MockImas
            with mock.patch.dict(sys.modules, {'nc2imas': fake_nc2imas}):
                sys.modules.pop('imas_pf_passive', None)
                try:
                    writer = importlib.import_module('imas_pf_passive')
                    fake = MockImas()
                    writer.imas_pf_passive(fake, geometry_path=str(EXAMPLE_GEOMETRY))
                    ids = fake.put_result
                    self.assertIsNotNone(ids)
                    self.assertEqual(ids.ids_properties.homogeneous_time, 0)
                    np.testing.assert_array_equal(ids.time, [0.1, 0.2, 0.3])
                    np.testing.assert_array_equal(ids.loop[0].time, [0.1, 0.2, 0.3])
                    np.testing.assert_array_equal(ids.loop[1].time, [0.1, 0.3])
                    np.testing.assert_array_equal(ids.loop[1].current, [-1000.0, -2000.0])
                    self.assertIn('sha256=', ids.ids_properties.comment)
                    self.assertIn('pf_passive_geometry', ids.code.parameters)
                    self.assertEqual(len(ids.loop[2].element), 1)
                finally:
                    sys.modules.pop('imas_pf_passive', None)
        self._with_cdf(run)

    def test_si_prefix_case_and_manifest_typo(self):
        self.assertEqual(_source_scale('A'), 1.0)
        self.assertEqual(_source_scale('kA'), 1000.0)
        self.assertEqual(_source_scale('MA'), 1.0e6)
        self.assertEqual(_source_scale('mA'), 1.0e-3)
        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            _source_scale('ma')
        def typo(d):
            d['loops'][0]['resistence_ohm'] = 1.0
        path = self._manifest(typo)
        with self.assertRaisesRegex(ValueError, 'unsupported keys'):
            load_geometry_manifest(path)

    def test_manifest_strict_matching(self):
        path = self._manifest(lambda d: d['loops'].pop(1))
        with self.assertRaisesRegex(ValueError, 'does not map present CDF channels'):
            self._with_cdf(lambda nc: build_passive_plan(nc, geometry_path=path, machine='SYNTHETIC'))
        with self.assertRaisesRegex(ValueError, 'differs from CDF device'):
            self._with_cdf(lambda nc: build_passive_plan(nc, geometry_path=EXAMPLE_GEOMETRY, machine='NSTX'))
        def bad_sign(d):
            d['loops'][0]['cdf_current_sign'] = 0
        path = self._manifest(bad_sign)
        with self.assertRaisesRegex(ValueError, 'must explicitly be'):
            load_geometry_manifest(path)

    def test_geometry_only_when_cdf_has_no_ps(self):
        filename = Path(self.tempdir.name) / 'empty.CDF'
        with netcdf_file(filename, 'w') as nc:
            nc.createDimension('TIME', 2)
            t = nc.createVariable('TIME', 'd', ('TIME',))
            t[:] = [0.0, 1.0]
        with netcdf_file(filename, 'r', mmap=False) as nc:
            self.assertEqual(build_passive_plan(nc).loops, ())
            plan = build_passive_plan(nc, machine='SYNTHETIC', geometry_path=EXAMPLE_GEOMETRY)
            self.assertEqual(len(plan.loops), 3)
            self.assertTrue(all(loop.channel is None for loop in plan.loops))
            ids = ids_node()
            self.assertEqual(populate_passive_ids(ids, plan, 4), 0)

    def test_six_geometry_types(self):
        def add_geometries(d):
            originals = [
                ('oblique', {'type':'oblique', 'r':1.0, 'z':0.0, 'length_alpha':0.1,
                             'length_beta':0.05, 'alpha':0.2, 'beta':0.1}),
                ('arc', {'type':'arcs_of_circle', 'r':[1.0, 1.1], 'z':[0.0, 0.1],
                         'curvature_radii':[0.1, 0.1]}),
                ('thick', {'type':'thick_line', 'first_point':{'r':1.0,'z':0.0},
                           'second_point':{'r':1.1,'z':0.2}, 'thickness':0.02}),
            ]
            for name, geo in originals:
                d['loops'].append({'name': name, 'elements': [
                    {'name':name, 'turns_with_sign':1, 'geometry':geo}]})
        path = self._manifest(add_geometries)
        self._with_cdf(lambda nc: populate_passive_ids(
            ids_node(), build_passive_plan(nc, geometry_path=path, machine='SYNTHETIC'), 4))

    def test_real_fixtures_do_not_fabricate_currents(self):
        roots = Path(__file__).parent / 'data' / 'nstx'
        fixtures = [roots / case / f'{case}.CDF'
                    for case in ('141711Y01', '141711Y02')]
        for path in fixtures:
            with self.subTest(path=path.name), netcdf_file(path, 'r', mmap=False) as nc:
                self.assertEqual(len(build_passive_plan(nc).loops), 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
