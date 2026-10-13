#!/usr/bin/python
"""Pure source-side pf_active regression (no installed IMAS/netCDF4 needed).

The fixtures are synthetic: the early non-ISolver CDF examples do not carry
KK_/VV_/VS_ coil/supply channels.  These tests establish code and schema-node
shape behavior, NOT machine winding orientation or coil-current provenance.
"""

from __future__ import annotations

import importlib
import json
import logging
from pathlib import Path
import sys
import tempfile
import types
from types import SimpleNamespace
from unittest import mock
import unittest

import numpy as np
from scipy.io import netcdf_file

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pf_active_mapping import (  # noqa: E402
    build_active_plan, read_active_channels, populate_active_ids,
    load_active_manifest, _voltage_scale,
)

GEOMETRY = Path(__file__).with_name("pf_active_geometry.synthetic.json")


class Node(SimpleNamespace):
    """Mock structured IMAS node."""


class StructureArray(list):
    def __init__(self, factory):
        super().__init__()
        self.factory = factory

    def resize(self, n):
        self.clear()
        self.extend(self.factory() for _ in range(n))


def trace_node():
    return Node(time=None, data=None)


def element_node():
    return Node(name=None, description=None, identifier=None, turns_with_sign=None,
                area=None, geometry=Node(
                    geometry_type=None, outline=Node(), rectangle=Node(),
                    oblique=Node(), arcs_of_circle=Node(), annulus=Node(),
                    thick_line=Node(first_point=Node(), second_point=Node())))


def coil_node():
    return Node(name=None, identifier=None, description=None, resistance=None,
                current=trace_node(), voltage=trace_node(),
                element=StructureArray(element_node))


def supply_node():
    return Node(name=None, identifier=None, description=None,
                current=trace_node(), voltage=trace_node())


def circuit_node():
    return Node(name=None, identifier=None, description=None, connections=None,
                current=trace_node(), voltage=trace_node())


def ids_node():
    return Node(coil=StructureArray(coil_node), supply=StructureArray(supply_node),
                circuit=StructureArray(circuit_node),
                ids_properties=Node(comment="test", homogeneous_time=None),
                code=Node(parameters="<parameters/>"), time=None)


class TestActiveMapping(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.cdf = Path(self.tmp.name) / "test.CDF"
        with netcdf_file(self.cdf, "w") as ds:
            ds.createDimension("TIME", 3)
            ds.createDimension("TIME3", 2)
            t = ds.createVariable("TIME", "d", ("TIME",))
            t[:] = [0.1, 0.2, 0.3]
            t3 = ds.createVariable("TIME3", "d", ("TIME3",))
            t3[:] = [0.1, 0.3]
            for name, dimension, values, units in (
                ("KK_PF1AU", "TIME", [1.0, 2.0, 3.0], "kA"),
                ("KK_CS1", "TIME", [100, 200, 300], "mA"),
                ("VV_PF1AU", "TIME3", [10, 20], "V"),
                ("VS_PF1AU", "TIME", [0.01, 0.02, 0.03], "kV"),
                ("CC_EQ1", "TIME", [90, 80, 70], "A"),
                ("KT_PF1AU", "TIME", [5, 10, 15], "kA"),
            ):
                v = ds.createVariable(name, "d", (dimension,))
                v[:] = values
                v.units = units
                v.long_name = name

    def with_nc(self, action):
        with netcdf_file(self.cdf, "r", mmap=False) as nc:
            return action(nc)

    def manifest(self, mutate):
        data = json.loads(GEOMETRY.read_text(encoding="utf-8"))
        mutate(data)
        path = Path(self.tmp.name) / "machine.json"
        path.write_text(json.dumps(data), encoding="utf-8")
        return path

    def test_channel_detection_units_and_native_times(self):
        def run(nc):
            out = read_active_channels(nc)
            self.assertIn("KK_PF1AU", out)
            self.assertNotIn("KT_PF1AU", out)
            self.assertNotIn("CC_EQ1", out)
            np.testing.assert_array_equal(out["KK_PF1AU"].values, [1000, 2000, 3000])
            np.testing.assert_array_equal(out["KK_CS1"].values, [0.1, 0.2, 0.3])
            np.testing.assert_array_equal(out["VV_PF1AU"].time, [0.1, 0.3])
            np.testing.assert_array_equal(out["VS_PF1AU"].values, [10, 20, 30])
            self.assertIn("CC_EQ1", read_active_channels(nc, include_equilibrium=True))
        self.with_nc(run)

    def test_cdf_only_no_target_as_actual(self):
        def run(nc):
            plan = build_active_plan(nc)
            self.assertEqual([c.name for c in plan.coils], ["CS1", "PF1AU"])
            self.assertEqual(len(plan.supplies), 1)
            self.assertIsNone(plan.supplies[0].current)
            self.assertEqual(plan.ignored_targets, ("KT_PF1AU",))
            self.assertEqual(plan.ignored_equilibrium, ("CC_EQ1",))
            ids = ids_node()
            counts = populate_active_ids(ids, plan, 4)
            self.assertEqual(counts["coil_current"], 2)
            self.assertEqual(counts["coil_voltage"], 1)
            self.assertEqual(counts["supply_current"], 0)
            self.assertEqual(counts["supply_voltage"], 1)
            self.assertIsNone(ids.supply[0].current.data)
            np.testing.assert_array_equal(ids.coil[1].voltage.time, [0.1, 0.3])
            np.testing.assert_array_equal(ids.supply[0].voltage.data, [10, 20, 30])
        self.with_nc(run)

    def test_authoritative_geometry_signs_and_dd3_dd4_circuits(self):
        def run(nc):
            plan = build_active_plan(nc, geometry_path=GEOMETRY, machine="SYNTHETIC")
            self.assertEqual(len(plan.coils), 3)
            self.assertEqual(len(plan.supplies), 1)
            self.assertEqual(len(plan.circuits), 1)
            self.assertIn("CC_EQ1", [c.current.channel.variable for c in plan.coils])
            ids4 = ids_node()
            counts = populate_active_ids(ids4, plan, 4)
            self.assertEqual(counts["supply_current"], 1)
            np.testing.assert_allclose(ids4.coil[2].current.data, [30., 80./3., 70./3.])
            np.testing.assert_array_equal(ids4.coil[0].current.data, [-1000, -2000, -3000])
            np.testing.assert_array_equal(ids4.supply[0].current.data, [1000, 2000, 3000])
            np.testing.assert_array_equal(ids4.supply[0].voltage.data, [-10, -20, -30])
            np.testing.assert_array_equal(ids4.coil[0].voltage.data, [10, 20])
            self.assertEqual(ids4.coil[0].element[0].geometry.geometry_type, 2)
            self.assertEqual(len(ids4.coil[1].element[0].geometry.outline.r), 5)
            self.assertEqual(ids4.circuit[0].connections.shape, (3, 4))
            ids3 = ids_node()
            populate_active_ids(ids3, plan, 3)
            self.assertEqual(ids3.coil[0].identifier, "PF1_AU")
            self.assertEqual(len(ids3.coil[1].element[0].geometry.outline.r), 4)
            self.assertEqual(ids3.circuit[0].connections.shape, (4, 8))
        self.with_nc(run)

    def test_writer_end_to_end_mock_ids(self):
        def run(nc):
            class MockImas:
                device = "SYNTHETIC"
                dd_major_version = "4"
                logger = logging.getLogger("active-regression")

                def __init__(self):
                    self.nc = nc
                    self.record = None

                def new_ids(self, name):
                    self_name = name
                    assert self_name == "pf_active"
                    return ids_node()

                def fill_properties(self, obj):
                    return obj

                def fill_code_properties(self, obj):
                    return obj

                def get_ncvar(self, var):
                    return np.asarray(self.nc.variables[var][:]) if var in self.nc.variables else np.array([])

                def putdb_imas(self, obj):
                    self.record = obj

            fake = types.ModuleType("nc2imas")
            fake.nc2imas = MockImas
            with mock.patch.dict(sys.modules, {"nc2imas": fake}):
                sys.modules.pop("imas_pf_active", None)
                try:
                    module = importlib.import_module("imas_pf_active")
                    converter = MockImas()
                    module.imas_pf_active(converter, geometry_path=str(GEOMETRY))
                    ids = converter.record
                    self.assertIsNotNone(ids)
                    self.assertEqual(ids.ids_properties.homogeneous_time, 0)
                    np.testing.assert_array_equal(ids.time, [0.1, 0.2, 0.3])
                    self.assertIn("pf_active_geometry", ids.code.parameters)
                    self.assertEqual(len(ids.coil), 3)
                    self.assertEqual(len(ids.circuit), 1)
                finally:
                    sys.modules.pop("imas_pf_active", None)
        self.with_nc(run)

    def test_manifest_rejects_typos_bad_sign_and_unmapped_sources(self):
        path = self.manifest(lambda d: d["coils"][0].update(resistence_ohm=0.1))
        with self.assertRaisesRegex(ValueError, "unsupported keys"):
            load_active_manifest(path)
        path = self.manifest(lambda d: d["coils"][0].update(current_sign=0))
        with self.assertRaisesRegex(ValueError, "must explicitly be"):
            self.with_nc(lambda nc: build_active_plan(nc, geometry_path=path))
        path = self.manifest(lambda d: d["coils"][0].pop("cdf_current"))
        with self.assertRaisesRegex(ValueError, "without cdf_current"):
            self.with_nc(lambda nc: build_active_plan(nc, geometry_path=path))
        path = self.manifest(lambda d: d["coils"][0].update(cdf_current="KT_PF1AU"))
        with self.assertRaisesRegex(ValueError, "expected one of"):
            self.with_nc(lambda nc: build_active_plan(nc, geometry_path=path))
        path = self.manifest(lambda d: d["coils"][2].pop("cc_effective_turns"))
        with self.assertRaisesRegex(ValueError, "requires authoritative cc_effective_turns"):
            self.with_nc(lambda nc: build_active_plan(nc, geometry_path=path))
        path = self.manifest(lambda d: d["coils"][2].update(cc_effective_turns=0))
        with self.assertRaisesRegex(ValueError, "must be nonzero"):
            self.with_nc(lambda nc: build_active_plan(nc, geometry_path=path))
        path = self.manifest(lambda d: d["coils"][0].update(cc_effective_turns=4))
        with self.assertRaisesRegex(ValueError, "only for CC_"):
            self.with_nc(lambda nc: build_active_plan(nc, geometry_path=path))
        path = self.manifest(lambda d: d["coils"][2].pop("cdf_current"))
        with self.assertRaisesRegex(ValueError, "without cdf_current"):
            self.with_nc(lambda nc: build_active_plan(nc, geometry_path=path))
        with self.assertRaisesRegex(ValueError, "differs from CDF device"):
            self.with_nc(lambda nc: build_active_plan(nc, geometry_path=GEOMETRY, machine="NSTX"))

    def test_missing_partial_data_has_no_dummy_or_faked_current(self):
        empty = Path(self.tmp.name) / "empty.CDF"
        with netcdf_file(empty, "w") as nc:
            nc.createDimension("TIME", 2)
            var = nc.createVariable("TIME", "d", ("TIME",))
            var[:] = [0, 1]
            var = nc.createVariable("KT_PF1AU", "d", ("TIME",))
            var[:] = [1, 2]
        with netcdf_file(empty, "r", mmap=False) as nc:
            plan = build_active_plan(nc)
            self.assertEqual(plan.coils, ())
            self.assertEqual(plan.supplies, ())
            self.assertEqual(plan.ignored_targets, ("KT_PF1AU",))
        # Geometry-only machine specification still creates legitimate coil descriptions.
        with netcdf_file(empty, "r", mmap=False) as nc:
            plan = build_active_plan(nc, geometry_path=GEOMETRY, machine="SYNTHETIC")
            self.assertEqual(len(plan.coils), 3)
            self.assertTrue(all(c.current is None for c in plan.coils))
        self.assertEqual(_voltage_scale("mV"), 1e-3)
        with self.assertRaises(ValueError):
            _voltage_scale("mv")

    def test_invalid_source_units_or_time_shape_raise(self):
        wrong = Path(self.tmp.name) / "wrong_units.CDF"
        with netcdf_file(wrong, "w") as nc:
            nc.createDimension("TIME", 2)
            t = nc.createVariable("TIME", "d", ("TIME",))
            t[:] = [0.0, 1.0]
            v = nc.createVariable("VS_SUPPLY", "d", ("TIME",))
            v[:] = [1.0, 2.0]
            v.units = "mv"  # ambiguous mV/MV prefix
        with netcdf_file(wrong, "r", mmap=False) as nc:
            with self.assertRaisesRegex(ValueError, "Unsupported/ambiguous"):
                read_active_channels(nc)
        wrong_time = Path(self.tmp.name) / "wrong_time.CDF"
        with netcdf_file(wrong_time, "w") as nc:
            nc.createDimension("TIME", 2)
            t = nc.createVariable("TIME", "d", ("TIME",))
            t[:] = [0.3, 0.1]
            v = nc.createVariable("KK_COIL", "d", ("TIME",))
            v[:] = [1.0, 2.0]
            v.units = "A"
        with netcdf_file(wrong_time, "r", mmap=False) as nc:
            with self.assertRaisesRegex(ValueError, "strictly increasing"):
                read_active_channels(nc)

    def test_circuit_dd_specific_and_duplicate_machine_connections(self):
        def run(nc):
            path = self.manifest(lambda d: d["circuits"][0].pop("connections_dd3"))
            plan = build_active_plan(nc, geometry_path=path, machine="SYNTHETIC")
            ids3 = ids_node()
            counts = populate_active_ids(ids3, plan, 3)
            self.assertEqual(counts["circuits"], 0)
            ids4 = ids_node()
            counts = populate_active_ids(ids4, plan, 4)
            self.assertEqual(counts["circuits"], 1)
            path = self.manifest(lambda d: d["coils"][1].update(name="PF1_AU"))
            with self.assertRaisesRegex(ValueError, "duplicate machine identifiers"):
                build_active_plan(nc, geometry_path=path, machine="SYNTHETIC")
            path = self.manifest(lambda d: d["circuits"][0].update(connections_dd4=[[1, 0, 0]] * 2))
            with self.assertRaisesRegex(ValueError, "expected 4 columns"):
                build_active_plan(nc, geometry_path=path, machine="SYNTHETIC")
        self.with_nc(run)

    def test_writer_skips_no_coil_without_dummy_time(self):
        empty = Path(self.tmp.name) / "only_target.CDF"
        with netcdf_file(empty, "w") as nc:
            nc.createDimension("TIME", 2)
            t = nc.createVariable("TIME", "d", ("TIME",))
            t[:] = [0.0, 1.0]
            target = nc.createVariable("KT_CS1", "d", ("TIME",))
            target[:] = [0.1, 0.2]
        with netcdf_file(empty, "r", mmap=False) as nc:
            class Fake:
                device = "SYNTHETIC"
                logger = logging.getLogger("pfactive-no-current")
                def __init__(self):
                    self.nc = nc
                def new_ids(self, name):
                    raise AssertionError("No IDS may be manufactured when source is only a target")
            fake = types.ModuleType("nc2imas")
            fake.nc2imas = Fake
            with mock.patch.dict(sys.modules, {"nc2imas": fake}):
                sys.modules.pop("imas_pf_active", None)
                try:
                    module = importlib.import_module("imas_pf_active")
                    module.imas_pf_active(Fake())
                finally:
                    sys.modules.pop("imas_pf_active", None)

    def test_real_regression_cases_have_no_active_channels(self):
        # Only known non-ISolver fixtures. The new ITER/NSTX cases DO have
        # coil channels and are audited separately by make verify-isolver.
        roots = Path(__file__).parent / "data" / "nstx"
        paths = [roots / case / f"{case}.CDF"
                 for case in ("141711Y01", "141711Y02")]
        for p in paths:
            with self.subTest(p=p.name), netcdf_file(p, "r", mmap=False) as nc:
                plan = build_active_plan(nc)
                self.assertEqual(len(plan.coils), 0)
                self.assertEqual(len(plan.supplies), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
