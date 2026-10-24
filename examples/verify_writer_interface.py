#!/usr/bin/python
"""Verify the writer interface: construction does not write, ``build()`` returns the IDS, ``run()`` submits.

Pure tests (no IMAS) of ``writer_base.IdsWriter`` and ``run_writer``, plus a static check that every
migrated writer module

* defines its writer as an ``IdsWriter`` subclass with a ``build`` method,
* never submits to the database itself (no ``putdb_imas*`` call, no ``create_ids``),
* is wrapped in ``run_writer`` by ``convert.make_writers`` (as is every other writer, so migrated and
  legacy writers share one runner).

    python examples/verify_writer_interface.py
"""

import ast
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from writer_base import IdsWriter, run_writer

MIGRATED = ('imas_tf', 'imas_wall', 'imas_lh_antennas', 'imas_ic_antennas', 'imas_nbi', 'imas_ec_launchers',
            'imas_sawteeth', 'imas_dataset_fair', 'imas_workflow', 'imas_summary',
            'imas_transport_solver_numerics', 'imas_pf_active', 'imas_pf_passive',
            'imas_radiation', 'imas_distribution_sources', 'imas_edge_sources', 'imas_waves',
            'imas_core_sources', 'imas_core_profiles', 'imas_equilibrium', 'imas_distributions',
            'imas_core_transport')
# writers whose constructor only stores options (a geometry path, a flag)
OPTIONS_ONLY_INIT = {'imas_ec_launchers', 'imas_pf_active', 'imas_pf_passive', 'imas_equilibrium'}


class FakeConverter:
    def __init__(self):
        self.logger = None
        self.calls = []

    def putdb_imas(self, ids):
        self.calls.append(('put', ids))

    def putdb_imas_profile(self, ids, occurrence):
        self.calls.append(('profile', ids, occurrence))


class Pair(IdsWriter):
    built = 0

    def build(self):
        type(self).built += 1
        return ['A', ('B', 3), 'C']


class Nothing(IdsWriter):
    def build(self):
        return []


class InterfaceTests(unittest.TestCase):
    def test_constructing_a_writer_neither_builds_nor_writes(self):
        converter = FakeConverter()
        Pair.built = 0
        Pair(converter)
        self.assertEqual((Pair.built, converter.calls), (0, []))

    def test_build_returns_the_ids_without_writing(self):
        converter = FakeConverter()
        self.assertEqual(Pair(converter).build(), ['A', ('B', 3), 'C'])
        self.assertEqual(converter.calls, [])

    def test_run_submits_in_order_with_and_without_an_explicit_occurrence(self):
        converter = FakeConverter()
        count = Pair(converter).run()
        self.assertEqual(count, 3)
        self.assertEqual(converter.calls, [('put', 'A'), ('profile', 'B', 3), ('put', 'C')])

    def test_an_empty_build_submits_nothing(self):
        converter = FakeConverter()
        self.assertEqual(Nothing(converter).run(), 0)
        self.assertEqual(converter.calls, [])

    def test_build_must_be_implemented(self):
        with self.assertRaisesRegex(NotImplementedError, 'build'):
            IdsWriter(FakeConverter()).build()

    def test_run_writer_submits_migrated_writers_and_leaves_legacy_ones_alone(self):
        converter = FakeConverter()
        self.assertEqual(run_writer(Pair(converter)), 3)

        class Legacy:                       # did all its work (and its puts) in the constructor
            run = 'not a method of the writer protocol'

        self.assertIsNone(run_writer(Legacy()))
        self.assertEqual(len(converter.calls), 3)

    def test_a_build_that_fails_submits_nothing(self):
        class Broken(IdsWriter):
            def build(self):
                ids = ['first']
                raise RuntimeError(f'mapping failed after building {ids}')

        converter = FakeConverter()
        with self.assertRaisesRegex(RuntimeError, 'mapping failed'):
            Broken(converter).run()
        self.assertEqual(converter.calls, [])


class MigratedWriterTests(unittest.TestCase):
    def tree(self, module):
        return ast.parse((ROOT / f'{module}.py').read_text(encoding='utf-8'))

    def test_migrated_writers_are_idswriters_with_build_and_never_submit_themselves(self):
        for module in MIGRATED:
            with self.subTest(module=module):
                tree = self.tree(module)
                cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == module)
                self.assertEqual([getattr(b, 'id', None) for b in cls.bases], ['IdsWriter'])
                methods = {n.name for n in cls.body if isinstance(n, ast.FunctionDef)}
                self.assertIn('build', methods)
                self.assertNotIn('create_ids', methods)
                if module in OPTIONS_ONLY_INIT:
                    init = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == '__init__')
                    init_calls = {n.func.attr for n in ast.walk(init) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)}
                    self.assertEqual(init_calls, {'__init__'})   # only super().__init__: construction does nothing
                else:
                    self.assertNotIn('__init__', methods)      # construction does nothing
                called = {n.func.attr for n in ast.walk(cls) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)}
                self.assertFalse({c for c in called if c.startswith('putdb_imas')}, called)
                build = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == 'build')
                returns = [n for n in ast.walk(build) if isinstance(n, ast.Return)]
                self.assertTrue(returns and all(r.value is not None for r in returns),
                                'every path of build() returns the IDS list (possibly empty)')

    def test_every_writer_that_convert_imports_is_migrated(self):
        imported = {alias.name for node in self.tree('convert').body if isinstance(node, ast.ImportFrom)
                    and (node.module or '').startswith('imas_') for alias in node.names}
        self.assertEqual(len(imported), 22, sorted(imported))
        self.assertEqual(imported, set(MIGRATED), 'a writer is imported by convert.py but not checked here (or vice versa)')

    def test_every_writer_factory_goes_through_run_writer(self):
        tree = self.tree('convert')
        factory = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == 'make_writers')
        lambdas = [n for n in ast.walk(factory) if isinstance(n, ast.Lambda)]
        self.assertGreaterEqual(len(lambdas), 22)
        for node in lambdas:
            self.assertIsInstance(node.body, ast.Call)
            self.assertEqual(node.body.func.id, 'run_writer')


if __name__ == '__main__':
    unittest.main(verbosity=1)
