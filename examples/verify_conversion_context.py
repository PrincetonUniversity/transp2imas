#!/usr/bin/python
"""Verify the conversion context that IDS writers use instead of the whole converter.

Behaviour of ``conversion_context`` (pure tests, no IMAS) and static checks tying it to the code:

* a view forwards its listed names to the owner at call time and refuses any other name, with a message that
  names the view; it is read-only;
* no name is listed by two views;
* every listed name exists on ``nc2imas`` (a method, class attribute, or an attribute its constructors set),
  so a list cannot name something that was renamed away;
* every ``self.<view>.<name>`` in a writer module uses a listed name and sits in an ``IdsWriter`` subclass;
* writers do not use ``self.im`` at all any more.

    python examples/verify_conversion_context.py
"""

import ast
from pathlib import Path
import re
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import conversion_context as cc
from writer_base import IdsWriter

WRITER_FILES = sorted(p for p in ROOT.glob('imas_*.py') if p.name not in ('imas_class.py', 'imas_store.py'))
VIEW_ATTRIBUTES = set(cc.VIEWS) | {'ctx'}
# Names set on the converter outside the class bodies: the COCOS sign fields are attached by ``_init_cocos_signs``.
SET_ELSEWHERE = {'ip_sign', 'bt_sign', 'psi_sign', 'q_sign', 'phi_sign', 'f_sign', 'jparallel_sign', 'vloop_sign',
                 'jphi_sign', 'toroidal_rotation_sign', 'poloidal_field_sign', 'theta_raw_from_target'}


def class_members(path, class_name):
    """Methods, class attributes and ``self.X =`` assignments of ``class_name`` in ``path``."""
    tree = ast.parse(path.read_text(encoding='utf-8'))
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == class_name)
    names = set()
    for node in cls.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            names.update(t.id for t in node.targets if isinstance(t, ast.Name))
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
    for node in ast.walk(cls):
        if isinstance(node, ast.Attribute) and isinstance(node.ctx, ast.Store) \
                and isinstance(node.value, ast.Name) and node.value.id == 'self':
            names.add(node.attr)
    return names


def converter_names():
    names = set(SET_ELSEWHERE)
    for file, cls in (('nc2imas.py', 'nc2imas'), ('imas_class.py', 'imas_connection'),
                      ('netcdf_class.py', 'netcdf_connection')):
        names |= class_members(ROOT / file, cls)
    return names


def view_uses(path):
    """[(class name, view, name, line)] of every ``self.<view>.<name>`` in ``path``."""
    tree = ast.parse(path.read_text(encoding='utf-8'))
    found = []
    for cls in (n for n in ast.walk(tree) if isinstance(n, ast.ClassDef)):
        for node in ast.walk(cls):
            if (isinstance(node, ast.Attribute) and isinstance(node.value, ast.Attribute)
                    and isinstance(node.value.value, ast.Name) and node.value.value.id == 'self'
                    and node.value.attr in VIEW_ATTRIBUTES):
                found.append((cls, node.value.attr, node.attr, node.lineno))
    return found


class Owner:
    """Stand-in for the converter: only what the test reads."""
    ip_sign = -1
    mapping = 'ledger'
    logger = 'log'
    debug_max = False

    def get_ncvar(self, name, scale=1.0):
        return (name, scale)


class ViewBehaviour(unittest.TestCase):
    def test_a_listed_name_is_forwarded_at_call_time(self):
        owner = Owner()
        view = cc.SignsView(owner)
        self.assertEqual(view.ip_sign, -1)
        owner.ip_sign = 1                       # live state, not a copy
        self.assertEqual(view.ip_sign, 1)
        self.assertEqual(cc.CdfView(owner).get_ncvar('TE', 2.0), ('TE', 2.0))

    def test_an_unlisted_name_is_refused_and_names_the_view(self):
        view = cc.SignsView(Owner())
        with self.assertRaises(AttributeError) as caught:
            view.get_ncvar
        self.assertIn("'signs'", str(caught.exception))
        self.assertIn('ip_sign', str(caught.exception))

    def test_a_listed_name_the_owner_lacks_is_an_attribute_error_so_getattr_defaults_work(self):
        self.assertEqual(getattr(cc.CdfView(Owner()), 'vars_', 'absent'), 'absent')

    def test_views_are_read_only(self):
        with self.assertRaises(AttributeError):
            cc.SignsView(Owner()).ip_sign = 5

    def test_no_name_is_listed_by_two_views(self):
        seen = {}
        for view, cls in cc.VIEWS.items():
            for name in cls.names:
                self.assertNotIn(name, seen, f'{name!r} is listed by {seen.get(name)!r} and {view!r}')
                seen[name] = view

    def test_the_context_offers_every_view_mapping_and_logger(self):
        context = cc.ConversionContext(Owner())
        for view in cc.VIEWS:
            self.assertIsInstance(getattr(context, view), cc.ContextView)
        self.assertEqual((context.mapping, context.logger, context.debug_max), ('ledger', 'log', False))

    def test_of_reuses_the_owner_context_and_wraps_a_double(self):
        owner = Owner()
        owner.context = cc.ConversionContext(owner)
        self.assertIs(cc.ConversionContext.of(owner), owner.context)
        wrapped = cc.ConversionContext.of(Owner())
        self.assertIsInstance(wrapped, cc.ConversionContext)
        self.assertIsInstance(cc.ConversionContext.of(object()).cdf, cc.CdfView)

    def test_a_writer_gets_views_from_whatever_im_it_holds(self):
        class Writer(IdsWriter):
            pass
        writer = Writer.__new__(Writer)         # built without __init__, as several tests do
        writer.im = Owner()
        self.assertEqual(writer.signs.ip_sign, -1)
        first = writer.ctx
        self.assertIs(writer.ctx, first)        # cached per converter
        writer.im = Owner()
        self.assertIsNot(writer.ctx, first)     # a different converter gets its own views


class SurfaceChecks(unittest.TestCase):
    def test_every_listed_name_exists_on_the_converter(self):
        available = converter_names()
        for view, cls in cc.VIEWS.items():
            missing = [n for n in cls.names if n not in available]
            self.assertEqual(missing, [], f'view {view!r} lists names the converter does not have')

    def test_writers_use_only_listed_names_and_only_inside_writer_classes(self):
        problems = []
        for path in WRITER_FILES:
            for cls, view, name, line in view_uses(path):
                if not any('IdsWriter' in ast.unparse(b) for b in cls.bases):
                    problems.append(f'{path.name}:{line} {cls.name} is not an IdsWriter but uses self.{view}.{name}')
                    continue
                allowed = {'mapping', 'logger', 'debug_max', 'supports_occurrences'} | set(cc.VIEWS) if view == 'ctx' else set(cc.VIEWS[view].names)
                if name not in allowed:
                    problems.append(f'{path.name}:{line} self.{view}.{name} is not offered by the {view!r} view')
        self.assertEqual(problems, [])

    def test_writers_do_not_use_self_im_at_all(self):
        """Everything a writer needs comes through a view, ``self.ctx`` or ``self.logger``."""
        pattern = re.compile(r'\bself\.im\b')
        problems = [f'{path.name}:{n}: {line.strip()}' for path in WRITER_FILES
                    for n, line in enumerate(path.read_text(encoding='utf-8').splitlines(), 1)
                    if pattern.search(line) and not line.lstrip().startswith('#')]
        self.assertEqual(problems, [], 'use a view (self.cdf, self.signs, ...) or add the name to one')


if __name__ == '__main__':
    unittest.main(verbosity=1)
