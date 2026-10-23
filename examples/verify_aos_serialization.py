#!/usr/bin/python
"""Regression checks for IDS arrays of structures (AoS) that get erased on write.

IMAS-Python's ``IDSStructArray.resize(n)`` defaults to ``keep=False`` and recreates every
element, so "allocate, fill, then ``resize(used)``" silently erases the filled entries
while the writer still reports success.  Mock-based tests miss it when the fake ``resize``
keeps entries, so this script uses real IDS objects:

1. static policy: no ``resize(used)``-style shrink after filling (``keep=True`` or
   ``append_successful_entries`` must be used),
2. ``append_successful_entries`` with real ``IDSFactory`` objects (DD 3.39.0 and 4.1.1),
   including a builder that fails half way,
3. (not with ``--quick``) real conversions of a bundled NSTX case through the radiation,
   distribution_sources and edge_sources writers, read back from HDF5: every entry the
   writer reports must hold data.

    python examples/verify_aos_serialization.py [--quick]
"""

import argparse
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from utilities import append_successful_entries

try:
    import imas
except ImportError:      # pragma: no cover - environment dependent
    imas = None

SKIP_PARTS = {'.git', '.venv', 'build', '__pycache__', 'transp2imas-review'}
# ``x.resize(used)`` / ``x.resize(n_used)`` ... without keep=True: the destructive shrink pattern.
SHRINK = re.compile(r'\.resize\(\s*(?:used|n_used|nused|count|kept|n_kept)\s*\)')
CASE = ROOT / 'examples' / 'data' / 'nstx' / '141711Y05' / '141711Y05.CDF'
RUN = {'quick': False}


class StaticPolicyTests(unittest.TestCase):
    def test_no_destructive_shrink_after_fill(self):
        offenders = []
        for path in sorted(ROOT.rglob('*.py')):
            if any(part in SKIP_PARTS for part in path.parts):
                continue
            for number, line in enumerate(path.read_text(encoding='utf-8').splitlines(), 1):
                if SHRINK.search(line) and 'keep=True' not in line and not line.lstrip().startswith('#'):
                    offenders.append(f'{path.relative_to(ROOT)}:{number}: {line.strip()}')
        self.assertEqual(offenders, [], 'resize(n) erases populated AoS entries; use '
                         'append_successful_entries or resize(n, keep=True)')


@unittest.skipIf(imas is None, 'IMAS-Python is required')
class HelperTests(unittest.TestCase):
    DD_VERSIONS = ('3.39.0', '4.1.1')

    def test_only_successful_entries_are_kept_and_complete(self):
        for dd in self.DD_VERSIONS:
            with self.subTest(dd=dd):
                factory = imas.IDSFactory(dd)
                ids = factory.new('distribution_sources')

                def ok(name, index):
                    def build(source):
                        source.species.type.name = name
                        source.species.type.index = index
                        source.process.resize(1)
                        source.process[0].type.name = 'NBI'
                        return True
                    return build

                def fails_halfway(source):
                    source.species.type.name = 'stale'           # partial data of a failed builder
                    source.process.resize(1)
                    source.process[0].type.index = 99
                    return False

                used = append_successful_entries(lambda: factory.new('distribution_sources'), ids, 'source',
                                                 [ok('a', 1), fails_halfway, ok('b', 2), fails_halfway])
                self.assertEqual(used, 2)
                self.assertEqual(len(ids.source), 2)
                self.assertEqual([str(s.species.type.name) for s in ids.source], ['a', 'b'])
                self.assertEqual([int(s.species.type.index) for s in ids.source], [1, 2])
                self.assertEqual([str(s.process[0].type.name) for s in ids.source], ['NBI', 'NBI'])
                # nothing the failed builders wrote may leak into the appended entries
                self.assertFalse(any(s.process[0].type.index.has_value for s in ids.source))

    def test_no_entries_leaves_the_array_empty(self):
        factory = imas.IDSFactory('4.1.1')
        ids = factory.new('radiation')
        self.assertEqual(append_successful_entries(lambda: factory.new('radiation'), ids, 'process',
                                                   [lambda p: False, lambda p: False]), 0)
        self.assertEqual(len(ids.process), 0)


def convert_and_read_back(cdf, ids_names, dd_version, outdir):
    """Run convert.py for ``ids_names`` and return the opened HDF5 DBEntry."""
    result = subprocess.run(
        [sys.executable, str(ROOT / 'convert.py'), '-i', str(cdf), '-b', 'HDF5', '-u', str(outdir),
         '--dd-version', dd_version, '-new', '--ids', ','.join(ids_names)],
        cwd=ROOT, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stdout[-1500:] + result.stderr[-1500:])
    match = re.search(r'ids stored in (\S+)', result.stdout)
    if not match:
        raise RuntimeError('convert.py did not report the stored entry:\n' + result.stdout[-1500:])
    path = Path(match.group(1))
    database, major, shot, run = path.parts[-4:]
    entry = imas.DBEntry(imas.ids_defs.HDF5_BACKEND, database, int(shot), int(run), str(path.parents[3]),
                         dd_version=dd_version)
    entry.open()
    return entry


@unittest.skipIf(imas is None, 'IMAS-Python is required')
class ConvertedEntryTests(unittest.TestCase):
    """Real HDF5 write/readback of the three writers that used the destructive resize."""

    @classmethod
    def setUpClass(cls):
        if RUN['quick'] or not CASE.is_file():
            raise unittest.SkipTest('end-to-end readback skipped (--quick or case not present)')
        cls.tmp = tempfile.TemporaryDirectory()
        cls.entry = convert_and_read_back(CASE, ['radiation', 'distribution_sources', 'edge_sources'],
                                          '4.1.1', cls.tmp.name)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'entry'):
            cls.entry.close()
            cls.tmp.cleanup()

    def test_radiation_processes_hold_data(self):
        ids = self.entry.get('radiation')
        self.assertGreater(len(ids.process), 0)
        for process in ids.process:
            self.assertTrue(str(process.identifier.name), 'process identifier erased')
            self.assertEqual(len(process.profiles_1d), len(ids.time), 'process profiles erased')

    def test_distribution_sources_hold_data(self):
        ids = self.entry.get('distribution_sources')
        self.assertGreater(len(ids.source), 0)
        for source in ids.source:
            self.assertEqual(str(source.species.type.name), 'ion')
            self.assertEqual(str(source.process[0].type.name), 'NBI')
            self.assertEqual(len(source.profiles_1d), len(ids.time))

    def test_edge_sources_hold_data(self):
        ids = self.entry.get('edge_sources')
        self.assertGreater(len(ids.source), 0)
        for source in ids.source:
            self.assertEqual(len(source.ggd), len(ids.time), 'source ggd erased')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='skip the HDF5 end-to-end readback')
    args, rest = parser.parse_known_args(argv)
    RUN['quick'] = args.quick
    result = unittest.main(argv=[sys.argv[0]] + rest, exit=False, verbosity=1).result
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
