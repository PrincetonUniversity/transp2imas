#!/usr/bin/python
"""Verify the write queue, the conversion runner and the Data Dictionary adapter (no IMAS needed).

The review asks that "no database writes occur merely as a side effect of constructing a writer".
Writers still build and submit their IDS in their constructor, so the runner turns deferral on and
flushes after each writer returned:

* constructing a writer writes nothing; the runner's flush does, in submission order,
* a writer that raises writes none of its own output, is journaled FAIL with its name, and later
  writers do not run; earlier writers' output stays,
* the readback manifest names each IDS by the writer that produced it, and only accepted puts count,
* a backend failure during the flush is attributed to the writer being flushed,
* ``DDAdapter`` answers the version questions writers ask.

    python examples/verify_conversion_runner.py
"""

import importlib.util
from pathlib import Path
import sys
import types
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from conversion_plan import ConversionJournal, ConversionRunner
from dd_adapter import DDAdapter, parse_dd_version, resolve_dd_version
from imas_store import ImasStore, WriteRequest
from mapping_status import MappingLedger


_CONNECTION_CLASS = []


def load_connection_class():
    """imas_class.py with a stand-in ``imas`` module, so no IMAS installation is needed.

    Loaded once, with numpy imported beforehand: a module first imported inside ``patch.dict`` on
    ``sys.modules`` is evicted again afterwards, and numpy cannot be imported twice per process.
    """
    if _CONNECTION_CLASS:
        return _CONNECTION_CLASS[0]
    import numpy  # noqa: F401
    spec = importlib.util.spec_from_file_location('_imas_class_for_runner_test', ROOT / 'imas_class.py')
    module = importlib.util.module_from_spec(spec)
    with mock.patch.dict(sys.modules, {'imas': types.ModuleType('imas')}):
        spec.loader.exec_module(module)
    _CONNECTION_CLASS.append(module.imas_connection)
    return module.imas_connection


def make_connection(put_log, fail_on=None):
    """A connection object with a recording backend (``ids.put``) and no database."""
    conn = object.__new__(load_connection_class())
    conn.occur = 0
    conn._written_ids = []
    conn._active_writer = None
    conn.mapping = MappingLedger()

    def put(ids, occurrence):
        if ids == fail_on:
            raise ValueError(f'IMAS invalid coordinates in {ids}')
        put_log.append((ids, occurrence))

    conn.ids = mock.Mock()
    conn.ids.put.side_effect = put
    return conn


class DDAdapterTests(unittest.TestCase):
    def test_version_questions(self):
        dd3, dd4 = DDAdapter('3.39.0'), DDAdapter('4.1.1')
        self.assertTrue(dd3.is_dd3 and not dd3.is_dd4 and dd3.major == '3')
        self.assertTrue(dd4.is_dd4 and not dd4.is_dd3 and dd4.major == '4')
        self.assertTrue(dd3.before(3, 42) and not dd3.at_least(3, 42))
        self.assertTrue(dd4.at_least(3, 42) and dd4.at_least(4, 1, 0) and not dd4.at_least(4, 1, 2))

    def test_parsing_and_selection(self):
        self.assertEqual(parse_dd_version('4.1'), (4, 1, 0))
        self.assertEqual(parse_dd_version('3.42.0-dev'), (3, 42, 0))
        default = lambda: '4.1.1'
        self.assertEqual(resolve_dd_version('3.39.0', default, {}), '3.39.0')
        self.assertEqual(resolve_dd_version('auto', default, {'IMAS_VERSION': '3.41.0'}), '3.41.0')
        self.assertEqual(resolve_dd_version(None, default, {}), '4.1.1')

    def test_connection_delegates_to_the_adapter(self):
        conn = object.__new__(load_connection_class())
        conn._set_dd('3.39.0')
        self.assertEqual((conn.dd_version, conn.dd_major_version, conn.dd_version_tuple), ('3.39.0', '3', (3, 39, 0)))
        self.assertTrue(conn.dd_before(4) and conn.is_dd3 and not conn.is_dd4 and not conn.dd_at_least(3, 42))


class StoreTests(unittest.TestCase):
    def test_immediate_and_deferred_submission(self):
        written = []
        store = ImasStore()
        store.submit(WriteRequest('a', 0, 'w'), lambda r: written.append(r.ids))
        self.assertEqual(written, ['a'])
        store.deferred = True
        store.submit(WriteRequest('b', 0, 'w'), lambda r: written.append(r.ids))
        store.submit(WriteRequest('c', 1, 'w'), lambda r: written.append(r.ids))
        self.assertEqual((written, store.pending), (['a'], 2))
        store.flush(lambda r: written.append(r.ids))
        self.assertEqual((written, store.pending), (['a', 'b', 'c'], 0))

    def test_failing_flush_drops_the_rest_and_discard_empties_the_queue(self):
        written = []

        def write(request):
            if request.ids == 'bad':
                raise ValueError('rejected')
            written.append(request.ids)

        store = ImasStore()
        store.deferred = True
        for name in ('ok', 'bad', 'never'):
            store.submit(WriteRequest(name, 0), write)
        with self.assertRaises(ValueError):
            store.flush(write)
        self.assertEqual((written, store.pending), (['ok'], 0))
        store.submit(WriteRequest('x', 0), write)
        self.assertEqual(store.discard(), 1)
        self.assertEqual(store.pending, 0)


class RunnerTests(unittest.TestCase):
    def runner(self, conn, writers):
        journal = ConversionJournal('case.CDF', tuple(writers))
        return ConversionRunner(conn, writers, journal), journal

    def test_construction_writes_nothing_and_the_flush_writes_in_order_with_attribution(self):
        puts = []
        conn = make_connection(puts)
        seen_during_construction = {}

        def writer(name, ids, occurrences=(0,)):
            def build():
                for occ in occurrences:
                    conn.putdb_imas_profile(f'{ids}{occ}', occ) if occ else conn.putdb_imas(ids)
                seen_during_construction[name] = list(puts)       # what the backend had when construction ended
            return build

        runner, journal = self.runner(conn, {'summary': writer('summary', 'S'),
                                             'core_transport': writer('core_transport', 'CT', (0, 1))})
        runner.run(('summary', 'core_transport'))
        self.assertEqual(seen_during_construction, {'summary': [], 'core_transport': [('S', 0)]})
        self.assertEqual(puts, [('S', 0), ('CT', 0), ('CT1', 1)])
        self.assertEqual(conn._written_ids, [{'ids': 'summary', 'occurrence': 0},
                                             {'ids': 'core_transport', 'occurrence': 0},
                                             {'ids': 'core_transport', 'occurrence': 1}])
        self.assertEqual([(w['ids'], w['status'], w['write_count']) for w in journal.data['writers']],
                         [('summary', 'WROTE', 1), ('core_transport', 'WROTE', 2)])
        self.assertFalse(conn.store.deferred)

    def test_a_failing_writer_writes_none_of_its_output_and_stops_the_run(self):
        puts, ran = [], []
        conn = make_connection(puts)

        def good():
            ran.append('good')
            conn.putdb_imas('G')

        def broken():
            ran.append('broken')
            conn.putdb_imas('B1')
            raise RuntimeError('mapping failed')

        def later():
            ran.append('later')

        runner, journal = self.runner(conn, {'good': good, 'broken': broken, 'later': later})
        with self.assertRaisesRegex(RuntimeError, 'mapping failed'):
            runner.run(('good', 'broken', 'later'))
        self.assertEqual(ran, ['good', 'broken'])
        self.assertEqual(puts, [('G', 0)])
        self.assertEqual(conn._written_ids, [{'ids': 'good', 'occurrence': 0}])
        rows = journal.data['writers']
        self.assertEqual([(r['ids'], r['status']) for r in rows], [('good', 'WROTE'), ('broken', 'FAIL')])
        self.assertIn('mapping failed', rows[1]['error'])
        self.assertEqual(conn.store.pending, 0)
        self.assertFalse(conn.store.deferred)

    def test_a_backend_rejection_at_flush_is_attributed_to_the_writer_and_not_recorded(self):
        puts = []
        conn = make_connection(puts, fail_on='INVALID')
        runner, journal = self.runner(conn, {'equilibrium': lambda: conn.putdb_imas('INVALID')})
        with self.assertRaisesRegex(ValueError, 'invalid coordinates'):
            runner.run(('equilibrium',))
        self.assertEqual((puts, conn._written_ids), ([], []))
        row = journal.data['writers'][0]
        self.assertEqual((row['ids'], row['status'], row['write_count']), ('equilibrium', 'FAIL', 0))
        self.assertIn('invalid coordinates', row['error'])

    def test_a_writer_that_writes_nothing_gets_a_recorded_reason(self):
        conn = make_connection([])
        runner, journal = self.runner(conn, {'wall': lambda: None})
        runner.run(('wall',))
        row = journal.data['writers'][0]
        self.assertEqual((row['status'], row['write_count']), ('SKIP', 0))
        self.assertIn('skipped', row['mapping'])

    def test_put_without_a_runner_still_writes_immediately(self):
        puts = []
        conn = make_connection(puts)
        conn.putdb_imas('X')
        self.assertEqual(puts, [('X', 0)])
        # no active writer and no IDS metadata to name it: the manifest falls back to 'unknown'
        self.assertEqual(conn._written_ids, [{'ids': 'unknown', 'occurrence': 0}])


class PublicationReportTests(unittest.TestCase):
    """A backend failure after earlier puts is a *partial publication*, and it is reported as such.

    The writes are sequential, not a transaction: nothing here rolls back, deletes or rewrites what was accepted.
    """

    @staticmethod
    def identity(ids, occurrence, writer='w'):
        return {'writer': writer, 'ids': ids, 'occurrence': occurrence}

    def store(self, requests):
        store = ImasStore()
        store.deferred = True
        for ids, occurrence in requests:
            store.submit(WriteRequest(ids, occurrence, 'w'), lambda r: None)
        return store

    def test_a_complete_flush_reports_every_request_as_accepted(self):
        store = self.store([('a', 0), ('b', 0), ('b', 1)])
        report = store.flush(lambda r: None)
        self.assertTrue(report.complete)
        self.assertEqual(report.as_dict(), {
            'complete': True, 'failed': None, 'unattempted': [],
            'requested': [self.identity('a', 0), self.identity('b', 0), self.identity('b', 1)],
            'accepted': [self.identity('a', 0), self.identity('b', 0), self.identity('b', 1)]})

    def test_a_failure_reports_requested_accepted_failed_and_unattempted_with_identities(self):
        written = []

        def write(request):
            if (request.ids, request.occurrence) == ('b', 1):
                raise ValueError('rejected coordinates')
            written.append((request.ids, request.occurrence))

        store = self.store([('a', 0), ('b', 0), ('b', 1), ('c', 0), ('c', 1)])
        with self.assertRaises(ValueError) as caught:                         # the original exception type is kept
            store.flush(write)
        report = caught.exception.publication_report.as_dict()
        self.assertFalse(report['complete'])
        self.assertEqual(report['requested'], [self.identity(i, o) for i, o in
                                               (('a', 0), ('b', 0), ('b', 1), ('c', 0), ('c', 1))])
        self.assertEqual(report['accepted'], [self.identity('a', 0), self.identity('b', 0)])
        self.assertEqual(report['failed'], {**self.identity('b', 1), 'error': 'ValueError: rejected coordinates'})
        self.assertEqual(report['unattempted'], [self.identity('c', 0), self.identity('c', 1)])
        self.assertEqual(written, [('a', 0), ('b', 0)])                        # accepted data stays; nothing was undone
        self.assertEqual(store.pending, 0)

    def test_the_identity_names_the_ids_by_its_metadata_when_the_writer_is_unknown(self):
        ids = types.SimpleNamespace(metadata=types.SimpleNamespace(name='core_transport'))
        store = ImasStore()
        store.deferred = True
        store.submit(WriteRequest(ids, 2, None), lambda r: None)
        self.assertEqual(store.flush(lambda r: None).as_dict()['requested'],
                         [{'writer': None, 'ids': 'core_transport', 'occurrence': 2}])

    def runner(self, conn, writers):
        journal = ConversionJournal('case.CDF', tuple(writers))
        return ConversionRunner(conn, writers, journal), journal

    def test_the_journal_records_a_backend_failure_after_earlier_writes_as_partial(self):
        conn = make_connection([], fail_on='CT1')

        def core_transport():
            conn.putdb_imas('CT0')
            conn.putdb_imas_profile('CT1', 1)
            conn.putdb_imas_profile('CT2', 2)

        runner, journal = self.runner(conn, {'summary': lambda: conn.putdb_imas('S'), 'core_transport': core_transport})
        with self.assertRaisesRegex(ValueError, 'invalid coordinates in CT1'):
            runner.run(('summary', 'core_transport'))
        publication = journal.data['writers'][1]['publication']
        self.assertEqual([r['occurrence'] for r in publication['accepted']], [0])
        self.assertEqual((publication['failed']['occurrence'], publication['failed']['error'].split(':')[0]), (1, 'ValueError'))
        self.assertEqual([r['occurrence'] for r in publication['unattempted']], [2])
        summary = journal.data['publication']
        self.assertEqual((summary['status'], summary['kind']), ('partial', 'backend_failure'))
        self.assertEqual([(r['ids'], r['occurrence']) for r in summary['accepted']],
                         [('summary', 0), ('core_transport', 0)])         # earlier writers' output is accepted too
        self.assertEqual(summary['failed']['ids'], 'core_transport')

    def test_a_build_failure_is_not_a_backend_failure(self):
        conn = make_connection([])

        def broken():
            raise RuntimeError('mapping failed')

        runner, journal = self.runner(conn, {'summary': lambda: conn.putdb_imas('S'), 'wall': broken})
        with self.assertRaisesRegex(RuntimeError, 'mapping failed'):
            runner.run(('summary', 'wall'))
        summary = journal.data['publication']
        self.assertEqual((summary['status'], summary['kind']), ('partial', 'build_failure'))
        self.assertEqual([(r['ids'], r['occurrence']) for r in summary['accepted']], [('summary', 0)])
        self.assertIsNone(summary['failed'])                                    # nothing was attempted for the failed writer

    def test_a_failure_before_any_write_is_nothing_published(self):
        conn = make_connection([])
        runner, journal = self.runner(conn, {'wall': lambda: (_ for _ in ()).throw(RuntimeError('mapping failed'))})
        with self.assertRaises(RuntimeError):
            runner.run(('wall',))
        self.assertEqual(journal.data['publication']['status'], 'none')

    def test_a_clean_run_is_complete(self):
        conn = make_connection([])
        runner, journal = self.runner(conn, {'summary': lambda: conn.putdb_imas('S')})
        runner.run(('summary',))
        journal.finish()
        self.assertEqual(journal.data['publication']['status'], 'complete')
        self.assertEqual(journal.data['status'], 'PASS')


class CompletionMarkerTests(unittest.TestCase):
    """``ids stored in`` is the line scripts read as "the conversion finished"; a failed run must not print it."""

    def closed(self, completed):
        import contextlib
        import io
        conn = make_connection([])
        conn.outdir = '/entry'
        conn.ids.close = mock.Mock()
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            conn.close_imas_database(completed=completed)
        return out.getvalue()

    def test_the_marker_is_printed_only_for_a_completed_conversion(self):
        self.assertIn('ids stored in /entry', self.closed(True))
        self.assertNotIn('ids stored in', self.closed(False))

    def test_the_default_is_unchanged_for_existing_callers(self):
        self.assertIn('ids stored in', self.closed(True))


if __name__ == '__main__':
    unittest.main(verbosity=1)
