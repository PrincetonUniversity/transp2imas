#!/usr/bin/python
"""Declarative IDS writer selection and durable conversion-write journaling.

This module deliberately has no IMAS or NetCDF dependency, allowing CI to
verify selection, journaling, and failure semantics without a site IMAS setup.
"""

import json
from pathlib import Path
import time


def _json_default(value):
    """Encode NumPy metadata from CDF/IMAS without discarding numeric types.

    The journal intentionally does not depend on NumPy itself.  NumPy scalars
    and arrays expose ``tolist()``, returning native JSON numbers and lists.
    Unexpected objects remain errors rather than silently becoming strings.
    """
    if isinstance(value, Path):
        return str(value)
    tolist = getattr(value, 'tolist', None)
    if callable(tolist):
        native = tolist()
        if native is not value:
            return native
    raise TypeError(f'Cannot serialize conversion journal value of type '
                    f'{type(value).__name__}')


WRITER_ORDER = (
    'workflow', 'wall', 'tf', 'lh_antennas', 'ic_antennas', 'ec_launchers',
    'nbi', 'distribution_sources', 'distributions', 'pf_active', 'core_profiles',
    'pf_passive', 'radiation', 'equilibrium', 'summary', 'core_sources',
    'waves', 'core_transport', 'transport_solver_numerics', 'sawteeth', 'edge_sources', 'dataset_fair',
)


# Backward-compatible aliases for legacy converter progress labels.
IDS_ALIASES = {'lh_antenna': 'lh_antennas', 'ic_antenna': 'ic_antennas'}

def select_ids(argument):
    """Return IDS names in dependency-safe writer order; fail on typos."""
    if argument is None or str(argument).strip().lower() == 'all':
        return WRITER_ORDER
    parts = [IDS_ALIASES.get(name.strip().lower(), name.strip().lower())
             for name in str(argument).split(',')]
    if not parts or any(not name for name in parts):
        raise ValueError('--ids requires a comma-separated list of IDS names')
    unknown = sorted(set(parts) - set(WRITER_ORDER))
    if unknown:
        raise ValueError('Unknown --ids name(s): ' + ', '.join(unknown) +
                         '; available: ' + ', '.join(WRITER_ORDER))
    return tuple(name for name in WRITER_ORDER if name in set(parts))


class ConversionRunner:
    """Run the selected writers in order, flush each writer's IDSs, and journal the outcome.

    The converter is expected to expose ``defer_writes``, ``flush_writes``, ``discard_writes``,
    ``mapping`` (a ``MappingLedger``) and the ``_written_ids`` manifest and ``_active_writer``
    attribution of the IMAS connection.  Writers are zero-argument factories that build their IDS
    objects and submit them; with deferral on, nothing reaches the database until the runner flushes
    after that writer returned. A build failure submits none of that writer's output. A backend
    failure during flush may leave already accepted IDSs in the entry: no rollback is provided.
    Accepted writes and the failure are journaled with the writer's name; remaining writers do not run.
    """

    def __init__(self, converter, writers, journal, clock=time.monotonic):
        self.converter = converter
        self.writers = writers
        self.journal = journal
        self.clock = clock

    def run(self, selected):
        nci = self.converter
        nci.defer_writes(True)
        try:
            for name in selected:
                print(name)
                before = len(nci._written_ids)
                started = self.clock()
                nci._active_writer = name
                nci.mapping.writer = name
                write_seconds = 0.0
                try:
                    self.writers[name]()
                    flush_started = self.clock()
                    try:
                        nci.flush_writes()
                    finally:
                        write_seconds = self.clock() - flush_started
                except Exception as exc:
                    nci.discard_writes()
                    report = getattr(exc, 'publication_report', None)      # set by ImasStore.flush on a backend failure
                    self.journal.record_writer(name, nci._written_ids[before:],
                                               self.clock() - started, error=exc,
                                               mapping=nci.mapping.summary(name),
                                               consumed=nci.mapping.consumed(name),
                                               reads=nci.mapping.read_summary(name),
                                               write_seconds=write_seconds,
                                               publication=None if report is None else report.as_dict())
                    self.journal.note_failure(exc)
                    raise
                else:
                    if len(nci._written_ids) == before and not nci.mapping.entries(name):
                        # A writer that wrote nothing and recorded no reason of its own.
                        nci.mapping.skipped('no IDS written; the writer logged its reason in the conversion log',
                                            target=name)
                    self.journal.record_writer(name, nci._written_ids[before:],
                                               self.clock() - started,
                                               mapping=nci.mapping.summary(name),
                                               consumed=nci.mapping.consumed(name),
                                               reads=nci.mapping.read_summary(name),
                                               write_seconds=write_seconds)
                finally:
                    nci._active_writer = None
        finally:
            nci.defer_writes(False)


class ConversionJournal:
    """Record what IMAS actually accepted, not merely which writers were called."""

    def __init__(self, source, selected):
        self.data = {
            'format_version': 1,
            'source_cdf': str(source),
            'selected_ids': list(selected),
            'machine': None,
            'shot': None,
            'run': None,
            'dd_version': None,
            'backend': None,
            'occurrence': None,
            'output_directory': None,
            'mapping_policy': None,
            'writers': [],
            'written_ids': [],
            'status': 'RUNNING',
            'error': None,
        }

    def attach_converter(self, converter):
        for field, attr in (
                ('machine', 'imas_dbname'), ('shot', 'shot'), ('run', 'run'),
                ('dd_version', 'dd_version'), ('occurrence', 'occur'),
                ('output_directory', 'outdir')):
            value = getattr(converter, attr)
            if value is not None:
                # NetCDF exposes shot/run/occurrence as NumPy int32/int64.
                # Convert them before adding them to the public JSON journal.
                value = (int(value) if field in ('shot', 'run', 'occurrence')
                         else str(value))
            self.data[field] = value
        backend = getattr(converter, 'backend', None)
        # Internal IMAS enum values are integers; normalize the display name.
        self.data['backend'] = (None if backend is None
                                else str(getattr(backend, 'name', backend)))

    def record_writer(self, name, writes, duration, error=None, mapping=None, consumed=None,
                      write_seconds=None, reads=None, publication=None):
        # Preserve only the stable readback contract, with native JSON types.
        # IMAS and NetCDF may supply NumPy scalar IDs and occurrence numbers.
        events = [{'ids': str(event['ids']),
                   'occurrence': int(event['occurrence'])}
                  for event in writes]
        status = 'FAIL' if error else ('WROTE' if events else 'SKIP')
        record = {'ids': name, 'status': status, 'write_count': len(events),
                  'seconds': round(float(duration), 6),
                  'occurrences': [event['occurrence'] for event in events]}
        if error:
            record['error'] = f'{type(error).__name__}: {error}'
        if mapping:
            # Absent optional inputs and skipped fields with reasons (see mapping_status.MappingLedger).
            record['mapping'] = mapping
        if write_seconds is not None:
            # Time in the backend (IMAS validation, put and compression) as part of 'seconds'; the rest
            # is reading, transforming and assigning.
            record['write_seconds'] = round(float(write_seconds), 6)
        if consumed:
            # CDF variables this writer actually read (measured, not inferred from the source text).
            record['consumed'] = sorted(consumed)
        if reads:
            # How often, not only which: total reads, distinct variables and the most repeated ones.
            record['reads'] = reads
        if publication is not None:
            # A backend failure during this writer's flush: requested, accepted, failed and unattempted writes.
            record['publication'] = publication
        self.data['writers'].append(record)
        self.data['written_ids'].extend(events)
        return status

    def _accepted(self):
        return [dict(event) for event in self.data['written_ids']]

    def note_failure(self, error):
        """Record what was published when ``error`` ended the run (see :meth:`fail`), without changing the status.

        ``publication.status`` is ``none`` (nothing was accepted), or ``partial``: writes are sequential, not a
        transaction, so IDSs accepted before the failure remain in the database.  ``kind`` separates a
        ``backend_failure`` (a put failed during a flush; the failing and unattempted requests are named) from a
        ``build_failure`` (a writer failed before its flush; none of its output was submitted).  Nothing is rolled
        back or deleted.
        """
        accepted = self._accepted()
        report = getattr(error, 'publication_report', None)
        failed = unattempted = None
        if report is not None:
            failed = {'ids': report.failed.get('writer') or report.failed['ids'],
                      'occurrence': report.failed['occurrence'], 'error': report.failed['error']}
            unattempted = [{'ids': item.get('writer') or item['ids'], 'occurrence': item['occurrence']}
                           for item in report.unattempted]
        self.data['publication'] = {'status': 'partial' if accepted else 'none',
                                    'kind': 'backend_failure' if report is not None else 'build_failure',
                                    'accepted': accepted, 'failed': failed, 'unattempted': unattempted or []}

    def fail(self, error):
        self.data['status'] = 'FAIL'
        self.data['error'] = f'{type(error).__name__}: {error}'
        self.note_failure(error)

    def finish(self):
        if self.data['status'] == 'RUNNING':
            self.data['status'] = 'PASS'
            self.data['publication'] = {'status': 'complete', 'kind': None, 'accepted': self._accepted(),
                                        'failed': None, 'unattempted': []}

    def save(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + '.tmp')
        try:
            tmp.write_text(json.dumps(self.data, indent=2, sort_keys=True,
                                      default=_json_default, allow_nan=False) + '\n',
                           encoding='utf-8')
            tmp.replace(path)
        finally:
            tmp.unlink(missing_ok=True)
