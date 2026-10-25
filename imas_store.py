#!/usr/bin/python
"""Write requests and the store that queues and flushes them.

Writers build IDS objects and *submit* them.  Whether a submission is written at once or held until
the conversion runner flushes it is decided here, not by the writer, so constructing a writer need
not write to the database.  The store has no IMAS dependency: the actual write is a callable the
connection supplies, and a request is recorded as accepted only after that callable returned.
"""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class WriteRequest:
    """One IDS to write: the object, its occurrence and the writer that produced it."""

    ids: object
    occurrence: int
    writer: object = None          # writer name at submission time, or None


def identity(request):
    """``{'writer', 'ids', 'occurrence'}`` of a request: the writer that produced it, the IDS name, the occurrence."""
    ids = request.ids
    name = getattr(getattr(ids, 'metadata', None), 'name', None) or str(ids)
    return {'writer': request.writer, 'ids': str(name), 'occurrence': int(request.occurrence)}


@dataclass
class FlushReport:
    """What one flush attempted.  Writes are sequential, not a transaction: ``accepted`` stay accepted.

    ``requested`` is every queued request in submission order; ``accepted`` the ones the backend took before any
    failure; ``failed`` the request that raised (with the error) or ``None``; ``unattempted`` the requests after it,
    which were dropped from the queue and never tried.  Identities are ``identity(request)``.
    """

    requested: list = field(default_factory=list)
    accepted: list = field(default_factory=list)
    failed: dict = None
    unattempted: list = field(default_factory=list)

    @property
    def complete(self):
        return self.failed is None and not self.unattempted

    def as_dict(self):
        return {'complete': self.complete, 'requested': list(self.requested), 'accepted': list(self.accepted),
                'failed': self.failed, 'unattempted': list(self.unattempted)}


class ImasStore:
    """Immediate or deferred submission of ``WriteRequest`` records."""

    def __init__(self):
        self.deferred = False
        self._queue = []

    @property
    def pending(self):
        return len(self._queue)

    def submit(self, request, write):
        """Write now, or queue when deferred.  ``write(request)`` performs and records the write."""
        if self.deferred:
            self._queue.append(request)
        else:
            write(request)

    def flush(self, write):
        """Write every queued request in submission order and return the :class:`FlushReport`.

        A failing write stops the flush: the failing and the remaining requests are dropped, so nothing is written
        twice, and the original exception reaches the caller (the runner), which names the writer.  The exception
        carries ``publication_report``, the :class:`FlushReport`, so a partial publication is reported as requested,
        accepted, failed and unattempted writes instead of being inferred from a log.  Requests accepted before the
        failure stay accepted: there is no rollback, deletion or rewrite.
        """
        queue, self._queue = self._queue, []
        report = FlushReport(requested=[identity(request) for request in queue])
        for position, request in enumerate(queue):
            try:
                write(request)
            except Exception as exc:
                report.failed = {**identity(request), 'error': f'{type(exc).__name__}: {exc}'}
                report.unattempted = [identity(rest) for rest in queue[position + 1:]]
                exc.publication_report = report
                raise
            report.accepted.append(identity(request))
        return report

    def discard(self):
        """Drop queued requests (a writer failed before its output was complete)."""
        dropped = len(self._queue)
        self._queue = []
        return dropped
