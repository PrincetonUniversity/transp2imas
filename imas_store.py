#!/usr/bin/python
"""Write requests and the store that queues and flushes them.

Writers build IDS objects and *submit* them.  Whether a submission is written at once or held until
the conversion runner flushes it is decided here, not by the writer, so constructing a writer need
not write to the database.  The store has no IMAS dependency: the actual write is a callable the
connection supplies, and a request is recorded as accepted only after that callable returned.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class WriteRequest:
    """One IDS to write: the object, its occurrence and the writer that produced it."""

    ids: object
    occurrence: int
    writer: object = None          # writer name at submission time, or None


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
        """Write every queued request in submission order.

        A failing write stops the flush: the failing and the remaining requests are dropped, so
        nothing is written twice and the exception reaches the caller (the runner), which names
        the writer.  Requests accepted before the failure stay accepted.
        """
        queue, self._queue = self._queue, []
        for request in queue:
            write(request)

    def discard(self):
        """Drop queued requests (a writer failed before its output was complete)."""
        dropped = len(self._queue)
        self._queue = []
        return dropped
