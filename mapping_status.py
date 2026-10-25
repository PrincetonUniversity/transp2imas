#!/usr/bin/python
"""Structured mapping status: separate absent optional data from broken required mappings.

A TRANSP CDF legitimately lacks many variables, and the converter must then skip the matching IMAS
fields.  That is a different situation from a variable that *is* in the CDF but cannot be mapped
(wrong shape, wrong coordinate, unreadable, schema mismatch): that is a defect in the mapping or the
data and must not look like "optional data absent".  This module provides

* ``MappingError``: raised for a broken mapping (required unit, shape, coordinate or schema error),
* ``MappingLedger``: records every optional input that was absent and every field that was skipped,
  with the reason, per writer, so the conversion journal can report them.

It has no IMAS or NetCDF dependency.
"""

MISSING_OPTIONAL = 'missing_optional'     # a source variable that is not in the CDF
SKIPPED = 'skipped'                       # a field not written, with the reason
UNAVAILABLE_IN_DD = 'unavailable_in_dd'   # an IDS or node that the target Data Dictionary lacks
ADJUSTED = 'adjusted'                     # a field written after a documented adjustment of the source
KINDS = (MISSING_OPTIONAL, SKIPPED, UNAVAILABLE_IN_DD, ADJUSTED)


class MappingError(RuntimeError):
    """A source quantity that exists cannot be mapped correctly (never an absent optional input)."""

    def __init__(self, reason, *, source=None, target=None, writer=None):
        self.reason = str(reason)
        self.source = source
        self.target = target
        self.writer = writer
        parts = [self.reason]
        if source:
            parts.append(f'source {source}')
        if target:
            parts.append(f'target {target}')
        if writer:
            parts.append(f'writer {writer}')
        super().__init__('; '.join(parts))


_POLICY = {'strict': True}


def set_strict(strict):
    """Strict (default): a broken mapping raises.  Lenient: it is recorded as a skipped field."""
    _POLICY['strict'] = bool(strict)


def is_strict():
    return _POLICY['strict']


class RecordingVariables:
    """Read-only mapping view of ``Dataset.variables`` that records every variable actually read.

    Only item access and ``get`` count as a read; membership tests, iteration, ``keys`` and
    ``items`` (used to list or probe names) do not.
    """

    def __init__(self, variables, ledger):
        self._variables = variables
        self._ledger = ledger

    def __getitem__(self, name):
        variable = self._variables[name]
        self._ledger.consume(name)
        return variable

    def get(self, name, default=None):
        if name in self._variables:
            return self[name]
        return default

    def peek(self, name):
        """The variable object without recording a read: for its declaration (dimensions, attributes), not its data."""
        return self._variables[name]

    def __contains__(self, name):
        return name in self._variables

    def __iter__(self):
        return iter(self._variables)

    def __len__(self):
        return len(self._variables)

    def keys(self):
        return self._variables.keys()

    def items(self):
        return self._variables.items()

    def values(self):
        return self._variables.values()


class RecordingDataset:
    """Delegates to a NetCDF ``Dataset`` but exposes ``variables`` through ``RecordingVariables``."""

    def __init__(self, dataset, ledger):
        self._dataset = dataset
        self.variables = RecordingVariables(dataset.variables, ledger)

    def __getattr__(self, name):
        return getattr(self._dataset, name)


class MappingLedger:
    """Per-writer record of optional inputs that were absent and fields that were skipped."""

    def __init__(self):
        self.writer = None            # name of the writer currently running, set by the runner
        self._entries = {}            # (writer, kind, source, target, reason) -> count
        self._consumed = {}           # writer (None = before any writer) -> {CDF variable names read}
        self._reads = {}              # writer (None = before any writer) -> {CDF variable name: number of reads}

    def record(self, kind, *, source=None, target=None, reason='', writer=None):
        if kind not in KINDS:
            raise ValueError(f'unknown mapping status kind {kind!r}; expected one of {KINDS}')
        key = (writer if writer is not None else self.writer, kind, source, target, str(reason))
        self._entries[key] = self._entries.get(key, 0) + 1

    def missing_optional(self, source, *, target=None):
        self.record(MISSING_OPTIONAL, source=source, target=target, reason='not in the CDF')

    def skipped(self, reason, *, source=None, target=None):
        self.record(SKIPPED, source=source, target=target, reason=reason)

    def unavailable_in_dd(self, target, reason):
        self.record(UNAVAILABLE_IN_DD, target=target, reason=reason)

    def adjusted(self, reason, *, source=None, target=None):
        """The field was written, but only after the source was adjusted (for example a level clamped into range)."""
        self.record(ADJUSTED, source=source, target=target, reason=reason)

    def consume(self, source):
        """Record that the CDF variable ``source`` was actually read (by the current writer)."""
        name = str(source)
        self._consumed.setdefault(self.writer, set()).add(name)
        counts = self._reads.setdefault(self.writer, {})
        counts[name] = counts.get(name, 0) + 1

    def consumed(self, writer=None):
        """Sorted names of CDF variables read by ``writer`` (``None``: before any writer ran)."""
        return sorted(self._consumed.get(writer, ()))

    def read_summary(self, writer=None, *, repeated_over=10, listed=5):
        """``{'total', 'variables', 'most_repeated'}`` of the CDF reads of ``writer``, or ``{}`` without reads.

        ``most_repeated`` lists up to ``listed`` variables read more than ``repeated_over`` times, most
        first: where a read sits inside a loop, which is what a caching decision needs to see.
        """
        counts = self._reads.get(writer, {})
        if not counts:
            return {}
        repeated = sorted(((n, name) for name, n in counts.items() if n > repeated_over), reverse=True)
        return {'total': sum(counts.values()), 'variables': len(counts),
                'most_repeated': {name: n for n, name in repeated[:listed]}}

    def entries(self, writer=None):
        """``[{'kind', 'source', 'target', 'reason', 'count'}]`` for one writer (all if ``None``)."""
        rows = []
        for (owner, kind, source, target, reason), count in self._entries.items():
            if writer is None or owner == writer:
                rows.append({'kind': kind, 'source': source, 'target': target, 'reason': reason, 'count': count})
        rows.sort(key=lambda r: (r['kind'], str(r['source']), str(r['target']), r['reason']))
        return rows

    def summary(self, writer, *, max_listed=200):
        """Compact JSON-ready view for the conversion journal.

        Absent optional variables are listed by name (capped), skipped fields and DD gaps with
        their reasons.  Counts are always complete even when the listing is capped.
        """
        rows = self.entries(writer)
        result = {}
        for kind in KINDS:
            of_kind = [r for r in rows if r['kind'] == kind]
            if not of_kind:
                continue
            if kind == MISSING_OPTIONAL:
                names = sorted({r['source'] for r in of_kind})
                result[kind] = {'count': len(names), 'sources': names[:max_listed],
                                'truncated': len(names) > max_listed}
            else:
                listed = [{k: r[k] for k in ('source', 'target', 'reason', 'count') if r[k] is not None}
                          for r in of_kind[:max_listed]]
                result[kind] = {'count': len(of_kind), 'fields': listed,
                                'truncated': len(of_kind) > max_listed}
        return result

    def broken(self, error, *, source=None, target=None):
        """Handle a broken mapping according to the policy.

        Strict: raise ``MappingError``.  Lenient: record it as a skipped field *with the reason*
        (never silently) and let the caller carry on without the field.  ``error`` is an exception
        or a message.
        """
        reason = f'{type(error).__name__}: {error}' if isinstance(error, BaseException) else str(error)
        if is_strict():
            raise MappingError(reason, source=source, target=target, writer=self.writer) from (
                error if isinstance(error, BaseException) else None)
        self.skipped(reason, source=source, target=target)
