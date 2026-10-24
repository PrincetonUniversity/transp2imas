#!/usr/bin/python
"""Data Dictionary version logic, separate from the database connection.

Writers ask version questions (``at_least(3, 42)``, ``is_dd4``) to choose renamed nodes.  Those
questions need only a version string, so they live here without an IMAS dependency and can be
tested, and later passed to writers, without a database.
"""

import os


def parse_dd_version(version):
    """Convert a DD version string such as 3.39.0 or 4.1.1 to a three-integer tuple."""
    values = []
    for token in str(version).split('.')[:3]:
        digits = ''.join(ch for ch in token if ch.isdigit())
        values.append(int(digits) if digits else 0)
    while len(values) < 3:
        values.append(0)
    return tuple(values)


def resolve_dd_version(requested, default_version, environ=None):
    """Pick the target DD: an explicit request, else IMAS_VERSION/IMAS_DD_VERSION, else the default.

    ``default_version`` is a callable returning the IMAS-Python default so it is only evaluated
    when nothing else selects a version.
    """
    environ = os.environ if environ is None else environ
    explicit = None if requested in (None, '', 'auto') else str(requested)
    if explicit:
        return explicit
    from_env = environ.get('IMAS_VERSION') or environ.get('IMAS_DD_VERSION')
    if from_env:
        return from_env
    return str(default_version())


class DDAdapter:
    """Version predicates for one Data Dictionary version."""

    def __init__(self, version):
        self.version = str(version)
        self.tuple = parse_dd_version(self.version)

    @property
    def major(self):
        return str(self.tuple[0])

    def at_least(self, major, minor=0, patch=0):
        return self.tuple >= (int(major), int(minor), int(patch))

    def before(self, major, minor=0, patch=0):
        return self.tuple < (int(major), int(minor), int(patch))

    @property
    def is_dd4(self):
        return self.tuple[0] >= 4

    @property
    def is_dd3(self):
        return self.tuple[0] == 3

    def __repr__(self):
        return f'DDAdapter({self.version})'
