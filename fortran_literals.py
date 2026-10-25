#!/usr/bin/python
"""Fortran scalar literals, as one contract for the converter's namelist reader and the curation package.

Pure Python with no dependencies, so standalone curation does not need the converter stack.  What is supported
is exactly one scalar per value:

* integers and reals, signed, with ``D``/``d`` or ``E``/``e`` exponents (``1.0D-3``, ``-2.5d+2``, ``.5``, ``-4.``),
* logicals (``.TRUE.``, ``.T.``, ``T``, ``.FALSE.``, ``.F.``, ``F`` and the spelled-out forms), any case,
* quoted strings with doubled quotes.

Everything else (arrays, expressions, repeated assignments, includes, ``nan``/``inf``, a number with trailing text)
is *not a supported literal*: the parsers return ``None`` for it.  ``None`` means "unknown", and callers keep the raw
text next to it.  A value that was supplied but cannot be read must never become ``0``, ``False`` or the value of an
earlier assignment.  There is no ``eval``.  This is not a Fortran parser: it does not resolve defaults, includes or
runtime updates.
"""

import math
import re

REAL = re.compile(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[dDeE][+-]?\d+)?\Z')
INTEGER = re.compile(r'[+-]?\d+\Z')
_D_EXPONENT = re.compile(r'(?<=[0-9.])[dD](?=[+-]?\d)')
TRUE = ('.TRUE.', '.T.', 'T', 'TRUE')
FALSE = ('.FALSE.', '.F.', 'F', 'FALSE')


def normalize_exponents(text):
    """``1.0D-3`` -> ``1.0E-3``: Fortran ``D`` exponents as the ``E`` that Python reads (numeric tokens only)."""
    return _D_EXPONENT.sub('E', text)


def parse_number(text):
    """The value of one numeric literal as ``int`` or ``float``, or ``None`` (not a literal, or not finite)."""
    token = text.strip()
    if INTEGER.match(token):
        return int(token)
    if REAL.match(token):
        value = float(normalize_exponents(token))
        return value if math.isfinite(value) else None
    return None


def parse_logical(text):
    """``True``/``False`` for a Fortran logical literal, or ``None``."""
    token = text.strip().upper()
    if token in TRUE:
        return True
    if token in FALSE:
        return False
    return None


def parse_string(text):
    """The content of a quoted string (doubled quotes undone), or ``None``."""
    token = text.strip()
    if len(token) >= 2 and token[0] in '\'"' and token[-1] == token[0]:
        return token[1:-1].replace(token[0] * 2, token[0])
    return None


def parse_scalar(text):
    """One supported scalar (``bool``, ``int``, ``float`` or ``str``), or ``None`` for anything else.

    A trailing comma is the list separator of a one-element list and is ignored; any other trailing text makes the
    value unsupported.
    """
    token = text.strip().rstrip(',').strip()
    for parse in (parse_logical, parse_number, parse_string):
        value = parse(token)
        if value is not None:
            return value
    return None
