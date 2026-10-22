#!/usr/bin/python
"""Pure-Python reader for TRANSP "AC" (ASCII compressed) files such as ``DATA<n>``.

The format is defined by the TRANSP ``aclib`` and ``i3encod`` libraries.  Each
item starts with a header line ``*<T><n> <NAME>`` (type letter, number of
dimensions), followed, for arrays, by a dimension line and encoded data lines.
Numbers use a base-64 alphabet ``0-9a-zA-Z<>``:

* ``I`` integers: sign-magnitude, ``' '`` or ``'-'`` plus base-64 digits.
* ``R`` real*4: 6 characters, IEEE-style (sign, 8-bit exponent, 27-bit mantissa).
* ``D`` real*8: 12 characters (exponent, two 30-bit decimal-mantissa halves).
* ``L`` logicals: one character per value.
* ``_0:xxx`` in real arrays encodes a run of zeros.

No TRANSP libraries are needed.
"""

import re
from pathlib import Path

import numpy as np

_ALPHABET = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ<>"
_CODE = {c: i for i, c in enumerate(_ALPHABET)}
_HEADER = re.compile(r"^\*([A-Z])(\d) (\S+)(.*)$")
_INT_TOKEN = re.compile(r"[ -][0-9a-zA-Z<>]+")


def _b64(text):
    value = 0
    for char in text:
        value = 64 * value + _CODE[char]
    return value


def _decode_int_token(token):
    value = _b64(token[1:])
    return -value if token[0] == "-" else value


_LOOKUP = np.full(256, -1, dtype=np.int64)
for _char, _value in _CODE.items():
    _LOOKUP[ord(_char)] = _value
_POW10 = {k: float(f"1e{k}") for k in range(-320, 321)}


def _decode_real4(field):
    if field[:2] == "0 ":
        return 0.0
    return float(_decode_fields(field, 6, "R")[0])


def _decode_real8(field):
    return float(_decode_fields(field, 12, "D")[0])


def _field_digits(fields, start, stop):
    """Base-64 value of characters ``start:stop`` of every fixed-width field."""
    value = np.zeros(fields.shape[0], dtype=np.int64)
    for col in range(start, stop):
        digit = _LOOKUP[fields[:, col]]
        if (digit < 0).any():
            raise ValueError("illegal character in encoded real data")
        value = value * 64 + digit
    return value


def _decode_fields(text, width, kind):
    """Vectorised decode of concatenated ``width``-character real fields.

    Zero runs (``_0:nnn``) are expanded; blank padding fields are dropped.
    """
    raw = np.frombuffer(text.encode("latin-1"), dtype=np.uint8)
    fields = raw[:raw.size // width * width].reshape(-1, width)
    blank = (fields == 32).all(axis=1)
    fields = fields[~blank]
    zero_run = (fields[:, 0] == ord("_")) & (fields[:, 1] == ord("0")) & (fields[:, 2] == ord(":"))
    plain_zero = (fields[:, 0] == ord("0")) & (fields[:, 1] == 32)
    counts = np.ones(fields.shape[0], dtype=np.int64)
    values = np.zeros(fields.shape[0])
    ok = ~(zero_run | plain_zero)
    if zero_run.any():
        counts[zero_run] = _field_digits(fields[zero_run], 3, 6)
    if ok.any():
        sel = fields[ok]
        if kind == "R":
            bits = _field_digits(sel, 0, 6)
            sign = np.where((bits >> 35) & 1, 1.0, -1.0)
            exponent = ((bits >> 27) & 0xFF).astype(np.int64)
            mantissa = (bits & ((1 << 27) - 1)).astype(np.float64)
            values[ok] = sign * np.ldexp(1.0 + mantissa / 2.0 ** 27, exponent - 127)
        else:
            exponent = _field_digits(sel, 0, 2)
            ic1 = _field_digits(sel, 2, 7).astype(np.float64)
            ic2 = _field_digits(sel, 7, 12).astype(np.float64)
            sign = np.where(exponent > 1000, 1.0, -1.0)
            exponent = np.where(exponent > 1000, exponent - 1000, exponent) - 500
            scale = np.array([_POW10[int(e)] for e in exponent])
            values[ok] = np.where(ic1 == 0, 0.0, sign * scale * (ic1 * 1.0e-9 + ic2 * 1.0e-18))
    return np.repeat(values, counts)


def _decode_reals(lines, count, width, kind):
    """Decode ``count`` reals from fixed-width fields, expanding zero runs."""
    padded = []
    for line in lines:
        if line:
            padded.append(line.ljust(-(-len(line) // width) * width))
    flat = _decode_fields("".join(padded), width, kind)
    if flat.size < count:
        raise ValueError(f"AC data ended after {flat.size} of {count} values")
    return flat[:count]


class ACItem:
    __slots__ = ("name", "kind", "ndim", "scalar_text", "lines")

    def __init__(self, name, kind, ndim, scalar_text):
        self.name = name
        self.kind = kind
        self.ndim = ndim
        self.scalar_text = scalar_text
        self.lines = []


class ACFile:
    """Lazy index of an AC file; items are decoded on first access."""

    def __init__(self, path):
        self.path = Path(path)
        self.items = {}
        self.header = ""
        current = None
        with open(self.path, "r", errors="replace") as handle:
            self.header = handle.readline().rstrip("\n")
            for raw in handle:
                line = raw.rstrip("\n")
                match = _HEADER.match(line)
                if match:
                    kind, ndim, name, rest = match.groups()
                    current = ACItem(name, kind, int(ndim), rest)
                    self.items[name] = current
                elif current is not None:
                    current.lines.append(line)

    def __contains__(self, name):
        return name in self.items

    def names(self):
        return list(self.items)

    def shape(self, name):
        """Fortran-order dimensions of an array item (empty tuple for scalars)."""
        item = self.items[name]
        if item.ndim == 0:
            return ()
        tokens = [t for t in re.split(r"[ +]+", item.lines[0].strip()) if t]
        return tuple(_b64(t) for t in tokens[:item.ndim])

    def get(self, name):
        """Return a scalar or an array reshaped Fortran-style, or ``None`` for text."""
        item = self.items[name]
        kind = item.kind
        if kind == "C":
            return None
        if item.ndim == 0:
            text = item.scalar_text
            if kind == "I":
                return _decode_int_token(" " + text.strip()) if text.strip()[0] != "-" else _decode_int_token(text.strip())
            if kind == "L":
                return text.strip().startswith("T")
            if kind == "D":
                return _decode_real8(text.strip()[:12])
            if kind == "R":
                return _decode_real4(text.strip()[:6])
            return None
        dims = self.shape(name)
        count = int(np.prod(dims))
        body = item.lines[1:]
        if kind == "D":
            flat = _decode_reals(body, count, 12, "D")
        elif kind == "R":
            flat = _decode_reals(body, count, 6, "R")
        elif kind == "I":
            values = [_decode_int_token(t) for line in body for t in _INT_TOKEN.findall(line)]
            flat = np.asarray(values[:count], dtype=np.int64)
        elif kind == "L":
            flat = np.array([c == "T" for line in body for c in line], dtype=bool)[:count]
        else:
            raise ValueError(f"Unsupported AC type {kind} for {name}")
        if flat.size != count:
            raise ValueError(f"{name}: decoded {flat.size} of {count} values")
        return flat.reshape(dims, order="F")
