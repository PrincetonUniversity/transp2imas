"""Reader for TRANSP's formatted TF/YF/XF plotting files.

XF data is indexed as byte ranges so large native profile files are not copied
into memory until a writer asks for a particular quantity.
"""

from pathlib import Path
import re

import numpy as np
from netCDF4 import Dataset

from mapping_status import MappingLedger, RecordingDataset

_RECORD = re.compile(rb"^\s*(\d+)\s+([A-Za-z][A-Za-z0-9_]*)\s*$")
_NUMBER = re.compile(rb"[+-]?(?:(?:\d+\.?\d*)|(?:\.\d+))(?:[EeDd][+-]?\d+)?|[+-]?(?:nan|inf(?:inity)?)", re.I)


def _pln_values(raw):
    """Parse whitespace-delimited or adjacent Fortran E/D formatted numeric fields."""
    tokens = _NUMBER.findall(raw)
    residue = _NUMBER.sub(b"", raw)
    if residue.strip():
        raise ValueError(f"Unrecognized text in PLN numeric data: {residue.strip()[:40]!r}")
    # TRANSP's formatted PLN records are single-precision reals. Round-trip
    # through float32 to match the CDF reader's promotion of float32 -> float64.
    return np.asarray([float(token.replace(b"D", b"E").replace(b"d", b"e")) for token in tokens],
                      dtype=np.float32).astype(np.float64)


class _PlnVariable:
    def __init__(self, source, name, dimensions, attributes=None):
        self.source = source
        self.name = name
        self.dimensions = dimensions
        self.dtype = np.dtype("float64")
        self.attributes = dict(attributes or {})

    @property
    def shape(self):
        # Match the NetCDF variable API used by diagnostics and writers. The
        # PLN reader caches the decoded array, so later indexing is unchanged.
        return self.source.read_variable(self.name).shape

    def __getitem__(self, key):
        return self.source.read_variable(self.name)[key]

    def set_auto_mask(self, enabled):
        # Match the small part of the netCDF4.Variable API used by get_ncvar.
        return None

    def ncattrs(self):
        return list(self.attributes)

    def getncattr(self, name):
        return self.attributes[name]

    def __getattr__(self, name):
        try:
            return self.attributes[name]
        except KeyError as exc:
            raise AttributeError(name) from exc


class _Variables(dict):
    def peek(self, name):
        return dict.__getitem__(self, name)


class _Dataset:
    def __init__(self, source):
        self.source = source
        self.variables = _Variables()
        self.groups = {}

    def getncattr(self, name):
        return self.source.attributes[name]

    def ncattrs(self):
        return list(self.source.attributes)

    def close(self):
        self.source.close()

    def isopen(self):
        return self.source.file is not None


class PlnInput:
    """NetCDF-like view of a matching TF.PLN, YF.PLN and XF.PLN set."""

    def __init__(self, xf_path, machine, time_selection="all"):
        self.xf_path = Path(xf_path)
        if not self.xf_path.name.upper().endswith("XF.PLN"):
            raise ValueError("Native TRANSP input must be the run's XF.PLN file")
        self.base = self.xf_path.with_name(self.xf_path.name[:-6])
        self.tf_path = self.base.with_name(self.base.name + "TF.PLN")
        self.yf_path = self.base.with_name(self.base.name + "YF.PLN")
        for path in (self.yf_path, self.xf_path):
            if not path.is_file():
                raise FileNotFoundError(f"Required native TRANSP file not found: {path}")
        candidate_cdf = self.base.with_name(self.base.name + ".CDF")
        self.cdf_path = candidate_cdf if candidate_cdf.is_file() else None
        if not self.tf_path.is_file():
            if self.cdf_path is None:
                raise FileNotFoundError(f"Native input needs {self.tf_path.name} or a matching CDF for variable labels")
        if time_selection not in ("last", "all"):
            raise ValueError("PLN time selection must be 'last' or 'all'")
        self.time_selection = time_selection
        self.file = None
        self.mapping = MappingLedger()
        self.ncfile = str(self.xf_path)
        self.device = machine
        self.runid, self.shot = self._identity()
        self.attributes = {"device": machine, "shot": self.shot, "Runid": self.runid}
        if self.tf_path.is_file():
            self.profiles, self.axes, self.scalar_names = self._read_tf()
            self.variable_schema = {}
            if self.cdf_path is not None:
                _, _, _, self.variable_schema = self._read_cdf_schema()
        else:
            self.profiles, self.axes, self.scalar_names, self.variable_schema = self._read_cdf_schema()
        self.file = self.xf_path.open("rb")
        self.profile_index, profile_times = self._index_xf()
        self.scalar_names, scalar_values, scalar_times = self._read_yf()
        if not len(profile_times) or not len(scalar_times):
            raise ValueError("XF.PLN and YF.PLN must each contain at least one time sample")
        if time_selection == "all" and len(profile_times) != len(scalar_times):
            raise ValueError(f"XF/YF time counts differ ({len(profile_times)} and {len(scalar_times)})")
        self._selected_indices = list(range(len(profile_times))) if time_selection == "all" else [-1]
        self.timep = profile_times if time_selection == "all" else profile_times[-1:]
        self.timesc = scalar_times if time_selection == "all" else scalar_times[-1:]
        scalar_values = scalar_values if time_selection == "all" else scalar_values[-1:]
        self._cache = {"TIME": self.timesc, "TIME3": self.timep}
        self._scalar_data = {name: scalar_values[:, i] for i, name in enumerate(self.scalar_names)}
        self.variables = {
            "TIME": _PlnVariable(self, "TIME", ("TIME",)),
            "TIME3": _PlnVariable(self, "TIME3", ("TIME3",)),
        }
        for name in self.scalar_names:
            schema = self.variable_schema.get(name, (("TIME",), None, {}))
            self.variables[name] = _PlnVariable(self, name, schema[0], schema[2])
        for name, records in self.profile_index.items():
            if len({records[i][2] for i in self._selected_indices}) != 1:
                continue
            schema = self.variable_schema.get(name)
            if schema is not None and schema[0] and schema[0][0] == "TIME3":
                shape = schema[1]
                expected = int(np.prod(shape[1:]))
                if records[self._selected_indices[0]][2] != expected:
                    continue
                dimensions = schema[0]
            else:
                axis_type = self.profiles.get(name)
                if axis_type not in self.axes:
                    continue
                dimensions = ("TIME3", self.axes[axis_type])
            attributes = schema[2] if schema is not None else {}
            self.variables[name] = _PlnVariable(self, name, dimensions, attributes)
        dataset = _Dataset(self)
        dataset.variables.update(self.variables)
        self.nc = RecordingDataset(dataset, self.mapping)
        self.vars_ = list(self.variables)

    def _identity(self):
        runid = self.xf_path.name[:-6]
        match = re.match(r"^(\d+)(.*)$", runid)
        if not match:
            raise ValueError(f"Cannot derive shot and run id from {self.xf_path.name}")
        return runid, int(match.group(1))

    def _read_cdf_schema(self):
        """Use a matching CDF's declarations to label YF columns when TF.PLN is absent."""
        profiles = {}
        scalar_names = []
        schema = {}
        with Dataset(self.cdf_path, "r") as source:
            for name, variable in source.variables.items():
                upper_name = name.upper()
                dims = tuple(variable.dimensions)
                shape = tuple(variable.shape)
                attributes = {key: variable.getncattr(key) for key in variable.ncattrs()}
                schema[upper_name] = (dims, shape, attributes)
                if dims == ("TIME",):
                    scalar_names.append(upper_name)
                elif len(dims) >= 2 and dims[0] == "TIME3":
                    if dims[1] == "X":
                        profiles[upper_name] = 1
                    elif dims[1] == "XB":
                        profiles[upper_name] = 2
        if not scalar_names or scalar_names[0] != "TIME":
            raise ValueError(f"Matching CDF {self.cdf_path} has no leading TIME variable")
        return profiles, {1: "X", 2: "XB"}, scalar_names[1:], schema

    def _read_tf(self):
        lines = self.tf_path.read_text(errors="replace").splitlines()
        header = lines[0] if lines else ""
        if len(header) < 44 or header[40:].split()[-1] != "F09":
            raise ValueError("Only the current formatted F09 native TRANSP PLN layout is supported")
        nscalar, nprofile = int(header[15:20]), int(header[20:25])
        nxr = abs(int(header[25:30]))
        index = 1
        if int(header[25:30]) < 0:
            nprofile = int(lines[index].strip())
            index += 1
        index += 1
        scalars = []
        for _ in range(nscalar):
            if index + 2 >= len(lines):
                raise ValueError("Truncated TF.PLN scalar labels")
            scalars.append(lines[index].split()[0].upper())
            index += 3
        profiles = {}
        for _ in range(nprofile):
            if index + 2 >= len(lines):
                raise ValueError("Truncated TF.PLN profile labels")
            parts = lines[index].split()
            if len(parts) < 2:
                raise ValueError("Malformed TF.PLN profile descriptor")
            profiles[parts[0].upper()] = int(parts[1])
            index += 3
        # The current IMAS writers accept only TRANSP's zone-centre and zone-boundary
        # radial grids. TF.PLN profile type 1/2 identifies those grids directly.
        axes = {1: "X", 2: "XB"}
        return profiles, axes, scalars

    def _index_xf(self):
        records = {}
        times = []
        current = None
        with self.xf_path.open("rb") as stream:
            while True:
                offset = stream.tell()
                line = stream.readline()
                if not line:
                    break
                header = _RECORD.match(line.rstrip(b"\r\n"))
                if not header:
                    if current is None and line.strip():
                        raise ValueError(f"Malformed XF.PLN record at byte {offset}")
                    continue
                count = int(header.group(1))
                name = header.group(2).decode("ascii").upper()
                data_start = stream.tell()
                consumed = 0
                while consumed < count:
                    value_line = stream.readline()
                    if not value_line:
                        raise ValueError(f"Truncated XF.PLN array {name}")
                    consumed += len(_NUMBER.findall(value_line))
                if consumed != count:
                    raise ValueError(f"XF.PLN array {name} has {consumed} values; expected {count}")
                end = stream.tell()
                if name == "TIME":
                    if count != 1:
                        raise ValueError("XF.PLN TIME record must contain one value")
                    stream.seek(data_start)
                    time_values = _pln_values(stream.read(end - data_start))
                    if len(time_values) != 1:
                        raise ValueError("Malformed XF.PLN TIME record")
                    times.append(time_values[0])
                else:
                    records.setdefault(name, []).append((data_start, end, count))
                current = name
        if not times:
            raise ValueError("No TIME records found in XF.PLN")
        for name, entries in records.items():
            if len(entries) != len(times):
                raise ValueError(f"XF.PLN quantity {name} has {len(entries)} records for {len(times)} times")
        return records, np.asarray(times, dtype=np.float64)

    def _read_yf(self):
        raw = _pln_values(self.yf_path.read_bytes())
        count = len(self.scalar_names)
        width = count + 1
        if raw.size % width:
            raise ValueError(f"YF.PLN contains {raw.size} values; expected complete rows of {width}")
        values = raw.reshape((-1, width))
        return self.scalar_names, values[:, 1:], values[:, 0]

    def read_variable(self, name):
        if name in self._cache:
            return self._cache[name]
        if name in self._scalar_data:
            data = self._scalar_data[name]
        elif name in self.profile_index:
            parts = []
            for i in self._selected_indices:
                start, end, count = self.profile_index[name][i]
                self.file.seek(start)
                raw = self.file.read(end - start)
                values = _pln_values(raw)
                if values.size != count:
                    raise ValueError(f"Could not parse {name} in XF.PLN ({values.size} values; expected {count})")
                parts.append(values)
            data = np.asarray(parts, dtype=np.float64)
            schema = self.variable_schema.get(name)
            if schema is not None and schema[1] and len(schema[1]) > 2:
                data = data.reshape((len(parts),) + tuple(schema[1][1:]))
        else:
            raise KeyError(name)
        self._cache[name] = data
        return data

    def close(self):
        if self.file is not None:
            self.file.close()
            self.file = None
