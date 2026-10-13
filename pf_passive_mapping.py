#!/usr/bin/python
"""Pure-Python/NumPy translation layer for TRANSP passive conductor outputs.

No IMAS or netCDF4 import is required. The module separates source reading,
geometry validation and IDS population so the physics and shape conventions can
be tested on normal GitHub CI runners. All geometries and electrical properties
MUST originate in a user-provided authoritative device model; none are guessed
from CDF currents, circuit labels, or plasma-current direction.

Mapping rules:
  * PS_<positive integer> is a passive-loop current, never a PF active coil.
  * CDF current values are normalized to amperes; no Ip/COCOS sign change.
  * Each loop gets its own source time base (`loop.time`).
  * The optional geometry manifest names exactly which PS_n channel drives each
    loop and explicitly records the channel-to-loop current sign.
  * Polygon outlines must be open under DD3 and closed under DD4.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
import re
from typing import Any

import numpy as np


_PS_PATTERN = re.compile(r"PS_([1-9][0-9]*)$")
_GEOMETRY_CODES = {
    "outline": 1,
    "rectangle": 2,
    "oblique": 3,
    "arcs_of_circle": 4,
    "annulus": 5,
    "thick_line": 6,
}


@dataclass(frozen=True)
class PassiveChannel:
    """One physical source-current time series as saved in a TRANSP CDF."""

    variable: str
    index: int
    name: str
    description: str
    time: np.ndarray
    current_ampere: np.ndarray
    source_units: str


@dataclass(frozen=True)
class PassiveLoop:
    """One validated IMAS loop, with optional source current and geometry."""

    name: str
    description: str
    channel: PassiveChannel | None
    current_sign: int
    resistance_ohm: float | None
    resistivity_ohm_m: float | None
    elements: tuple[dict[str, Any], ...]


@dataclass(frozen=True)
class PassivePlan:
    """Entire input-to-IDS mapping, safe to build without the IMAS runtime."""

    loops: tuple[PassiveLoop, ...]
    geometry_path: Path | None
    geometry_source: str | None
    geometry_json: str | None


def _text(value: Any) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return bytes(value).decode("ISO-8859-1", errors="replace").strip()
    return "" if value is None else str(value).strip()


def _variable_attribute(variable: Any, name: str) -> str:
    """Support both netCDF4.Variable and scipy.io.netcdf_file_variable."""
    try:
        value = variable.getncattr(name)
    except AttributeError:
        value = getattr(variable, name, None)
    return _text(value)


def _floating_array(var: Any, name: str) -> np.ndarray:
    values = np.ma.asarray(var[:], dtype=np.float64)
    if np.ma.isMaskedArray(values):
        values = values.filled(np.nan)
    try:
        return np.asarray(values, dtype=np.float64)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{name}: expected a numeric array") from exc


def _source_scale(units: str) -> float:
    """Normalize current units, preserving SI prefix case (mA != MA).

    An unqualified missing CDF unit attribute follows the established TRANSP
    PS_n amperes convention. Lowercase ``ma`` is ambiguous and is rejected.
    """
    compact = units.strip().replace(" ", "")
    lower = compact.lower()
    if lower in ("", "a", "amp", "amps", "ampere", "amperes"):
        return 1.0
    if compact in ("kA", "KA") or lower in ("kiloampere", "kiloamperes"):
        return 1.0e3
    if compact == "MA" or lower in ("megaampere", "megaamperes"):
        return 1.0e6
    if compact == "mA" or lower in ("milliampere", "milliamperes"):
        return 1.0e-3
    raise ValueError(f"Unsupported/ambiguous PS_n current units {units!r}; expected A, kA, MA or mA")


def read_passive_channels(nc: Any) -> dict[str, PassiveChannel]:
    """Read PS_n currents and their exact time bases from a CDF-like object.

    Raises on an ambiguous time base, duplicate numeric channel identifiers,
    invalid units, nonfinite data or mismatched shapes. Missing PS_n is normal
    for most TRANSP CDFs and returns an empty dictionary.
    """
    variables = nc.variables
    candidates = []
    for name in variables:
        match = _PS_PATTERN.fullmatch(name)
        if match:
            candidates.append((int(match.group(1)), name))
    candidates.sort()
    if not candidates:
        return {}

    seen_indices: set[int] = set()
    time_arrays = {}
    for name in ("TIME", "TIME3"):
        if name in variables:
            t = _floating_array(variables[name], name)
            if t.ndim != 1 or t.size == 0 or not np.all(np.isfinite(t)):
                raise ValueError(f"{name}: time base must be a finite 1D array")
            if np.any(np.diff(t) <= 0.0):
                raise ValueError(f"{name}: time base is not strictly increasing")
            time_arrays[name] = t

    out = {}
    for index, var_name in candidates:
        if index in seen_indices:
            raise ValueError(f"Multiple CDF variables identify the same PS_{index} channel")
        seen_indices.add(index)
        variable = variables[var_name]
        values = _floating_array(variable, var_name)
        if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
            raise ValueError(f"{var_name}: current must be a finite 1D array")
        dimensions = tuple(_text(x).upper() for x in getattr(variable, "dimensions", ()))
        if len(dimensions) != 1:
            raise ValueError(f"{var_name}: expected one time dimension, got {dimensions}")
        dim = dimensions[0]
        if dim in time_arrays:
            time = time_arrays[dim]
        else:
            # Old TRANSP CDF variants can give a nonstandard time dimension
            # name. Matching is permitted only when it is unambiguous.
            matches = [t for t in time_arrays.values() if len(t) == len(values)]
            if len(matches) != 1:
                raise ValueError(f"{var_name}: cannot determine time base from {dimensions}")
            time = matches[0]
        if len(time) != len(values):
            raise ValueError(
                f"{var_name}: {len(values)} currents do not match {len(time)} times ({dim})"
            )
        units = _variable_attribute(variable, "units")
        factor = _source_scale(units)
        description = _variable_attribute(variable, "long_name")
        label = re.search(r"<([^<>]+)>", description)
        display_name = label.group(1).strip() if label else var_name
        out[var_name] = PassiveChannel(
            variable=var_name,
            index=index,
            name=display_name or var_name,
            description=description,
            time=np.array(time, copy=True),
            current_ampere=np.array(values * factor, copy=True),
            source_units=units or "A (TRANSP convention; attribute absent)",
        )
    return out


def _object(value: Any, where: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{where}: expected JSON object")
    return value


def _keys(value: dict[str, Any], accepted: set[str], where: str) -> None:
    """Reject dangerous misspellings rather than silently losing machine data."""
    unexpected = set(value) - accepted
    if unexpected:
        raise ValueError(f"{where}: unsupported keys {sorted(unexpected)}")


def _positive(value: Any, where: str, *, allow_zero: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{where}: expected finite numeric value")
    number = float(value)
    if not math.isfinite(number) or (number < 0 if allow_zero else number <= 0):
        raise ValueError(f"{where}: expected {'nonnegative' if allow_zero else 'positive'} finite value")
    return number


def _finite(value: Any, where: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{where}: expected finite numeric value")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{where}: expected finite numeric value")
    return number


def _vector(values: Any, where: str) -> list[float]:
    if not isinstance(values, list) or not values:
        raise ValueError(f"{where}: expected nonempty numeric list")
    return [_finite(x, f"{where}[{i}]") for i, x in enumerate(values)]


def _point(value: Any, where: str) -> dict[str, float]:
    point = _object(value, where)
    _keys(point, {"r", "z"}, where)
    return {"r": _positive(point.get("r"), f"{where}.r"),
            "z": _finite(point.get("z"), f"{where}.z")}


def _geometry(value: Any, where: str) -> dict[str, Any]:
    geo = _object(value, where)
    shape = _text(geo.get("type"))
    allowed = {
        "rectangle": {"r", "z", "width", "height"},
        "outline": {"r", "z"},
        "oblique": {"r", "z", "length_alpha", "length_beta", "alpha", "beta"},
        "arcs_of_circle": {"r", "z", "curvature_radii"},
        "annulus": {"r", "z", "radius_inner", "radius_outer"},
        "thick_line": {"first_point", "second_point", "thickness"},
    }
    if shape in allowed:
        _keys(geo, allowed[shape] | {"type"}, where)
    if shape not in _GEOMETRY_CODES:
        raise ValueError(f"{where}.type: must be one of {', '.join(_GEOMETRY_CODES)}")
    result: dict[str, Any] = {"type": shape}
    if shape == "rectangle":
        for field in ("r", "width", "height"):
            result[field] = _positive(geo.get(field), f"{where}.{field}")
        result["z"] = _finite(geo.get("z"), f"{where}.z")
    elif shape == "outline":
        rr = _vector(geo.get("r"), f"{where}.r")
        zz = _vector(geo.get("z"), f"{where}.z")
        if len(rr) != len(zz):
            raise ValueError(f"{where}: outline R/Z arrays must have the same length")
        if len(rr) >= 2 and (rr[0], zz[0]) == (rr[-1], zz[-1]):
            rr, zz = rr[:-1], zz[:-1]
        if len(rr) < 3 or len(set(zip(rr, zz))) < 3 or any(r <= 0 for r in rr):
            raise ValueError(f"{where}: outline requires >=3 distinct points with R>0")
        area2 = sum(rr[i] * zz[(i + 1) % len(rr)] - rr[(i + 1) % len(rr)] * zz[i]
                    for i in range(len(rr)))
        if abs(area2) < 1.0e-16:
            raise ValueError(f"{where}: outline polygon is degenerate")
        result.update(r=rr, z=zz)
    elif shape == "annulus":
        result["r"] = _positive(geo.get("r"), f"{where}.r")
        result["z"] = _finite(geo.get("z"), f"{where}.z")
        result["radius_inner"] = _positive(geo.get("radius_inner"), f"{where}.radius_inner", allow_zero=True)
        result["radius_outer"] = _positive(geo.get("radius_outer"), f"{where}.radius_outer")
        if result["radius_outer"] <= result["radius_inner"]:
            raise ValueError(f"{where}: annulus outer radius must exceed inner radius")
    elif shape == "oblique":
        result["r"] = _positive(geo.get("r"), f"{where}.r")
        result["z"] = _finite(geo.get("z"), f"{where}.z")
        for field in ("length_alpha", "length_beta"):
            result[field] = _positive(geo.get(field), f"{where}.{field}")
        for field in ("alpha", "beta"):
            result[field] = _finite(geo.get(field), f"{where}.{field}")
    elif shape == "arcs_of_circle":
        for field in ("r", "z", "curvature_radii"):
            result[field] = _vector(geo.get(field), f"{where}.{field}")
        if len({len(result[k]) for k in ("r", "z", "curvature_radii")}) != 1 or len(result["r"]) < 2:
            raise ValueError(f"{where}: arcs arrays must have same length >=2")
        if any(r <= 0 for r in result["r"]) or any(x == 0 for x in result["curvature_radii"]):
            raise ValueError(f"{where}: arcs require positive R and nonzero curvature radii")
    elif shape == "thick_line":
        result["first_point"] = _point(geo.get("first_point"), f"{where}.first_point")
        result["second_point"] = _point(geo.get("second_point"), f"{where}.second_point")
        result["thickness"] = _positive(geo.get("thickness"), f"{where}.thickness")
        if result["first_point"] == result["second_point"]:
            raise ValueError(f"{where}: endpoints must differ")
    return result


def load_geometry_manifest(path: str | Path, machine: str | None = None) -> tuple[dict[str, Any], str]:
    """Load and strictly validate *authoritative* machine passive geometry.

    The manifest uses metres/radians/amperes/ohms. Device-specific geometry
    cannot be reconstructed from a current channel's TRANSP long_name.
    """
    source_path = Path(path).expanduser()
    content = source_path.read_text(encoding="utf-8")
    try:
        raw = _object(json.loads(content), "geometry manifest")
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid passive geometry JSON: {exc}") from exc
    _keys(raw, {"schema_version", "machine", "source_reference", "units", "loops"},
          "geometry manifest")
    if raw.get("schema_version") != 1 or isinstance(raw.get("schema_version"), bool):
        raise ValueError("geometry manifest.schema_version must be integer 1")
    device = _text(raw.get("machine"))
    if not device:
        raise ValueError("geometry manifest.machine is required")
    if machine and device.casefold() != _text(machine).casefold():
        raise ValueError(f"Geometry machine {device!r} differs from CDF device {machine!r}")
    source = _text(raw.get("source_reference"))
    if not source:
        raise ValueError("geometry manifest.source_reference must identify an authoritative device source")
    if raw.get("units") != "m":
        raise ValueError("geometry manifest.units must be 'm' (angles in radians)")
    loops = raw.get("loops")
    if not isinstance(loops, list) or not loops:
        raise ValueError("geometry manifest.loops must be a nonempty list")

    used_names = set()
    used_channels = set()
    parsed_loops = []
    for index, item in enumerate(loops):
        where = f"geometry manifest.loops[{index}]"
        item = _object(item, where)
        _keys(item, {"name", "description", "cdf_variable", "cdf_current_sign",
                     "resistance_ohm", "resistivity_ohm_m", "elements"}, where)
        name = _text(item.get("name"))
        if not name or name in used_names:
            raise ValueError(f"{where}.name must be unique and nonempty")
        used_names.add(name)
        cdf_var = item.get("cdf_variable")
        if cdf_var is not None:
            if not isinstance(cdf_var, str) or not _PS_PATTERN.fullmatch(cdf_var):
                raise ValueError(f"{where}.cdf_variable must be PS_<positive index>")
            if cdf_var in used_channels:
                raise ValueError(f"{where}: duplicate mapping for {cdf_var}")
            used_channels.add(cdf_var)
            current_sign = item.get("cdf_current_sign")
            if type(current_sign) is not int or current_sign not in (-1, 1):
                raise ValueError(f"{where}.cdf_current_sign must explicitly be +1 or -1")
        else:
            current_sign = 1  # no current to transform
            if "cdf_current_sign" in item:
                raise ValueError(f"{where}: cdf_current_sign without cdf_variable")
        elements = item.get("elements")
        if not isinstance(elements, list) or not elements:
            raise ValueError(f"{where}.elements must be nonempty")
        parsed_elements = []
        for j, element in enumerate(elements):
            ewhere = f"{where}.elements[{j}]"
            element = _object(element, ewhere)
            _keys(element, {"name", "description", "identifier", "turns_with_sign",
                            "area_m2", "geometry"}, ewhere)
            e_name = _text(element.get("name"))
            if not e_name:
                raise ValueError(f"{ewhere}.name is required")
            turns = _finite(element.get("turns_with_sign"), f"{ewhere}.turns_with_sign")
            if turns == 0:
                raise ValueError(f"{ewhere}: turns_with_sign must be nonzero")
            record = {
                "name": e_name,
                "description": _text(element.get("description")),
                "identifier": _text(element.get("identifier")),
                "turns_with_sign": turns,
                "area_m2": None,
                "geometry": _geometry(element.get("geometry"), f"{ewhere}.geometry"),
            }
            if element.get("area_m2") is not None:
                record["area_m2"] = _positive(element["area_m2"], f"{ewhere}.area_m2")
            parsed_elements.append(record)
        resistance = item.get("resistance_ohm")
        resistivity = item.get("resistivity_ohm_m")
        parsed_loops.append({
            "name": name,
            "description": _text(item.get("description")),
            "cdf_variable": cdf_var,
            "cdf_current_sign": current_sign,
            "resistance_ohm": None if resistance is None else _positive(resistance, f"{where}.resistance_ohm", allow_zero=True),
            "resistivity_ohm_m": None if resistivity is None else _positive(resistivity, f"{where}.resistivity_ohm_m", allow_zero=True),
            "elements": tuple(parsed_elements),
        })
    return {"machine": device, "source_reference": source, "loops": parsed_loops}, content


def build_passive_plan(
    nc: Any, *, geometry_path: str | Path | None = None, machine: str | None = None
) -> PassivePlan:
    """Build a verified loop/channel mapping; missing channels stay current-empty.

    Without geometry input, retain the legacy current-only export, but do not
    fabricate geometry, resistance, turns or connectivity. With geometry input,
    every present PS_n channel must have an explicit loop mapping. The manifest
    may additionally provide machine loops for which the CDF has no current.
    """
    channels = read_passive_channels(nc)
    if geometry_path is None:
        used_names = set()
        loops = []
        for channel in channels.values():
            name = channel.name
            if name in used_names:
                name = channel.variable
            used_names.add(name)
            loops.append(PassiveLoop(name=name, description=channel.description,
                                     channel=channel, current_sign=1,
                                     resistance_ohm=None, resistivity_ohm_m=None,
                                     elements=()))
        return PassivePlan(tuple(loops), None, None, None)

    geometry, content = load_geometry_manifest(geometry_path, machine=machine)
    listed = {loop["cdf_variable"] for loop in geometry["loops"] if loop["cdf_variable"]}
    unmatched = set(channels) - listed
    if unmatched:
        raise ValueError(
            "The passive geometry manifest does not map present CDF channels: "
            + ", ".join(sorted(unmatched))
        )
    loops = []
    for item in geometry["loops"]:
        channel = channels.get(item["cdf_variable"])
        loops.append(PassiveLoop(
            name=item["name"], description=item["description"],
            channel=channel, current_sign=item["cdf_current_sign"],
            resistance_ohm=item["resistance_ohm"],
            resistivity_ohm_m=item["resistivity_ohm_m"],
            elements=item["elements"],
        ))
    return PassivePlan(tuple(loops), Path(geometry_path).expanduser(),
                       geometry["source_reference"], content)


def populate_element(element: Any, spec: dict[str, Any], dd_major: int) -> None:
    """Write a validated physical conductor element into a DD3/DD4 IDS node."""
    element.name = spec["name"]
    if hasattr(element, "description") and spec["description"]:
        element.description = spec["description"]
    if hasattr(element, "identifier") and spec["identifier"]:
        element.identifier = spec["identifier"]
    element.turns_with_sign = spec["turns_with_sign"]
    if spec["area_m2"] is not None:
        element.area = spec["area_m2"]
    geometry = spec["geometry"]
    shape = geometry["type"]
    container = element.geometry
    container.geometry_type = _GEOMETRY_CODES[shape]
    if shape == "outline":
        rr = list(geometry["r"])
        zz = list(geometry["z"])
        if int(dd_major) >= 4:
            rr.append(rr[0])
            zz.append(zz[0])
        container.outline.r = np.asarray(rr, dtype=np.float64)
        container.outline.z = np.asarray(zz, dtype=np.float64)
    elif shape == "thick_line":
        target = container.thick_line
        target.first_point.r = geometry["first_point"]["r"]
        target.first_point.z = geometry["first_point"]["z"]
        target.second_point.r = geometry["second_point"]["r"]
        target.second_point.z = geometry["second_point"]["z"]
        target.thickness = geometry["thickness"]
    else:
        target = getattr(container, shape)
        for key, value in geometry.items():
            if key == "type":
                continue
            if isinstance(value, list):
                value = np.asarray(value, dtype=np.float64)
            setattr(target, key, value)


def populate_passive_ids(ids: Any, plan: PassivePlan, dd_major: int) -> int:
    """Fill passive loop nodes only; return the number of current channels.

    Create the full loop AoS *before* populating any element to avoid the
    IMAS-Python resize-recreates-AoS-elements pitfall.
    """
    ids.loop.resize(len(plan.loops))
    ncurrents = 0
    for index, spec in enumerate(plan.loops):
        loop = ids.loop[index]
        loop.name = spec.name
        if spec.description and hasattr(loop, "description"):
            loop.description = spec.description
        if spec.resistance_ohm is not None:
            loop.resistance = spec.resistance_ohm
        if spec.resistivity_ohm_m is not None:
            loop.resistivity = spec.resistivity_ohm_m
        loop.element.resize(len(spec.elements))
        for j, element in enumerate(spec.elements):
            populate_element(loop.element[j], element, dd_major)
        if spec.channel is not None:
            loop.time = np.array(spec.channel.time, dtype=np.float64, copy=True)
            loop.current = np.asarray(
                spec.current_sign * spec.channel.current_ampere, dtype=np.float64
            )
            ncurrents += 1
    return ncurrents
