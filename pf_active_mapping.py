#!/usr/bin/python
"""Validated TRANSP PF coil/supply channels and optional machine configuration.

This source-only module intentionally does not import IMAS.  The CDF provides
waveforms (KK_, VV_, VS_, and optionally CC_), but not reliable coil cross-
sections, effective turns, circuit connections, or current polarity relative
to the IMAS winding convention.  Those machine properties come exclusively
from an authoritative JSON manifest.

KT_ is a *target/reference* circuit current, NOT a measured supply current.
It is inventoried but deliberately not written to pf_active.supply.current.
In verified ISolver cases CC_<coil> is the corresponding coil ampere-turns
(CC/KK equals the effective turns), NOT another physical circuit current.
Binding CC_ to an IMAS coil current therefore requires an explicit source-backed
``cc_effective_turns`` divisor in the coil geometry manifest.
The future destination for target currents is pulse_schedule.

The mapping preserves native per-channel TIME/TIME3 coordinates, converts
SI units, and supports DD3/DD4 cross-section and circuit representations.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any

import numpy as np

# Reuse the geometry validation and DD-specific writing that is already
# regression-tested for pf_passive; both IDSs use the same element shapes.
from pf_passive_mapping import (
    _finite, _geometry, _keys, _object, _positive, _source_scale, _text,
    _variable_attribute, _floating_array, populate_element,
)

_CHANNEL_PATTERN = re.compile(r"^(KK|VV|KT|VS|CC)_([A-Za-z0-9][A-Za-z0-9_]*)$")
_CONTROLLER_PATTERN = re.compile(r"^(KD|KS)_([A-Za-z0-9][A-Za-z0-9_]*)$")
_CURRENT_PREFIXES = frozenset({"KK", "KT", "CC"})
_VOLTAGE_PREFIXES = frozenset({"VV", "VS"})


@dataclass(frozen=True)
class ActiveChannel:
    """Single time-dependent TRANSP CDF channel normalized into A or V."""

    variable: str
    prefix: str
    suffix: str
    time: np.ndarray
    values: np.ndarray
    description: str
    original_units: str


@dataclass(frozen=True)
class BoundChannel:
    """A CDF waveform and explicit (or unverified legacy) orientation sign."""

    channel: ActiveChannel
    sign: int
    factor: float = 1.0


@dataclass(frozen=True)
class CoilSpec:
    name: str
    description: str
    current: BoundChannel | None
    voltage: BoundChannel | None
    resistance_ohm: float | None
    elements: tuple[dict[str, Any], ...]


@dataclass(frozen=True)
class SupplySpec:
    name: str
    description: str
    current: BoundChannel | None
    voltage: BoundChannel | None


@dataclass(frozen=True)
class CircuitSpec:
    name: str
    description: str
    connections_dd3: np.ndarray | None
    connections_dd4: np.ndarray | None
    current: BoundChannel | None
    voltage: BoundChannel | None


@dataclass(frozen=True)
class ActivePlan:
    coils: tuple[CoilSpec, ...]
    supplies: tuple[SupplySpec, ...]
    circuits: tuple[CircuitSpec, ...]
    ignored_targets: tuple[str, ...]
    ignored_equilibrium: tuple[str, ...]
    manifest_path: Path | None
    manifest_source: str | None
    manifest_text: str | None
    ignored_controller: tuple[str, ...] = ()  # KD desired-current / KS controller-sigma


def _voltage_scale(units: str) -> float:
    compact = units.strip().replace(" ", "")
    lower = compact.lower()
    if lower in ("", "v", "volt", "volts"):
        return 1.0
    if compact in ("kV", "KV") or lower in ("kilovolt", "kilovolts"):
        return 1.0e3
    if compact == "MV" or lower in ("megavolt", "megavolts"):
        return 1.0e6
    if compact == "mV" or lower in ("millivolt", "millivolts"):
        return 1.0e-3
    raise ValueError(f"Unsupported/ambiguous PF voltage unit {units!r}; expected V, kV, MV or mV")


def read_active_channels(nc: Any, *, include_equilibrium: bool = False) -> dict[str, ActiveChannel]:
    """Find KK/VV/KT/VS/CC channels with strict unit/grid validation.

    A channel must be a finite 1D array on TIME, TIME3 or a uniquely
    length-matched legacy time coordinate.  We do *not* time-interpolate
    unrelated coil and supply traces.  Missing unit attributes follow the
    legacy TRANSP assumption of amperes for currents and volts for voltages.
    """
    variables = nc.variables
    selected = [(name, _CHANNEL_PATTERN.fullmatch(name)) for name in variables]
    selected = [(name, match) for name, match in selected if match and
                match.group(1) != "KT" and
                (include_equilibrium or match.group(1) != "CC")]
    if not selected:
        return {}
    times: dict[str, np.ndarray] = {}
    for name in ("TIME", "TIME3"):
        if name in variables:
            arr = _floating_array(variables[name], name)
            if arr.ndim != 1 or not len(arr) or not np.all(np.isfinite(arr)):
                raise ValueError(f"{name}: expected finite nonempty one-dimensional time")
            if np.any(np.diff(arr) <= 0):
                raise ValueError(f"{name}: time must be strictly increasing")
            times[name] = np.array(arr, copy=True)

    channels: dict[str, ActiveChannel] = {}
    for name, match in sorted(selected):
        prefix, suffix = match.groups()
        variable = variables[name]
        arr = _floating_array(variable, name)
        if arr.ndim != 1 or not len(arr) or not np.all(np.isfinite(arr)):
            raise ValueError(f"{name}: expected finite nonempty 1D waveform")
        dimensions = tuple(_text(d).upper() for d in getattr(variable, "dimensions", ()))
        if len(dimensions) != 1:
            raise ValueError(f"{name}: must have exactly one time dimension, got {dimensions}")
        if dimensions[0] in times:
            time = times[dimensions[0]]
        else:
            choices = [t for t in times.values() if len(t) == len(arr)]
            unique = []
            for t in choices:
                if not any(np.array_equal(t, prior) for prior in unique):
                    unique.append(t)
            if len(unique) != 1:
                raise ValueError(f"{name}: ambiguous/unavailable CDF time dimension {dimensions}")
            time = unique[0]
        if len(arr) != len(time):
            raise ValueError(f"{name}: {len(arr)} values but {len(time)} time coordinates")
        units = _variable_attribute(variable, "units")
        factor = _source_scale(units) if prefix in _CURRENT_PREFIXES else _voltage_scale(units)
        channels[name] = ActiveChannel(
            variable=name, prefix=prefix, suffix=suffix,
            time=np.asarray(time, dtype=np.float64).copy(),
            values=np.asarray(arr * factor, dtype=np.float64),
            description=_variable_attribute(variable, "long_name"),
            original_units=units or ("A (assumed)" if prefix in _CURRENT_PREFIXES else "V (assumed)"),
        )
    return channels


def _name_list(data: list[dict], where: str) -> None:
    names = [entry["name"] for entry in data]
    if len(names) != len(set(names)):
        raise ValueError(f"{where}: duplicate machine identifiers")


def _parse_elements(items: Any, where: str) -> tuple[dict[str, Any], ...]:
    if not isinstance(items, list) or not items:
        raise ValueError(f"{where}: coils require one or more sourced cross-section elements")
    out = []
    for index, item in enumerate(items):
        item_where = f"{where}[{index}]"
        item = _object(item, item_where)
        _keys(item, {"name", "description", "identifier", "turns_with_sign", "area_m2", "geometry"}, item_where)
        name = _text(item.get("name"))
        if not name:
            raise ValueError(f"{item_where}.name is required")
        turns = _finite(item.get("turns_with_sign"), f"{item_where}.turns_with_sign")
        if turns == 0:
            raise ValueError(f"{item_where}: turns_with_sign cannot be zero")
        area = item.get("area_m2")
        out.append({
            "name": name,
            "description": _text(item.get("description")),
            "identifier": _text(item.get("identifier")),
            "turns_with_sign": turns,
            "area_m2": None if area is None else _positive(area, f"{item_where}.area_m2"),
            "geometry": _geometry(item.get("geometry"), f"{item_where}.geometry"),
        })
    return tuple(out)


def _channel_binding(
    item: dict[str, Any], stem: str, allowed: set[str],
    channels: dict[str, ActiveChannel], where: str,
    *, cc_effective_turns: float | None = None,
) -> BoundChannel | None:
    var_name = item.get(f"cdf_{stem}")
    sign_name = f"{stem}_sign"
    if var_name is None:
        if item.get(sign_name) is not None:
            raise ValueError(f"{where}: {sign_name} provided without cdf_{stem}")
        if cc_effective_turns is not None:
            raise ValueError(f"{where}: cc_effective_turns provided without cdf_current")
        return None
    if not isinstance(var_name, str) or not _CHANNEL_PATTERN.fullmatch(var_name):
        raise ValueError(f"{where}.cdf_{stem}: invalid channel name {var_name!r}")
    prefix = var_name.split("_", 1)[0]
    if prefix not in allowed:
        raise ValueError(f"{where}.cdf_{stem}={var_name}: expected one of {sorted(allowed)}")
    sign = item.get(sign_name)
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError(f"{where}.{sign_name} must explicitly be +1 or -1")
    if prefix == "CC":
        # CC_ is the equilibrium coil ampere-turns, not KK_ circuit current.
        # Never infer physical turns from this particular CDF alone: use the
        # approved machine model to establish the conversion and sign.
        if cc_effective_turns is None:
            raise ValueError(f"{where}.cdf_{stem}={var_name}: "
                             "CC_ requires authoritative cc_effective_turns")
        turns = _finite(cc_effective_turns, f"{where}.cc_effective_turns")
        if turns == 0:
            raise ValueError(f"{where}.cc_effective_turns must be nonzero")
        factor = 1.0 / turns
    else:
        if cc_effective_turns is not None:
            raise ValueError(f"{where}: cc_effective_turns is only for CC_ coil current")
        factor = 1.0
    channel = channels.get(var_name)
    return BoundChannel(channel, sign, factor) if channel is not None else None


def _parse_connections(value: Any, ncols: int, dd: int, where: str) -> np.ndarray | None:
    if value is None:
        return None
    if not isinstance(value, list) or len(value) < 2 or not all(isinstance(row, list) for row in value):
        raise ValueError(f"{where}: expected at least two matrix rows")
    if any(len(row) != ncols for row in value):
        raise ValueError(f"{where}: expected {ncols} columns (DD{dd})")
    allowed = {0, 1} if dd == 3 else {-1, 0, 1}
    if any(type(v) is not int or v not in allowed for row in value for v in row):
        raise ValueError(f"{where}: invalid matrix entry for DD{dd}; allowed {sorted(allowed)}")
    return np.asarray(value, dtype=np.int32)


def load_active_manifest(path: str | Path, machine: str | None = None) -> tuple[dict[str, Any], str]:
    """Validate a sourced machine description before binding any CDF traces."""
    raw_text = Path(path).expanduser().read_text(encoding="utf-8")
    try:
        root = _object(json.loads(raw_text), "pf_active manifest")
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid pf_active geometry JSON: {exc}") from exc
    _keys(root, {"schema_version", "machine", "source_reference", "units", "coils", "supplies", "circuits"}, "pf_active manifest")
    if type(root.get("schema_version")) is not int or root["schema_version"] != 1:
        raise ValueError("pf_active manifest.schema_version must be 1")
    machine_name = _text(root.get("machine"))
    if not machine_name:
        raise ValueError("pf_active manifest.machine is required")
    if machine is not None and _text(machine).casefold() != machine_name.casefold():
        raise ValueError(f"PF active geometry machine {machine_name!r} differs from CDF device {machine!r}")
    source = _text(root.get("source_reference"))
    if not source:
        raise ValueError("pf_active manifest.source_reference must identify an authoritative machine source")
    if root.get("units") != "m":
        raise ValueError("pf_active manifest.units must be 'm' (angles in radians)")
    coil_items = root.get("coils", [])
    supply_items = root.get("supplies", [])
    circuit_items = root.get("circuits", [])
    if not all(isinstance(x, list) for x in (coil_items, supply_items, circuit_items)):
        raise ValueError("pf_active manifest: coils, supplies, circuits must be JSON arrays")
    if not coil_items and not supply_items:
        raise ValueError("pf_active manifest needs at least one coil or supply")
    coils = []
    for i, item in enumerate(coil_items):
        where = f"pf_active manifest.coils[{i}]"
        item = _object(item, where)
        _keys(item, {"name", "description", "cdf_current", "cdf_voltage", "current_sign", "voltage_sign", "resistance_ohm", "elements", "cc_effective_turns"}, where)
        name = _text(item.get("name"))
        if not name:
            raise ValueError(f"{where}.name is required")
        resistance = item.get("resistance_ohm")
        coils.append({
            "name": name, "description": _text(item.get("description")),
            "cdf_current": item.get("cdf_current"), "cdf_voltage": item.get("cdf_voltage"),
            "current_sign": item.get("current_sign"), "voltage_sign": item.get("voltage_sign"),
            "cc_effective_turns": item.get("cc_effective_turns"),
            "resistance_ohm": None if resistance is None else _positive(resistance, f"{where}.resistance_ohm", allow_zero=True),
            "elements": _parse_elements(item.get("elements"), f"{where}.elements"),
        })
    supplies = []
    for i, item in enumerate(supply_items):
        where = f"pf_active manifest.supplies[{i}]"
        item = _object(item, where)
        _keys(item, {"name", "description", "cdf_current", "cdf_voltage", "current_sign", "voltage_sign"}, where)
        name = _text(item.get("name"))
        if not name:
            raise ValueError(f"{where}.name is required")
        supplies.append({"name": name, "description": _text(item.get("description")),
                         "cdf_current": item.get("cdf_current"), "cdf_voltage": item.get("cdf_voltage"),
                         "current_sign": item.get("current_sign"), "voltage_sign": item.get("voltage_sign")})
    _name_list(coils, "coils")
    _name_list(supplies, "supplies")
    circuits = []
    ns, nc = len(supplies), len(coils)
    for i, item in enumerate(circuit_items):
        where = f"pf_active manifest.circuits[{i}]"
        item = _object(item, where)
        _keys(item, {"name", "description", "connections_dd3", "connections_dd4", "cdf_current", "cdf_voltage", "current_sign", "voltage_sign"}, where)
        name = _text(item.get("name"))
        if not name:
            raise ValueError(f"{where}.name is required")
        c3 = _parse_connections(item.get("connections_dd3"), 2 * (ns + nc), 3, f"{where}.connections_dd3")
        c4 = _parse_connections(item.get("connections_dd4"), ns + nc, 4, f"{where}.connections_dd4")
        if c3 is None and c4 is None:
            raise ValueError(f"{where}: circuits need sourced DD3 and/or DD4 connectivity")
        circuits.append({
            "name": name, "description": _text(item.get("description")),
            "connections_dd3": c3, "connections_dd4": c4,
            "cdf_current": item.get("cdf_current"), "cdf_voltage": item.get("cdf_voltage"),
            "current_sign": item.get("current_sign"), "voltage_sign": item.get("voltage_sign"),
        })
    _name_list(circuits, "circuits")
    return {"machine": machine_name, "source_reference": source,
            "coils": coils, "supplies": supplies, "circuits": circuits}, raw_text


def build_active_plan(
    nc: Any, *, geometry_path: str | Path | None = None,
    machine: str | None = None,
) -> ActivePlan:
    """Discover CDF coils/supplies or bind channels to authoritative geometry.

    In legacy mode, a KK/VV suffix creates a coil and a VS suffix creates a
    supply.  KT target currents and CC equilibrium-current channels are not
    taken as actual supply currents.  In manifest mode all present KK/VV/VS
    traces must be mapped, or an error is raised instead of silently losing
    physical outputs.  Geometry-only coils are permitted in a device manifest.
    """
    # Target (KT) and equilibrium (CC) identifiers are advisory. An odd
    # unused target array must not break conversion of sound KK/VV/VS data.
    targets = tuple(sorted(name for name in nc.variables if
                           (match := _CHANNEL_PATTERN.fullmatch(name)) and match.group(1) == "KT"))
    equilibrium = tuple(sorted(name for name in nc.variables if
                               (match := _CHANNEL_PATTERN.fullmatch(name)) and match.group(1) == "CC"))
    controller = tuple(sorted(name for name in nc.variables if
                              _CONTROLLER_PATTERN.fullmatch(name)))
    channels = read_active_channels(nc, include_equilibrium=geometry_path is not None)
    if geometry_path is None:
        coil_names = sorted({ch.suffix for ch in channels.values() if ch.prefix in {"KK", "VV"}})
        supply_names = sorted({ch.suffix for ch in channels.values() if ch.prefix == "VS"})
        coils = tuple(CoilSpec(
            name=name, description=f"TRANSP PF coil/circuit {name}; geometry and winding sign not verified",
            current=BoundChannel(channels[f"KK_{name}"], 1) if f"KK_{name}" in channels else None,
            voltage=BoundChannel(channels[f"VV_{name}"], 1) if f"VV_{name}" in channels else None,
            resistance_ohm=None, elements=(),
        ) for name in coil_names)
        supplies = tuple(SupplySpec(
            name=name, description=f"TRANSP PF source {name}; orientation unverified",
            current=None,
            voltage=BoundChannel(channels[f"VS_{name}"], 1),
        ) for name in supply_names)
        return ActivePlan(coils, supplies, (), targets, equilibrium, None, None, None, controller)

    manifest, contents = load_active_manifest(geometry_path, machine=machine)
    used_primary: set[str] = set()
    coils = []
    for i, item in enumerate(manifest["coils"]):
        where = f"coils[{i}]"
        current = _channel_binding(item, "current", {"KK", "CC"}, channels, where,
                                   cc_effective_turns=item["cc_effective_turns"])
        voltage = _channel_binding(item, "voltage", {"VV"}, channels, where)
        for source in (item.get("cdf_current"), item.get("cdf_voltage")):
            if source:
                if source in used_primary:
                    raise ValueError(f"{source}: CDF source assigned to multiple coils")
                used_primary.add(source)
        coils.append(CoilSpec(item["name"], item["description"], current, voltage,
                              item["resistance_ohm"], item["elements"]))
    supplies = []
    used_supply: set[str] = set()
    for i, item in enumerate(manifest["supplies"]):
        where = f"supplies[{i}]"
        current = _channel_binding(item, "current", {"KK"}, channels, where)
        voltage = _channel_binding(item, "voltage", {"VS"}, channels, where)
        for source in (item.get("cdf_current"), item.get("cdf_voltage")):
            if source:
                if source in used_supply:
                    raise ValueError(f"{source}: duplicate supply channel")
                used_supply.add(source)
        supplies.append(SupplySpec(item["name"], item["description"], current, voltage))
    circuits = []
    for i, item in enumerate(manifest["circuits"]):
        where = f"circuits[{i}]"
        circuits.append(CircuitSpec(
            item["name"], item["description"], item["connections_dd3"], item["connections_dd4"],
            _channel_binding(item, "current", {"KK"}, channels, where),
            _channel_binding(item, "voltage", {"VV", "VS"}, channels, where),
        ))
    required_coil = {name for name, ch in channels.items() if ch.prefix in {"KK", "VV"}}
    required_supply = {name for name, ch in channels.items() if ch.prefix == "VS"}
    missing = (required_coil - used_primary) | (required_supply - used_supply)
    if missing:
        raise ValueError("pf_active manifest does not map present CDF channels: "
                         + ", ".join(sorted(missing)))
    unused_equilibrium = tuple(name for name in equilibrium if name not in used_primary
                               and name not in used_supply and
                               all(name not in (entry["cdf_current"], entry["cdf_voltage"])
                                       for entry in manifest["circuits"]))
    return ActivePlan(tuple(coils), tuple(supplies), tuple(circuits), targets,
                      unused_equilibrium, Path(geometry_path).expanduser(),
                      manifest["source_reference"], contents, controller)


def all_time_coordinates(plan: ActivePlan, dd_major: int) -> tuple[np.ndarray, ...]:
    """Return all time arrays that are actually bound to IDS data nodes."""
    traces: list[np.ndarray] = []
    for spec in (*plan.coils, *plan.supplies):
        for binding in (spec.current, spec.voltage):
            if binding is not None:
                traces.append(binding.channel.time)
    for spec in plan.circuits:
        connections = spec.connections_dd4 if int(dd_major) >= 4 else spec.connections_dd3
        if connections is not None:
            for binding in (spec.current, spec.voltage):
                if binding is not None:
                    traces.append(binding.channel.time)
    return tuple(traces)


def _fill_trace(target: Any, binding: BoundChannel | None) -> bool:
    if binding is None:
        return False
    target.time = np.array(binding.channel.time, dtype=np.float64, copy=True)
    target.data = np.asarray(binding.sign * binding.factor * binding.channel.values, dtype=np.float64)
    return True


def _identity(obj: Any, name: str, description: str, dd_major: int) -> None:
    if int(dd_major) < 4 and hasattr(obj, "identifier"):
        obj.identifier = name
        obj.name = description or name
    else:
        obj.name = name
        if description and hasattr(obj, "description"):
            obj.description = description


def populate_active_ids(ids: Any, plan: ActivePlan, dd_major: int) -> dict[str, int]:
    """Populate exactly sourced nodes, avoiding the IMAS AoS resize pitfall.

    DD3 circuit matrices have 2*(Nsupply+Ncoil) side columns with 0/1;
    DD4 has (Nsupply+Ncoil) signed columns with 0/+1/-1.  We never
    convert one convention to the other without machine connectivity data.
    """
    counts = {"coil_current": 0, "coil_voltage": 0,
              "supply_current": 0, "supply_voltage": 0, "circuits": 0}
    ids.coil.resize(len(plan.coils))
    ids.supply.resize(len(plan.supplies))
    applicable_circuits = [c for c in plan.circuits if (
        c.connections_dd4 if int(dd_major) >= 4 else c.connections_dd3
    ) is not None]
    ids.circuit.resize(len(applicable_circuits))

    for idx, spec in enumerate(plan.coils):
        coil = ids.coil[idx]
        _identity(coil, spec.name, spec.description, dd_major)
        if spec.resistance_ohm is not None:
            coil.resistance = spec.resistance_ohm
        coil.element.resize(len(spec.elements))
        for el_idx, element_spec in enumerate(spec.elements):
            populate_element(coil.element[el_idx], element_spec, dd_major)
        counts["coil_current"] += bool(_fill_trace(coil.current, spec.current))
        counts["coil_voltage"] += bool(_fill_trace(coil.voltage, spec.voltage))

    for idx, spec in enumerate(plan.supplies):
        supply = ids.supply[idx]
        _identity(supply, spec.name, spec.description, dd_major)
        counts["supply_current"] += bool(_fill_trace(supply.current, spec.current))
        counts["supply_voltage"] += bool(_fill_trace(supply.voltage, spec.voltage))

    for idx, spec in enumerate(applicable_circuits):
        circuit = ids.circuit[idx]
        _identity(circuit, spec.name, spec.description, dd_major)
        circuit.connections = np.array(spec.connections_dd4 if int(dd_major) >= 4
                                       else spec.connections_dd3, dtype=np.int32, copy=True)
        _fill_trace(circuit.current, spec.current)
        _fill_trace(circuit.voltage, spec.voltage)
        counts["circuits"] += 1
    return counts
