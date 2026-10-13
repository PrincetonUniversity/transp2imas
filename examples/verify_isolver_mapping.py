#!/usr/bin/python
"""Verify real TRANSP ISolver PF outputs without the IMAS runtime.

This regression uses ITER 38530P94 and NSTX 142202U01. It checks that the
converter consumes actual circuit currents and voltages (KK/VV/VS), preserves
all native time points, identifies passive PS_n currents, and never confuses
CC_ equilibrium ampere-turns or KT_/KD_/KS_ controller traces with measured
PF coil/supply currents.

For each coil, the independently saved CC_ and KK_ traces satisfy
CC_(t) = N_effective * KK_(t), with an integer-like device-specific factor.
That is a numerical consistency check, not a license to invent a machine
geometry, winding polarity, or an IMAS effective-turns value from a CDF.

Run ``make verify-isolver``.  Results can be inspected in build/pf-isolver/.
Large source files are an optional local regression dataset, not part of CI.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import re
import sys

import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = ROOT / "examples"
sys.path.insert(0, str(ROOT))

from pf_active_mapping import build_active_plan, populate_active_ids  # noqa: E402
from pf_passive_mapping import build_passive_plan, populate_passive_ids  # noqa: E402
from verify_pf_active_mapping import ids_node as active_ids_node  # noqa: E402
from verify_pf_passive_mapping import ids_node as passive_ids_node  # noqa: E402

CASES = (
    ("ITER", "38530P94", EXAMPLES / "data/iter/38530P94.CDF",
     EXAMPLES / "data/iter/38530P94TR.DAT", 11, 119),
    ("NSTX", "142202U01", EXAMPLES / "data/nstx/142202U01/142202U01.CDF",
     EXAMPLES / "data/nstx/142202U01/142202U01TR.DAT", 9, 82),
)


def verify(condition: bool, message: str) -> None:
    """Raise an actionable assertion without printing huge source arrays."""
    if not condition:
        raise AssertionError(message)


def as_array(variable) -> np.ndarray:
    """Load only one selected PF channel, in double precision."""
    return np.asarray(variable[:], dtype=np.float64)


def isolver_coils_from_trdat(trdat: Path) -> set[str]:
    """Recover declared PFC_NAMES values, including indexed continuations.

    This intentionally is not a complete TR.DAT parser; it reads only string
    assignments to PFC_NAMES in the two example inputs.
    """
    names: set[str] = set()
    for raw in trdat.read_text(encoding="ISO-8859-1").splitlines():
        line = raw.split("!", 1)[0].strip()
        if re.match(r"^PFC_NAMES(?:\s*\(\s*\d+\s*\))?\s*=", line, re.I):
            names.update(x.upper() for x in re.findall(r"'([^']+)'", line))
    verify(bool(names), f"{trdat}: missing PFC_NAMES declarations")
    return names


def verify_run(cdf: Path, trdat: Path, machine: str, ncoil: int, npassive: int):
    """Check the CDF source against real and mock-IDS PF data structures."""
    verify(cdf.is_file(), f"missing ISolver fixture: {cdf}")
    verify(trdat.is_file(), f"missing ISolver namelist: {trdat}")
    trdat_text = trdat.read_text(encoding="ISO-8859-1")
    verify(bool(re.search(r"(?im)^\s*LEVGEO\s*=\s*12\b", trdat_text)),
           f"{trdat.name}: expected LEVGEO=12")
    expected_names = isolver_coils_from_trdat(trdat)
    verify(len(expected_names) == ncoil, f"{trdat.name}: PFC_NAMES count mismatch")

    active_rows, passive_rows = [], []
    with netcdf_file(cdf, "r", mmap=False) as nc:
        active = build_active_plan(nc)
        passive = build_passive_plan(nc)
        verify(len(active.coils) == ncoil, f"{cdf.name}: expected {ncoil} coils, got {len(active.coils)}")
        verify(len(active.supplies) == ncoil, f"{cdf.name}: expected {ncoil} supplies, got {len(active.supplies)}")
        verify(len(passive.loops) == npassive,
               f"{cdf.name}: expected {npassive} passive loops, got {len(passive.loops)}")
        verify({x.name for x in active.coils} == expected_names,
               f"{cdf.name}: CDF KK_ names differ from namelist PFC_NAMES")
        verify(len(active.ignored_targets) == ncoil, f"{cdf.name}: KT_ channels must stay targets")
        verify(len(active.ignored_equilibrium) == ncoil, f"{cdf.name}: CC_ channels must stay separate")
        verify(len(active.ignored_controller) == 2 * ncoil,
               f"{cdf.name}: KD_ desired and KS_ uncertainty channels must stay separate")
        verify(all(s.current is None for s in active.supplies),
               f"{cdf.name}: no actual supply current is directly established by VS_ or KT_")

        active_ids = active_ids_node()
        counts = populate_active_ids(active_ids, active, 4)
        verify(counts['coil_current'] == ncoil and counts['coil_voltage'] == ncoil,
               f"{cdf.name}: coil waveform count changed")
        verify(counts['supply_voltage'] == ncoil and counts['supply_current'] == 0,
               f"{cdf.name}: supply waveform count changed")
        for idx, spec in enumerate(active.coils):
            suffix = spec.name
            kk_name, cc_name, vv_name = (f"KK_{suffix}", f"CC_{suffix}", f"VV_{suffix}")
            kk, cc = as_array(nc.variables[kk_name]), as_array(nc.variables[cc_name])
            verify(np.max(np.abs(kk)) > 0.0, f"{kk_name}: zero-current waveform cannot establish CC/KK")
            iref = int(np.argmax(np.abs(kk)))
            effective_turns = int(round(cc[iref] / kk[iref]))
            verify(effective_turns != 0, f"{cc_name}: invalid effective turns")
            error = float(np.max(np.abs(cc - effective_turns * kk)) /
                          max(1.0, float(np.max(np.abs(cc)))))
            verify(error < 2e-6,
                   f"{cdf.name}/{cc_name}: CC is not an integer multiple of KK; error {error:g}")
            verify(spec.current.channel.variable == kk_name and spec.voltage.channel.variable == vv_name,
                   f"{cdf.name}: incorrect current/voltage channel binding")
            current_node = active_ids.coil[idx].current
            voltage_node = active_ids.coil[idx].voltage
            np.testing.assert_array_equal(current_node.data, kk, err_msg=kk_name)
            np.testing.assert_array_equal(current_node.time, spec.current.channel.time, err_msg=kk_name)
            np.testing.assert_array_equal(voltage_node.data, as_array(nc.variables[vv_name]), err_msg=vv_name)
            np.testing.assert_array_equal(voltage_node.time, spec.voltage.channel.time, err_msg=vv_name)
            active_rows.append({
                "coil": suffix, "cdf_current": kk_name, "cdf_ampere_turns": cc_name,
                "cc_per_kk_integer": effective_turns, "max_rel_cc_closure_error": f"{error:.7g}",
                "cdf_terminal_voltage": vv_name, "cdf_supply_voltage": f"VS_{suffix}",
                "n_time": len(current_node.time),
            })
        for idx, spec in enumerate(active.supplies):
            expected_voltage = f"VS_{spec.name}"
            verify(spec.voltage.channel.variable == expected_voltage, f"{cdf.name}: bad supply {spec.name}")
            np.testing.assert_array_equal(active_ids.supply[idx].voltage.data,
                                          as_array(nc.variables[expected_voltage]))
            np.testing.assert_array_equal(active_ids.supply[idx].voltage.time,
                                          spec.voltage.channel.time)

        passive_ids = passive_ids_node()
        written = populate_passive_ids(passive_ids, passive, 4)
        verify(written == npassive, f"{cdf.name}: expected {npassive} passive current traces")
        indices = [spec.channel.index for spec in passive.loops]
        verify(indices == list(range(1, npassive + 1)),
               f"{cdf.name}: missing or out-of-order PS_n channels")
        for idx, spec in enumerate(passive.loops):
            ch = spec.channel
            verify(ch is not None, f"{cdf.name}: unexpected empty current at loop {idx}")
            np.testing.assert_array_equal(passive_ids.loop[idx].current, as_array(nc.variables[ch.variable]))
            np.testing.assert_array_equal(passive_ids.loop[idx].time, ch.time)
            verify(bool(re.fullmatch(r"PS_[1-9][0-9]*", ch.variable)),
                   f"{cdf.name}: invalid passive variable {ch.variable}")
            passive_rows.append({"cdf_variable": ch.variable,
                                 "transp_label": ch.name,
                                 "description": ch.description,
                                 "n_time": len(ch.time),
                                 "units": ch.source_units})
        verify(len({x['transp_label'] for x in passive_rows}) == npassive,
               f"{cdf.name}: duplicated passive current device labels")
    return active_rows, passive_rows


def save_csv(filename: Path, rows: list[dict]) -> None:
    """Write small, auditable channel-index tables (not waveform copies)."""
    with filename.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report-dir', type=Path, default=None,
                        help='Directory for per-device PF source-channel CSV reports')
    args = parser.parse_args()
    if args.report_dir:
        args.report_dir.mkdir(parents=True, exist_ok=True)

    for machine, runid, cdf, trdat, ncoil, npassive in CASES:
        active, passive = verify_run(cdf, trdat, machine, ncoil, npassive)
        max_error = max(float(x['max_rel_cc_closure_error']) for x in active)
        print(f"{machine} {runid}: PASS: {len(active)} coil currents/voltages/supplies, "
              f"{len(passive)} passive loops; CC/KK max relative error={max_error:.3g}")
        if args.report_dir:
            save_csv(args.report_dir / f"{runid}-pf-active.csv", active)
            save_csv(args.report_dir / f"{runid}-pf-passive.csv", passive)
    print('ISolver PF source mapping: PASS (real ITER + NSTX fixtures; no IMAS runtime)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
