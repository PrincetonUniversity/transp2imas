#!/usr/bin/python
"""Inventory data needed by waves, distributions, sawteeth, radiation detail, and PF passive.

This check deliberately depends only on scipy so it can be run without an IMAS
installation.  It reports whether the CDF/TR.DAT contain enough information to
exercise the next converter IDSs; it does not write an IMAS entry.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from scipy.io import netcdf_file
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sawteeth_mapping import build_sawteeth_plan


def _decode(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace").strip()
    return str(value).strip()


def _read_tr_value(text: str, key: str):
    # TRANSP namelist syntax used by the supplied inputs: one assignment per
    # line, comments introduced by '!'.  Keep this intentionally narrow for an
    # inventory/reporting tool rather than trying to replace the converter's
    # full namelist parser.
    pat = re.compile(rf"^\s*{re.escape(key)}\s*=\s*(.*?)\s*$", re.I)
    for raw in text.splitlines():
        line = raw.split("!", 1)[0]
        m = pat.match(line)
        if m:
            return m.group(1).strip().rstrip(",")
    return None


def _truth(value):
    if value is None:
        return None
    token = value.strip().strip("'\"").lower().replace(".", "")
    if token in {"true", "t"}:
        return True
    if token in {"false", "f"}:
        return False
    return None


def _matches(names, pattern):
    rx = re.compile(pattern, re.I)
    return sorted(name for name in names if rx.fullmatch(name))


def _contains(names, token):
    token = token.upper()
    return sorted(name for name in names if token in name.upper())


def _describe(nc, names):
    lines = []
    for name in names:
        var = nc.variables[name]
        units = _decode(getattr(var, "units", ""))
        long_name = _decode(getattr(var, "long_name", ""))
        dims = ",".join(var.dimensions)
        lines.append(f"    {name}: ({dims}) {units}; {long_name}")
    return lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cdf", type=Path)
    ap.add_argument("trdat", type=Path)
    args = ap.parse_args()

    tr_text = args.trdat.read_text(errors="replace")
    ec_model = _read_tr_value(tr_text, "EC_MODEL")
    nllh = _truth(_read_tr_value(tr_text, "NLLH"))
    nlsaw = _truth(_read_tr_value(tr_text, "NLSAW"))
    nlfbm = _truth(_read_tr_value(tr_text, "NLFBM"))
    aimp = _read_tr_value(tr_text, "AIMP")
    xzimp = _read_tr_value(tr_text, "XZIMP")

    with netcdf_file(str(args.cdf), "r", mmap=False) as nc:
        names = list(nc.variables)

        ec_deposition = _matches(names, r"(?:PEECH|ECCUR)\d*")
        ec_launch = _matches(names, r"PECIN\d*")
        ic_antenna = _matches(names, r"(?:PICHA|FREQA|RFPWRAB)\d+")
        ic_resolved_profiles = _matches(names, r"PIC_F\d+")
        ic_profiles = [n for n in ("QICHA", "ICCUR", "XGRID_NPHI", "PIC_NPHI",
                                    "PEICH", "PIICH", "QRFBTOT") if n in names]
        rf_summaries = [n for n in ("PRFB_D", "PRFFI") if n in names]

        fbm = _matches(names, r"FBM(?:_.*)?")
        bmvol = _matches(names, r"BMVOL(?:_.*)?")
        bdens2 = _matches(names, r"BDENS2.*")
        eba2pl = _matches(names, r"EBA2PL.*")
        eba2pp = _matches(names, r"EBA2PP.*")

        saw_names = sorted(set(_contains(names, "SAW")))
        # XIQ1 is useful context (q=1-surface location) but is not an event
        # time or a crash/mixing diagnostic by itself.
        saw_context = [n for n in ("XIQ1",) if n in names]
        crash_names = sorted(set(_contains(names, "CRASH") + _contains(names, "MIX")))

        ec_supported = bool(ec_deposition)
        ic_supported = bool(ic_resolved_profiles or _matches(names, r"RFPWRAB\d+")) or (
            len(set(int(re.search(r"(\d+)$", n).group(1)) for n in ic_antenna if re.search(r"(\d+)$", n))) == 1
            and "QICHA" in names
        )
        waves_supported = ec_supported or ic_supported

        # Full FBM output needs a phase-space distribution plus its geometric
        # and energy/pitch descriptors.  BDENS2/BMVOL alone are reduced moments.
        distribution_complete = bool(fbm and bmvol and eba2pl and eba2pp)
        saw_plan = build_sawteeth_plan(nc)
        saw_events_available = bool(saw_plan is not None and saw_plan.events > 0)

        rad_generic = [n for n in ("PRLS_TOK", "PRLX_SINGL") if n in names]
        passive_currents = _matches(names, r"PS_\d+")
        edge_state = [n for n in ("T0RECYC", "T0GASFL", "OM0RECYC", "OM0GASFL") if n in names]
        edge_boundary = [n for n in ("PEEDG", "PIEDG") if n in names]
        edge_integrals = [n for n in ("GASD", "RCYD", "GASIZ", "GASIN") if n in names]

        print("TRANSP next-IDS inventory")
        print(f"CDF:    {args.cdf}")
        print(f"TR.DAT: {args.trdat}")
        print()
        print("Namelist state")
        print(f"  EC_MODEL = {ec_model}")
        print(f"  NLLH     = {nllh}")
        print(f"  NLSAW    = {nlsaw}")
        print(f"  NLFBM    = {nlfbm}")
        print(f"  AIMP     = {aimp}")
        print(f"  XZIMP    = {xzimp}")
        print()

        print("waves")
        print(f"  EC deposition/current profiles: {ec_deposition or 'none'}")
        print(f"  EC launcher-only profiles:      {ec_launch or 'none'}")
        print(f"  IC antenna traces:              {ic_antenna or 'none'}")
        print(f"  IC resolved absorption profiles:{ic_resolved_profiles or 'none'}")
        print(f"  IC spatial profiles/grid:       {ic_profiles or 'none'}")
        print(f"  RF summary quantities:          {rf_summaries or 'none'}")
        print(f"  Supported waves IDS payload:    {'YES' if waves_supported else 'NO'}")
        for line in _describe(nc, rf_summaries):
            print(line)
        print()

        print("fast-ion distributions")
        print(f"  FBM:     {fbm or 'none'}")
        print(f"  BMVOL:   {bmvol or 'none'}")
        print(f"  BDENS2:  {bdens2 or 'none'}")
        print(f"  EBA2PL:  {eba2pl or 'none'}")
        print(f"  EBA2PP:  {eba2pp or 'none'}")
        print(f"  Complete FBM distribution payload: {'YES' if distribution_complete else 'NO'}")
        for line in _describe(nc, bmvol + bdens2):
            print(line)
        print()

        print("sawteeth")
        print(f"  Sawtooth-named variables: {saw_names or 'none'}")
        print(f"  q=1 context:              {saw_context or 'none'}")
        print(f"  Crash/mixing variables:   {crash_names or 'none'}")
        print(f"  TLASTSAW new updates:     {saw_plan.events if saw_plan else 0}")
        print(f"  Exportable crash events:  {'YES' if saw_events_available else 'NO'}")
        for line in _describe(nc, saw_names + saw_context + crash_names):
            print(line)
        print()

        print("radiation species detail")
        print(f"  Generic single-impurity aliases: {rad_generic or 'none'}")
        print(f"  Namelist impurity A/Z:           AIMP={aimp}, XZIMP={xzimp}")
        for line in _describe(nc, rad_generic):
            print(line)
        print()

        print("legacy edge-source comments / disposition")
        print(f"  Neutral source-state descriptors: {edge_state or 'none'}")
        print(f"  Boundary energy terms:            {edge_boundary or 'none'}")
        print(f"  Integrated neutral/source rates:  {edge_integrals or 'none'}")
        for line in _describe(nc, edge_state + edge_boundary + edge_integrals):
            print(line)
        print("  Note: state descriptors and boundary transport terms are not volumetric edge_sources terms.")
        print()

        print("pf_passive")
        print(f"  Passive current channels PS_n: {passive_currents or 'none'}")
        for line in _describe(nc, passive_currents):
            print(line)
        print()

        if not waves_supported:
            print("RESULT waves: skip for this run (no supported absorbed/deposited wave payload).")
        else:
            print("RESULT waves: payload present; exercise the waves writer in an IMAS runtime.")
        if not distribution_complete:
            print("RESULT distributions: incomplete; do not synthesize a phase-space distribution from moments.")
        else:
            print("RESULT distributions: full FBM-style payload appears present.")
        if not saw_events_available:
            print("RESULT sawteeth: no crash-event payload to export for this run.")
        else:
            print("RESULT sawteeth: native TLASTSAW crash history is exportable; trigger mechanism and mixing/inversion radii remain unavailable.")
        if rad_generic:
            print("RESULT radiation detail: generic TOK/SINGL aliases present; resolve species from AIMP/XZIMP before mapping.")
        if not passive_currents:
            print("RESULT pf_passive: no PS_n current payload in this run; geometry/circuit audit cannot be validated here.")
        else:
            print("RESULT pf_passive: passive-current payload present; geometry/circuit audit can proceed.")


if __name__ == "__main__":
    main()
