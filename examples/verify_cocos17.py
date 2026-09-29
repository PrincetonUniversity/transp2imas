#!/usr/bin/python
"""Print the expected DD4/COCOS17 signs for a TRANSP CDF/TR.DAT pair.

This is a source-side sanity check. It does not open the generated IMAS database.
"""
import argparse
import re
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file


def read_bool(text, key):
    # strip full-line comments before searching
    active = "\n".join(line for line in text.splitlines()
                       if not line.lstrip().startswith(("!", "#")))
    m = re.search(rf"(?im)^\s*{re.escape(key)}\s*=\s*\.?\s*(true|false|t|f)\s*\.?", active)
    if not m:
        return None
    return m.group(1).lower() in ("true", "t")


def rng(a):
    a = np.asarray(a)
    return float(np.nanmin(a)), float(np.nanmax(a))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cdf")
    ap.add_argument("--dat", help="TR.DAT path (default: <runid>TR.DAT next to CDF)")
    ap.add_argument("--ip-direction", choices=("auto", "ccw", "cw"), default="auto",
                    help="override plasma-current direction when NLJCCW is absent")
    ap.add_argument("--bt-direction", choices=("auto", "ccw", "cw"), default="auto",
                    help="override toroidal-field direction when NLBCCW is absent")
    ns = ap.parse_args()
    cdf = Path(ns.cdf)
    if cdf.suffix.upper() != ".CDF":
        cdf = cdf.with_suffix(".CDF")
    dat = Path(ns.dat) if ns.dat else cdf.with_name(cdf.stem + "TR.DAT")
    text = dat.read_text(errors="replace")
    jccw = read_bool(text, "NLJCCW") if ns.ip_direction == "auto" else ns.ip_direction == "ccw"
    bccw = read_bool(text, "NLBCCW") if ns.bt_direction == "auto" else ns.bt_direction == "ccw"
    missing = []
    if jccw is None:
        missing.append("NLJCCW")
        jccw = True
    if bccw is None:
        missing.append("NLBCCW")
        bccw = True
    if missing:
        print("WARNING: %s absent; using converter fallback CCW direction. "
              "Use --ip-direction/--bt-direction to override." % ", ".join(missing))
    si = 1.0 if jccw else -1.0
    sb = 1.0 if bccw else -1.0
    sq = si * sb
    spsi = -si
    print(f"target COCOS=17; NLJCCW={jccw}; NLBCCW={bccw}")
    print(f"expected multipliers: Ip={si:+.0f} Bphi/Phi/F={sb:+.0f} q={sq:+.0f} psi={spsi:+.0f}")

    nc = netcdf_file(str(cdf), "r", mmap=True)
    try:
        for name, fac in (("PCURC", si), ("BZ", sb), ("Q", sq),
                          ("PLFLX2PI", spsi), ("TRFLX", sb),
                          ("CUR", si), ("CUROH", si), ("CURBS", si)):
            if name not in nc.variables:
                continue
            raw = np.asarray(nc.variables[name].data, dtype=float).copy()
            out = fac * raw
            print(f"{name:10s} raw={rng(raw)!s:>28s}  expected IMAS={rng(out)}")
            del raw, out

        # Core-profile current semantic check (DD4): CUR -> j_phi, while
        # j_total is <J.B>/B0 from PLJB/BZ. Both are converted to A/m^2.
        if all(name in nc.variables for name in ("CUR", "PLJB", "BZ")):
            cur = np.asarray(nc.variables["CUR"].data, dtype=float).copy()
            pljb = np.asarray(nc.variables["PLJB"].data, dtype=float).copy()
            bz = np.asarray(nc.variables["BZ"].data, dtype=float).copy()
            j_phi = si * cur * 1.0e4
            j_total = si * (pljb * 1.0e4) / bz[:, None]
            print(f"j_phi      expected A/m2 range={rng(j_phi)}")
            print(f"j_total    expected A/m2 range={rng(j_total)}  [PLJB/BZ]")
            del cur, pljb, bz, j_phi, j_total
    finally:
        nc.close()


if __name__ == "__main__":
    main()
