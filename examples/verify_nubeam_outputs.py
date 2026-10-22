#!/usr/bin/python
"""Check only TR.DAT-selected DATA snapshots and FI_TAR archives."""

import argparse
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from nubeam_outputs import discover_nubeam_outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir", nargs="?", type=Path,
                        default=ROOT / "examples" / "data" / "nstx" / "120967E03_send")
    parser.add_argument("--run-id", default="120967E03")
    args = parser.parse_args()
    data, fi = discover_nubeam_outputs(args.case_dir, args.run_id)
    for output in data:
        header = output.data_file.read_text(errors="replace").splitlines()[0]
        match = re.match(
            rf"\s*{re.escape(args.run_id)}\s+(\d+)\s+([+-]?[\d.]+(?:[EeDd][+-]?\d+)?)",
            header,
        )
        assert match, f"Unexpected DATA header: {header!r}"
        assert int(match.group(1)) == output.index
        saved_time = float(match.group(2).replace("D", "E").replace("d", "e"))
        assert abs(saved_time - output.time) < 1.0e-5, (
            f"{output.data_file.name} time {saved_time} != TR.DAT time {output.time}"
        )
        print(f"DATA{output.index}: {saved_time:.7g} s ({output.data_file.stat().st_size} bytes)")
    for output in fi:
        print(f"FI_TAR.GZ{output.index}: TR.DAT time {output.time:.7g} s "
              f"({output.fi_archive.stat().st_size} bytes)")
    print("PASS: output lists follow TR.DAT outtim/fi_outtim; no particle-state files read")


if __name__ == "__main__":
    main()
