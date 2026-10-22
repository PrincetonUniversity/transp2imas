#!/usr/bin/python
"""Locate NUBEAM outputs selected by the TRANSP TR.DAT time lists.

The NUBEAM conversion path uses only TRANSP DATA snapshots and FI_TAR
archives.  CPU particle-state CDFs are intentionally not discovered or read.
"""

from dataclasses import dataclass
from pathlib import Path
import re


@dataclass(frozen=True)
class NubeamOutput:
    index: int
    time: float
    data_file: Path | None = None
    fi_archive: Path | None = None


def _time_list(trdat_text, variable):
    """Read one comma-separated Fortran real list from TR.DAT."""
    text = re.sub(r"!.*", "", trdat_text)
    match = re.search(
        rf"(?im)^\s*{re.escape(variable)}\s*(?:\(\s*\d+\s*\))?\s*=\s*([^\n/]*)",
        text,
    )
    if not match:
        return []
    values = []
    for token in match.group(1).split(","):
        token = token.strip().rstrip("&").strip()
        if not token:
            continue
        repeat = re.fullmatch(r"(\d+)\s*\*\s*([+-]?[\d.]+(?:[EeDd][+-]?\d+)?)", token)
        if repeat:
            values.extend([float(repeat.group(2).replace("D", "E").replace("d", "e"))]
                          * int(repeat.group(1)))
            continue
        try:
            values.append(float(token.replace("D", "E").replace("d", "e")))
        except ValueError as exc:
            raise ValueError(f"Cannot parse {variable} time value {token!r}") from exc
    # TRANSP uses large positive sentinels for unused output slots.
    return [value for value in values if value < 1000.0]


def discover_nubeam_outputs(case_dir, run_id, trdat_file=None):
    """Return DATA and FI_TAR files for each time requested by TR.DAT.

    ``outtim`` controls numbered ``DATA<n>`` files; ``fi_outtim`` controls
    numbered ``FI_TAR.GZ<n>`` archives. The two lists are independent and
    can have different lengths. Any missing requested output is an error.
    """
    case_dir = Path(case_dir)
    trdat_file = Path(trdat_file) if trdat_file else case_dir / f"{run_id}TR.DAT"
    text = trdat_file.read_text(errors="replace")
    data_times = _time_list(text, "outtim")
    fi_times = _time_list(text, "fi_outtim")
    if not data_times:
        raise ValueError(f"TR.DAT contains no active outtim values: {trdat_file}")
    if not fi_times:
        raise ValueError(f"TR.DAT contains no active fi_outtim values: {trdat_file}")

    data_outputs = []
    for index, time in enumerate(data_times, start=1):
        path = case_dir / f"{run_id}.DATA{index}"
        if not path.is_file():
            raise FileNotFoundError(f"TR.DAT requests {path.name}, but it is missing")
        data_outputs.append(NubeamOutput(index=index, time=time, data_file=path))

    fi_outputs = []
    for index, time in enumerate(fi_times, start=1):
        path = case_dir / f"{run_id}_FI_TAR.GZ{index}"
        if not path.is_file():
            raise FileNotFoundError(f"TR.DAT requests {path.name}, but it is missing")
        fi_outputs.append(NubeamOutput(index=index, time=time, fi_archive=path))
    return data_outputs, fi_outputs
