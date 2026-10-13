#!/usr/bin/python
"""Inspect a TRANSP CDF and optionally plot its final safety-factor profiles."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cdf", type=Path, help="TRANSP CDF file")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("q-final6.png"),
        help="plot destination (default: %(default)s)",
    )
    parser.add_argument("--show", action="store_true", help="display the plot interactively")
    args = parser.parse_args()

    with Dataset(args.cdf, "r") as dataset:
        names = list(dataset.variables)
        print(f"file: {args.cdf}")
        print(f"data model: {dataset.data_model}")
        print(f"variables: {len(names)}")

        for required in ("XB", "Q"):
            if required not in dataset.variables:
                raise SystemExit(f"required variable not found: {required}")

        xb = np.asarray(dataset.variables["XB"][:])
        q = np.asarray(dataset.variables["Q"][:])
        print(f"XB shape: {xb.shape}")
        print(f"Q shape: {q.shape}")

    figure, axis = plt.subplots(figsize=(8, 6))
    for index in range(max(0, q.shape[0] - 6), q.shape[0]):
        axis.plot(xb[index], q[index], label=f"time index {index}")
    axis.set(xlabel="rho boundary (XB)", ylabel="safety factor (Q)", xlim=(0, 1))
    axis.grid(True)
    axis.legend()
    figure.tight_layout()
    figure.savefig(args.output, dpi=150)
    print(f"wrote: {args.output}")
    if args.show:
        plt.show()
    plt.close(figure)


if __name__ == "__main__":
    main()
