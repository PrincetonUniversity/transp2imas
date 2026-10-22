#!/usr/bin/python
"""Read the standalone FLR fast-ion distribution CDF written by NUBEAM.

This module intentionally depends only on NetCDF/SciPy and NumPy.  It does
not import TRANSP or NUBEAM Python/Fortran libraries.  The NUBEAM
``*_fbm_ptcl.cdf`` format contains the distribution and energy grid; its
radial coordinate and time are supplied by the companion TRANSP CDF.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import netcdf_file


@dataclass(frozen=True)
class FastIonDistribution:
    """A selected NUBEAM FLR distribution and its native axes.

    ``values`` is ordered ``(radial_zone, pitch, energy)`` in
    particles m^-3 keV^-1. ``pitch`` is v_parallel/v on NUBEAM's native
    uniform grid. ``radial_coordinate`` is the native TRANSP MC-grid ``x``
    coordinate and is deliberately not relabelled as rho_tor_norm.
    """

    species_label: str
    time: float
    radial_coordinate: np.ndarray
    pitch: np.ndarray
    energy_kev: np.ndarray
    values: np.ndarray

    def density(self):
        """Integrate over NUBEAM energy and pitch bins (particles m^-3)."""
        energy_edges = _center_edges(self.energy_kev)
        energy_width = np.diff(energy_edges)
        pitch_width = 2.0 / self.pitch.size
        # NUBEAM's delta-solid-angle convention assigns dSA=1/Npitch to
        # each pitch bin (see the NUBEAM FBM writer comments).
        return np.sum(self.values * energy_width[None, None, :], axis=(1, 2)) * pitch_width / 2.0


def _center_edges(centers):
    """Infer finite-volume edges from monotonic center coordinates."""
    centers = np.asarray(centers, dtype=np.float64)
    if centers.ndim != 1 or centers.size < 2 or np.any(np.diff(centers) <= 0.0):
        raise ValueError("Bin centers must be a strictly increasing 1D array")
    edges = np.empty(centers.size + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return edges


def _read_variable(path, name):
    with netcdf_file(str(path), "r", mmap=False) as dataset:
        if name not in dataset.variables:
            raise KeyError(f"{name!r} is absent from {path}")
        return np.array(dataset.variables[name][:], dtype=np.float64, copy=True)


def read_fbm_ptcl(distribution_file, transport_cdf, species_label,
                  requested_time=None):
    """Read one species from ``*_fbm_ptcl.cdf`` at the nearest TRANSP time.

    The CDF's NetCDF dimension order is ``(energy, pitch, radius, gyro)``;
    the returned array is transposed to ``(radius, pitch, energy)``.  The
    writer emits one time-independent FLR file, so ``requested_time`` selects
    the closest time slice of the companion TRANSP ``TIME3``/``MCINDX`` grid.
    """
    distribution_file = Path(distribution_file)
    transport_cdf = Path(transport_cdf)
    label = str(species_label).strip()
    e_name = f"E_{label}"
    f_name = f"FBM_PTCL_{label}"

    with netcdf_file(str(distribution_file), "r", mmap=False) as dataset:
        if e_name not in dataset.variables or f_name not in dataset.variables:
            raise KeyError(f"Species {label!r} is not present in {distribution_file}")
        energy = np.array(dataset.variables[e_name][:], dtype=np.float64, copy=True)
        raw = np.array(dataset.variables[f_name][:], dtype=np.float64, copy=True)

    times = _read_variable(transport_cdf, "TIME3").reshape(-1)
    radial = _read_variable(transport_cdf, "MCINDX")
    if radial.ndim != 2 or radial.shape[0] != times.size:
        raise ValueError("MCINDX must have shape (TIME3, radial_zone)")
    if raw.ndim != 4:
        raise ValueError(f"{f_name} must have four axes (gyro, radius, pitch, energy)")
    # scipy's NetCDF reader exposes dimensions in the file declaration order;
    # NUBEAM writes (energy, pitch, radius, gyro).
    if raw.shape[0] == 1 and raw.shape[-1] == energy.size:
        # The Fortran writer declares (energy, pitch, radius, gyro), which
        # SciPy presents in reverse order: (gyro, radius, pitch, energy).
        values = raw[0]
    elif raw.shape[-1] == 1 and raw.shape[0] == energy.size:
        values = raw[:, :, :, 0].transpose(2, 1, 0)
    else:
        raise ValueError(f"Unexpected shape for {f_name}: {raw.shape}")

    if values.shape[0] != radial.shape[1] or values.shape[2] != energy.size:
        raise ValueError(
            f"Distribution axes {values.shape} do not match radial grid "
            f"{radial.shape[1]} and energy grid {energy.size}"
        )
    if requested_time is None:
        time_index = times.size - 1
    else:
        time_index = int(np.argmin(np.abs(times - float(requested_time))))
    pitch = -1.0 + (np.arange(values.shape[1], dtype=np.float64) + 0.5) * (
        2.0 / values.shape[1]
    )
    return FastIonDistribution(
        species_label=label,
        time=float(times[time_index]),
        radial_coordinate=np.array(radial[time_index], copy=True),
        pitch=pitch,
        energy_kev=energy,
        values=values,
    )
