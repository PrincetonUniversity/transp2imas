#!/usr/bin/python
"""Source model for NUBEAM fast-ion distributions (TRANSP ``DATA<n>`` + ``FI_TAR.GZ<n>``).

Inputs are limited to the numbered TRANSP ``DATA<n>`` snapshot (distribution,
energy/pitch axes, MC-zone volumes, zone layout) and the permitted equilibrium
members of the matching ``FI_TAR.GZ<n>`` archive (plasma state flux-surface
geometry and the orientation flags NUBEAM itself used).  NUBEAM CPU
particle-state members are never extracted or parsed.

NUBEAM conventions that matter here (cf. NUBEAM ``depsub.f90``,
``fbm_wnetcdf_r8.f90``, ``nbi_fbm_info.f90``):

* ``FBM`` is the guiding-centre distribution, ``FBM_PTCL`` the particle
  (finite-Larmor-radius) distribution, both in ``particles cm^-3 eV^-1`` with
  the pitch normalised so that ``n = sum f dE dxi / 2``.
* Pitch ``xi = v_parallel / v`` is measured along the direction of the toroidal
  plasma current, *not* along B: ``xi_NUBEAM = nsjdotb * xi_B`` where
  ``nsjdotb = +1`` when ``Jphi`` is parallel to ``Bphi`` and ``-1`` otherwise.
* The MC zone grid is ``NZONE_FB`` rows of equal width in ``x`` (square root of
  normalised toroidal flux); row ``i`` has ``NTHZSM[i+lcentr-1] -
  NTHZSM[i+lcentr-2]`` zones, uniform in the geometric poloidal angle between
  ``THBDY0`` and ``THBDY1``.  Zone volumes (``BMVOL``) are in ``cm^3``.
"""

from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
import re
import tarfile

import numpy as np

from nubeam_outputs import discover_nubeam_outputs
from transp_acfile import ACFile

CM3_TO_M3 = 1.0e-6
# Members of FI_TAR archives that may be read.  Anything else (in particular
# ``*ptcl_state*`` CPU particle files) is skipped without being extracted.
PERMITTED_MEMBERS = {
    "ps_ts1_state": re.compile(r"_ps_ts1_state\.cdf$"),
    "xplasma": re.compile(r"_debug_xplasma_nubeam\.cdf$"),
}
FORBIDDEN_MEMBER = re.compile(r"ptcl_state")


class NubeamInputError(RuntimeError):
    """Required NUBEAM output is missing or inconsistent."""


@dataclass
class FluxSurfaceMap:
    """Map NUBEAM ``(x, theta)`` to ``(R, Z)`` in metres via spline interpolation."""

    rho: np.ndarray
    theta: np.ndarray
    r_geo: np.ndarray  # (ntheta, nrho)
    z_geo: np.ndarray
    r_axis: float
    z_axis: float

    def __post_init__(self):
        from scipy.interpolate import RectBivariateSpline
        theta = np.asarray(self.theta, dtype=float)
        if abs((theta[-1] - theta[0]) - 2.0 * np.pi) > 1.0e-6:
            raise NubeamInputError("equilibrium theta grid does not span one full turn")
        # Pad one point each side so cubic splines see the periodicity.
        pad = 3
        ext_theta = np.concatenate([theta[-pad - 1:-1] - 2.0 * np.pi, theta, theta[1:pad + 1] + 2.0 * np.pi])

        def padded(array):
            return np.concatenate([array[-pad - 1:-1], array, array[1:pad + 1]], axis=0)

        self._r = RectBivariateSpline(self.rho, ext_theta, padded(self.r_geo).T, kx=3, ky=3)
        self._z = RectBivariateSpline(self.rho, ext_theta, padded(self.z_geo).T, kx=3, ky=3)

    def rz(self, x, theta):
        x, theta = np.broadcast_arrays(np.asarray(x, float), np.asarray(theta, float))
        theta = (theta + np.pi) % (2.0 * np.pi) - np.pi
        return self._r.ev(x, theta), self._z.ev(x, theta)

    def jacobian(self, x, theta):
        """Return R, ``|dR dZ / dx dtheta|`` at the points."""
        x, theta = np.broadcast_arrays(np.asarray(x, float), np.asarray(theta, float))
        theta = (theta + np.pi) % (2.0 * np.pi) - np.pi
        r = self._r.ev(x, theta)
        det = (self._r.ev(x, theta, dx=1) * self._z.ev(x, theta, dy=1)
               - self._r.ev(x, theta, dy=1) * self._z.ev(x, theta, dx=1))
        return r, np.abs(det)


@dataclass
class McGrid:
    """NUBEAM irregular 2D Monte-Carlo zone grid."""

    nrows: int
    zones_per_row: np.ndarray  # (nrows,)
    theta0: float
    theta1: float

    @property
    def nzones(self):
        return int(self.zones_per_row.sum())

    @property
    def x_edges(self):
        return np.arange(self.nrows + 1) / self.nrows

    def zone_table(self):
        """Per-zone ``(row, index_in_row, x_in, x_out, theta_lo, theta_hi)``."""
        rows = []
        edges = self.x_edges
        span = self.theta1 - self.theta0
        for irow, count in enumerate(self.zones_per_row):
            for j in range(int(count)):
                rows.append((irow, j, edges[irow], edges[irow + 1],
                             self.theta0 + span * j / count, self.theta0 + span * (j + 1) / count))
        return np.array(rows)

    def zone_centres(self):
        table = self.zone_table()
        return 0.5 * (table[:, 2] + table[:, 3]), 0.5 * (table[:, 4] + table[:, 5])

    def zone_volumes(self, fmap, order=6):
        """Toroidal volumes (m^3) of every zone by Gauss-Legendre quadrature."""
        nodes, weights = np.polynomial.legendre.leggauss(order)
        volumes = []
        for _irow, _j, x0, x1, t0, t1 in self.zone_table():
            xs = 0.5 * (x1 - x0) * nodes + 0.5 * (x1 + x0)
            ts = 0.5 * (t1 - t0) * nodes + 0.5 * (t1 + t0)
            xx, tt = np.meshgrid(xs, ts, indexing="ij")
            r, det = fmap.jacobian(xx, tt)
            ww = np.outer(weights, weights) * 0.25 * (x1 - x0) * (t1 - t0)
            volumes.append(2.0 * np.pi * np.sum(ww * r * det))
        return np.array(volumes)

    def polygon_area(self, fmap, samples=12):
        """Poloidal-plane areas (m^2) of zones from curved-edge polygons."""
        areas = []
        for _irow, _j, x0, x1, t0, t1 in self.zone_table():
            s = np.linspace(0.0, 1.0, samples + 1)
            x_path = np.concatenate([x0 + (x1 - x0) * s, np.full(samples, x1),
                                     x1 + (x0 - x1) * s, np.full(samples, x0)])
            t_path = np.concatenate([np.full(samples + 1, t0), t0 + (t1 - t0) * s[1:],
                                     np.full(samples + 1, t1), t1 + (t0 - t1) * s[1:]])
            # Closing segment is implied; drop duplicate end points.
            r, z = fmap.rz(x_path, t_path)
            areas.append(0.5 * abs(np.sum(r * np.roll(z, -1) - np.roll(r, -1) * z)))
        return np.array(areas)

    def connectivity(self):
        """Shared-node connectivity of the zone grid (all indices 0-based).

        Returns a dict with

        * ``nodes``: ``(x_index, Fraction_of_theta_turn)`` keys; the magnetic
          axis is a single node and ``theta0``/``theta1`` are the same point,
        * ``faces``: node-index lists in counter-clockwise order (triangles for
          the central row, quadrilaterals elsewhere),
        * ``edges``: node-index pairs,
        * ``face_edges``: per-face edge indices in the order of the face nodes,
        * ``neighbours``: per-face, per-edge lists of adjacent face indices
          (several where a row boundary joins zones of a finer neighbouring
          row, none on the plasma boundary and for the axis point).
        """
        keys, index = [], {}

        def node(ix, frac):
            key = (0, Fraction(0)) if ix == 0 else (ix, frac % 1)
            if key not in index:
                index[key] = len(keys)
                keys.append(key)
            return index[key]

        faces, spans, row_faces = [], [], []
        for irow, count in enumerate(self.zones_per_row):
            count = int(count)
            row_faces.append([])
            for j in range(count):
                f0, f1 = Fraction(j, count), Fraction(j + 1, count)
                inner0, inner1 = node(irow, f0), node(irow, f1)
                outer0, outer1 = node(irow + 1, f0), node(irow + 1, f1)
                faces.append([inner0, outer0, outer1] if irow == 0 else [inner0, outer0, outer1, inner1])
                spans.append((f0, f1))
                row_faces[-1].append(len(faces) - 1)

        edges, edge_index = [], {}
        face_edges = []
        for nodes in faces:
            ids = []
            for k, a in enumerate(nodes):
                b = nodes[(k + 1) % len(nodes)]
                key = frozenset((a, b))
                if key not in edge_index:
                    edge_index[key] = len(edges)
                    edges.append((a, b))
                ids.append(edge_index[key])
            face_edges.append(ids)

        def overlapping(row, lo, hi):
            return [f for f in row_faces[row] if spans[f][0] < hi and spans[f][1] > lo]

        owners = {}
        for iface, ids in enumerate(face_edges):
            for e in ids:
                owners.setdefault(e, []).append(iface)
        neighbours = []
        rows_of = np.repeat(np.arange(self.nrows), self.zones_per_row)
        for iface, nodes in enumerate(faces):
            row = int(rows_of[iface])
            lo, hi = spans[iface]
            per_edge = []
            for k, e in enumerate(face_edges[iface]):
                a, b = nodes[k], nodes[(k + 1) % len(nodes)]
                xa, xb = keys[a][0], keys[b][0]
                if xa == xb and xa > 0:                       # arc of constant x
                    other = row + 1 if xa == row + 1 else row - 1
                    per_edge.append(overlapping(other, lo, hi) if 0 <= other < self.nrows else [])
                elif xa == xb == 0:                           # degenerate (not produced)
                    per_edge.append([])
                else:                                         # radial edge of constant theta
                    per_edge.append([f for f in owners[e] if f != iface])
            neighbours.append(per_edge)
        return dict(nodes=keys, faces=faces, edges=edges, face_edges=face_edges, neighbours=neighbours)


@dataclass
class FastIonSpecies:
    label: str
    a: float
    z: float
    fusion_product: bool
    energy_centres_ev: np.ndarray
    energy_edges_ev: np.ndarray
    gc: np.ndarray               # (nE, npitch, nzone) particles cm^-3 eV^-1, NUBEAM pitch
    particle: np.ndarray | None  # same, finite-Larmor-radius distribution
    index: int = 1               # 1-based NUBEAM fast-species index (matches SPECIES_LOST)

    @property
    def cdf_suffix(self):
        return cdf_species_suffix(self.a, self.z, self.fusion_product)

    def cdf_names(self, stem, allow_plain=False):
        """Candidate TRANSP CDF names of per-species quantity ``stem``, best first.

        Fusion products only have the quantities of ``FUSION_STEMS`` (under the F family).
        The plain (unsuffixed) name is the sum over beam species, so it is offered only
        when ``allow_plain`` (exactly one beam species) and never for fusion products.
        """
        suffix = self.cdf_suffix
        if self.fusion_product:
            family = FUSION_STEMS.get(stem)
            return [] if family is None or suffix is None else [f"{family}_{suffix}"]
        return ([f"{stem}_{suffix}"] if suffix else []) + ([stem] if allow_plain else [])


# TRANSP names per-species CDF variables with a suffix.  Beam (and minority) ions use
# H, D, T, 3 (He3), 4 (He4) and N/A/K/X (Ne, Ar, Kr, Xe); fusion products use a mirrored
# "F" variable family (FDENS_4, FTRAP_T, ...) with _3, _4, _T and _P suffixes.  Only some
# beam quantities have a fusion-product counterpart (below).
BEAM_SUFFIX = {(1, 1): "H", (1, 2): "D", (1, 3): "T", (2, 3): "3", (2, 4): "4"}
HEAVY_BEAM_SUFFIX = {10: "N", 18: "A", 36: "K", 54: "X"}
FUSION_SUFFIX = {(1, 1): "P", (1, 3): "T", (2, 3): "3", (2, 4): "4"}
FUSION_STEMS = {"BDENS": "FDENS", "BDENS2": "FDENS2", "BTRAP": "FTRAP", "BSTH": "FSTH", "BPHXB": "FPHXB"}


def cdf_species_suffix(a, z, fusion_product=False):
    """TRANSP CDF variable-name suffix of a fast species (``None`` if unknown)."""
    key = (int(round(z)), int(round(a)))
    if fusion_product:
        return FUSION_SUFFIX.get(key)
    return BEAM_SUFFIX.get(key) or HEAVY_BEAM_SUFFIX.get(key[0], "Z")


@dataclass
class LostIons:
    """NUBEAM lost-orbit records gathered over ``[t_start, t_end]`` (SI units).

    ``pitch`` keeps the NUBEAM sign convention (along the plasma current).
    ``prompt`` marks prompt (first-orbit) losses (NUBEAM ``GOOSE_LOST`` = 1).
    """

    t_start: float
    t_end: float
    time: np.ndarray
    nbi_unit: np.ndarray
    energy_fraction: np.ndarray
    weight: np.ndarray        # real particles represented by each record
    r: np.ndarray             # m
    z: np.ndarray             # m
    energy_ev: np.ndarray
    pitch: np.ndarray
    prompt: np.ndarray
    species_index: np.ndarray

    @property
    def count(self):
        return int(self.weight.size)

    def loss_power(self):
        """Mean power (W) carried by the records over the gathering window."""
        span = self.t_end - self.t_start
        return float(np.sum(self.weight * self.energy_ev * 1.602176634e-19) / span) if span > 0 else float("nan")


@dataclass
class FastIonSnapshot:
    time: float
    run_id: str
    data_file: Path
    fi_archive: Path
    species: list
    grid: McGrid
    zone_volume_m3: np.ndarray            # NUBEAM BMVOL converted to m^3
    pitch_centres: np.ndarray             # NUBEAM convention
    pitch_edges: np.ndarray
    flux_map: FluxSurfaceMap
    nsjdotb: int                          # +1 when Jphi parallel to Bphi
    nsjdotb_source: str
    b_axis_vac: float | None = None
    lost: LostIons | None = None
    notes: list = field(default_factory=list)


_BOOL = re.compile(r"(?im)^\s*{name}\s*=\s*\.?([TtFf])")


def read_orientation_flags(trdat_text):
    """Return ``(jphi_ccw, bphi_ccw)`` from TR.DAT, ``None`` for missing flags."""
    out = []
    for name in ("NLJCCW", "NLBCCW"):
        match = re.search(_BOOL.pattern.format(name=name), re.sub(r"!.*", "", trdat_text), re.M | re.I)
        out.append(None if match is None else match.group(1).upper() == "T")
    return tuple(out)


def read_run_nsjdotb(trdat_file):
    """Sign relation ``+1`` (Jphi || Bphi) / ``-1`` from TR.DAT, or ``None``."""
    jccw, bccw = read_orientation_flags(Path(trdat_file).read_text(errors="replace"))
    if jccw is None or bccw is None:
        return None
    return 1 if jccw == bccw else -1


def pitch_edges_from_centres(centres):
    """Pitch bin boundaries as built by the NUBEAM kernel writer (``-1 ... +1``)."""
    centres = np.asarray(centres, float)
    return np.concatenate([[-1.0], 0.5 * (centres[1:] + centres[:-1]), [1.0]])


def _ids_of_species(a, z):
    z_int = int(round(z))
    a_int = int(round(a))
    names = {(1, 1): "H", (1, 2): "D", (1, 3): "T", (2, 3): "He3", (2, 4): "He4"}
    return names.get((z_int, a_int), f"Z{z_int}A{a_int}")


def read_data_snapshot(data_file):
    """Read the distribution-related content of one ``DATA<n>`` file."""
    ac = ACFile(data_file)
    for name in ("TA", "NSFAST", "NZNBMEA", "EFBMB", "XXKSID", "FBM", "BMVOL",
                 "NTHZSM", "LCENTR", "NFBZNSI", "THBDY0", "THBDY1", "ABEAMS", "XZBEAMS"):
        if name not in ac:
            raise NubeamInputError(f"{Path(data_file).name} has no {name}; "
                                   "cannot build the NUBEAM fast-ion distribution")
    nzones = int(ac.get("NFBZNSI"))
    cumulative = np.asarray(ac.get("NTHZSM"), dtype=int)
    lcentr = int(ac.get("LCENTR"))
    counts = None
    for nrows in range(1, cumulative.size - lcentr + 2):
        if cumulative[lcentr + nrows - 2] == nzones:
            counts = np.diff(cumulative[lcentr - 2:lcentr + nrows - 1])
            break
    if counts is None or counts.sum() != nzones or (counts <= 0).any():
        raise NubeamInputError("NTHZSM/LCENTR do not describe NFBZNSI MC zones")
    if "NZONE_FB" in ac and int(ac.get("NZONE_FB")) != counts.size:
        raise NubeamInputError("NZONE_FB disagrees with the NTHZSM zone-row layout")
    grid = McGrid(counts.size, counts, float(ac.get("THBDY0")), float(ac.get("THBDY1")))

    nspecies = int(ac.get("NSFAST"))
    fbm = ac.get("FBM")
    flr = ac.get("FBM_PTCL") if "FBM_PTCL" in ac else None
    energy = np.atleast_2d(ac.get("EFBM"))
    energy_edges = np.atleast_2d(ac.get("EFBMB"))
    nebins = np.atleast_1d(ac.get("NZNBMEA"))
    abeams = np.atleast_1d(ac.get("ABEAMS"))
    zbeams = np.atleast_1d(ac.get("XZBEAMS"))
    fprod = np.atleast_1d(ac.get("NLFPROD")) if "NLFPROD" in ac else np.zeros(nspecies, bool)
    pitch = np.asarray(ac.get("XXKSID"), dtype=float)
    npitch = int(ac.get("NZNBMA")) if "NZNBMA" in ac else pitch.size
    pitch = pitch[:npitch]
    species = []
    for ns in range(nspecies):
        ne = int(nebins[ns])
        if ne <= 0 or not np.any(fbm[:ne, :npitch, :, ns]):
            continue
        particle = None
        if flr is not None:
            particle = flr[:ne, :npitch, :, ns, 0]
        species.append(FastIonSpecies(
            label=_ids_of_species(abeams[ns], zbeams[ns]), a=float(abeams[ns]), z=float(zbeams[ns]),
            fusion_product=bool(fprod[ns]) if ns < fprod.size else False,
            energy_centres_ev=energy[:ne, ns].copy(), energy_edges_ev=energy_edges[:ne + 1, ns].copy(),
            gc=fbm[:ne, :npitch, :, ns].copy(), particle=particle, index=ns + 1))
    if not species:
        raise NubeamInputError(f"{Path(data_file).name} contains no non-empty fast-ion distribution")
    volumes = np.asarray(ac.get("BMVOL"), dtype=float)[:nzones] * CM3_TO_M3
    return dict(time=float(ac.get("TA")), header=ac.header, species=species, grid=grid,
                volumes=volumes, pitch=pitch, lost=_read_lost_ions(ac))


_LOST_FIELDS = ("GTIME_LOST", "YNBSCE_LOST", "YNBIEN_LOST", "WGHT_LOST", "RMJION_LOST",
                "YYION_LOST", "YZE_LOST", "XKSID_LOST", "GOOSE_LOST", "SPECIES_LOST")


def _read_lost_ions(ac):
    """Lost-orbit records of a DATA file (``None`` if absent); unused slots (weight 0) dropped."""
    if not all(name in ac for name in _LOST_FIELDS + ("GTIME1_LOST", "GTIME2_LOST")):
        return None
    arrays = {name: np.atleast_1d(np.asarray(ac.get(name), dtype=float)) for name in _LOST_FIELDS}
    used = arrays["WGHT_LOST"] > 0.0
    if not used.any():
        return None
    pick = {name: values[used] for name, values in arrays.items()}
    return LostIons(
        t_start=float(ac.get("GTIME1_LOST")), t_end=float(ac.get("GTIME2_LOST")),
        time=pick["GTIME_LOST"], nbi_unit=np.rint(pick["YNBSCE_LOST"]).astype(int),
        energy_fraction=np.rint(pick["YNBIEN_LOST"]).astype(int), weight=pick["WGHT_LOST"],
        r=pick["RMJION_LOST"] * 0.01, z=pick["YYION_LOST"] * 0.01,
        energy_ev=pick["YZE_LOST"], pitch=pick["XKSID_LOST"],
        prompt=pick["GOOSE_LOST"] <= 1.0 + 1.0e-3, species_index=np.rint(pick["SPECIES_LOST"]).astype(int))


def read_archive_members(fi_archive):
    """Stream permitted members of an ``FI_TAR`` archive into memory.

    The archive is read sequentially.  Members that are not explicitly
    permitted (including every CPU particle-state file) are skipped without
    being extracted or parsed.
    """
    found = {}
    wanted = set(PERMITTED_MEMBERS)
    with tarfile.open(fi_archive, mode="r|gz") as tar:
        for member in tar:
            if not member.isfile():
                continue
            base = Path(member.name).name
            if FORBIDDEN_MEMBER.search(base):
                continue
            for key, pattern in PERMITTED_MEMBERS.items():
                if key in wanted and pattern.search(base):
                    found[key] = tar.extractfile(member).read()
                    wanted.discard(key)
                    break
            if "ps_ts1_state" in found and not wanted:
                break
    if "ps_ts1_state" not in found:
        raise NubeamInputError(f"{Path(fi_archive).name} has no *_ps_ts1_state.cdf equilibrium member")
    return found


def _read_equilibrium(ps_bytes):
    import netCDF4
    with netCDF4.Dataset("ps_ts1_state.cdf", mode="r", memory=ps_bytes) as nc:
        def var(name):
            return np.asarray(nc[name][:], dtype=float) if name in nc.variables else None
        rho, theta = var("rho_eq"), var("th_eq")
        r_geo, z_geo = var("R_geo"), var("Z_geo")
        if any(v is None for v in (rho, theta, r_geo, z_geo)):
            raise NubeamInputError("plasma state lacks rho_eq/th_eq/R_geo/Z_geo")
        fmap = FluxSurfaceMap(rho, theta, r_geo, z_geo,
                              float(var("R_axis")), float(var("Z_axis")))
        b_vac = var("B_axis_vac")
    return fmap, None if b_vac is None else float(b_vac)


def _read_xplasma_flags(xp_bytes):
    import netCDF4
    with netCDF4.Dataset("xplasma.cdf", mode="r", memory=xp_bytes) as nc:
        if "__Bphi_CCW" not in nc.variables or "__Jphi_CCW" not in nc.variables:
            return None
        bphi, jphi = int(nc["__Bphi_CCW"][...]), int(nc["__Jphi_CCW"][...])
    return 1 if bphi * jphi > 0 else -1


def match_outputs(data_outputs, fi_outputs, tol=1.0e-5):
    """Pair DATA snapshots with FI archives of the same TR.DAT output time.

    The two lists have independent lengths.  Returns ``(pairs, unmatched_data,
    unmatched_fi)``; only matched pairs can produce a distribution.
    """
    pairs, unmatched_data = [], []
    used = set()
    for data in data_outputs:
        match = next((fi for fi in fi_outputs if abs(fi.time - data.time) <= tol and fi.index not in used), None)
        if match is None:
            unmatched_data.append(data)
        else:
            used.add(match.index)
            pairs.append((data, match))
    unmatched_fi = [fi for fi in fi_outputs if fi.index not in used]
    return pairs, unmatched_data, unmatched_fi


def read_fast_ion_snapshots(case_dir, run_id, trdat_file=None, logger=None):
    """Read every matching DATA/FI_TAR pair selected by the run's TR.DAT lists."""
    case_dir = Path(case_dir)
    trdat = Path(trdat_file) if trdat_file else case_dir / f"{run_id}TR.DAT"
    data_outputs, fi_outputs = discover_nubeam_outputs(case_dir, run_id, trdat)
    pairs, lost_data, lost_fi = match_outputs(data_outputs, fi_outputs)
    note = logger.warning if logger else print
    for out in lost_data:
        note(f"distributions: DATA{out.index} (t={out.time:g} s) has no FI_TAR archive at that time; skipped")
    for out in lost_fi:
        note(f"distributions: FI_TAR.GZ{out.index} (t={out.time:g} s) has no DATA snapshot at that time; skipped")
    run_nsjdotb = read_run_nsjdotb(trdat)
    snapshots = []
    for data, fi in pairs:
        src = read_data_snapshot(data.data_file)
        if abs(src["time"] - data.time) > 1.0e-5:
            raise NubeamInputError(f"{data.data_file.name} time {src['time']:g} != TR.DAT {data.time:g}")
        members = read_archive_members(fi.fi_archive)
        fmap, b_vac = _read_equilibrium(members["ps_ts1_state"])
        flag, source = run_nsjdotb, "TR.DAT NLJCCW/NLBCCW"
        if "xplasma" in members:
            xp_flag = _read_xplasma_flags(members["xplasma"])
            if xp_flag is not None:
                if flag is not None and xp_flag != flag:
                    note(f"distributions: TR.DAT and NUBEAM xplasma disagree on Jphi.Bphi sign; "
                         f"using NUBEAM value {xp_flag:+d}")
                flag, source = xp_flag, "NUBEAM xplasma __Jphi_CCW/__Bphi_CCW"
        if flag is None:
            raise NubeamInputError("cannot determine the Jphi/Bphi orientation (NLJCCW/NLBCCW) "
                                   "needed to define the pitch sign")
        snapshots.append(FastIonSnapshot(
            time=src["time"], run_id=run_id, data_file=data.data_file, fi_archive=fi.fi_archive,
            species=src["species"], grid=src["grid"], zone_volume_m3=src["volumes"],
            pitch_centres=src["pitch"], pitch_edges=pitch_edges_from_centres(src["pitch"]),
            flux_map=fmap, nsjdotb=flag, nsjdotb_source=source, b_axis_vac=b_vac,
            lost=src["lost"]))
    return snapshots


def species_moments(snapshot, species, distribution="gc"):
    """Zone particle numbers, energies and parallel energies in SI units.

    Returns a dict with per-zone densities ``n`` (m^-3), ``energy`` (J m^-3),
    ``energy_parallel`` (J m^-3) and totals ``N``, ``E``, ``E_par``.
    ``energy_parallel`` uses ``v_parallel^2 = xi^2 v^2`` of the NUBEAM pitch (the
    sign convention does not matter for a square).
    """
    f = species.gc if distribution == "gc" else species.particle
    if f is None:
        raise NubeamInputError("particle (FLR) distribution is not available")
    de = np.diff(species.energy_edges_ev)
    dxi = np.diff(snapshot.pitch_edges)
    cell = f * 1.0e6 * de[:, None, None] * dxi[None, :, None] / 2.0   # m^-3 per cell
    energy_j = species.energy_centres_ev * 1.602176634e-19
    n = cell.sum(axis=(0, 1))
    e_dens = np.einsum("epz,e->z", cell, energy_j)
    epar = np.einsum("epz,e,p->z", cell, energy_j, snapshot.pitch_centres ** 2)
    vol = snapshot.zone_volume_m3
    return dict(n=n, energy=e_dens, energy_parallel=epar,
                N=float(np.sum(n * vol)), E=float(np.sum(e_dens * vol)), E_par=float(np.sum(epar * vol)))
