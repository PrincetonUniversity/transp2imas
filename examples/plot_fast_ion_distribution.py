#!/usr/bin/python
"""Plot NUBEAM fast-ion distributions from a converted IDS or from TRANSP outputs.

Sources
  --case-dir DIR --run-id ID   read DATA<n> / FI_TAR.GZ<n> directly (TR.DAT times)
  --database DB --shot N --run N [--user U] [--backend HDF5] [--dd-version V]
                               read the ``distributions`` IDS written by convert.py
  --uri URI                    read an IMAS URI instead

Figures per selected time and distribution (all in IMAS conventions: pitch is
v_parallel/v along B, energy in keV, densities in m^-3):
  zones      poloidal map of the MC zones coloured by fast-ion density
  phase      f(E, pitch) summed over the plasma and in the core/mid/edge zone rows
  ep         energy-pitch map of particles per cell summed over the plasma
  moments    energy spectrum, pitch distribution, radial density/pressure profile

Example
  python examples/plot_fast_ion_distribution.py --case-dir examples/data/nstx/120967E03_send \\
      --run-id 120967E03 --save-dir build/fi-plots
"""

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

EV = 1.602176634e-19
# viridis for magnitude maps, neutral ink for text, a fixed order for series.
SEQUENTIAL = "viridis"
SERIES = ("#2a6fb0", "#c8751e", "#3a8f5f")
INK, MUTED, GRID = "#1f2933", "#5f6b76", "#d5dae0"


@dataclass
class PhaseSpace:
    """One distribution at one time, in IMAS conventions."""

    label: str
    gyro_type: int
    time: float
    polygons: list           # per zone: (n, 2) R, Z corner array
    zone_rows: np.ndarray    # row index of each zone
    zone_volume: np.ndarray  # m^3
    energy_edges_ev: np.ndarray
    pitch_edges: np.ndarray  # ascending, along B
    fd3v: np.ndarray         # (nzone, nE, npitch)  m^-3
    r_axis: float
    z_axis: float
    mass_kg: float

    @property
    def title(self):
        kind = "particle position (FLR)" if self.gyro_type == 1 else "guiding centre"
        return f"{self.label} fast ions, {kind}, t = {self.time:g} s"

    def n_zone(self):
        return self.fd3v.sum(axis=(1, 2))

    def total(self):
        return float(np.sum(self.n_zone() * self.zone_volume))

    def row_average(self, zone_values):
        rows = self.zone_rows
        vol = np.bincount(rows, weights=self.zone_volume)
        return np.bincount(rows, weights=zone_values * self.zone_volume) / vol

    def kinetic_energy_density(self):
        centre = 0.5 * (self.energy_edges_ev[1:] + self.energy_edges_ev[:-1]) * EV
        return np.einsum("zep,e->z", self.fd3v, centre)

    def dn_de(self):
        """dN/dE in particles per keV, volume integrated."""
        weighted = np.einsum("zep,z->e", self.fd3v, self.zone_volume)
        return weighted / (np.diff(self.energy_edges_ev) / 1e3)

    def dn_dpitch(self):
        weighted = np.einsum("zep,z->p", self.fd3v, self.zone_volume)
        return weighted / np.diff(self.pitch_edges)


def smooth_polygons(coords, faces, rows, axis, samples=14):
    """Zone outlines with curved constant-x edges.

    ``coords`` are the grid nodes, ``faces`` their node lists (counter-clockwise;
    ``[inner0, outer0, outer1, inner1]`` or the central ``[axis, outer0, outer1]``)
    and ``rows`` the row of each face.  Nodes on one flux surface are joined by a
    periodic spline of the polar radius about the axis; constant-theta edges stay
    straight.  Falls back to chords when a surface has too few nodes.
    """
    from scipy.interpolate import CubicSpline
    coords = np.asarray(coords, float)
    level = {}                                     # node -> flux-surface index
    for face, row in zip(faces, rows):
        outer = face[1:3]
        for node in outer:
            level[node] = int(row) + 1
        if len(face) == 4:
            for node in (face[0], face[3]):
                level[node] = int(row)
    surfaces = {}
    for ix in sorted(set(level.values())):
        nodes = [n for n, v in level.items() if v == ix]
        rel = coords[nodes] - np.asarray(axis)
        phi = np.arctan2(rel[:, 1], rel[:, 0])
        order = np.argsort(phi)
        phi, rad = phi[order], np.hypot(rel[order, 0], rel[order, 1])
        if phi.size >= 4 and np.all(np.diff(phi) > 1e-9):
            spline = CubicSpline(np.r_[phi, phi[0] + 2 * np.pi], np.r_[rad, rad[0]], bc_type="periodic")
            surfaces[ix] = spline
    def arc(a, b, ix):
        ra, rb = coords[a] - axis, coords[b] - axis
        pa, pb = np.arctan2(ra[1], ra[0]), np.arctan2(rb[1], rb[0])
        if ix not in surfaces:
            return np.array([coords[a], coords[b]])
        sweep = (pb - pa + np.pi) % (2 * np.pi) - np.pi     # signed, shorter way round
        ang = pa + np.linspace(0.0, sweep, samples)
        start = surfaces[ix].x[0]
        wrapped = start + (ang - start) % (2 * np.pi)
        rad = surfaces[ix](wrapped)
        return np.column_stack([axis[0] + rad * np.cos(ang), axis[1] + rad * np.sin(ang)])
    polygons = []
    for face, row in zip(faces, rows):
        row = int(row)
        if len(face) == 3:
            path = [coords[face[0]][None], arc(face[1], face[2], row + 1)]
        else:
            path = [coords[face[0]][None], arc(face[1], face[2], row + 1),
                    arc(face[3], face[0], row)]
        polygons.append(np.vstack(path))
    return polygons


def phase_spaces_from_ids(ids):
    """Yield :class:`PhaseSpace` objects for every distribution and time of an IDS."""
    from imas_distributions import AMU
    times = np.asarray(ids.time)
    for dist in ids.distribution:
        element = dist.species.ion.element[0]
        label = str(getattr(dist.species.ion, "name", "") or getattr(dist.species.ion, "label", "") or "ion")
        process = str(dist.process[0].type.name) if len(dist.process) else ""
        if process and process != "NBI":          # e.g. fusion T vs beam T share the species name
            label = f"{label} ({process})"
        mass = float(element.a) * AMU
        for it, tval in enumerate(times):
            ggd = dist.ggd[it]
            space = ggd.grid.space
            nodes = np.array([o.geometry for o in space[0].objects_per_dimension[0].object])
            faces = space[0].objects_per_dimension[2].object
            rows = np.array([int(f.geometry[0]) - 1 for f in faces])
            face_nodes = [list(np.asarray(f.nodes, dtype=int) - 1) for f in faces]
            volume = np.array([float(f.geometry[2]) for f in faces])
            e_edges = np.array([o.geometry[0] for o in space[1].objects_per_dimension[0].object])
            p_edges = np.array([o.geometry[0] for o in space[2].objects_per_dimension[0].object])
            values = np.asarray(ggd.expansion_fd3v[0].grid_subset[0].values)
            fd3v = values.reshape((len(faces), e_edges.size - 1, p_edges.size - 1), order="F")
            axis = nodes[face_nodes[0][0]]   # the central row is made of triangles meeting at the axis
            if len(ids.magnetic_axis.r) > it and len(ids.magnetic_axis.z) > it:
                axis = (float(ids.magnetic_axis.r[it]), float(ids.magnetic_axis.z[it]))
            polygons = smooth_polygons(nodes, face_nodes, rows, axis)
            yield PhaseSpace(label, int(dist.gyro_type), float(tval), polygons, rows, volume,
                             e_edges, p_edges, fd3v, float(axis[0]), float(axis[1]), mass)


def phase_spaces_from_source(case_dir, run_id):
    """Same content straight from DATA/FI_TAR (applies the pitch/B convention)."""
    import imas_distributions as idm
    import nubeam_fi as nf
    out = []
    for snap in nf.read_fast_ion_snapshots(case_dir, run_id):
        conn = snap.grid.connectivity()
        theta0, span = snap.grid.theta0, snap.grid.theta1 - snap.grid.theta0
        coords = []
        for ix, frac in conn["nodes"]:
            if ix == 0:
                coords.append((snap.flux_map.r_axis, snap.flux_map.z_axis))
            else:
                r, z = snap.flux_map.rz(ix / snap.grid.nrows, theta0 + float(frac) * span)
                coords.append((float(r), float(z)))
        coords = np.array(coords)
        rows = np.repeat(np.arange(snap.grid.nrows), snap.grid.zones_per_row)
        polygons = smooth_polygons(coords, conn["faces"], rows,
                                   (snap.flux_map.r_axis, snap.flux_map.z_axis))
        edges, perm = idm.pitch_to_field_direction(snap.pitch_edges, snap.nsjdotb)
        for sp in snap.species:
            for gyro, f in ((2, sp.gc), (1, sp.particle)):
                if f is None:
                    continue
                fd3v, _f = idm.phase_space_cells(f, sp, snap, perm)
                shape = (snap.grid.nzones, sp.energy_edges_ev.size - 1, edges.size - 1)
                label = f"{sp.label} (fusion)" if sp.fusion_product else sp.label
                out.append(PhaseSpace(label, gyro, snap.time, polygons, rows, snap.zone_volume_m3,
                                      sp.energy_edges_ev, edges, fd3v.reshape(shape, order="F"),
                                      snap.flux_map.r_axis, snap.flux_map.z_axis, sp.a * idm.AMU))
    return out


def _style(ax):
    ax.grid(True, color=GRID, lw=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)


def _density_f(ps, mask):
    """f in (E, pitch) averaged over the zones selected by ``mask`` (m^-6 s^3 equivalent shape)."""
    vol = ps.zone_volume[mask]
    cell = np.einsum("zep,z->ep", ps.fd3v[mask], vol) / vol.sum()
    de = np.diff(ps.energy_edges_ev)[:, None] / 1e3
    dp = np.diff(ps.pitch_edges)[None, :]
    return cell / (de * dp)       # m^-3 keV^-1 per unit pitch


def figure_zones(ps, plt):
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    fig, ax = plt.subplots(figsize=(5.2, 6.4), constrained_layout=True)
    n = ps.n_zone()
    coll = PolyCollection(ps.polygons, array=n, cmap=SEQUENTIAL, norm=Normalize(0, n.max()),
                          edgecolors="face", linewidths=0.4)
    ax.add_collection(coll)
    ax.plot([ps.r_axis], [ps.z_axis], marker="+", color=INK, ms=8, lw=0)
    allp = np.vstack(ps.polygons)
    ax.set_xlim(allp[:, 0].min() - 0.05, allp[:, 0].max() + 0.05)
    ax.set_ylim(allp[:, 1].min() - 0.05, allp[:, 1].max() + 0.05)
    ax.set_aspect("equal")
    ax.set_xlabel("R (m)", color=INK)
    ax.set_ylabel("Z (m)", color=INK)
    cbar = fig.colorbar(coll, ax=ax, shrink=0.8)
    cbar.set_label("fast-ion density (m$^{-3}$)", color=INK)
    ax.set_title(f"{ps.title}\nN = {ps.total():.3e}, {len(ps.polygons)} MC zones", color=INK, fontsize=9)
    _style(ax)
    ax.grid(False)
    return fig


def figure_phase(ps, plt):
    from matplotlib.colors import LogNorm
    nrows = int(ps.zone_rows.max()) + 1
    selections = [("all zones", np.ones(ps.zone_rows.size, bool)),
                  (f"core row {1}", ps.zone_rows == 0),
                  (f"mid row {nrows // 2 + 1}", ps.zone_rows == nrows // 2),
                  (f"edge row {nrows}", ps.zone_rows == nrows - 1)]
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.8), constrained_layout=True, sharey=True)
    energy_kev = ps.energy_edges_ev / 1e3
    maps = [_density_f(ps, mask) for _name, mask in selections]
    vmax = max(m.max() for m in maps)
    for ax, (name, _mask), f in zip(axes, selections, maps):
        mesh = ax.pcolormesh(ps.pitch_edges, energy_kev, np.where(f > 0, f, np.nan), cmap=SEQUENTIAL,
                             norm=LogNorm(vmax * 1e-4, vmax), shading="flat", rasterized=True)
        ax.set_title(name, color=INK, fontsize=9)
        ax.set_xlabel("pitch  v$_\\parallel$/v  (along B)", color=INK)
        ax.axvline(0.0, color=MUTED, lw=0.6, ls=":")
        _style(ax)
        ax.grid(False)
    axes[0].set_ylabel("energy (keV)", color=INK)
    cbar = fig.colorbar(mesh, ax=axes, shrink=0.9)
    cbar.set_label("f  (m$^{-3}$ keV$^{-1}$ per unit pitch)", color=INK)
    fig.suptitle(ps.title, color=INK, fontsize=10)
    return fig


def figure_moments(ps, plt):
    fig, (ax_e, ax_p, ax_n, ax_w) = plt.subplots(1, 4, figsize=(15, 3.6), constrained_layout=True)
    centre = 0.5 * (ps.energy_edges_ev[1:] + ps.energy_edges_ev[:-1]) / 1e3
    ax_e.plot(centre, ps.dn_de(), color=SERIES[0], lw=1.6)
    ax_e.set_yscale("log")
    ax_e.set_ylim(ps.dn_de().max() * 1e-5, ps.dn_de().max() * 2)
    ax_e.set_xlabel("energy (keV)", color=INK)
    ax_e.set_ylabel("dN/dE (keV$^{-1}$)", color=INK)
    pc = 0.5 * (ps.pitch_edges[1:] + ps.pitch_edges[:-1])
    ax_p.plot(pc, ps.dn_dpitch(), color=SERIES[0], lw=1.6)
    ax_p.axvline(0.0, color=MUTED, lw=0.6, ls=":")
    ax_p.set_xlabel("pitch  v$_\\parallel$/v  (along B)", color=INK)
    ax_p.set_ylabel("dN/d(pitch)", color=INK)
    nrows = int(ps.zone_rows.max()) + 1
    rho = (np.arange(nrows) + 0.5) / nrows
    ax_n.plot(rho, ps.row_average(ps.n_zone()), color=SERIES[0], lw=1.6)
    ax_n.set_ylabel("fast-ion density (m$^{-3}$)", color=INK)
    pressure = 2.0 / 3.0 * ps.row_average(ps.kinetic_energy_density()) / 1e3
    ax_w.plot(rho, pressure, color=SERIES[1], lw=1.6)
    ax_w.set_ylabel("fast-ion pressure (kPa)", color=INK)
    for ax in (ax_n, ax_w):
        ax.set_xlabel(r"$\rho_{tor,norm}$ (zone-row centre)", color=INK)
        ax.set_ylim(bottom=0)
    for ax in (ax_e, ax_p, ax_n, ax_w):
        _style(ax)
    fig.suptitle(ps.title, color=INK, fontsize=10)
    return fig


def figure_ep(ps, plt):
    """Particles per (energy, pitch) cell summed over the plasma, energy on x, pitch on y."""
    particles = np.einsum("zep,z->ep", ps.fd3v, ps.zone_volume).T       # (pitch, energy)
    energy_kev = 0.5 * (ps.energy_edges_ev[1:] + ps.energy_edges_ev[:-1]) / 1e3
    pitch = 0.5 * (ps.pitch_edges[1:] + ps.pitch_edges[:-1])
    peak = particles.max()
    exponent = int(np.floor(np.log10(peak)))
    scale = 10.0 ** exponent
    data = np.ma.masked_where(particles <= 0.0, particles / scale)
    from matplotlib.ticker import MaxNLocator
    levels = MaxNLocator(nbins=9, steps=[1, 2, 2.5, 5, 10]).tick_values(0.0, peak / scale)
    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    cs = ax.contourf(energy_kev, pitch, data, levels=levels, cmap=SEQUENTIAL, extend="max")
    ax.axhline(0.0, color="#2a8f3c", lw=1.4)
    ax.set_xlim(energy_kev[0] - 0.5 * np.diff(ps.energy_edges_ev[:2]) / 1e3, ps.energy_edges_ev[-1] / 1e3)
    ax.set_ylim(-1.0, 1.0)
    ax.set_xlabel("Energy [keV]", color=INK)
    ax.set_ylabel("Pitch  v$_\\parallel$/v  (along B)", color=INK)
    cbar = fig.colorbar(cs, ax=ax)
    cbar.set_label(f"particles per cell  [$\\times$10$^{{{exponent}}}$]", color=INK)
    ax.set_title(ps.title, color=INK, fontsize=9)
    _style(ax)
    ax.grid(False)
    return fig


FIGURES = {"zones": figure_zones, "phase": figure_phase, "moments": figure_moments, "ep": figure_ep}


def make_figures(ps, which=None):
    import matplotlib.pyplot as plt
    return {name: FIGURES[name](ps, plt) for name in (which or FIGURES)}


def public_imasdb_uri(args):
    """URI of an HDF5 entry in ``~/public/imasdb`` (what convert.py writes by default).

    Used when ``--user`` is ``public`` or omitted, because the legacy
    ``user='public'`` form needs ``IMAS_HOME`` while the URI form does not.
    """
    if args.backend.upper() != "HDF5" or args.user not in (None, "public") or args.shot is None:
        return None
    major = (args.dd_version or "4").split(".")[0]
    path = Path.home() / "public" / "imasdb" / args.database / major / str(args.shot) / str(args.run)
    return f"imas:hdf5?path={path}" if path.is_dir() else None


def load(args):
    if args.case_dir:
        return phase_spaces_from_source(args.case_dir, args.run_id)
    import imas
    uri = args.uri or public_imasdb_uri(args)
    if uri:
        print("Opening", uri)
        dd_version = args.dd_version
        if not dd_version and "/imasdb/" in uri:     # directory name after the database is the DD major
            dd_version = "3.39.0" if f"/{args.database}/3/" in uri else "4.1.1"
        entry = imas.DBEntry(uri, "r", dd_version=dd_version) if dd_version else imas.DBEntry(uri, "r")
    else:
        backend = imas.ids_defs.HDF5_BACKEND if args.backend.upper() == "HDF5" else imas.ids_defs.MDSPLUS_BACKEND
        opts = {"dd_version": args.dd_version} if args.dd_version else {}
        entry = imas.DBEntry(backend, args.database, args.shot, args.run, args.user, **opts)
        entry.open()
    try:
        ids = entry.get("distributions", occurrence=args.occurrence)
        return list(phase_spaces_from_ids(ids))
    finally:
        entry.close()


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case-dir", type=Path)
    p.add_argument("--run-id")
    p.add_argument("--uri")
    p.add_argument("--database", default="NSTX")
    p.add_argument("--shot", type=int)
    p.add_argument("--run", type=int)
    p.add_argument("--user")
    p.add_argument("--backend", default="HDF5")
    p.add_argument("--dd-version")
    p.add_argument("--occurrence", type=int, default=0)
    p.add_argument("--time", type=float, help="plot only the time slice closest to this value (s)")
    p.add_argument("--gyro-type", type=int, choices=(1, 2), help="1 particle position, 2 guiding centre")
    p.add_argument("--figures", nargs="+", choices=sorted(FIGURES), default=sorted(FIGURES))
    p.add_argument("--save-dir", type=Path, help="write PNG files here instead of opening windows")
    ns = p.parse_args(argv)
    if bool(ns.case_dir) == bool(ns.uri or ns.shot is not None):
        p.error("choose either --case-dir/--run-id or an IDS source (--uri or --shot/--run)")
    if ns.case_dir and not ns.run_id:
        p.error("--run-id is required with --case-dir")
    return ns


def main(argv=None):
    ns = parse_args(argv)
    if ns.save_dir:
        import matplotlib
        matplotlib.use("Agg")
    spaces = load(ns)
    if ns.gyro_type:
        spaces = [ps for ps in spaces if ps.gyro_type == ns.gyro_type]
    if ns.time is not None and spaces:
        nearest = min({ps.time for ps in spaces}, key=lambda t: abs(t - ns.time))
        spaces = [ps for ps in spaces if ps.time == nearest]
    if not spaces:
        print("No distributions selected")
        return 1
    import matplotlib.pyplot as plt
    for ps in spaces:
        print(f"{ps.title}: N = {ps.total():.4e}, zones = {len(ps.polygons)}, "
              f"E bins = {ps.energy_edges_ev.size - 1}, pitch bins = {ps.pitch_edges.size - 1}")
        for name, fig in make_figures(ps, ns.figures).items():
            if ns.save_dir:
                ns.save_dir.mkdir(parents=True, exist_ok=True)
                safe = re.sub(r"[^A-Za-z0-9]+", "-", ps.label).strip("-")
                path = ns.save_dir / f"fi_{safe}_g{ps.gyro_type}_t{ps.time:.3f}_{name}.png"
                fig.savefig(path, dpi=150)
                print("  saved", path)
                plt.close(fig)
    if not ns.save_dir:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
