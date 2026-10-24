#!/usr/bin/python
"""Populate the IMAS ``distributions`` IDS from NUBEAM fast-ion distributions.

The distribution is written on NUBEAM's own phase-space grid, built from three
GGD spaces:

1. space 1, ``primary_standard``: the irregular poloidal Monte-Carlo zone grid
   (``R``, ``Z`` nodes; triangular central zones and quadrilateral zones),
2. space 2, ``secondary_structured``: kinetic-energy bins (``eV``),
3. space 3, ``secondary_structured``: pitch bins (``v_parallel/v`` along B).

``expansion_fd3v`` holds ``f d3v`` (``m^-3``: particles per m^3 in each
phase-space cell) on the implicit ``4d_hypervolumes`` grid subset, ordered with
the zone index fastest, then energy, then pitch (Fortran order).  ``expansion``
holds ``f`` itself (``m^-6 s^3``) using the exact cell velocity volume
``d3v = 2 pi (v_hi^3 - v_lo^3)/3 * dxi``.

Pitch and COCOS.  NUBEAM measures the pitch along the toroidal plasma current,
``xi_NUBEAM = nsjdotb * xi_B``.  IMAS defines ``pitch`` along the magnetic
field, so ``xi_IMAS = nsjdotb * xi_NUBEAM`` (``nsjdotb = +1`` for Jphi parallel
to Bphi).  This is a physical-direction relation independent of the COCOS
index (11 for DD3, 17 for DD4), which only changes psi/theta/q sign
conventions; the grid here is expressed in ``R, Z`` so no further COCOS
transformation applies.  When ``nsjdotb = -1`` the pitch bins are reversed so
that the axis stays ascending.
"""

from pathlib import Path
import re

import numpy as np
from scipy import constants

from nubeam_fi import (NubeamInputError, McGrid, read_fast_ion_snapshots,
                       species_moments)
from nubeam_outputs import discover_nubeam_outputs
from ids_identity import IdentityHelper
from writer_base import IdsWriter

# GGD / coordinate identifier indices (DD utilities identifiers).
COORD_R, COORD_Z, COORD_ENERGY, COORD_PITCH = 4, 3, 301, 403
SUBSET_CELLS, SUBSET_4D = 5, 200
COORD_TIME_LOST, COORD_NBI_UNIT, COORD_PROMPT = -1, -2, -3   # private coordinate identifiers for markers
SUBSET_ENERGY_BINS, SUBSET_PITCH_BINS = -1, -2
GRID_INDEX = -1
NBI_PROCESS_INDEX = 1
EV = constants.e
AMU = constants.atomic_mass


def regrid_to_rows(values, zone_edges, weights, nrows):
    """Weighted average of zone-centred ``values`` over ``nrows`` equal-width rows in x.

    ``zone_edges`` (nzone + 1) are the source zone boundaries in x, ``weights`` the
    per-zone volume (for volume densities) or cross-section area (for current
    densities).  Zones are split in proportion to their overlap with a row, so the
    weighted integral is conserved.
    """
    values, weights = np.asarray(values, float), np.asarray(weights, float)
    row_edges = np.arange(nrows + 1) / nrows
    lo, hi = np.asarray(zone_edges[:-1], float), np.asarray(zone_edges[1:], float)
    out = np.empty(nrows)
    for k in range(nrows):
        frac = np.clip(np.minimum(row_edges[k + 1], hi) - np.maximum(row_edges[k], lo), 0.0, None) / (hi - lo)
        norm = np.sum(frac * weights)
        out[k] = np.sum(frac * weights * values) / norm if norm > 0 else 0.0
    return out


class StandaloneHelper(IdentityHelper):
    """Minimal DD-version helper used when no ``nc2imas`` instance is available."""

    def __init__(self, dd_version):
        parts = [int(re.sub(r"\D", "", token) or 0) for token in str(dd_version).split(".")[:3]]
        self.dd_version_tuple = tuple(parts + [0] * (3 - len(parts)))
        super().__init__(lambda major, minor=0, patch=0: self.dd_version_tuple >= (major, minor, patch))


def pitch_to_field_direction(pitch_edges, nsjdotb):
    """Return ascending IMAS pitch edges and the NUBEAM->IMAS bin permutation."""
    edges = np.asarray(pitch_edges, float)
    nbins = edges.size - 1
    if nsjdotb >= 0:
        return edges.copy(), np.arange(nbins)
    return -edges[::-1], np.arange(nbins)[::-1]


def phase_space_cells(f_nubeam, species, snapshot, permutation):
    """``f d3v`` (m^-3) and ``f`` (m^-6 s^3) flattened with zone fastest."""
    de = np.diff(species.energy_edges_ev)
    dxi = np.diff(snapshot.pitch_edges)
    f = f_nubeam[:, permutation, :]
    dxi = dxi[permutation]
    fd3v = f * 1.0e6 * de[:, None, None] * dxi[None, :, None] / 2.0
    mass = species.a * AMU
    v_edges = np.sqrt(2.0 * species.energy_edges_ev * EV / mass)
    vol_v = 2.0 * np.pi * (v_edges[1:] ** 3 - v_edges[:-1] ** 3) / 3.0
    d3v = vol_v[:, None, None] * dxi[None, :, None]
    f_si = fd3v / d3v
    order = (2, 0, 1)   # (zone, energy, pitch) with zone fastest in Fortran order
    return (np.ascontiguousarray(np.transpose(fd3v, order)).ravel(order="F"),
            np.ascontiguousarray(np.transpose(f_si, order)).ravel(order="F"))


# Fusion-product source processes (distribution_source_identifier); He4 comes from several
# reactions (D+T, D+He3, ...) that the data cannot tell apart, so it uses the generic nuclear id.
FUSION_PROCESS = {(1, 1): ("D_D_to_T_H", 105), (1, 3): ("D_D_to_T_H", 105),
                  (2, 3): ("D_D_to_3He_n", 106), (2, 4): ("nuclear", 100)}


class DistributionsBuilder:
    """Fill a ``distributions`` IDS from :class:`nubeam_fi.FastIonSnapshot` objects.

    After :meth:`fill`, ``meta`` lists for each ``ids.distribution`` entry the NUBEAM
    species (first snapshot), the distribution kind (``gc``/``particle``) and gyro type.
    """

    def __init__(self, helper):
        self.h = helper
        self.meta = []

    @staticmethod
    def _identifier(node, name, index, description):
        node.name = name
        node.index = int(index)
        node.description = description

    def _coordinate_types(self, space, coords):
        space.coordinates_type.resize(len(coords))
        for k, (name, index, description) in enumerate(coords):
            if self.h.dd_at_least(4, 0, 0):
                self._identifier(space.coordinates_type[k], name, index, description)
            else:
                space.coordinates_type[k] = index

    def _axis_space(self, space, edges, coord, description):
        self._identifier(space.identifier, "secondary_structured", 3, description)
        self._identifier(space.geometry_type, "standard", 0, "Standard geometry")
        self._coordinate_types(space, [coord])
        space.objects_per_dimension.resize(2)
        nodes, bins = space.objects_per_dimension[0], space.objects_per_dimension[1]
        self._identifier(nodes.geometry_content, "node_coordinates", 1, "Node coordinates")
        nodes.object.resize(edges.size)
        for k, value in enumerate(edges):
            nodes.object[k].geometry = np.array([value], dtype=np.float64)
        bins.object.resize(edges.size - 1)
        for k in range(edges.size - 1):
            bins.object[k].nodes = np.array([k + 1, k + 2], dtype=np.int32)

    def _poloidal_space(self, space, snapshot):
        grid, fmap = snapshot.grid, snapshot.flux_map
        conn = grid.connectivity()
        theta0, span = grid.theta0, grid.theta1 - grid.theta0
        self._identifier(space.identifier, "primary_standard", 1,
                         "NUBEAM irregular poloidal Monte-Carlo zone grid")
        self._identifier(space.geometry_type, "standard", 0, "Standard geometry")
        self._coordinate_types(space, [("r", COORD_R, "Major radius"), ("z", COORD_Z, "Vertical coordinate z")])
        space.objects_per_dimension.resize(3)
        nodes, edges, faces = space.objects_per_dimension
        self._identifier(nodes.geometry_content, "node_coordinates", 1, "Node (R, Z) coordinates in m")
        node_r, node_z = [], []
        for ix, frac in conn["nodes"]:
            if ix == 0:
                node_r.append(fmap.r_axis)
                node_z.append(fmap.z_axis)
            else:
                r, z = fmap.rz(ix / grid.nrows, theta0 + float(frac) * span)
                node_r.append(float(r))
                node_z.append(float(z))
        nodes.object.resize(len(conn["nodes"]))
        for k in range(len(conn["nodes"])):
            nodes.object[k].geometry = np.array([node_r[k], node_z[k]], dtype=np.float64)
        edges.object.resize(len(conn["edges"]))
        for k, (a, b) in enumerate(conn["edges"]):
            edges.object[k].nodes = np.array([a + 1, b + 1], dtype=np.int32)
        self._identifier(faces.geometry_content, "face_indices_volume", 31,
                         "(row, zone in row, toroidal zone volume in m^3) of NUBEAM MC zone")
        areas = grid.polygon_area(fmap)
        table = grid.zone_table()
        faces.object.resize(grid.nzones)
        for k in range(grid.nzones):
            obj = faces.object[k]
            obj.nodes = np.array([n + 1 for n in conn["faces"][k]], dtype=np.int32)
            obj.measure = float(areas[k])
            obj.geometry = np.array([table[k, 0] + 1, table[k, 1] + 1, snapshot.zone_volume_m3[k]],
                                    dtype=np.float64)
            obj.boundary.resize(len(conn["face_edges"][k]))
            for m, e in enumerate(conn["face_edges"][k]):
                obj.boundary[m].index = e + 1
                obj.boundary[m].neighbours = np.array([n + 1 for n in conn["neighbours"][k][m]],
                                                      dtype=np.int32)

    def _subsets(self, grid_node, nzones, nenergy, npitch):
        grid_node.grid_subset.resize(4)
        cells, energy, pitch, phase = grid_node.grid_subset
        self._identifier(cells.identifier, "cells", SUBSET_CELLS, "All NUBEAM poloidal MC zones (implicit)")
        cells.dimension = 3
        self._identifier(energy.identifier, "energy_bins", SUBSET_ENERGY_BINS, "All energy bins of space 2")
        energy.dimension = 2
        energy.element.resize(nenergy)
        for k in range(nenergy):
            energy.element[k].object.resize(1)
            obj = energy.element[k].object[0]
            obj.space, obj.dimension, obj.index = 2, 2, k + 1
        self._identifier(pitch.identifier, "pitch_bins", SUBSET_PITCH_BINS, "All pitch bins of space 3")
        pitch.dimension = 2
        pitch.element.resize(npitch)
        for k in range(npitch):
            pitch.element[k].object.resize(1)
            obj = pitch.element[k].object[0]
            obj.space, obj.dimension, obj.index = 3, 2, k + 1
        self._identifier(phase.identifier, "4d_hypervolumes", SUBSET_4D,
                         "All (zone, energy, pitch) phase-space cells; zone index fastest, then energy, then pitch")
        phase.dimension = 4

    def _fill_species(self, dist, species):
        dist.species.type.name = "ion"
        dist.species.type.index = 2
        dist.species.type.description = "Ion species in a single/average state"
        ion = dist.species.ion
        ion.element.resize(1)
        ion.element[0].a = float(species.a)
        self.h.set_nuclear_charge(ion.element[0], species.z)
        ion.element[0].atoms_n = 1
        if int(round(species.z)) == 1:
            ion.z_ion = 1.0
        self.h.set_species_name(ion, species.label)

    def fill(self, ids, snapshots):
        """Fill ``ids`` (homogeneous time) from time-ordered snapshots."""
        keys = [sp.index for sp in snapshots[0].species]
        for snap in snapshots[1:]:
            keys = [key for key in keys if key in {sp.index for sp in snap.species}]
        if not keys:
            raise NubeamInputError("fast-ion species differ between the selected output times")
        pick = lambda snap, key: next(sp for sp in snap.species if sp.index == key)
        ids.ids_properties.homogeneous_time = 1
        ids.time = np.array([snap.time for snap in snapshots], dtype=np.float64)

        kinds = []   # (species index, kind, gyro_type)
        for key in keys:
            kinds.append((key, "gc", 2))
            if all(pick(snap, key).particle is not None for snap in snapshots):
                kinds.append((key, "particle", 1))
        ids.distribution.resize(len(kinds))
        nt = len(snapshots)
        # Lost-orbit markers hold ion positions, so they belong to the particle-position entry.
        marker_kind = {key: ("particle" if (key, "particle", 1) in kinds else "gc") for key in keys}
        self.meta = []
        for dist, (key, kind, gyro) in zip(ids.distribution, kinds):
            first = pick(snapshots[0], key)
            self.meta.append(dict(index=key, kind=kind, gyro=gyro, species=first))
            self._fill_species(dist, first)
            dist.gyro_type = gyro
            dist.is_delta_f = 0
            self._fill_process(dist, first)
            dist.ggd.resize(nt)
            dist.global_quantities.resize(nt)
            dist.profiles_1d.resize(nt)
            for it, snap in enumerate(snapshots):
                sp = pick(snap, key)
                self._fill_time_slice(dist, it, snap, sp, kind)
                if kind == marker_kind[key]:
                    self._fill_lost_markers(dist, it, snap, sp, nt)

    def _fill_process(self, dist, species):
        process = dist.process
        process.resize(1)
        if species.fusion_product:
            name, index = FUSION_PROCESS.get((int(round(species.z)), int(round(species.a))), ("nuclear", 100))
            self._identifier(process[0].type, name, index, "Source from nuclear reaction (fusion product)")
            return
        self._identifier(process[0].type, "NBI", NBI_PROCESS_INDEX, "Source from neutral beam injection")
        self._identifier(process[0].nbi_energy, "total", 0, "Sum over all NBI energy components")

    def _fill_lost_markers(self, dist, it, snap, sp, nt):
        """NUBEAM lost-orbit records as markers (weights = real particles per record).

        Coordinates: R (m), Z (m), kinetic energy (eV), pitch along B (IMAS sign),
        and three private ones: loss time (s), NBI unit (0 = not from a beam) and a
        prompt (first-orbit) loss flag (1/0).
        """
        lost = snap.lost
        if lost is None:
            return
        select = lost.species_index == sp.index
        if not select.any():
            return
        if len(dist.markers) != nt:
            dist.markers.resize(nt)
        marker = dist.markers[it]
        coords = [("r", COORD_R, "Major radius, m"), ("z", COORD_Z, "Vertical coordinate z, m"),
                  ("energy_kinetic", COORD_ENERGY, "Kinetic energy, eV"),
                  ("pitch", COORD_PITCH, "v_parallel/v along the magnetic field"),
                  ("time_lost", COORD_TIME_LOST, "Time at which the ion was lost, s"),
                  ("nbi_unit", COORD_NBI_UNIT, "NBI unit the ion came from (0: not from a beam)"),
                  ("prompt_loss", COORD_PROMPT, "1: prompt (first-orbit) loss, 0: other")]
        marker.coordinate_identifier.resize(len(coords))
        for k, (name, index, description) in enumerate(coords):
            self._identifier(marker.coordinate_identifier[k], name, index, description)
        marker.weights = lost.weight[select]
        marker.positions = np.column_stack([
            lost.r[select], lost.z[select], lost.energy_ev[select], snap.nsjdotb * lost.pitch[select],
            lost.time[select], lost.nbi_unit[select].astype(float), lost.prompt[select].astype(float)])

    def _fill_time_slice(self, dist, it, snap, sp, kind):
        f_n = sp.gc if kind == "gc" else sp.particle
        edges, perm = pitch_to_field_direction(snap.pitch_edges, snap.nsjdotb)
        fd3v, f_si = phase_space_cells(f_n, sp, snap, perm)
        ggd = dist.ggd[it]
        grid = ggd.grid
        self._identifier(grid.identifier, "nubeam_mc_phase_space", GRID_INDEX,
                         "NUBEAM MC zones (R,Z) x energy x pitch(v_parallel/v along B)")
        grid.space.resize(3)
        self._poloidal_space(grid.space[0], snap)
        self._axis_space(grid.space[1], sp.energy_edges_ev, ("energy_kinetic", COORD_ENERGY, "Kinetic energy in eV"),
                         "Kinetic energy bins")
        self._axis_space(grid.space[2], edges, ("pitch", COORD_PITCH, "v_parallel/v along the magnetic field"),
                         "Pitch bins, v_parallel/v with parallel along B (IMAS convention)")
        self._subsets(grid, snap.grid.nzones, sp.energy_edges_ev.size - 1, edges.size - 1)
        for aos, values in ((ggd.expansion_fd3v, fd3v), (ggd.expansion, f_si)):
            aos.resize(1)
            aos[0].grid_subset.resize(1)
            node = aos[0].grid_subset[0]
            node.grid_index = GRID_INDEX
            node.grid_subset_index = SUBSET_4D
            node.values = values

        moments = species_moments(snap, sp, "gc" if kind == "gc" else "particle")
        gq = dist.global_quantities[it]
        gq.particles_n = moments["N"]
        gq.particles_fast_n = moments["N"]
        gq.energy = moments["E"]
        gq.energy_fast = moments["E"]
        gq.energy_fast_parallel = moments["E_par"]

        # Row-averaged (flux-zone) moments on x = rho_tor_norm zone-row centres.
        prof = dist.profiles_1d[it]
        rows = np.repeat(np.arange(snap.grid.nrows), snap.grid.zones_per_row)
        vol = snap.zone_volume_m3
        row_vol = np.bincount(rows, weights=vol)

        def row_average(zone_values):
            return np.bincount(rows, weights=zone_values * vol) / row_vol

        prof.grid.rho_tor_norm = (np.arange(snap.grid.nrows) + 0.5) / snap.grid.nrows
        dens = row_average(moments["n"])
        prof.density = dens
        prof.density_fast = dens
        w_total = row_average(moments["energy"])
        w_par = row_average(moments["energy_parallel"])
        prof.pressure = 2.0 / 3.0 * w_total
        prof.pressure_fast = 2.0 / 3.0 * w_total
        prof.pressure_fast_parallel = 2.0 * w_par
        if self.h.dd_at_least(4, 0, 0):
            # DD3.39 has no perpendicular pressure node.
            prof.pressure_fast_perpendicular = w_total - w_par


class imas_distributions(IdsWriter):
    """Writer registered as the ``distributions`` IDS in ``convert.py``."""

    def _locate_case(self):
        cdf = Path(self.cdf.ncfile)
        run_id = re.sub(r"\.cdf$", "", cdf.name, flags=re.I)
        return cdf.parent, run_id

    def build(self):
        self.logger.info("\ndistributions")
        case_dir, run_id = self._locate_case()
        trdat = case_dir / f"{run_id}TR.DAT"
        if not trdat.is_file():
            self.logger.info(f"distributions: {trdat.name} not found next to the CDF; skipping")
            self.ctx.mapping.skipped(f'{trdat.name} not found next to the CDF', target='distributions')
            return []
        try:
            discover_nubeam_outputs(case_dir, run_id, trdat)
        except FileNotFoundError as exc:
            self.logger.warning(f"distributions: requested NUBEAM output unavailable ({exc}); skipping")
            self.ctx.mapping.skipped(f'requested NUBEAM output unavailable: {exc}', target='distributions')
            return []
        except ValueError as exc:
            self.logger.info(f"distributions: no NUBEAM fast-ion output requested ({exc}); skipping")
            self.ctx.mapping.skipped(f'no NUBEAM fast-ion output requested: {exc}', target='distributions')
            return []
        snapshots = read_fast_ion_snapshots(case_dir, run_id, trdat, logger=self.logger)
        if not snapshots:
            self.logger.warning("distributions: no DATA/FI_TAR pair shares an output time; skipping")
            self.ctx.mapping.skipped('no DATA/FI_TAR pair shares an output time', target='distributions')
            return []
        snapshots.sort(key=lambda snap: snap.time)

        ids = self.dd.new_ids("distributions")
        ids = self.fill.fill_properties(ids)
        self.fill.fill_code_properties(ids)
        builder = DistributionsBuilder(IdentityHelper(self.dd.dd_at_least))
        builder.fill(ids, snapshots)
        self._fill_machine_fields(ids, snapshots)
        self._fill_cdf_power_losses(ids, snapshots, builder.meta)
        self._fill_cdf_torque_current(ids, snapshots, builder.meta)
        self._fill_cdf_balances(ids, snapshots, builder.meta)
        sign = snapshots[0].nsjdotb
        ids.ids_properties.comment = (
            f"{ids.ids_properties.comment}; NUBEAM fast-ion distribution from DATA/FI_TAR at "
            f"{', '.join(f'{snap.time:g}' for snap in snapshots)} s; pitch = v_parallel/v along B "
            f"(NUBEAM pitch is along Ip; nsjdotb={sign:+d} from {snapshots[0].nsjdotb_source}); "
            "expansion_fd3v = f*d3v in m^-3 on the 4d_hypervolumes subset (zone fastest, then energy, pitch)")
        self.logger.info(f"distributions: built {len(ids.distribution)} distribution(s) "
                         f"at {len(snapshots)} time(s)")
        return [ids]

    def _cdf_at(self, name, times):
        """CDF ``TIME``-dependent scalar interpolated to ``times``; ``None`` if absent."""
        values = self.cdf.get_ncvar(name)
        if not self.fill.is_valid(values):
            return None
        return np.interp(times, self.grid.timesc, np.asarray(values, dtype=np.float64))

    @staticmethod
    def _sole_beam(meta):
        """True when exactly one beam (non-fusion) species is written: only then are the
        unsuffixed TRANSP beam totals unambiguous."""
        return len({m["species"].index for m in meta if not m["species"].fusion_product}) == 1

    def _series_for(self, stem, species, sole, times):
        """Per-species ``TIME`` series via the species' candidate CDF names (``None`` if absent)."""
        for name in species.cdf_names(stem, allow_plain=sole):
            values = self._cdf_at(name, times)
            if values is not None:
                return values
        return None

    def _profile_for(self, stem, species, sole, times):
        for name in species.cdf_names(stem, allow_plain=sole):
            values = self._cdf_profile_at(name, times)
            if values is not None:
                return values
        return None

    def _beam_total_profile(self, stem, species, sole, times):
        """Unsuffixed all-beam profile (``CURB``, ``TQBE``, ...): beam species, one beam only."""
        if species.fusion_product or not sole:
            return None
        return self._cdf_profile_at(stem, times)

    def _fill_cdf_power_losses(self, ids, snapshots, meta):
        """Orbit-loss and charge-exchange loss powers (W, positive = loss) of each beam species.

        TRANSP ``BPLIM_<sp>`` is the fast-ion orbit loss power; the CX loss is the
        sum of ``BPCXI_<sp>`` (inside the plasma) and ``BPCXX_<sp>`` (outside).
        They are time-dependent CDF scalars interpolated to the FI output times.
        Fusion products have no power counterpart in TRANSP (only particle losses).
        Skipped for DD versions without these nodes.
        """
        times = np.array([snap.time for snap in snapshots])
        sole = self._sole_beam(meta)
        for dist, m in zip(ids.distribution, meta):
            spec = m["species"]
            orbit = self._series_for("BPLIM", spec, sole, times)
            cx_in, cx_out = (self._series_for("BPCXI", spec, sole, times),
                             self._series_for("BPCXX", spec, sole, times))
            for it, gq in enumerate(dist.global_quantities):
                if orbit is not None and hasattr(gq, "power_first_orbit"):
                    gq.power_first_orbit = float(orbit[it])
                if cx_in is not None and cx_out is not None and hasattr(gq, "power_charge_exchange"):
                    gq.power_charge_exchange = float(cx_in[it] + cx_out[it])

    def _cdf_profile_at(self, name, times):
        """CDF ``(TIME3, X)`` profile linearly interpolated in time to ``times``."""
        values = self.cdf.get_ncvar(name)
        if not self.fill.is_valid(values):
            return None
        values = np.asarray(values, dtype=np.float64)
        time3 = np.asarray(self.cdf.get_ncvar("TIME3"), dtype=np.float64)
        if values.ndim != 2 or values.shape[0] != time3.size:
            return None
        return np.array([[np.interp(t, time3, values[:, ix]) for ix in range(values.shape[1])] for t in times])

    def _cdf_grid(self, times):
        """``(zone_edges, DVOL, DAREA)`` of the TRANSP grid at ``times``, or ``None``."""
        dvol, darea = self._cdf_profile_at("DVOL", times), self._cdf_profile_at("DAREA", times)
        xb = self.cdf.get_ncvar("XB")
        if dvol is None or darea is None or not self.fill.is_valid(xb):
            return None
        xb = np.asarray(xb, dtype=np.float64)
        return np.concatenate([[0.0], np.asarray(xb[0] if xb.ndim == 2 else xb)]), dvol, darea

    @staticmethod
    def _set_first(node, names, value):
        for name in names:
            if hasattr(node, name):
                setattr(node, name, value)
                return True
        return False

    def _fill_cdf_torque_current(self, ids, snapshots, meta):
        """Fast-ion J x B torque and beam-driven current from TRANSP CDF profiles.

        * ``TQJXB_<sp>`` (N m cm^-3) is the toroidal torque of the fast-ion radial
          current: global ``torque_tor_j_radial`` = sum(TQJXB * DVOL) in N m and the
          profile ``torque_phi_j_radial`` in N m^-2, regridded to the MC rows.
        * ``CURB`` (A cm^-2) is the shielded beam-driven current density (includes the
          electron back-current).  It is the total over all beams, so it is attached
          only for a sole beam species: global ``current_phi`` =
          sum(CURB * DAREA) in A and the profile in A m^-2.
        The radial fast-ion currents ``CURBRORB``/``CURBRABD`` have no counterpart
        node in ``distributions``; their effect is carried by the J x B torque.
        """
        times = np.array([snap.time for snap in snapshots])
        grid = self._cdf_grid(times)
        if grid is None:
            return
        edges, dvol, darea = grid
        sole = self._sole_beam(meta)
        for dist, m in zip(ids.distribution, meta):
            spec = m["species"]
            current = self._beam_total_profile("CURB", spec, sole, times)
            unshielded = self._beam_total_profile("UCURB", spec, sole, times)
            torque = self._profile_for("TQJXB", spec, sole, times)
            torque_scalar = None if torque is not None else self._series_for("BPHXB", spec, sole, times)
            for it, snap in enumerate(snapshots):
                gq, prof = dist.global_quantities[it], dist.profiles_1d[it]
                nrows = snap.grid.nrows
                if torque is not None:
                    self._set_first(gq, ("torque_tor_j_radial",), float(np.sum(torque[it] * dvol[it])))
                    self._set_first(prof, ("torque_phi_j_radial", "torque_tor_j_radial"),
                                    regrid_to_rows(torque[it], edges, dvol[it], nrows) * 1.0e6)
                elif torque_scalar is not None:
                    # No profile (fusion products): use TRANSP's scalar J x B torque.
                    self._set_first(gq, ("torque_tor_j_radial",), float(torque_scalar[it]))
                if current is not None:
                    self._set_first(gq, ("current_phi", "current_tor"), float(np.sum(current[it] * darea[it])))
                    self._set_first(prof, ("current_phi", "current_tor"),
                                    regrid_to_rows(current[it], edges, darea[it], nrows) * 1.0e4)
                if unshielded is not None:
                    # Fast-ion current without the electron back-current (TRANSP UCURB).
                    self._set_first(prof, ("current_fast_phi", "current_fast_tor"),
                                    regrid_to_rows(unshielded[it], edges, darea[it], nrows) * 1.0e4)

    @staticmethod
    def _pick(node, *names):
        for name in names:
            if hasattr(node, name):
                return getattr(node, name)
        return None

    def _fill_cdf_balances(self, ids, snapshots, meta):
        """Collisional exchange, thermalization, NBI source and trapped fraction.

        Global values come from TRANSP time series (W, s^-1, or N m from a volume
        integral of a torque density); profiles are TRANSP zone profiles regridded to
        the MC zone rows (densities in SI).  Beam power to electrons/ions
        (``BPTE``/``BPTI``, ``PBE``/``PBI``) and the matching torques
        (``TQBE``/``TQBI``) fill ``collisions``; ``BPTH``/``BSTH``/``TQBTH`` fill
        ``thermalization``; the MC deposition (``SBDEPMC``, ``PBDEPMC``, ``BDEP``, and
        ``BPHDEP0`` on the particle-position entry only) fills one NBI ``source``.
        ``BTRAP`` (banana fraction) times the fast-ion density gives the trapped density.
        The thermal-ion exchange is one aggregate entry (TRANSP does not split it).

        Names follow TRANSP's species suffixes (``nubeam_fi.cdf_species_suffix``).  Fusion
        products only get the quantities that exist in TRANSP's F family (thermalization
        particles, trapped fraction, density).  Unsuffixed beam totals (``TQBE``, ``CURB``,
        ...) are attached only when exactly one beam species is written.
        """
        times = np.array([snap.time for snap in snapshots])
        grid = self._cdf_grid(times)
        if grid is None:
            return
        edges, dvol, _darea = grid
        sole = self._sole_beam(meta)
        for dist, m in zip(ids.distribution, meta):
            spec = m["species"]
            ser = lambda stem: self._series_for(stem, spec, sole, times)
            prof = lambda stem: self._profile_for(stem, spec, sole, times)
            total = lambda stem: self._beam_total_profile(stem, spec, sole, times)
            integral = lambda q: None if q is None else np.sum(q * dvol, axis=1)
            p_e, p_i, p_th = ser("BPTE"), ser("BPTI"), ser("BPTH")
            n_th = ser("BSTH")
            tq_e, tq_i, tq_th = total("TQBE"), total("TQBI"), total("TQBTH")
            t_e, t_i, t_th = integral(tq_e), integral(tq_i), integral(tq_th)
            # Deposited torque is only taken at the ion position (BPHDEP0, = BPHDP).  The
            # guiding-centre value BPHDEPGC (and BPHINJS) is ~1e13 in TRANSP, evidently not
            # in N m, so it is not used.
            src_n, src_p = ser("SBDEPMC"), ser("PBDEPMC")
            src_t = ser("BPHDEP0") if int(dist.gyro_type) == 1 else None
            q_e, q_i, q_th = prof("PBE"), prof("PBI"), prof("PBTH")
            dep, bdens, trap = prof("BDEP"), prof("BDENS"), prof("BTRAP")
            for it, snap in enumerate(snapshots):
                gq, pr = dist.global_quantities[it], dist.profiles_1d[it]
                nrows = snap.grid.nrows
                row = lambda values, scale: regrid_to_rows(values[it], edges, dvol[it], nrows) * scale
                # collisions: electrons (global) and one aggregate thermal-ion entry
                g_el, g_coll = self._pick(gq.collisions, "electrons"), gq.collisions
                if p_e is not None:
                    self._set_first(g_el, ("power_thermal",), float(p_e[it]))
                if t_e is not None:
                    self._set_first(g_el, ("torque_thermal_phi", "torque_thermal_tor"), float(t_e[it]))
                if p_i is not None or t_i is not None:
                    g_coll.ion.resize(1)
                    self.fill.set_species_name(g_coll.ion[0], "thermal ions (all species)")
                    if p_i is not None:
                        g_coll.ion[0].power_thermal = float(p_i[it])
                    if t_i is not None:
                        self._set_first(g_coll.ion[0], ("torque_thermal_phi", "torque_thermal_tor"), float(t_i[it]))
                p_el = self._pick(pr.collisions, "electrons") if hasattr(pr.collisions, "electrons") else pr.collisions
                if q_e is not None:
                    self._set_first(p_el, ("power_thermal",), row(q_e, 1.0e6))
                if tq_e is not None:
                    self._set_first(p_el, ("torque_thermal_phi", "torque_thermal_tor"), row(tq_e, 1.0e6))
                if q_i is not None or tq_i is not None:
                    pr.collisions.ion.resize(1)
                    self.fill.set_species_name(pr.collisions.ion[0], "thermal ions (all species)")
                    if q_i is not None:
                        pr.collisions.ion[0].power_thermal = row(q_i, 1.0e6)
                    if tq_i is not None:
                        self._set_first(pr.collisions.ion[0], ("torque_thermal_phi", "torque_thermal_tor"),
                                        row(tq_i, 1.0e6))
                # thermalization
                g_th = self._pick(gq, "thermalization", "thermalisation")
                if g_th is not None:
                    for value, names in ((n_th, ("particles",)), (p_th, ("power",))):
                        if value is not None:
                            self._set_first(g_th, names, float(value[it]))
                    if t_th is not None:
                        self._set_first(g_th, ("torque",), float(t_th[it]))
                p_thm = self._pick(pr, "thermalization", "thermalisation")
                if p_thm is not None:
                    if q_th is not None:
                        self._set_first(p_thm, ("energy",), row(q_th, 1.0e6))
                    if tq_th is not None:
                        self._set_first(p_thm, ("momentum_phi", "momentum_tor"), row(tq_th, 1.0e6))
                # NBI source term
                if any(v is not None for v in (src_n, src_p, src_t)):
                    gq.source.resize(1)
                    node = gq.source[0]
                    DistributionsBuilder._identifier(node.identifier.type, "particle_source_process", 2,
                                     "Particle source process (distribution/process)")
                    node.identifier.process_index = 1
                    if src_n is not None:
                        node.particles = float(src_n[it])
                    if src_p is not None:
                        node.power = float(src_p[it])
                    if src_t is not None:
                        self._set_first(node, ("torque_phi", "torque_tor"), float(src_t[it]))
                if dep is not None:
                    pr.source.resize(1)
                    DistributionsBuilder._identifier(pr.source[0].identifier.type, "particle_source_process", 2,
                                     "Particle source process (distribution/process)")
                    pr.source[0].identifier.process_index = 1
                    pr.source[0].particles = row(dep, 1.0e6)
                # Trapped fast-ion density = TRANSP banana fraction (weighted by TRANSP's fast-ion
                # density) x the density of this distribution.  TRANSP's 1-D BDENS is the density
                # at the ion position (it matches the FLR distribution to 1-3% in every row), so
                # it differs strongly from the guiding-centre density in the outer rows and is
                # not used as the density itself.
                if trap is not None and hasattr(pr, "trapped") and len(pr.density_fast):
                    weights = dvol[it] * (bdens[it] if bdens is not None else 1.0)
                    fraction = regrid_to_rows(trap[it], edges, weights, nrows)
                    pr.trapped.density_fast = fraction * np.asarray(pr.density_fast)

    def _fill_machine_fields(self, ids, snapshots):
        times = np.array([snap.time for snap in snapshots])
        ids.magnetic_axis.r = np.array([snap.flux_map.r_axis for snap in snapshots])
        ids.magnetic_axis.z = np.array([snap.flux_map.z_axis for snap in snapshots])
        reference = self.signs.vacuum_field_reference(times)
        if reference is not None:
            ids.vacuum_toroidal_field.r0, ids.vacuum_toroidal_field.b0 = reference


if __name__ == "__main__":
    print("Please use: python convert.py")
