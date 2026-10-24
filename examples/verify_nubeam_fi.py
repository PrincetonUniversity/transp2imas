#!/usr/bin/python
"""Verify the NUBEAM fast-ion distribution -> IMAS ``distributions`` conversion.

Everything except ``LocalSampleTests`` runs on small synthetic TRANSP outputs
(AC-format ``DATA<n>``, ``FI_TAR.GZ<n>`` archives, ``TR.DAT``) generated here,
so no TRANSP/NUBEAM libraries and no large sample are needed.  The local NSTX
case (``examples/data/nstx/120967E03_send``, never committed) is checked too
when it is present.

    python examples/verify_nubeam_fi.py [--quick] [--sample DIR]
"""

import argparse
import io
import math
import sys
import os
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import imas_distributions as idm
import nubeam_fi as nf
from nubeam_outputs import discover_nubeam_outputs
from transp_acfile import ACFile
from mapping_status import MappingLedger
from utilities import vacuum_field_reference

try:
    import imas
except ImportError:      # pragma: no cover - environment dependent
    imas = None
try:
    import netCDF4
except ImportError:      # pragma: no cover - environment dependent
    netCDF4 = None

ALPHABET = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ<>"
SAMPLE_DIR = ROOT / "examples" / "data" / "nstx" / "120967E03_send"
SAMPLE_ARGS = {"dir": SAMPLE_DIR}


# ---------------------------------------------------------------- AC encoders
# Registry entries exercised here: dist_orbit_cx_losses, dist_collisions, dist_thermalization, dist_nbi_source,
# dist_trapped_density, dist_torque_current
def b64(value, width):
    out = ""
    for _ in range(width):
        out = ALPHABET[value % 64] + out
        value //= 64
    return out


def enc_real8(z):
    if z == 0.0:
        return "ns0000000000"
    exp10 = math.floor(math.log10(abs(z))) + 1
    mant = abs(z) / 10.0 ** exp10
    ic1 = int(mant * 1e9)
    ic2 = max(0, min(999999999, int((mant - ic1 * 1e-9) * 1e18)))   # residual can round negative
    return b64(exp10 + 500 + (1000 if z > 0 else 0), 2) + b64(ic1, 5) + b64(ic2, 5)


def enc_real4(z):
    if z == 0.0:
        return "0     "
    frac, exp2 = math.frexp(abs(z))          # abs(z) = frac * 2**exp2, 0.5 <= frac < 1
    mantissa = int(round((2.0 * frac - 1.0) * 2 ** 27))
    bits = ((1 if z > 0 else 0) << 35) | ((exp2 - 1 + 127) << 27) | mantissa
    return b64(bits, 6)


def enc_array(values, width, encode):
    fields, zeros = [], 0

    def flush():
        nonlocal zeros
        if zeros > 1:
            fields.append(("_0:" + b64(zeros, 3)).ljust(width, "0"))
        elif zeros == 1:
            fields.append(encode(0.0))
        zeros = 0
    for value in values:
        if value == 0.0:
            zeros += 1
        else:
            flush()
            fields.append(encode(float(value)))
    flush()
    per_line = 72 // width
    return ["".join(fields[i:i + per_line]) for i in range(0, len(fields), per_line)]


def ac_item(name, value, kind=None):
    """Return AC-file lines for a scalar or Fortran-ordered numpy array."""
    value = np.asarray(value)
    kind = kind or ("I" if value.dtype.kind in "iu" else "L" if value.dtype.kind == "b" else "D")
    if value.ndim == 0:
        if kind == "I":
            text = (" " if value >= 0 else "-") + b64(abs(int(value)), 1)
        elif kind == "L":
            text = "T" if value else "F"
        else:
            text = enc_real8(float(value))
        return [f"*{kind}0 {name:<8} {text}"]
    lines = [f"*{kind}{value.ndim} {name}", "".join(" " + b64(d, 1 if d < 64 else 2) for d in value.shape)]
    flat = value.ravel(order="F")
    if kind == "L":
        text = "".join("T" if v else "F" for v in flat)
        return lines + [text[i:i + 80] for i in range(0, len(text), 80)]
    if kind == "D":
        lines += enc_array(flat, 12, enc_real8)
    elif kind == "R":
        lines += enc_array(flat, 6, enc_real4)
    elif kind == "I":
        lines += ["".join((" " if v >= 0 else "-") + b64(abs(int(v)), 2) for v in flat[i:i + 20])
                  for i in range(0, flat.size, 20)]
    return lines


# ----------------------------------------------------------- synthetic case
R0, MINOR, KAPPA = 1.0, 0.5, 1.6
ROWS = 3
ZONES = np.array([4, 8, 12])
NE, NP = 5, 6
TIMES = (0.1, 0.2, 0.3)


def circular_equilibrium():
    rho = np.linspace(0.0, 1.0, 11)
    theta = np.linspace(-np.pi, np.pi, 33)
    r_geo = R0 + MINOR * rho[None, :] * np.cos(theta)[:, None]
    z_geo = KAPPA * MINOR * rho[None, :] * np.sin(theta)[:, None]
    return rho, theta, r_geo, z_geo


def analytic_row_volumes():
    edges = np.arange(ROWS + 1) / ROWS
    return 2.0 * np.pi * R0 * np.pi * MINOR ** 2 * KAPPA * np.diff(edges ** 2)


def synthetic_distribution(seed=1):
    rng = np.random.default_rng(seed)
    nz = int(ZONES.sum())
    f = rng.random((NE, NP, nz)) * 1.0e8
    # Make the distribution visibly asymmetric in pitch.
    f *= np.linspace(0.2, 1.8, NP)[None, :, None]
    return f


LOST_N = 5


def lost_items(time):
    """Synthetic lost-orbit records: slot 0 is the unused padding NUBEAM leaves."""
    pad = lambda values: np.concatenate([[0.0], values])
    k = np.arange(1, LOST_N + 1, dtype=float)
    return [ac_item("GTIME1_LOST", time - 1.0e-3), ac_item("GTIME2_LOST", time),
            ac_item("GTIME_LOST", pad(time - 1.0e-3 + 1.0e-4 * k)), ac_item("YNBSCE_LOST", pad(1.0 + k % 2)),
            ac_item("YNBIEN_LOST", pad(1.0 + k % 3)), ac_item("WGHT_LOST", pad(1.0e12 * k)),
            ac_item("RMJION_LOST", pad(140.0 + k)), ac_item("YYION_LOST", pad(-50.0 + 10.0 * k)),
            ac_item("YZE_LOST", pad(1.0e4 * k)), ac_item("XKSID_LOST", pad(0.1 * k)),
            ac_item("GOOSE_LOST", pad(np.where(k % 2 == 0, 1.0001, 7.5))), ac_item("SPECIES_LOST", pad(np.ones(LOST_N)))]


def data_file_text(time, f, volumes_m3, flr=True, species=None):
    """DATA text; ``species`` is a list of ``(A, Z, fusion_product, f)`` (default: one D beam)."""
    species = species or [(2.0, 1.0, False, f)]
    ns = len(species)
    nbeam = sum(1 for sp in species if not sp[2])
    cumulative = np.concatenate([[0], np.cumsum(ZONES), np.full(4, int(ZONES.sum()) + 7)])
    centres = (np.arange(NE) + 0.5) * 20000.0
    edges = np.arange(NE + 1) * 20000.0
    pitch = np.linspace(-1.0, 1.0, NP + 1)
    pitch = 0.5 * (pitch[1:] + pitch[:-1])
    fbm = np.stack([sp[3] for sp in species], axis=-1)                       # (NE, NP, nz, ns)
    items = [
        ac_item("TA", time), ac_item("NSFAST", np.array(ns)), ac_item("NSBEAM", np.array(nbeam)),
        ac_item("ABEAMS", np.array([sp[0] for sp in species])), ac_item("XZBEAMS", np.array([sp[1] for sp in species])),
        ac_item("NLFPROD", np.array([sp[2] for sp in species])),
        ac_item("NZNBMEA", np.array([NE] * ns + [0])), ac_item("NZNBMA", np.array(NP)),
        ac_item("EFBM", np.tile(centres[:, None], (1, ns))), ac_item("EFBMB", np.tile(edges[:, None], (1, ns))),
        ac_item("XXKSID", pitch), ac_item("NFBZNSI", np.array(int(ZONES.sum()))),
        ac_item("NTHZSM", cumulative), ac_item("LCENTR", np.array(2)), ac_item("NZONE_FB", np.array(ROWS)),
        ac_item("THBDY0", -np.pi), ac_item("THBDY1", np.pi),
        ac_item("BMVOL", np.concatenate([volumes_m3 * 1.0e6, np.zeros(9)])),
        ac_item("FBM", fbm, kind="R"),
    ]
    if flr:
        items.append(ac_item("FBM_PTCL", fbm[..., None] * 0.9, kind="R"))
    items += lost_items(time)
    lines = [f"RUNX{1:5d}{time:16.7E}     1     5"] + [line for item in items for line in item]
    return "\n".join(lines) + "\n"


def make_plasma_state_bytes(netcdf_module):
    rho, theta, r_geo, z_geo = circular_equilibrium()
    path = Path(tempfile.mkdtemp()) / "ps.cdf"
    with netcdf_module.Dataset(path, "w", format="NETCDF3_CLASSIC") as nc:
        nc.createDimension("nrho", rho.size)
        nc.createDimension("nth", theta.size)
        for name, array, dims in (("rho_eq", rho, ("nrho",)), ("th_eq", theta, ("nth",)),
                                  ("R_geo", r_geo, ("nth", "nrho")), ("Z_geo", z_geo, ("nth", "nrho"))):
            nc.createVariable(name, "f8", dims)[:] = array
        for name, value in (("R_axis", R0), ("Z_axis", 0.0), ("B_axis_vac", 0.5)):
            nc.createVariable(name, "f8", ())[...] = value
    return path.read_bytes()


def make_xplasma_bytes(netcdf_module, bphi_ccw, jphi_ccw):
    path = Path(tempfile.mkdtemp()) / "xp.cdf"
    with netcdf_module.Dataset(path, "w", format="NETCDF3_CLASSIC") as nc:
        nc.createVariable("__Bphi_CCW", "i4", ())[...] = bphi_ccw
        nc.createVariable("__Jphi_CCW", "i4", ())[...] = jphi_ccw
    return path.read_bytes()


def write_archive(path, run_id, ps_bytes, xp_bytes=None):
    members = {f"./{run_id}_ps_ts1_state.cdf": ps_bytes,
               f"./{run_id}_nbi_ptcl_state_cpu0_32.cdf": b"NOT A NETCDF FILE - must never be parsed",
               f"./{run_id}_debug_nbi_ptcl_state_cpu1_32.cdf": b"ALSO NOT NETCDF"}
    if xp_bytes is not None:
        members[f"./{run_id}_debug_xplasma_nubeam.cdf"] = xp_bytes
    with tarfile.open(path, "w:gz") as tar:
        for name, payload in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            tar.addfile(info, io.BytesIO(payload))


def make_case(directory, outtim=TIMES, fi_outtim=TIMES, jccw=True, bccw=False, with_xplasma=True,
              species=None):
    """Create a synthetic TRANSP case directory; return the DATA distributions."""
    run = "RUNX"
    directory = Path(directory)
    fmt = lambda values: ",".join(f"{v:g}" for v in values)
    (directory / f"{run}TR.DAT").write_text(
        f"NLJCCW={'.TRUE.' if jccw else '.FALSE.'}\nNLBCCW={'.TRUE.' if bccw else '.FALSE.'}\n"
        f"outtim={fmt(outtim)}\nfi_outtim={fmt(fi_outtim)}\n")
    rho, theta, r_geo, z_geo = circular_equilibrium()
    volumes = nf.McGrid(ROWS, ZONES, -np.pi, np.pi).zone_volumes(
        nf.FluxSurfaceMap(rho, theta, r_geo, z_geo, R0, 0.0))
    dists = {}
    for index, time in enumerate(outtim, start=1):
        f = synthetic_distribution(seed=index)
        dists[time] = f
        spec = None if species is None else [(a, z, fp, synthetic_distribution(seed=10 * index + k))
                                             for k, (a, z, fp) in enumerate(species)]
        (directory / f"{run}.DATA{index}").write_text(data_file_text(time, f, volumes, species=spec))
    if netCDF4 is not None:
        ps = make_plasma_state_bytes(netCDF4)
        xp = make_xplasma_bytes(netCDF4, 1 if bccw else -1, 1 if jccw else -1) if with_xplasma else None
    for index, _time in enumerate(fi_outtim, start=1):
        if netCDF4 is not None:
            write_archive(directory / f"{run}_FI_TAR.GZ{index}", run, ps, xp)
        else:
            (directory / f"{run}_FI_TAR.GZ{index}").write_bytes(b"")
    return run, dists, volumes


class AcCodecTests(unittest.TestCase):
    def test_known_encodings_from_transp_files(self):
        from transp_acfile import _decode_real8
        self.assertEqual(_decode_real8("ns0000000000"), 0.0)
        self.assertAlmostEqual(_decode_real8("nu79MU000000"), 12.0, places=12)

    def test_round_trip_all_types_with_zero_runs(self):
        rng = np.random.default_rng(0)
        real = rng.normal(size=(4, 5)) * 10.0 ** rng.integers(-8, 9, size=(4, 5))
        real[1, 1:4] = 0.0
        real[3, :] = 0.0
        integers = np.array([[0, 5, -7], [123, -4000, 63]])
        lines = ["RUN 1 0.5 1 5"] + ac_item("SC", 3.25) + ac_item("NI", np.array(11)) + \
            ac_item("FLAG", np.array(True)) + ac_item("A8", real) + ac_item("A4", real, kind="R") + \
            ac_item("IA", integers)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.DATA1"
            path.write_text("\n".join(lines) + "\n")
            ac = ACFile(path)
            self.assertAlmostEqual(ac.get("SC"), 3.25, places=12)
            self.assertEqual(ac.get("NI"), 11)
            self.assertTrue(ac.get("FLAG"))
            np.testing.assert_allclose(ac.get("A8"), real, rtol=1e-12, atol=0)
            np.testing.assert_allclose(ac.get("A4"), real, rtol=2e-8, atol=0)
            np.testing.assert_array_equal(ac.get("IA"), integers)
            self.assertEqual(ac.shape("A8"), (4, 5))


class GridTests(unittest.TestCase):
    def setUp(self):
        rho, theta, r_geo, z_geo = circular_equilibrium()
        self.fmap = nf.FluxSurfaceMap(rho, theta, r_geo, z_geo, R0, 0.0)
        self.grid = nf.McGrid(ROWS, ZONES, -np.pi, np.pi)

    def test_flux_surface_map_matches_analytic_surfaces(self):
        r, z = self.fmap.rz(0.7, np.array([-3.0, 0.3, 3.1]))
        theta = np.array([-3.0, 0.3, 3.1])
        np.testing.assert_allclose(r, R0 + MINOR * 0.7 * np.cos(theta), atol=2e-4)
        np.testing.assert_allclose(z, KAPPA * MINOR * 0.7 * np.sin(theta), atol=2e-4)

    def test_zone_volumes_match_analytic_torus_rows(self):
        volumes = self.grid.zone_volumes(self.fmap)
        rows = np.repeat(np.arange(ROWS), ZONES)
        np.testing.assert_allclose(np.bincount(rows, weights=volumes), analytic_row_volumes(), rtol=2e-4)

    def test_polygon_areas_sum_to_plasma_cross_section(self):
        total = self.grid.polygon_area(self.fmap).sum()
        self.assertAlmostEqual(total, np.pi * MINOR ** 2 * KAPPA, delta=2e-3 * total)

    def test_connectivity_is_consistent_and_symmetric(self):
        conn = self.grid.connectivity()
        self.assertEqual(len(conn["faces"]), self.grid.nzones)
        self.assertEqual(sum(1 for face in conn["faces"] if len(face) == 3), ZONES[0])
        self.assertEqual(len([k for k in conn["nodes"] if k[0] == 0]), 1)   # single axis node
        for iface, per_edge in enumerate(conn["neighbours"]):
            self.assertEqual(len(per_edge), len(conn["faces"][iface]))
            for adjacent in per_edge:
                for other in adjacent:
                    self.assertTrue(any(iface in lst for lst in conn["neighbours"][other]))
        # Hanging nodes: a zone of the 4-zone row borders two zones of the 8-zone row.
        self.assertEqual(sorted(len(lst) for lst in conn["neighbours"][ZONES[0]]), [1, 1, 1, 2])

    def test_zone_centres_follow_nubeam_kernel_convention(self):
        x, th = self.grid.zone_centres()
        self.assertAlmostEqual(x[0], 0.5 / ROWS)
        self.assertAlmostEqual(th[0], -np.pi + 0.5 * 2 * np.pi / ZONES[0])


class ConventionTests(unittest.TestCase):
    def test_orientation_flags_and_sign(self):
        text = "NLJCCW=.TRUE.   ! comment\n NLBCCW = F\n"
        self.assertEqual(nf.read_orientation_flags(text), (True, False))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "t.DAT"
            path.write_text(text)
            self.assertEqual(nf.read_run_nsjdotb(path), -1)
            path.write_text("NLJCCW=T\nNLBCCW=T\n")
            self.assertEqual(nf.read_run_nsjdotb(path), 1)
            path.write_text("NLJCCW=T\n")
            self.assertIsNone(nf.read_run_nsjdotb(path))

    def test_pitch_is_flipped_only_when_j_and_b_antiparallel(self):
        edges = nf.pitch_edges_from_centres(np.array([-0.5, 0.0, 0.5]))
        same, perm = idm.pitch_to_field_direction(edges, +1)
        np.testing.assert_allclose(same, edges)
        flipped, perm = idm.pitch_to_field_direction(edges, -1)
        np.testing.assert_allclose(flipped, -edges[::-1])
        self.assertTrue(np.all(np.diff(flipped) > 0))
        np.testing.assert_array_equal(perm, [2, 1, 0])

    def test_phase_space_cells_units_ordering_and_closure(self):
        species = nf.FastIonSpecies("D", 2.0, 1.0, False, (np.arange(NE) + 0.5) * 20000.0,
                                    np.arange(NE + 1) * 20000.0, synthetic_distribution(), None)
        centres = np.linspace(-1, 1, NP + 1)
        centres = 0.5 * (centres[1:] + centres[:-1])
        snap = nf.FastIonSnapshot(0.1, "R", Path("d"), Path("a"), [species], nf.McGrid(ROWS, ZONES, -np.pi, np.pi),
                                  np.ones(int(ZONES.sum())), centres, nf.pitch_edges_from_centres(centres),
                                  None, -1, "test")
        _edges, perm = idm.pitch_to_field_direction(snap.pitch_edges, -1)
        fd3v, f_si = idm.phase_space_cells(species.gc, species, snap, perm)
        nz = int(ZONES.sum())
        cells = fd3v.reshape((nz, NE, NP), order="F")
        # zone index fastest; pitch axis reversed relative to NUBEAM.
        expected = species.gc[2, 5, 7] * 1e6 * 20000.0 * (2.0 / NP) / 2.0
        self.assertAlmostEqual(cells[7, 2, 0], expected, delta=1e-9 * expected)
        # zone density equals the NUBEAM integral n = sum f dE dxi / 2 (converted to m^-3)
        density = species.gc.sum(axis=(0, 1)) * 1e6 * 20000.0 * (2.0 / NP) / 2.0
        np.testing.assert_allclose(cells.sum(axis=(1, 2)), density, rtol=1e-12)
        # f * d3v == fd3v, with d3v = 2 pi (v_hi^3 - v_lo^3)/3 dxi
        v = np.sqrt(2 * species.energy_edges_ev * idm.EV / (2.0 * idm.AMU))
        d3v = 2 * np.pi * (v[3] ** 3 - v[2] ** 3) / 3.0 * (2.0 / NP)
        self.assertAlmostEqual(f_si.reshape((nz, NE, NP), order="F")[4, 2, 3] * d3v,
                               cells[4, 2, 3], delta=1e-9 * cells[4, 2, 3])


class RegridTests(unittest.TestCase):
    def test_regrid_conserves_weighted_integral_and_pairs_zones(self):
        rng = np.random.default_rng(3)
        values, weights = rng.random(20), rng.random(20) + 0.5
        edges = np.arange(21) / 20
        rows = idm.regrid_to_rows(values, edges, weights, 10)
        np.testing.assert_allclose(rows, (values[0::2] * weights[0::2] + values[1::2] * weights[1::2])
                                   / (weights[0::2] + weights[1::2]))
        # Non-commensurate target: the weighted integral is still conserved.
        rows7 = idm.regrid_to_rows(values, edges, weights, 7)
        frac = np.array([[max(0.0, min((k + 1) / 7, edges[i + 1]) - max(k / 7, edges[i])) * 20 for i in range(20)]
                         for k in range(7)])
        np.testing.assert_allclose(np.sum(rows7 * (frac * weights).sum(axis=1)), np.sum(values * weights))


class OutputMatchingTests(unittest.TestCase):
    def test_independent_time_lists_are_paired_by_time(self):
        with tempfile.TemporaryDirectory() as tmp:
            make_case(tmp, outtim=(0.1, 0.2, 0.3), fi_outtim=(0.3, 0.2), with_xplasma=False)
            data, fi = discover_nubeam_outputs(tmp, "RUNX")
            self.assertEqual((len(data), len(fi)), (3, 2))
            pairs, lost_data, lost_fi = nf.match_outputs(data, fi)
            self.assertEqual([(d.index, f.index) for d, f in pairs], [(2, 2), (3, 1)])
            self.assertEqual([d.index for d in lost_data], [1])
            self.assertEqual(lost_fi, [])


@unittest.skipIf(netCDF4 is None, "netCDF4 is required to build synthetic FI archives")
class ArchiveReaderTests(unittest.TestCase):
    def test_particle_state_members_are_never_extracted(self):
        with tempfile.TemporaryDirectory() as tmp:
            make_case(tmp, outtim=(0.1,), fi_outtim=(0.1,))
            opened = []
            original = tarfile.TarFile.extractfile

            def spy(self, member):
                opened.append(member.name if hasattr(member, "name") else str(member))
                return original(self, member)
            with mock.patch.object(tarfile.TarFile, "extractfile", spy):
                members = nf.read_archive_members(Path(tmp) / "RUNX_FI_TAR.GZ1")
            self.assertEqual(set(members), {"ps_ts1_state", "xplasma"})
            self.assertFalse([name for name in opened if "ptcl_state" in name])

    def test_snapshots_from_synthetic_case(self):
        with tempfile.TemporaryDirectory() as tmp:
            run, dists, volumes = make_case(tmp, outtim=(0.1, 0.2), fi_outtim=(0.2, 0.1), bccw=False, jccw=True)
            snaps = nf.read_fast_ion_snapshots(tmp, run)
            self.assertEqual([round(s.time, 6) for s in snaps], [0.1, 0.2])
            for snap in snaps:
                self.assertEqual(snap.nsjdotb, -1)
                self.assertEqual(snap.grid.nrows, ROWS)
                np.testing.assert_array_equal(snap.grid.zones_per_row, ZONES)
                np.testing.assert_allclose(snap.zone_volume_m3, volumes, rtol=1e-6)
                np.testing.assert_allclose(snap.species[0].gc, dists[round(snap.time, 6)], rtol=3e-8)
                np.testing.assert_allclose(snap.species[0].particle, 0.9 * dists[round(snap.time, 6)], rtol=3e-8)
                np.testing.assert_allclose(snap.grid.zone_volumes(snap.flux_map), snap.zone_volume_m3, rtol=3e-4)

    def test_lost_orbit_records_are_read_in_si_units(self):
        with tempfile.TemporaryDirectory() as tmp:
            run, _d, _v = make_case(tmp, outtim=(0.1,), fi_outtim=(0.1,))
            lost = nf.read_fast_ion_snapshots(tmp, run)[0].lost
            self.assertEqual(lost.count, LOST_N)                       # padding slot dropped
            np.testing.assert_allclose(lost.r, 1.40 + 0.01 * np.arange(1, LOST_N + 1))
            np.testing.assert_allclose(lost.z, -0.5 + 0.1 * np.arange(1, LOST_N + 1))
            np.testing.assert_array_equal(lost.prompt, [False, True, False, True, False])
            self.assertAlmostEqual(lost.loss_power(), np.sum(1e12 * np.arange(1, 6) * 1e4 * np.arange(1, 6)) * 1.602176634e-19 / 1e-3,
                                   delta=1e-9 * lost.loss_power())

    def test_missing_orientation_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            run, _d, _v = make_case(tmp, outtim=(0.1,), fi_outtim=(0.1,), with_xplasma=False)
            trdat = Path(tmp) / f"{run}TR.DAT"
            trdat.write_text("outtim=0.1\nfi_outtim=0.1\n")
            with self.assertRaises(nf.NubeamInputError):
                nf.read_fast_ion_snapshots(tmp, run)


NX = 20


def make_cdf(path, multi=False):
    """Synthetic TRANSP CDF; every quantity is linear in time so interpolation is exact."""
    import logging  # noqa: F401  (kept local: only the CDF fixture needs it)
    time = np.array([0.0, 0.2, 0.4])
    x = (np.arange(NX) + 0.5) / NX
    xb = (np.arange(NX) + 1.0) / NX
    k = np.arange(NX)
    series = {"BPLIM_D": 1.0e5, "BPCXI_D": 2.0e4, "BPCXX_D": 3.0e5, "BPTE_D": 2.0e6, "BPTI_D": 1.0e6,
              "BPTH_D": 6.0e4, "BSTH_D": 3.0e20, "SBDEPMC_D": 4.0e20, "PBDEPMC_D": 4.0e6, "BPHDEP0_D": 2.0,
              "BPHDEPGC_D": 1.0e13, "RAXIS": 100.0, "BZXR": 40.0}
    profiles = {"PBE_D": 2.0 + 0 * k, "PBI_D": 1.0 + 0 * k, "PBTH_D": 0.1 + 0 * k, "TQBE": 1.0e-6 * (1 + k),
                "TQBI": 2.0e-6 * (1 + k), "TQBTH": 3.0e-7 * (1 + k), "TQJXB_D": 1.0e-8 * (1 + k),
                "BDEP_D": 5.0e12 + 0 * k, "BDENS_D": 1.0e12 + 0 * k, "BTRAP_D": 0.2 + 0 * k,
                "CURB": 10.0 + k, "UCURB": 20.0 + k}
    if multi:
        # Second beam species (T), an unsuffixed total that must NOT be used when two beam
        # species exist, and fusion-product (F family) variables for T and He4.
        series.update({"BPLIM_T": 4.0e4, "BPTE_T": 3.0e6, "BPTE": 7.0e7, "BSTH_T": 5.0e20,
                       "FSTH_T": 1.0e19, "FSTH_4": 2.0e19, "FPHXB_T": 0.07, "FPHXB_4": 0.05})
        profiles.update({"BTRAP_T": 0.3 + 0 * k, "BDENS_T": 2.0e12 + 0 * k, "FTRAP_T": 0.35 + 0 * k,
                         "FTRAP_4": 0.45 + 0 * k, "FDENS_T": 1.0e10 + 0 * k, "FDENS_4": 3.0e10 + 0 * k})
    with netCDF4.Dataset(path, "w", format="NETCDF3_CLASSIC") as nc:
        nc.createDimension("TIME", time.size)
        nc.createDimension("TIME3", time.size)
        nc.createDimension("X", NX)
        nc.createDimension("XB", NX)
        nc.createVariable("TIME", "f8", ("TIME",))[:] = time
        nc.createVariable("TIME3", "f8", ("TIME3",))[:] = time
        nc.createVariable("X", "f8", ("TIME3", "X"))[:] = np.tile(x, (time.size, 1))
        nc.createVariable("XB", "f8", ("TIME3", "XB"))[:] = np.tile(xb, (time.size, 1))
        nc.createVariable("DVOL", "f8", ("TIME3", "X"))[:] = np.tile(1000.0 * (1 + k), (time.size, 1))
        nc.createVariable("DAREA", "f8", ("TIME3", "X"))[:] = np.tile(500.0 * (1 + k), (time.size, 1))
        for name, base in series.items():
            nc.createVariable(name, "f8", ("TIME",))[:] = base * (1.0 + time)
        for name, base in profiles.items():
            nc.createVariable(name, "f8", ("TIME3", "X"))[:] = np.outer(1.0 + time, base)
    return series, profiles


class ArchivePrefetchTests(unittest.TestCase):
    """Several FI_TAR archives are read at once but consumed strictly in order."""

    def test_results_come_back_in_input_order_whatever_the_completion_order(self):
        import threading
        import time

        def reader(name):
            time.sleep({"a": 0.15, "b": 0.0, "c": 0.05, "d": 0.0}[name])     # a finishes last
            return name.upper()
        self.assertEqual(list(nf.read_archives_in_order(["a", "b", "c", "d"], workers=4, reader=reader)),
                         ["A", "B", "C", "D"])

    def test_an_error_is_raised_at_its_turn_after_the_earlier_results(self):
        def reader(name):
            if name == "bad":
                raise nf.NubeamInputError("no equilibrium member")
            return name
        got = []
        with self.assertRaisesRegex(nf.NubeamInputError, "no equilibrium"):
            for item in nf.read_archives_in_order(["x", "bad", "y"], workers=3, reader=reader):
                got.append(item)
        self.assertEqual(got, ["x"])

    def test_at_most_twice_the_workers_are_read_ahead(self):
        import threading
        started, lock = [], threading.Lock()
        gate = threading.Event()

        def reader(index):
            with lock:
                started.append(index)
            gate.wait(timeout=5)
            return index
        generator = nf.read_archives_in_order(list(range(20)), workers=2, reader=reader)
        gate.set()
        self.assertEqual(next(generator), 0)
        with lock:
            ahead = len(started)
        self.assertLessEqual(ahead, 2 * 2 + 1)            # window of 4 plus the one being consumed
        self.assertEqual(list(generator), list(range(1, 20)))

    def test_single_worker_and_single_archive_are_sequential(self):
        order = []

        def reader(index):
            order.append(index)
            return index
        self.assertEqual(list(nf.read_archives_in_order([1, 2, 3], workers=1, reader=reader)), [1, 2, 3])
        self.assertEqual(order, [1, 2, 3])
        self.assertEqual(list(nf.read_archives_in_order([9], workers=4, reader=reader)), [9])

    def test_worker_count_follows_the_environment_and_the_archive_count(self):
        with mock.patch.dict(os.environ, {"TRANSP2IMAS_FI_WORKERS": "1"}):
            self.assertEqual(nf.archive_workers(10), 1)
        with mock.patch.dict(os.environ, {"TRANSP2IMAS_FI_WORKERS": "16"}):
            self.assertEqual(nf.archive_workers(3), 3)
            self.assertEqual(nf.archive_workers(100), 16)
        env = {k: v for k, v in os.environ.items() if k != "TRANSP2IMAS_FI_WORKERS"}
        with mock.patch.dict(os.environ, env, clear=True):
            self.assertTrue(1 <= nf.archive_workers(100) <= 4)
            self.assertEqual(nf.archive_workers(0), 1)

    @unittest.skipIf(netCDF4 is None, 'synthetic archives need netCDF4')
    def test_parallel_and_sequential_snapshots_agree(self):
        with tempfile.TemporaryDirectory() as tmp:
            run, _dists, _volumes = make_case(tmp, outtim=(0.1, 0.2, 0.3), fi_outtim=(0.3, 0.1, 0.2))
            with mock.patch.dict(os.environ, {"TRANSP2IMAS_FI_WORKERS": "1"}):
                sequential = nf.read_fast_ion_snapshots(tmp, run)
            with mock.patch.dict(os.environ, {"TRANSP2IMAS_FI_WORKERS": "3"}):
                parallel = nf.read_fast_ion_snapshots(tmp, run)
        self.assertEqual([round(s.time, 6) for s in sequential], [round(s.time, 6) for s in parallel])
        for a, b in zip(sequential, parallel):
            np.testing.assert_array_equal(a.flux_map.r_geo, b.flux_map.r_geo)
            self.assertEqual((a.nsjdotb, a.b_axis_vac), (b.nsjdotb, b.b_axis_vac))


class CdfShim:
    """Just enough of ``nc2imas`` for the distributions writer, backed by a real CDF."""

    cm2m = 0.01
    bt_sign = -1.0

    def __init__(self, cdf_path, dd_version):
        import logging
        self.ncfile = str(cdf_path)
        self.logger = logging.getLogger("verify_nubeam_fi")
        self.dd_version = dd_version
        self._helper = idm.StandaloneHelper(dd_version)
        self._nc = netCDF4.Dataset(cdf_path)
        self.timesc = np.array(self._nc["TIME"][:])
        self.mapping = MappingLedger()
        self.written = None

    def get_ncvar(self, name, mult=None):
        if name not in self._nc.variables:
            return None
        values = np.array(self._nc[name][:], dtype=float)
        return values * mult if mult is not None else values

    @staticmethod
    def is_valid(*values, name=None):
        return all(v is not None for v in values)

    def vacuum_field_reference(self, times):
        raxis, bzxr = self.get_ncvar("RAXIS", self.cm2m), self.get_ncvar("BZXR", self.cm2m)
        if raxis is None or bzxr is None:
            return None
        return vacuum_field_reference(raxis, bzxr, self.bt_sign, times, self.timesc)

    def new_ids(self, name):
        return imas.IDSFactory(self.dd_version).new(name)

    def fill_properties(self, ids):
        ids.ids_properties.comment = "synthetic"
        return ids

    def fill_code_properties(self, ids):
        return ids

    def putdb_imas(self, ids):
        self.written = ids

    def dd_at_least(self, *version):
        return self._helper.dd_at_least(*version)

    def set_species_name(self, node, value):
        self._helper.set_species_name(node, value)

    def set_nuclear_charge(self, node, value):
        self._helper.set_nuclear_charge(node, value)


class SpeciesNamingTests(unittest.TestCase):
    def test_transp_cdf_suffixes(self):
        beam = {(1, 1.0): "H", (1, 1.9995): "D", (1, 3.0): "T", (2, 3.0): "3", (2, 4.0): "4",
                (10, 20.0): "N", (18, 40.0): "A", (36, 84.0): "K", (54, 131.0): "X", (6, 12.0): "Z"}
        for (z, a), suffix in beam.items():
            self.assertEqual(nf.cdf_species_suffix(a, z, False), suffix, (z, a))
        fusion = {(1, 1.0): "P", (1, 3.0): "T", (2, 3.0): "3", (2, 4.0): "4"}
        for (z, a), suffix in fusion.items():
            self.assertEqual(nf.cdf_species_suffix(a, z, True), suffix)
        self.assertIsNone(nf.cdf_species_suffix(2.0, 1.0, True))      # no fusion deuteron family

    def test_candidate_names_beam_versus_fusion(self):
        mk = lambda a, z, fp: nf.FastIonSpecies("x", a, z, fp, np.zeros(1), np.zeros(2), np.zeros((1, 1, 1)), None)
        beam_d, fus_he4 = mk(2.0, 1.0, False), mk(4.0, 2.0, True)
        self.assertEqual(beam_d.cdf_names("BPTE"), ["BPTE_D"])
        self.assertEqual(beam_d.cdf_names("BPTE", allow_plain=True), ["BPTE_D", "BPTE"])
        self.assertEqual(fus_he4.cdf_names("BTRAP", allow_plain=True), ["FTRAP_4"])    # never plain
        self.assertEqual(fus_he4.cdf_names("BSTH"), ["FSTH_4"])
        self.assertEqual(fus_he4.cdf_names("BPLIM"), [])                              # no power counterpart


class RealCdfNamingTests(unittest.TestCase):
    """The suffix rules reproduce the names found in a real CDF with D beams and fusion products."""

    ITER_CDF = ROOT / "examples" / "data" / "iter" / "38530P94.CDF"

    def test_suffix_rules_resolve_to_existing_variables(self):
        if not self.ITER_CDF.is_file() or netCDF4 is None:
            self.skipTest("ITER example CDF not present")
        with netCDF4.Dataset(self.ITER_CDF) as nc:
            names = set(nc.variables)
        mk = lambda a, z, fp: nf.FastIonSpecies("x", a, z, fp, np.zeros(1), np.zeros(2), np.zeros((1, 1, 1)), None)
        beam_d = mk(2.0, 1.0, False)
        for stem in ("BPLIM", "BPTE", "BPTH", "BSTH", "TQJXB", "PBE", "BDEP", "BDENS", "BTRAP", "SBDEPMC"):
            self.assertIn(beam_d.cdf_names(stem)[0], names, stem)
        for a, z in ((3.0, 1.0), (3.0, 2.0), (4.0, 2.0)):             # fusion T, He3, He4
            fusion = mk(a, z, True)
            for stem in ("BDENS", "BTRAP", "BSTH", "BPHXB"):
                self.assertIn(fusion.cdf_names(stem)[0], names, (a, z, stem))
            self.assertEqual(fusion.cdf_names("BPLIM"), [])
            self.assertNotIn(f"BPLIM_{fusion.cdf_suffix}", names)    # TRANSP indeed has no such variable


SPECIES_MIX = [(2.0, 1.0, False), (3.0, 1.0, False), (3.0, 1.0, True), (4.0, 2.0, True)]   # D, T, fusion T, fusion He4


@unittest.skipIf(imas is None or netCDF4 is None, "IMAS-Python and netCDF4 are required")
class MultiSpeciesTests(unittest.TestCase):
    def _case(self, species, dd_version="4.1.1", multi=True):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        run, _d, _v = make_case(tmp.name, outtim=(0.1, 0.2), fi_outtim=(0.1, 0.2), species=species)
        make_cdf(Path(tmp.name) / f"{run}.CDF", multi=multi)
        shim = CdfShim(Path(tmp.name) / f"{run}.CDF", dd_version)
        writer = idm.imas_distributions(shim)
        self.assertIsNone(shim.written, 'constructing the writer must not write')
        writer.run()
        self.assertIsNotNone(shim.written)
        shim.written.validate()
        return shim.written, nf.read_fast_ion_snapshots(tmp.name, run)

    def test_species_are_keyed_by_index_not_label(self):
        ids, snaps = self._case(SPECIES_MIX)
        self.assertEqual(len(ids.distribution), 8)              # 4 species x (GC, particle)
        names = [str(d.species.ion.name) for d in ids.distribution]
        self.assertEqual(names, ["D", "D", "T", "T", "T", "T", "He4", "He4"])   # beam T and fusion T share a label
        proc = [(str(d.process[0].type.name), int(d.process[0].type.index)) for d in ids.distribution]
        self.assertEqual(proc[0], ("NBI", 1))
        self.assertEqual(proc[2], ("NBI", 1))
        self.assertEqual(proc[4], ("D_D_to_T_H", 105))          # fusion triton
        self.assertEqual(proc[6], ("nuclear", 100))             # fusion He4: reaction not identifiable
        # every entry carries ITS species' distribution (distinct synthetic data per species)
        nz = int(ZONES.sum())
        for n, (dist, m) in enumerate(zip(ids.distribution, [(k, g) for k in range(4) for g in (2, 1)])):
            values = np.asarray(dist.ggd[0].expansion_fd3v[0].grid_subset[0].values).reshape((nz, NE, NP), order="F")
            snap = snaps[0]
            expected = synthetic_distribution(seed=10 * 1 + m[0]) * (1.0 if m[1] == 2 else 0.9)
            ratio = values.sum() / expected.sum()
            self.assertAlmostEqual(ratio / (1e6 * 20000.0 * (2.0 / NP) / 2.0), 1.0, delta=1e-6, msg=f"entry {n}")
        # lost-orbit markers (species 1) appear only on the D particle-position entry
        counts = [len(d.markers) for d in ids.distribution]
        self.assertEqual(counts, [0, 2, 0, 0, 0, 0, 0, 0])

    def test_two_beam_species_use_their_own_suffixed_variables(self):
        ids, _snaps = self._case(SPECIES_MIX)
        for it, t in enumerate((0.1, 0.2)):
            lin = 1.0 + t
            d, tb = ids.distribution[0].global_quantities[it], ids.distribution[2].global_quantities[it]
            self.assertAlmostEqual(float(d.collisions.electrons.power_thermal), 2.0e6 * lin, delta=1e-4)
            self.assertAlmostEqual(float(tb.collisions.electrons.power_thermal), 3.0e6 * lin, delta=1e-4)   # BPTE_T, not 7e7
            self.assertAlmostEqual(float(tb.power_first_orbit), 4.0e4 * lin, delta=1e-6)
            self.assertAlmostEqual(float(tb.thermalization.particles), 5.0e20 * lin, delta=1e8)
            # unsuffixed all-beam totals are ambiguous with two beam species and must be skipped
            for dist in (ids.distribution[0], ids.distribution[2]):
                gq, pr = dist.global_quantities[it], dist.profiles_1d[it]
                self.assertFalse(gq.current_phi.has_value)
                self.assertFalse(pr.current_fast_phi.has_value)
                self.assertFalse(gq.collisions.electrons.torque_thermal_phi.has_value)
            # trapped fraction of the T beam uses BTRAP_T
            pr = ids.distribution[2].profiles_1d[it]
            np.testing.assert_allclose(np.asarray(pr.trapped.density_fast),
                                       0.3 * lin * np.asarray(pr.density_fast), rtol=1e-12)

    def test_fusion_products_only_get_their_f_family_quantities(self):
        ids, _snaps = self._case(SPECIES_MIX)
        for it, t in enumerate((0.1, 0.2)):
            lin = 1.0 + t
            for entry, thermal, trapped, torque in ((4, 1.0e19, 0.35, 0.07), (6, 2.0e19, 0.45, 0.05)):
                dist = ids.distribution[entry]
                gq, pr = dist.global_quantities[it], dist.profiles_1d[it]
                self.assertAlmostEqual(float(gq.thermalization.particles), thermal * lin, delta=1e8)
                self.assertAlmostEqual(float(gq.torque_tor_j_radial), torque * lin, delta=1e-12)   # scalar FPHXB_*
                np.testing.assert_allclose(np.asarray(pr.trapped.density_fast),
                                           trapped * lin * np.asarray(pr.density_fast), rtol=1e-12)
                for absent in (gq.power_first_orbit, gq.power_charge_exchange, gq.collisions.electrons.power_thermal,
                               gq.thermalization.power, gq.current_phi):
                    self.assertFalse(absent.has_value)
                self.assertEqual(len(gq.source), 0)
                self.assertEqual(len(gq.collisions.ion), 0)

    def test_sole_beam_with_fusion_products_attaches_beam_totals_to_the_beam_only(self):
        ids, _snaps = self._case([(2.0, 1.0, False), (4.0, 2.0, True)], multi=False)
        beam_gc, fusion_gc = ids.distribution[0], ids.distribution[2]
        self.assertTrue(beam_gc.global_quantities[0].current_phi.has_value)             # CURB (beam total)
        self.assertTrue(beam_gc.profiles_1d[0].current_fast_phi.has_value)
        self.assertTrue(beam_gc.global_quantities[0].collisions.electrons.torque_thermal_phi.has_value)   # TQBE
        self.assertFalse(fusion_gc.global_quantities[0].current_phi.has_value)
        self.assertFalse(fusion_gc.global_quantities[0].collisions.electrons.torque_thermal_phi.has_value)


@unittest.skipIf(imas is None or netCDF4 is None, "IMAS-Python and netCDF4 are required")
class WriterCdfQuantityTests(unittest.TestCase):
    """The CDF-derived distribution quantities, through the real writer and a synthetic CDF."""

    def _write(self, dd_version):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        run, _dists, _vol = make_case(tmp.name, outtim=(0.1, 0.2), fi_outtim=(0.1, 0.2))
        series, profiles = make_cdf(Path(tmp.name) / f"{run}.CDF")
        shim = CdfShim(Path(tmp.name) / f"{run}.CDF", dd_version)
        idm.imas_distributions(shim).run()
        self.assertIsNotNone(shim.written, "writer produced no IDS")
        shim.written.validate()
        return shim.written, series, profiles

    @staticmethod
    def _first(node, *names):
        for name in names:
            if hasattr(node, name):
                return getattr(node, name)
        raise AssertionError(f"none of {names} in {node}")

    def test_vacuum_field_satisfies_r0_times_b0_equals_bt_sign_bzxr(self):
        ids, _series, _profiles = self._write("4.1.1")
        times = np.asarray(ids.time)
        expected = -np.interp(times, [0.0, 0.2, 0.4], 0.01 * 40.0 * (1.0 + np.array([0.0, 0.2, 0.4])))   # bt_sign * BZXR
        np.testing.assert_allclose(float(ids.vacuum_toroidal_field.r0) * np.asarray(ids.vacuum_toroidal_field.b0),
                                   expected, rtol=1e-12)

    def test_balances_torques_currents_and_trapped_in_each_dd(self):
        dvol = 1000.0 * (1 + np.arange(NX))
        darea = 500.0 * (1 + np.arange(NX))
        for dd in ("3.39.0", "4.1.1"):
            ids, series, profiles = self._write(dd)
            for it, t in enumerate((0.1, 0.2)):
                lin = 1.0 + t
                for dist in ids.distribution:
                    gq, pr = dist.global_quantities[it], dist.profiles_1d[it]
                    tag = f"{dd} t={t} gyro={int(dist.gyro_type)}"
                    if dd != "3.39.0":
                        self.assertAlmostEqual(float(gq.power_first_orbit), 1.0e5 * lin, delta=1e-6, msg=tag)
                        self.assertAlmostEqual(float(gq.power_charge_exchange), 3.2e5 * lin, delta=1e-6, msg=tag)
                    el = gq.collisions.electrons
                    self.assertAlmostEqual(float(el.power_thermal), 2.0e6 * lin, delta=1e-4, msg=tag)
                    tq_e = np.sum(1.0e-6 * (1 + np.arange(NX)) * lin * dvol)
                    self.assertAlmostEqual(float(self._first(el, "torque_thermal_phi", "torque_thermal_tor")), tq_e,
                                           delta=1e-12, msg=tag)
                    self.assertEqual(len(gq.collisions.ion), 1)
                    self.assertAlmostEqual(float(gq.collisions.ion[0].power_thermal), 1.0e6 * lin, delta=1e-4)
                    th = self._first(gq, "thermalization", "thermalisation")
                    self.assertAlmostEqual(float(th.particles), 3.0e20 * lin, delta=1e8)
                    self.assertAlmostEqual(float(th.power), 6.0e4 * lin, delta=1e-6)
                    self.assertAlmostEqual(float(th.torque), np.sum(3.0e-7 * (1 + np.arange(NX)) * lin * dvol), delta=1e-12)
                    src = gq.source[0]
                    self.assertEqual(int(src.identifier.type.index), 2)
                    self.assertAlmostEqual(float(src.particles), 4.0e20 * lin, delta=1e8)
                    self.assertAlmostEqual(float(src.power), 4.0e6 * lin, delta=1e-6)
                    torque = self._first(src, "torque_phi", "torque_tor")
                    if int(dist.gyro_type) == 1:             # ion position only; GC value is not in N m
                        self.assertAlmostEqual(float(torque), 2.0 * lin, delta=1e-9, msg=tag)
                    else:
                        self.assertNotAlmostEqual(float(torque), 1.0e13 * lin, delta=1.0)
                    # J x B torque and driven current (global integrals, A / N m)
                    self.assertAlmostEqual(float(self._first(gq, "torque_tor_j_radial")),
                                           np.sum(1.0e-8 * (1 + np.arange(NX)) * lin * dvol), delta=1e-12)
                    self.assertAlmostEqual(float(self._first(gq, "current_phi", "current_tor")),
                                           np.sum((10.0 + np.arange(NX)) * lin * darea), delta=1e-6)
                    # profiles on the 10 MC rows: constant CDF profiles stay constant, units to SI
                    pe = self._first(pr.collisions, "electrons").power_thermal if hasattr(pr.collisions, "electrons") \
                        else pr.collisions.power_thermal
                    np.testing.assert_allclose(np.asarray(pe), 2.0e6 * lin, rtol=1e-12)
                    np.testing.assert_allclose(np.asarray(pr.trapped.density_fast),
                                               0.2 * lin * np.asarray(pr.density_fast), rtol=1e-12)
                    np.testing.assert_allclose(np.asarray(pr.source[0].particles), 5.0e12 * lin * 1.0e6, rtol=1e-12)
                    # independent brute-force overlap of the 20 TRANSP zones with the ROWS MC rows
                    share = np.array([[max(0.0, min((r + 1) / ROWS, (z + 1) / NX) - max(r / ROWS, z / NX)) * NX
                                       for z in range(NX)] for r in range(ROWS)])
                    fast = np.asarray(self._first(pr, "current_fast_phi", "current_fast_tor"))
                    expected = (share * darea * (20.0 + np.arange(NX)) * lin).sum(axis=1) \
                        / (share * darea).sum(axis=1) * 1.0e4
                    np.testing.assert_allclose(fast, expected, rtol=1e-12)
                    # the profile integral over the MC rows equals the CDF integral (conservation)
                    tq_prof = np.asarray(self._first(pr, "torque_phi_j_radial", "torque_tor_j_radial"))
                    row_vol = (share * dvol).sum(axis=1)
                    self.assertAlmostEqual(float(np.sum(tq_prof * row_vol * 1.0e-6)),
                                           np.sum(1.0e-8 * (1 + np.arange(NX)) * lin * dvol), delta=1e-12)


@unittest.skipIf(imas is None or netCDF4 is None, "IMAS-Python and netCDF4 are required")
class IdsTests(unittest.TestCase):
    DD_VERSIONS = ("3.39.0", "4.1.1")

    def _fill(self, dd_version, **case):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        run, dists, volumes = make_case(tmp.name, **case)
        snaps = nf.read_fast_ion_snapshots(tmp.name, run)
        ids = imas.IDSFactory(dd_version).new("distributions")
        idm.DistributionsBuilder(idm.StandaloneHelper(dd_version)).fill(ids, snaps)
        return ids, snaps, dists

    def test_structure_values_and_closure_in_each_dd(self):
        for dd in self.DD_VERSIONS:
            with self.subTest(dd=dd):
                ids, snaps, _ = self._fill(dd, outtim=(0.1, 0.2), fi_outtim=(0.1, 0.2))
                ids.validate()
                self.assertEqual(len(ids.distribution), 2)           # GC and particle
                self.assertEqual([int(d.gyro_type) for d in ids.distribution], [2, 1])
                np.testing.assert_allclose(ids.time, [0.1, 0.2])
                nz = int(ZONES.sum())
                for dist in ids.distribution:
                    for it, snap in enumerate(snaps):
                        ggd = dist.ggd[it]
                        grid = ggd.grid
                        self.assertEqual(len(grid.space), 3)
                        faces = grid.space[0].objects_per_dimension[2].object
                        self.assertEqual(len(faces), nz)
                        self.assertEqual({int(s.identifier.index) for s in grid.grid_subset}, {5, -1, -2, 200})
                        values = np.asarray(ggd.expansion_fd3v[0].grid_subset[0].values)
                        self.assertEqual(values.size, nz * NE * NP)
                        vol = np.array([o.geometry[2] for o in faces])
                        total = float(np.sum(values.reshape((nz, NE, NP), order="F").sum(axis=(1, 2)) * vol))
                        self.assertAlmostEqual(total, float(dist.global_quantities[it].particles_fast_n),
                                               delta=1e-9 * total)
                        # profile closure: volume-weighted row densities integrate to N
                        rows = np.repeat(np.arange(ROWS), ZONES)
                        row_vol = np.bincount(rows, weights=vol)
                        self.assertAlmostEqual(float(np.sum(dist.profiles_1d[it].density * row_vol)), total,
                                               delta=1e-9 * total)

    def test_lost_markers_use_field_direction_pitch(self):
        for dd in self.DD_VERSIONS:
            for jccw, bccw, sign in ((True, True, 1), (True, False, -1)):
                with self.subTest(dd=dd, sign=sign):
                    ids, snaps, _ = self._fill(dd, outtim=(0.1,), fi_outtim=(0.1,), jccw=jccw, bccw=bccw)
                    ids.validate()
                    self.assertEqual(len(ids.distribution[0].markers), 0)          # GC entry has none
                    marker = ids.distribution[1].markers[0]                       # particle-position entry
                    self.assertEqual([str(c.name) for c in marker.coordinate_identifier],
                                     ["r", "z", "energy_kinetic", "pitch", "time_lost", "nbi_unit", "prompt_loss"])
                    pos = np.asarray(marker.positions)
                    self.assertEqual(pos.shape, (LOST_N, 7))
                    np.testing.assert_allclose(pos[:, 3], sign * 0.1 * np.arange(1, LOST_N + 1))
                    np.testing.assert_allclose(np.asarray(marker.weights), 1e12 * np.arange(1, LOST_N + 1))
                    np.testing.assert_allclose(pos[:, 0], snaps[0].lost.r)

    def test_pitch_axis_follows_field_direction(self):
        for jccw, bccw, sign in ((True, True, 1), (True, False, -1), (False, False, 1), (False, True, -1)):
            with self.subTest(jccw=jccw, bccw=bccw):
                ids, snaps, dists = self._fill("4.1.1", outtim=(0.1,), fi_outtim=(0.1,), jccw=jccw, bccw=bccw)
                self.assertEqual(snaps[0].nsjdotb, sign)
                gc = ids.distribution[0]
                edges = np.array([o.geometry[0] for o in
                                  gc.ggd[0].grid.space[2].objects_per_dimension[0].object])
                self.assertTrue(np.all(np.diff(edges) > 0))
                nz = int(ZONES.sum())
                cells = np.asarray(gc.ggd[0].expansion_fd3v[0].grid_subset[0].values).reshape(
                    (nz, NE, NP), order="F")
                by_pitch = cells.sum(axis=(0, 1))
                # Bins are equally wide, so only the order may change: reversed when nsjdotb = -1.
                source = dists[0.1].sum(axis=(0, 2))
                scale = by_pitch.sum() / source.sum()
                expected = source if sign > 0 else source[::-1]
                np.testing.assert_allclose(by_pitch, expected * scale, rtol=1e-9)

    def test_plot_script_renders_all_figures(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest("matplotlib is not installed")
        matplotlib.use("Agg")
        sys.path.insert(0, str(ROOT / "examples"))
        import plot_fast_ion_distribution as plotter
        ids, _snaps, _ = self._fill("4.1.1", outtim=(0.1,), fi_outtim=(0.1,))
        spaces = list(plotter.phase_spaces_from_ids(ids))
        self.assertEqual(len(spaces), 2)
        self.assertAlmostEqual(spaces[0].total(), float(ids.distribution[0].global_quantities[0].particles_fast_n),
                               delta=1e-9 * spaces[0].total())
        import matplotlib.pyplot as plt
        with tempfile.TemporaryDirectory() as out:
            for name, fig in plotter.make_figures(spaces[0]).items():
                path = Path(out) / f"{name}.png"
                fig.savefig(path)
                plt.close(fig)
                self.assertGreater(path.stat().st_size, 1000)

    def test_ids_round_trip_through_hdf5(self):
        for dd in self.DD_VERSIONS:
            with self.subTest(dd=dd):
                ids, _snaps, _ = self._fill(dd, outtim=(0.1,), fi_outtim=(0.1,))
                with tempfile.TemporaryDirectory() as out:
                    uri = f"imas:hdf5?path={out}/db"
                    with imas.DBEntry(uri, "w", dd_version=dd) as entry:
                        entry.put(ids)
                    with imas.DBEntry(uri, "r", dd_version=dd) as entry:
                        back = entry.get("distributions")
                    a = np.asarray(ids.distribution[0].ggd[0].expansion_fd3v[0].grid_subset[0].values)
                    b = np.asarray(back.distribution[0].ggd[0].expansion_fd3v[0].grid_subset[0].values)
                    np.testing.assert_array_equal(a, b)
                    self.assertEqual(int(back.distribution[0].ggd[0].grid.identifier.index), idm.GRID_INDEX)


class LocalSampleTests(unittest.TestCase):
    """Optional checks against the untracked NSTX case; skipped in CI."""

    @classmethod
    def setUpClass(cls):
        sample = SAMPLE_ARGS["dir"]
        if netCDF4 is None or sample is None or not Path(sample).is_dir():
            raise unittest.SkipTest("local NSTX sample is not selected or not present")
        cls.snaps = nf.read_fast_ion_snapshots(sample, "120967E03")

    def test_pairs_follow_trdat_times(self):
        self.assertEqual([round(s.time, 3) for s in self.snaps], [0.51, 0.745])

    def test_reconstructed_zones_reproduce_nubeam_volumes(self):
        for snap in self.snaps:
            mine = snap.grid.zone_volumes(snap.flux_map)
            np.testing.assert_allclose(mine, snap.zone_volume_m3, rtol=5e-5)

    def test_distribution_closure_and_sign_conventions(self):
        for snap in self.snaps:
            self.assertEqual(snap.nsjdotb, -1)       # NSTX: Ip CCW, Bphi CW
            sp = snap.species[0]
            gc, flr = nf.species_moments(snap, sp, "gc"), nf.species_moments(snap, sp, "particle")
            self.assertGreater(gc["N"], 1.0e18)
            self.assertAlmostEqual(gc["N"] / flr["N"], 1.0, places=4)   # FLR conserves particles
            mean_pitch = np.sum(sp.gc.sum(axis=(0, 2)) * snap.pitch_centres) / sp.gc.sum()
            self.assertGreater(mean_pitch, 0.2)       # co-current NBI is positive along Ip


    def test_lost_orbit_power_matches_transp_orbit_loss(self):
        """Sum(weight * energy)/window of the lost records reproduces BPLIM_D."""
        cdf = Path(SAMPLE_ARGS["dir"]) / "120967E03.CDF"
        if not cdf.is_file():
            self.skipTest("TRANSP CDF not present")
        with netCDF4.Dataset(cdf) as nc:
            time = np.array(nc["TIME"][:])
            for snap in self.snaps:
                power = float(nc["BPLIM_D"][int(np.argmin(abs(time - snap.time)))])
                self.assertAlmostEqual(snap.lost.loss_power() / power, 1.0, delta=2e-3)
                self.assertGreater(snap.lost.prompt.mean(), 0.5)

    def test_cdf_quantities_match_transp_scalars_on_real_data(self):
        """New CDF-derived nodes reproduce TRANSP's own scalars (and its profile integrals)."""
        cdf = Path(SAMPLE_ARGS["dir"]) / "120967E03.CDF"
        if not cdf.is_file():
            self.skipTest("TRANSP CDF not present")
        shim = CdfShim(cdf, "4.1.1")
        ids = imas.IDSFactory("4.1.1").new("distributions")
        builder = idm.DistributionsBuilder(shim)
        builder.fill(ids, self.snaps)
        writer = idm.imas_distributions.__new__(idm.imas_distributions)
        writer.im, writer.logger = shim, shim.logger
        for method in (writer._fill_cdf_power_losses, writer._fill_cdf_torque_current, writer._fill_cdf_balances):
            method(ids, self.snaps, builder.meta)
        ids.validate()
        t = shim.timesc
        nc = shim._nc
        t3 = np.array(nc["TIME3"][:])
        for it, snap in enumerate(self.snaps):
            i, j = int(np.argmin(abs(t - snap.time))), int(np.argmin(abs(t3 - snap.time)))
            dvol = np.array(nc["DVOL"][j])
            for dist in ids.distribution:
                gq, pr = dist.global_quantities[it], dist.profiles_1d[it]
                for node, name in ((gq.collisions.electrons.power_thermal, "BPTE_D"),
                                   (gq.collisions.ion[0].power_thermal, "BPTI_D"),
                                   (gq.thermalization.power, "BPTH_D"), (gq.source[0].power, "PBDEPMC_D"),
                                   (gq.source[0].particles, "SBDEPMC_D"), (gq.thermalization.particles, "BSTH_D")):
                    self.assertAlmostEqual(float(node) / float(nc[name][i]), 1.0, delta=2e-3, msg=name)
                # CDF-internal consistency: volume integral of the profile = TRANSP's scalar
                for prof_name, scalar in (("PBE_D", "BPTE_D"), ("PBI_D", "BPTI_D"), ("PBTH_D", "BPTH_D")):
                    integral = float(np.sum(np.array(nc[prof_name][j]) * dvol))
                    self.assertAlmostEqual(integral / float(nc[scalar][i]), 1.0, delta=2e-3, msg=prof_name)
                # MC-row profile reproduces the electron heating when integrated over MC zone volumes
                rows = np.repeat(np.arange(snap.grid.nrows), snap.grid.zones_per_row)
                row_vol = np.bincount(rows, weights=snap.zone_volume_m3)
                heating = float(np.sum(np.asarray(dist.profiles_1d[it].collisions.electrons.power_thermal) * row_vol))
                self.assertAlmostEqual(heating / float(nc["BPTE_D"][i]), 1.0, delta=5e-3)
                # trapped fraction in (0, 1) and unshielded current exceeds the shielded one
                trapped = np.asarray(pr.trapped.density_fast) / np.asarray(pr.density_fast)
                self.assertTrue(np.all((trapped > 0.1) & (trapped < 0.8)), trapped)   # BTRAP spans 0.15-0.55
                self.assertGreater(float(np.sum(pr.current_fast_phi)), float(np.sum(pr.current_phi)))
            for dist in ids.distribution:
                if int(dist.gyro_type) == 1:
                    self.assertAlmostEqual(float(dist.global_quantities[it].source[0].torque_phi)
                                           / float(nc["BPHDEP0_D"][i]), 1.0, delta=2e-3)

    def test_transp_1d_density_is_the_particle_position_density(self):
        """BDENS_D follows the FLR distribution (1-3% per row), not the guiding-centre one."""
        cdf = Path(SAMPLE_ARGS["dir"]) / "120967E03.CDF"
        if not cdf.is_file():
            self.skipTest("TRANSP CDF not present")
        with netCDF4.Dataset(cdf) as nc:
            t3 = np.array(nc["TIME3"][:])
            for snap in self.snaps:
                j = int(np.argmin(abs(t3 - snap.time)))
                per_zone = np.array(nc["BDENS_D"][j]) * np.array(nc["DVOL"][j])          # particles per TRANSP zone
                bdens = per_zone.reshape(-1, 2).sum(axis=1)                              # 20 zones -> 10 rows
                rows = np.repeat(np.arange(snap.grid.nrows), snap.grid.zones_per_row)
                sp = snap.species[0]
                gc = np.bincount(rows, weights=nf.species_moments(snap, sp, "gc")["n"] * snap.zone_volume_m3)
                flr = np.bincount(rows, weights=nf.species_moments(snap, sp, "particle")["n"] * snap.zone_volume_m3)
                np.testing.assert_allclose(bdens, flr, rtol=0.05)
                self.assertGreater(bdens[-1] / gc[-1], 5.0)         # the GC density is far below it at the edge

    def test_moments_match_transp_cdf_gc_moments(self):
        """Zone density and <E_par>, <E_perp> equal TRANSP's BDENS2_D/BEPLL2_D/BEPRP2_D."""
        cdf = Path(SAMPLE_ARGS["dir"]) / "120967E03.CDF"
        if not cdf.is_file():
            self.skipTest("TRANSP CDF not present")
        with netCDF4.Dataset(cdf) as nc:
            t3 = np.array(nc["TIME3"][:])
            for snap in self.snaps:
                it = int(np.argmin(abs(t3 - snap.time)))
                nz = snap.grid.nzones
                dens = np.array(nc["BDENS2_D"][it, :nz]) * 1.0e6
                epar, eperp = np.array(nc["BEPLL2_D"][it, :nz]), np.array(nc["BEPRP2_D"][it, :nz])
                sp = snap.species[0]
                cell = sp.gc * np.diff(sp.energy_edges_ev)[:, None, None] * np.diff(snap.pitch_edges)[None, :, None] / 2
                tot = cell.sum(axis=(0, 1))
                mine_par = np.einsum("epz,e,p->z", cell, sp.energy_centres_ev, snap.pitch_centres ** 2) / tot
                mine_perp = np.einsum("epz,e,p->z", cell, sp.energy_centres_ev, 1 - snap.pitch_centres ** 2) / tot
                ok = dens > 0
                np.testing.assert_allclose(nf.species_moments(snap, sp, "gc")["n"][ok], dens[ok], rtol=1e-4)
                np.testing.assert_allclose(mine_par[ok], epar[ok], rtol=1e-4)
                np.testing.assert_allclose(mine_perp[ok], eperp[ok], rtol=1e-4)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="skip the (slow) local NSTX sample checks")
    parser.add_argument("--sample", type=Path, default=SAMPLE_DIR, help="optional local NSTX case directory")
    args, rest = parser.parse_known_args(argv)
    SAMPLE_ARGS["dir"] = None if args.quick else args.sample
    result = unittest.main(argv=[sys.argv[0]] + rest, exit=False, verbosity=1).result
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
