# NUBEAM fast-ion distribution functions

`convert.py --ids distributions` (also part of `--ids all`) writes the NUBEAM
fast-ion distribution to the IMAS `distributions` IDS. The converter uses only
the numbered TRANSP `DATA<n>` snapshots and `FI_TAR.GZ<n>` archives selected by
the TR.DAT `outtim` and `fi_outtim` lists, and needs no TRANSP/NUBEAM library.
NUBEAM CPU particle-state members of the archives are skipped while streaming
and are never extracted or parsed. Runs without `fi_outtim`, or without the
requested files, are skipped with a log message.

## Inputs and matching

| Source | Used for |
|---|---|
| `DATA<n>` (TRANSP "AC" ASCII-compressed, read by `transp_acfile.py`) | `FBM` (guiding centre) and `FBM_PTCL` (particle/FLR) distributions, `EFBM/EFBMB`, `XXKSID`, `BMVOL`, `NTHZSM`, `THBDY0/1`, species |
| `FI_TAR.GZ<n>` member `*_ps_ts1_state.cdf` | flux-surface geometry `R_geo, Z_geo(theta, rho)` used by NUBEAM |
| `FI_TAR.GZ<n>` member `*_debug_xplasma_nubeam.cdf` | `__Bphi_CCW/__Jphi_CCW` orientation NUBEAM actually used |
| `TR.DAT` | `outtim`, `fi_outtim`, `NLJCCW`, `NLBCCW` |

`outtim` and `fi_outtim` have independent lengths. A DATA snapshot is
converted only when a FI archive exists at the same time (1e-5 s); unmatched
outputs are skipped with a warning. Code: `nubeam_outputs.py`, `nubeam_fi.py`.

## Grid (reconstructed from NUBEAM's MC grid)

`NZONE_FB` rows of equal width in `x = rho_tor_norm`; row `i` has
`NTHZSM` differences zones (4, 8, ..., 40 for the NSTX case), uniform in the
geometric poloidal angle over `THBDY0..THBDY1`. The central row is made of
triangles meeting at the axis. Corners are mapped to `(R, Z)` with a bicubic
spline of the plasma-state `R_geo/Z_geo`. Reconstructed zone volumes agree
with NUBEAM's `BMVOL` to <1e-5 for the NSTX case (a check in the tests).

IDS layout per time slice (`distribution(i)/ggd(it)`), three spaces:

1. `primary_standard`: poloidal MC zones. Nodes are `(R, Z)`; faces carry
   `measure` (poloidal area, m^2), `geometry = [row, zone in row, BMVOL in m^3]`
   and boundary edges with neighbours (several across row boundaries where
   zone counts differ). Edges are straight chords between nodes; the true
   edges are constant-`x` / constant-`theta` curves.
2. `secondary_structured`: kinetic-energy bins (`energy_kinetic`, eV).
3. `secondary_structured`: pitch bins (`pitch`).

Grid subsets: `cells` (5, implicit), `energy_bins` (-1), `pitch_bins` (-2) and
`4d_hypervolumes` (200, implicit). Values are on subset 200, flattened with the
zone index fastest, then energy, then pitch (Fortran order).

* `expansion_fd3v[0]` (`m^-3`): particles per m^3 in each phase-space cell,
  `FBM[cm^-3 eV^-1] * 1e6 * dE[eV] * dxi / 2`.
* `expansion[0]` (`m^-6 s^3`): `f = fd3v / d3v`, `d3v = 2 pi (v_hi^3 - v_lo^3)/3 * dxi`
  (non-relativistic, NUBEAM's ion mass `A * amu`).

One `distribution` is written per species and per kind: `gyro_type = 2`
(guiding centre, `FBM`) and `gyro_type = 1` (particle position, `FBM_PTCL`).
`global_quantities` (particles, energy, parallel energy) and zone-row
`profiles_1d` (density, pressures) are derived from the same data.

## Pitch sign and COCOS

NUBEAM's `depsub.f90` states that its `v_parallel` "is relative to the direction
of the toroidal current, not b": `xi_NUBEAM = nsjdotb * xi_B`, with
`nsjdotb = +1` when `Jphi` is parallel to `Bphi`. IMAS `pitch` is along B, so
`xi_IMAS = nsjdotb * xi_NUBEAM`; for `nsjdotb = -1` the pitch bins are reversed
to stay ascending. This is a relation between physical directions and is
independent of the COCOS index (11 for DD3, 17 for DD4); no extra `ip_sign` or
`bt_sign` factor is applied, and the grid is in `R, Z`, so no angle or flux-sign
conversion is involved. `nsjdotb` is taken from NUBEAM's own xplasma flags when
available and cross-checked against `NLJCCW/NLBCCW`.

For the NSTX case (Ip counter-clockwise, Bphi clockwise, `nsjdotb = -1`) the
co-current beam ions have mean NUBEAM pitch +0.46 and IMAS pitch -0.46.

## Species handling

Fast species are keyed by their NUBEAM index (1..`NSFAST`), not by their (Z, A) label, so a
beam triton and a fusion triton are separate entries. `NLFPROD` marks fusion products.
Each species gets a guiding-centre and a particle-position distribution; the `process` is
`NBI` for beam species and a nuclear reaction for fusion products (`D_D_to_T_H` for
protons and tritons, `D_D_to_3He_n` for He3, generic `nuclear` for He4 since D+T, D+He3
and other reactions cannot be told apart).

TRANSP names per-species CDF variables with a suffix, and `nubeam_fi.cdf_species_suffix`
reproduces it (checked against the ITER example CDF):

| Species | Suffix | Species | Suffix |
|---|---|---|---|
| beam H, D, T | `H`, `D`, `T` | fusion proton | `P` |
| beam He3, He4 | `3`, `4` | fusion T, He3, He4 | `T`, `3`, `4` |
| beam Ne, Ar, Kr, Xe | `N`, `A`, `K`, `X` | other heavy beam ions | `Z` |

Beam species use the beam variables (`BPTE_D`, `TQJXB_T`, ...). Fusion products live in a
mirrored `F` family that only has some of them: `FDENS`, `FTRAP` (trapped fraction), `FSTH`
(thermalization particles) and `FPHXB` (J x B torque, scalar), so fusion-product entries get
the thermalization particle rate, the trapped density and the J x B torque, and nothing
else (TRANSP has no per-species power, collision or source variables for them). Unsuffixed
beam totals (`TQBE`, `TQBI`, `TQBTH`, `CURB`, `UCURB`, and the unsuffixed fallbacks of the
per-species names) are attached only when exactly one beam species is written, and never to
fusion products. A beam J x B torque is the volume integral of `TQJXB`; TRANSP's own scalar
`BPHXB` is a different, noisy estimator of the same torque (about 10% apart away from sign
changes, for example 0.184 vs 0.201 N m at 0.745 s for NSTX), so the two are not identical.

Only deuterium beams have been run end to end on real FI output. Hydrogen, tritium and
helium beams, and fusion products, rely on the naming table above and the synthetic
multi-species tests; RF-minority species (index above the beams, not marked as fusion
products) are treated like beam species for naming, which is unverified.

## Additional NUBEAM quantities

Written next to the distribution (all at the FI output times):

| IDS node | Source | Notes |
|---|---|---|
| `global_quantities.power_first_orbit` | CDF `BPLIM_<sp>` | W, positive = loss |
| `global_quantities.power_charge_exchange` | CDF `BPCXI_<sp>` + `BPCXX_<sp>` | W, inside + outside the plasma |
| `global_quantities.torque_tor_j_radial` | CDF `TQJXB_<sp>` x `DVOL` | N m, torque of the fast-ion radial current |
| `profiles_1d.torque_phi_j_radial` | CDF `TQJXB_<sp>` | N m^-2, regridded to the MC rows |
| `global_quantities.current_phi`, `profiles_1d.current_phi` | CDF `CURB` x `DAREA` | A and A m^-2; shielded beam-driven current including electron back-current. `CURB` is the total over all beams, so it is attached only when one fast species is written |
| `profiles_1d.current_fast_phi` | CDF `UCURB` | A m^-2; unshielded fast-ion current (without electron back-current), single-species runs only |
| `global_quantities.collisions.electrons.{power,torque}_thermal*` and `profiles_1d` equivalents | CDF `BPTE_<sp>`, `TQBE` (global), `PBE_<sp>`, `TQBE` (profile) | W, N m, W m^-3, N m^-2; the profile integral equals the scalar (checked) |
| `global_quantities.collisions.ion[0]` and profile | CDF `BPTI_<sp>`, `TQBI`, `PBI_<sp>` | one aggregate entry named "thermal ions (all species)"; TRANSP does not split it per ion species |
| `global_quantities.thermalization`, `profiles_1d.thermalization` | CDF `BSTH_<sp>`, `BPTH_<sp>`, `TQBTH`, `PBTH_<sp>` | particles s^-1, W, N m; profile energy and momentum densities (no particle-density profile exists in TRANSP) |
| `global_quantities.source[0]`, `profiles_1d.source[0]` (NBI process) | CDF `SBDEPMC_<sp>`, `PBDEPMC_<sp>`, `BPHDEP0_<sp>`, `BDEP_<sp>` | the torque only on the particle-position entry; see the caveat below |
| `profiles_1d.trapped.density_fast` | CDF `BTRAP_<sp>` x the density of the same distribution | banana fraction weighted by TRANSP's fast-ion density, rows 0.13 (core) to 0.65 (edge) for NSTX |
| `distribution(i).markers` | `DATA<n>` lost-orbit records | one marker per lost ion record, see below |

CDF scalars and profiles are interpolated in time to the FI output times. Profiles are
regridded from the TRANSP zones to the 10 MC zone rows with a conservative
volume- (or area-) weighted average. The radial fast-ion currents
(`CURBRORB`, `CURBRABD`, ...) have no node in `distributions`; their effect is
carried by the J x B torque. DD 3.39 uses `current_tor` / `torque_tor_j_radial`
and has no `power_*` nodes, so those are skipped there.

Caveats found while checking against the CDF:

* `BPHDEPGC_<sp>` (guiding-centre deposited momentum) and `BPHINJS_<sp>` are of order
  1e13 where the other torques are of order 1 N m, so they are not in N m. They are not
  used: the source torque is taken from `BPHDEP0_<sp>` (= `BPHDP`) at the ion position,
  so the guiding-centre entry has no source torque.
* TRANSP's 1-D `BDENS_<sp>` is the density at the **ion position**, not at the guiding
  centre: it agrees with the particle-position (`FBM_PTCL`) distribution to 1-3% in every
  MC row, whereas the guiding-centre density (`FBM`, `BDENS2_<sp>`) falls much faster
  outward (a factor 22 below `BDENS` in the outermost row, 0.9 in the core; totals agree to
  0.1%). The ions' Larmor radius spreads them beyond the surfaces their guiding centres
  occupy. The trapped density therefore uses only the `BTRAP` fraction times the density of
  the distribution being written.
* The beam-driven current is the total over all beams, so current nodes are skipped when
  more than one fast species is written.
* Beam-target and beam-beam neutron rates (`BTNTS_DD`, `BBNTS_DD`) are already written
  by `imas_summary.py` (`summary.fusion.neutron_rates.dd`); for the NSTX case beam-target
  is 94.8% and beam-beam 5.1% of the DD rate, and the total equals the DD total.

### Lost ions (`markers`)

NUBEAM's lost-orbit records (`*_LOST` in `DATA<n>`, gathered over
`GTIME1_LOST..GTIME2_LOST`, the last time step) are written as markers of the
particle-position (`gyro_type = 1`) distribution of their species, as ion
positions are at the loss point. `weights` is the number of real particles per
record. Coordinates: `r`, `z` (m), `energy_kinetic` (eV), `pitch` (along B, same
sign conversion as the distribution), and private identifiers `time_lost` (-1,
s), `nbi_unit` (-2, 0 = not from a beam) and `prompt_loss` (-3, 1 = prompt
first-orbit loss). Check on NSTX: `sum(weight * energy) / window` equals the
CDF orbit-loss power `BPLIM_D` to 5e-4 (218.2 vs 218.3 kW at 0.51 s), and about
84% of the lost weight is prompt loss.

## Verification and plots

```sh
make -C examples quick                  # includes verify_nubeam_fi.py --quick (synthetic, ~2 s)
python examples/verify_nubeam_fi.py     # also checks the local NSTX sample (~40 s)
python examples/verify_nubeam_outputs.py [CASE_DIRECTORY]
python examples/plot_fast_ion_distribution.py --case-dir CASE --run-id ID --save-dir build/fi-plots
python examples/plot_fast_ion_distribution.py --database NSTX --shot 120967 --run 503 --user DIR
```

Tests cover the AC codec, grid geometry against an analytic torus, time-list
pairing, archive streaming (particle-state members never extracted), pitch
orientation for all Ip/Bphi combinations, unit/ordering/closure of the cell
values, IDS validity and HDF5 round trips for DD 3.39.0 and 4.1.1.

Checked on the NSTX case: total fast ions 5.58e18 (0.51 s) and 5.49e18
(0.745 s); `FBM` and `FBM_PTCL` conserve particles to better than 1e-4. The
integrated distribution is within about 7% of NUBEAM's Monte-Carlo weight
total `WNMCTOT` (5.99e18); this reflects NUBEAM's own binning and was not
reduced further. The exploratory `nubeam_fbm.py` reader is not part of the
conversion path.
