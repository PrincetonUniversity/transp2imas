# Automatically verify every bundled TRANSP example

The verifier discovers all `examples/data/**/*.CDF` files recursively. It does
**not** rely on an out-of-date list of NSTX or DIII-D paths. Adding a new CDF
will automatically add it to the next source-verification run.

From the repository root:

```bash
make check                 # quick policy, writer, provenance and plot tests
make verify-all            # all 10 example CDFs + both real ISolver PF cases
make verify-waves          # run EC/IC closures for every CDF; N/A when absent
make verify-core-profiles  # verify Te/Ti peaking on TIME/TIME3 in all CDFs
make verify-inputs         # input XML preservation for every CDF
make verify-deep           # OPTIONAL, slow: additional legacy numerical scripts on all CDFs
```

`make verify-all` creates `build/verification/summary.md`, `summary.json`, and
`RUNID.log` for every case. These report CDF grid validity, EC/IC absorption
and current closure, all active/passive channel mappings, input XML and
source-file provenance. No CDF is dropped if a previous case fails. Source
checks do not depend on IMAS-Python. The real ITER/NSTX detailed PF checks also
run via `verify_isolver_mapping.py` and generate reports in `build/pf-isolver/`.

`make verify-deep` adds optional legacy source scripts for every CDF. These may
read large arrays multiple times and are significantly slower. The aggregate
report distinguishes PASS, WARN, FAIL, and N/A. A missing TR.DAT is reported,
not silently substituted from another run. For a strict input/archive check,
use `python3 examples/verify_all_examples.py --quick --strict-inputs`; this is
expected to fail until `147634T16TR.DAT` is supplied. For a strict 0.5% EC
current check, use `--strict-current`; ITER 38530P94 currently has a 0.543%
max time-point residual. Default verification reports that residual as a WARN
(up to 0.6%) rather than calling it an exact numerical closure.

## Current supplied examples

| Run | Device | Particular test / input |
|---|---|---|
| 89846 | C-Mod | TORIC two-channel ICRF power absorption |
| 147634B13 | DIII-D | TORBEAM EC; adjacent `ech_inp.nml` |
| 147634G14 | DIII-D | GENRAY EC; adjacent `genray_ech.in` |
| 147634T13 | DIII-D | TORAY EC; `toray.in`; TIME != TIME3 |
| 147634T16 | DIII-D | PT_SOLVER/MMM; `t16_pt.nml`, `mmm.nml`; **TR.DAT missing** |
| 38530P94 | ITER | 11 active PF circuits, 119 passive loops, EC + IC waves |
| 141711Y01 | NSTX | Reference source/transport checks |
| 141711Y02 | NSTX | Reference source/transport checks |
| 141711Y05 | NSTX | PT_SOLVER `y05_pt.nml` referenced from TR.DAT |
| 142202U01 | NSTX | 9 active PF circuits, 82 passive loops |

The adjacent TORAY/GENRAY/TORBEAM model files and unreferenced T16 transport
namelists are preserved in input XML with an explicit
`role="ADJACENT_UNREFERENCED"` marker. This demonstrates presence of a source
file, **not** that the particular file was executed by TRANSP. Referenced
`PT_TEMPLATE` files have priority over same-directory guesses.

For IMAS Data Dictionary validation and actual writing/reading, use an IMAS
runtime separately: `make verify-dd` followed by conversion and readback.

## Equilibrium TIME/TIME3 regression

Run `make verify-equilibrium` to check the native CDF dimension/coordinate of
all equilibrium sources and to evaluate the actual LCFS moment functions on
**every** bundled example.  It explicitly covers DIII-D 147634T13, whose
`TIME` (4 points) and `TIME3` (7 points) do not have the same number of
samples.  A synthetic test also catches the case of equal lengths but unequal
sample times.  The new `equilibrium-time-grid` check runs within `make
verify-all` and appears in each run's JSON verification report.

See [TIME/TIME3 mapping](time-grids.md).

## Neoclassical and neutral source verification

The per-run `verify-all` report includes `neoclassical_neutrals` with source
channel names, explicit wall/volume/recombination neutral populations and
units/native TIME3 grid checks. The stand-alone target
`make verify-neoclassical-neutrals` produces a focused JSON/Markdown report
for all CDFs. It also executes actual writer methods against IMAS-shaped
objects, and does not imply a completed IMAS DD4.1.1 write/readback.


### MMM07 reduced anomalous transport

`make verify-mmm` inspects every available example. Of the ten current cases,
DIII-D `147634T16` has nonzero MMM07 (Weiland, DRIBM/DRBM, ETG, MTM)
outputs; other declared MMM source arrays are zero-filled placeholders.
Source component closure and model code XML are checked, but these source
checks are not a substitute for an actual IMAS write/readback.
See [MMM mapping](mmm-transport.md).
