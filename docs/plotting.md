# Plotting existing IMAS equilibrium data with plot_ids.py

`examples/plot_ids.py` **reads an already converted IMAS database**. It is not
a CDF viewer. Run it in an environment where IMAS-Python and Matplotlib are
installed, typically after activating `~/venvs/imas/bin/activate`.

From the root of `transp2imas`:

```bash
# Default: 2D poloidal flux at equilibrium time slice 0, interactive window
python3 examples/plot_ids.py 38530P94 ITER HDF5 psirz

# Plot three magnetic field components in separate subplots
python3 examples/plot_ids.py 38530P94 ITER HDF5 bfield --time-index 2

# Write a noninteractive safety-factor plot to disk
python3 examples/plot_ids.py 38530P94 ITER HDF5 q --time-index 0 \
    --dd-version 4.1.1 --save build/iter_q.png

# NSTX magnetic equilibrium
python3 examples/plot_ids.py 142202U01 NSTX HDF5 psirz --time-index 0
```

Syntax: `plot_ids.py RUNID MACHINE [BACKEND] [FIELD] [--time-index N]`
with backend `HDF5` (default), `HDF`, `MDSPLUS`, or `MDS`. Fields:
`psirz`, `bfield`, `q`, `psi1d`, `pressure`, and `j_tor`. The time-index is
**zero-based** and refers to `equilibrium.time_slice[]`, not seconds; a
`time_index` out of range gives an informative error.

The run number is derived from the TRANSP suffix by the same routine as
`convert.py`. For example `38530P94` corresponds to IMAS shot 38530 and run
1694. If conversion used `-r 1`, this derived run will not match; plot the
matching run with `--run 1`, e.g.:

```bash
python3 examples/plot_ids.py 38530P94 ITER HDF5 psirz --run 1 --time-index 0
```

`--save` accepts PNG or PDF output, avoids GUI interaction, and makes the
parent folder automatically. `--user` accepts the IMAS DBEntry username or
explicit database location, and `--occurrence` selects an IDS occurrence.
Use `python3 examples/plot_ids.py --help` for the actual CLI.

The plots are restricted to the equilibrium IDS. Two-dimensional plots need
`equilibrium.time_slice[i].profiles_2d[0]`; when unavailable, try one of the
1D profile fields instead. `plot_ids.py` has a source-independent mock-
plotting regression (`examples/verify_plot_ids.py`). An actual database read
and full DD 4.1.1 plot still require the user's IMAS environment.
