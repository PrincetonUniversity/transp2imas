# transp2imas

**TRANSP → IMAS converter** for transporting TRANSP NetCDF (`<runid>.CDF`) results into IMAS Interface Data Structures (IDSs). When available, the matching `<runid>TR.DAT`, PT_SOLVER namelists, and other referenced model inputs provide model configuration and provenance. The converter supports IMAS Data Dictionary 3.x/4.x families and HDF5 and MDSplus backends.

## Requirements and setup

Python 3.10+ with IMAS-Python, NumPy, SciPy, netCDF4 and f90nml. See [requirements.txt](requirements.txt) and [setup details](docs/legacy-usage.md).

```bash
git clone git@github.com:PrincetonUniversity/transp2imas.git
cd transp2imas
make install
source .venv/bin/activate
```

To install into an already active environment (including the command-line
entry points), use:

```bash
python -m pip install .
```

For development, use `python -m pip install -e .` so edits in the checkout are
immediately available in that environment.

Alternatively, activate an existing IMAS-Python environment. On systems with an old IMAS Access Layer in `PYTHONPATH`, remove the conflicting setting before running the converter:

```bash
unset PYTHONPATH
source "$HOME/venvs/imas/bin/activate"
python3 -c 'import imas; print(imas.__file__)'
```

## Quick examples

Convert a sample ITER discharge to HDF5, targeting IMAS DD 4.1.1:

```bash
python3 convert.py -i examples/data/iter/38530P94.CDF -b HDF5 -new --dd-version 4.1.1
```

Convert only the transport IDS (including PT_SOLVER model metadata when available):

```bash
python3 convert.py -i examples/data/iter/38530P94.CDF -b HDF5 -new --dd-version 4.1.1 --ids core_transport
```

`--report-json FILE` writes a conversion journal that lists, per writer, the optional CDF variables that were absent and the fields that were skipped with their reasons. A CDF variable that exists but cannot be mapped (for example with the wrong radial size) aborts the conversion; `--lenient-mapping` records it as a skipped field instead. See [required versus optional mappings](docs/mapping-status.md); the journal also lists the CDF variables each writer read, and the [mapping registry](docs/mapping-registry.md) states what is mapped, deferred and verified.


Convert an NSTX example with one active PT_SOLVER region:

```bash
python3 convert.py -i examples/data/nstx/141711Y05/141711Y05.CDF -b HDF5 -new --dd-version 4.1.1 --ids core_transport
```

`-new` creates the entry for the selected database, pulse and run (`DBEntry.create()`) and, on the HDF5 backend, discards what an existing entry at that location held (checked: converting `--ids summary -new` after `--ids tf -new` leaves only `summary`). It does not allocate a run number; pick one with `-r`. Without `-new` the entry must already exist (otherwise the run stops with "Make sure it exists, else use -new"), and each converted IDS is added to it or replaces the same IDS, leaving the others alone. So use `-new` for the first conversion of a run and omit it to add further IDSs. Run `python3 convert.py --help` for all options, including run selection and MDSplus storage.

`--dd-version` selects the target Data Dictionary. When it is omitted or `auto`, `IMAS_VERSION` / `IMAS_DD_VERSION` from the environment is used, else the IMAS-Python default (3.39.0 on many installations), so pass `--dd-version 4.1.1` explicitly for DD 4.

## What gets exported

`transp2imas` provides IDS mappings for plasma profiles and equilibrium, sources, transport, heating systems, magnetic coils and passive structures, sawteeth and simulation metadata, subject to source-data availability and Data Dictionary support. It also preserves provenance of TRANSP and external-model inputs when those files are available.

For predictive TRANSP runs, `core_transport` occurrence **0** carries the combined transport solution. Additional occurrences represent active PT_SOLVER regions (axial, confinement, edge), each with anomalous and neoclassical model selection. Whole-radius model profiles may be retained, with regional validity ranges and numerical settings in supplementary XML. Available PT_SOLVER runtime channels (`NPTRITER`, `PT_CPTIME`, `PT_WTIME`) are serialized losslessly as supplementary XML together with their native `TIME` coordinate (not the `TIME3` grid), and, for DD >= 3.42, also as named control parameters of `transport_solver_numerics` on that native `TIME` grid (timers keep their CDF hours; no residuals or convergence status are invented; DD 3.39.0 has unnamed control parameters, so nothing is written there). Missing quantities are not fabricated. See [PT_SOLVER region mapping](docs/pt-solver-regions.md) and [input provenance](docs/input-provenance.md).

## Verification

```bash
make check                 # fast source checks
make verify-pt-solver      # region and model provenance tests
make verify-all            # example mapping tests
make verify-end-to-end PYTHON=python3 DD_VERSION=4.1.1 BACKEND=HDF5
```

The end-to-end target requires a working IMAS-Python installation and tests real database writes/readback; source-only tests cannot establish successful serialization.

## Data curation (initial framework)

Screen input provenance, native grids, profile/gradient validity, declared composition and model
controls without modifying the CDF or writing an IMAS entry:

```bash
python3 curate.py -i RUN.CDF --trdat RUNTR.DAT --output-dir build/curation/RUN
make verify-curation
```

Optional `--evidence my-evidence.json` supplies diagnostic/uncertainty declarations and configuration
references. Reports include JSON, Markdown and local evaluated/flagged masks. Missing evidence stays
unassessed; a screening pass is not scientific certification. Curation needs no IMAS:
`pip install -r requirements-curation.txt`; check a sidecar without a CDF with
`python3 curate.py --validate-only --evidence my-evidence.json`, or run a rule subset with `--policy`.
See [curation guide](docs/curation/README.md).

A local web workbench reads curation and conversion reports and, with `--enable-run`, lets you choose a CDF, the IMAS
backend, Data Dictionary version and IDSs, preview and run the conversion with a live log, and plot the converted data
(`python3 -m gui --root build`); see [the workbench](docs/gui.md).

## Repository layout and documentation

- `convert.py` — command-line entry point.
- `nc2imas.py`, `netcdf_class.py`, `imas_*.py` — conversion and IDS writers.
- `examples/` — test cases, scripts and regression fixtures.
- `docs/` — architecture, mappings, verification and detailed usage.
- [`docs/README.md`](docs/README.md) — documentation index.
- [`docs/end-to-end.md`](docs/end-to-end.md) — real IMAS write/readback tests.
- [`docs/pt-solver-regions.md`](docs/pt-solver-regions.md) — PT_SOLVER mapping and occurrence semantics.
- [`CHANGELOG.md`](CHANGELOG.md) — version history and implementation changes.
- [`docs/legacy-usage.md`](docs/legacy-usage.md) — detailed historical setup and troubleshooting notes.

## License

BSD-3-Clause-LBNL (PPPL Revision). See [LICENSE](LICENSE).
