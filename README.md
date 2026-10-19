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

Convert an NSTX example with one active PT_SOLVER region:

```bash
python3 convert.py -i examples/data/nstx/141711Y05/141711Y05.CDF -b HDF5 -new --dd-version 4.1.1 --ids core_transport
```

The `-new` flag creates a fresh entry; omit it only when you intend to work with an existing entry. Run `python3 convert.py --help` for all options, including run selection and MDSplus storage.

## What gets exported

`transp2imas` provides IDS mappings for plasma profiles and equilibrium, sources, transport, heating systems, magnetic coils and passive structures, sawteeth and simulation metadata, subject to source-data availability and Data Dictionary support. It also preserves provenance of TRANSP and external-model inputs when those files are available.

For predictive TRANSP runs, `core_transport` occurrence **0** carries the combined transport solution. Additional occurrences represent active PT_SOLVER regions (axial, confinement, edge), each with anomalous and neoclassical model selection. Whole-radius model profiles may be retained, with regional validity ranges and numerical settings in supplementary XML. Available PT_SOLVER runtime channels such as `NPTRITER`, `PT_CPTIME`, and `PT_WTIME` may be serialized as supplementary XML, rather than standardized IMAS fields. Missing quantities are not fabricated. See [PT_SOLVER region mapping](docs/pt-solver-regions.md) and [input provenance](docs/input-provenance.md).

## Verification

```bash
make check                 # fast source checks
make verify-pt-solver      # region and model provenance tests
make verify-all            # example mapping tests
make verify-end-to-end PYTHON=python3 DD_VERSION=4.1.1 BACKEND=HDF5
```

The end-to-end target requires a working IMAS-Python installation and tests real database writes/readback; source-only tests cannot establish successful serialization.

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
