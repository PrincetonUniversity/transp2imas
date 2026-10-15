# transp2imas

`transp2imas` converts a TRANSP NetCDF output file and its optional `TR.DAT`
namelist into an IMAS database entry. It supports HDF5 and MDSplus storage and
schema-aware output for the IMAS DD3 and DD4 families.

## Quick start

Activate an environment containing IMAS-Python, NumPy, SciPy, netCDF4, and
f90nml. If a legacy IMAS Access Layer has populated `PYTHONPATH`, clear it so
that the active virtual environment supplies the `imas` package:

```bash
unset IMAS_PREFIX IMAS_VERSION IMAS_HOME IMAS_PATH PYTHONPATH
source "$HOME/venvs/imas/bin/activate"
make test PYTHON="$HOME/venvs/imas/bin/python"
```

Run `make help` for the supported developer workflows. In particular,
`make check` runs fast static checks, `make verify-all` automatically audits
all 10 currently bundled CDFs, `make verify-waves` examines every CDF for EC/IC
wave deposition (including the new TORAY case), and
`make verify-dd` inspects the installed Data Dictionary layouts.
`make verify-equilibrium` checks native CDF TIME/TIME3 coordinates and LCFS
boundaries for every bundled example; see [time-grid guide](docs/time-grids.md).
`make verify-end-to-end PYTHON=python3 DD_VERSION=4.1.1 BACKEND=HDF5`
runs **real** IMAS conversion and readback for every available CDF on machines
with IMAS-Python. It records per-IDS write/skip/failure and checks native
`TIME3` coordinates against stored `equilibrium`/`core_profiles`. See
[real IMAS readback testing](docs/end-to-end.md). This is deliberately separate
from GitHub's source-only CI.
`make verify-sawteeth` checks native `TIME` crash history and q=1 radii
on **all 10 CDFs**; only C-Mod 89846 has observed new crashes. See
[sawteeth mapping](docs/sawteeth.md).
`make verify-pf-active` validates the PF coil/supply mapping and the synthetic
DD3/DD4 machine geometry, without requiring IMAS. `make verify-isolver`
validates the two real ITER/NSTX ISolver current datasets (local large fixtures).


The static checks also reject duplicate Python dictionary literal keys and bare `except:` handlers. Converter resources are managed explicitly: `convert.py` runs `nc2imas` as a context manager so both the IMAS database entry and the source NetCDF file are closed even if a writer fails.

For a repository-managed environment instead, run `make install`; it creates
`.venv` and installs [requirements.txt](requirements.txt).

## Repository layout

- `convert.py` is the command-line converter.
- `nc2imas.py`, `netcdf_class.py`, and `imas_*.py` implement the conversion.
- `examples/` contains inspection, plotting, and verification programs with a
  dedicated Makefile.
- `examples/data/nstx/<runid>/` contains the active general NSTX regression
  cases, one directory per CDF/TR.DAT pair.
- `examples/data/iter/` contains the ITER 38530P94 ISolver CDF/TR.DAT and PT_SOLVER namelist.
- `examples/data/nstx/142202U01/` contains the NSTX 142202U01 ISolver CDF/TR.DAT.
  The >100 MB CDF is a local-only regression fixture unless stored via Git LFS.
- `examples/data/d3d/` contains the TORBEAM and GENRAY EC-wave regression
  cases, plus TORAY and PT_SOLVER/MMM examples.
- `examples/data/cmod/` contains the multi-frequency C-Mod ICRF-wave
  regression case.
- `docs/` contains architecture, development, testing, mapping, release, and
  detailed technical-reference documentation.
- `.github/workflows/` contains the committed GitHub Actions CI and release
  workflows. GitHub only discovers repository workflows from this directory.

## Documentation

Start with [docs/README.md](docs/README.md). The main references are:

- [all-example verification](docs/examples-verification.md);
- [TIME/TIME3 coordinate mapping](docs/time-grids.md);
- [plotting converted IMAS IDSs](docs/plotting.md);
- [real IMAS end-to-end conversion/readback](docs/end-to-end.md);
- [sawteeth crash history and q=1 mapping](docs/sawteeth.md);
- [input XML and code provenance](docs/input-provenance.md);
- [architecture](docs/architecture.md);
- [development guide](docs/development.md);
- [testing and validation](docs/testing.md);
- [mapping workflow and roadmap](docs/mapping-workflow.md);
- [release process](docs/release-process.md);
- [technical mapping reference](docs/technical-reference.md);
- [active PF coil, supply and circuit mapping](docs/pf-active.md);
- [real ITER/NSTX ISolver validation and fixture storage](docs/isolver-fixtures.md);
- [passive PF loop mapping](docs/pf-passive.md);
- [contribution guide](CONTRIBUTING.md).

## Repository Python policy

All tracked Python files, including examples and verification programs, begin
with exactly `#!/usr/bin/python` and include a module-level docstring. The
header convention is enforced by `make check` and GitHub CI. The project itself
requires Python 3; on systems where `/usr/bin/python` is not installed, invoke
programs explicitly with the active Python 3 interpreter.

## License

This project is licensed under the BSD-3-Clause-LBNL (PPPL Revision). See
[LICENSE](LICENSE) for the complete terms.

## IMAS-Python interface

This version uses the current **IMAS-Python** object API rather than the legacy
Access-Layer Python calling pattern. In particular:

- backend constants are taken from `imas.ids_defs`;
- empty IDS objects are created from `DBEntry.factory.new(...)`;
- IDSs are written with `DBEntry.put(ids, occurrence=...)`;
- `DBEntry.create()` and `DBEntry.open()` use exception-based error handling;
- `IMAS_VERSION` (or `IMAS_DD_VERSION`) is treated as an optional Data Dictionary
  version. If neither is set, IMAS-Python uses its configured default.

The code is intended for IMAS-Python 2.x. The pulse/run form of `DBEntry` is kept
because it maps directly to the existing `transp2imas` command-line interface and
remains supported by current IMAS-Python.

## Conversion corrections

The following TRANSP-to-IMAS corrections are included:

- `CUR`, `CUROH`, and `CURBS` are converted from A/cm^2 to A/m^2 using `1.0e4`.
  The previous code accidentally used the cm^2-to-m^2 area factor (`1.0e-4`),
  making current densities smaller by `1.0e8`.
- TRANSP `DVOL` and `DAREA` are per-zone quantities. IMAS `volume` and `area` grid
  fields are enclosed/cumulative quantities, so cumulative sums are written. Raw
  zone values are still used internally for volume and area integrations.
- The existing ITER direction convention (`sgn=-1`) is applied consistently to
  plasma current, toroidal magnetic field, current-density profiles, and signed
  magnetic flux quantities in `core_profiles` and `equilibrium`.

 
The application transp2imas converts a TRANSP NetCDF output file (\<runid\>.CDF), aided
by information and data from the namelist file, to populate an IMAS tree.
The namelist is used if available in the current directory in the format
`<runid>TR.DAT`.

There are two backends supported for IMAS IDS storage: MDSPLUS and HDF5.
The MDSPLUS storage has 3 types of files: .tree, .datafile, and .characteristics.
The default location is ~/public/imasdb/\<machine\>/\<imas\_version\>/\<occur\>/.

For HDF5, the default location is the current directory, and each IDS is stored separately
as an .h5 file within \<machine\>/\<version\>/\<shot\>/\<run\> directory.


The converter requires Python 3.10 or newer and has been exercised with Python
3.10, 3.11, and 3.12.

The environment modules valid on PPPL portal are listed here.

```
#gcc 9.3 (env_gcc.txt)
module load python/3.10.7
module load gcc/9.3.0
module load openmpi/4.0.3
module load netcdf-c/4.7.3
module load szip/2.1.1
module load imas/3.37.1-4.11.0

#intel ifort 2019  (env_intel.txt)
module load python/3.10.7
module load flexlm
module load intel/2019.u3
module load openmpi
module load netcdf-c
module load netcdf-fortran
module load imas/3.37.1-4.11.0

#intel ifort 2019 IMAS 3.39  (env_intel_3.39)
module load python/3.11.4
module load flexlm
module load intel/2019.u3
module load openmpi
module load netcdf-c
module load netcdf-fortran
module load imas/3.39.0-4.11.5

```

## Convert CDF to IDS
TRANSP run output CDF file or its symlink is required to be present in the current directory. Only one CDF file is converted at a time. The input \<runid\>.CDF file and corresponding namelist \<runid\>TR.DAT (if available) are used in some data structures. An IDS is written out only if relevant input is present as valid CDF file; otherwise a message is added in a log.out file.


A Python virtual environment with NumPy, SciPy, netCDF4, f90nml, and
IMAS-Python is required. The machine name (for example, ITER or D3D) and shot
number are read from the CDF to avoid mismatched command-line metadata.

The first step is to load modules.

source \<intel_or_gcc_imas_modules\>

If python virtual environment is already setup, skip this this section (of 6 lines).
```
cd ~/
python3 -m venv  python_venv
(Now we can cd  to anywhere we want)
source ~/python_venv/bin/activate
pip install --upgrade wheel setuptools pip
pip install flake8 matplotlib netCDF4 f90nml numpy  PyYAML scipy
#for imas3.39, need this too
pip install pyqt5 pyqtgraph

```

To run convert script,

```
git clone transp2imas from github (git@github.com:PrincetonUniversity/transp2imas.git).
cd transp2imas
dir=`pwd`
#Add the repo path to PYTHONPATH.
export PYTHONPATH=$PYTHONPATH:$dir
cd <work_directory>

source <python_venv>/bin/activate
.   env_gcc.txt  (or . env_intel.txt)

python convert.py  -i <runid> [ -b MDS/HDF -r <ids_run_num> -new ]
# runid is TRANSP runid; MDS/HDF for backend
#  <runid>.CDF is the CDF file required.
#  Also used if available <runid>TR.DAT
# machine name (eg. ITER, D3D) is retrieved from CDF file.

Remember to use -new if creating fresh, to avoid accumulating old copies making the database huge.
Excluding this option may causes error sometimes, which is being fixed.

All options:
usage: convert.py [-h] [-i INPUT_CDF] [-d INPUT_DIR] [-b {MDS,MDSPLUS,HDF,HDF5}] [-r RUN] [-new] [-o OCCUR]  \
       [-u USER] [-m MACHINE] [-mom {yes,no}] [-zones NZONES] [-ecnml {yes,no}]  \
       [--debug] [--debug-max] [-v] [-l {debug,info,warn,error,critical}]

Convert TRANSP output CDF to IMAS IDS (To run all CDF files in <dir_path> and its subdirs: -i all -d <dir_path>)

options:
  -h, --help            show this help message and exit
  -i INPUT_CDF, --input-cdf INPUT_CDF
                        CDF TRANSP output file, with or without .CDF
  -d INPUT_DIR, --input-dir INPUT_DIR
                        input directory (file can be in its subdir), if file not in pwd
  -b {MDS,MDSPLUS,HDF,HDF5}, --backend {MDS,MDSPLUS,HDF,HDF5}
                        IMAS backend
  -r RUN, --run RUN     run number of o/p IDS (overrides default converted number: eg: Z10 -> 2610)
  -new, --new-db        create new database
  -o OCCUR, --occur OCCUR
                        occurrence
  -u USER, --user USER  username or path of database
  -m MACHINE, --machine MACHINE
                        machine. eg: ITER
  -mom {yes,no}, --moments {yes,no}
                        convert moments to equilibrium, default is no
  -zones NZONES, --nzones NZONES
                        num poloidal grid points per flux surface = max(65,min(251,(2*nzones+1))). Default is NR
  -ecnml {yes,no}, --ec-nml-launch {yes,no}
                        get EC launch position from namelist (for testing)
  --debug               extra debug prints
  --debug-max           All possible debug output including plots
  -v, --verbose         increase output verbosity
  -l {debug,info,warn,error,critical}, --log {debug,info,warn,error,critical}
                        Provide logging level


To convert a CDF file in a different directory :
$ python convert.py -i </path/CDFname>
or, $ python convert.py -i <CDFname> -d <dir_path>

To convert all CDF files in a directory and its subdirectories:
$ python convert.py -i all -d <dir_path>

Examples
$ python convert.py -i 38530P87  -d  <path>/
$ python convert.py -i 38530P87_<text>.CDF  -d  <path>/

The option -u can be used as -u `pwd` or -u ~/work or -u <absolute_path>. If this is not used the default path ~/public/imasdb/<machine>/version/ is used.
The CDF file (extension is optional) can be a renamed one seperated by non-digit, e.g, an underscore.

The namelists in corresponding directories will be used if available, in all cases.

```
The above steps will create MDS+ backend if HDF is not opted for. A log file log.out is written including which variables were not found in the input CDF file, for each IDS written out. In the above command line arguments, the run number is the end part of MDS output file names, while occurrence can be selected only when the ids is accessed, or from VIZ. If an existing database is appended, the size of MDS+ ids_...datafile increases with the current run.

The fixed boundary equilibrium is calculated from moments using option -mom or --moments. This options is disabled by default, since it is not needed if the source CDF file has free-boundary equilibrium data, which is filled in equilibrium.time_slice.profiles_2d.psi. Also this method is very slow, and should be optimized. This method is needed for the fixed boundary equilibrium, since in wich case the data in time_slice.profiles_2d.psi causes error in TRANSP. This option can be enabled by

```
-mom yes  (OR, --moments yes)
```

## IMAS VIZ for visualization

### Use the PPPL portal module

```
source <python_venv>/bin/activate
.   env_gcc.txt (or . env_intel.txt)
module load imas-viz/2.6.1
python  $VIZ_HOME/imasviz/VizGUI/QtVIZ_GUI.py

```

After opening the VIZ GUI, the database selector shows the available IDS in ~/public/imasdb,
for HDF and MDS+ separately. According to user document, by changing username, ids from that user area can be seen and open. Opening IDS from other directories (not from public/imasdb) seem to be not possible.
Selecting and double clicking the desired IDS will just populate the
fields on the upper section. Then, click "open" in order to open the IDS. After opening, the filled IDS's will be shown in blue color. Double-click each of those in order to load, or right-click and select "get .. data for occurrence #". It seems that only the first loaded occurrence shows contents, even though all the occurrence can be loaded. There is no unload option in orer to see the contents of other occurrences. The only way is to re-open the VIZ. Upon double clicking the message window on the bottom part shows the status, and any errors in red color. This is one way for checking if IDS have valid entries.

For HDF5 backend, the VIZ shows the shot number and occurrence, but not the separate .h5 files.

To navigate, in general, right-click on headers to see options, and try if that works. Many of the options don't seem to work.

By default, the "occurrence 0" shows that only time slice 0, as seen by , and right-click on that shows option to group time slices.

Right-clicking on the equilibrium and core-pofiles filled-in IDS's will show a 'plugin' option to view "visualization of 1D nodes along coordinate axis". In that 1D view, select the "Refresh Plot" checkbox at the bottom, and change the time slice to view more times.

On the "Array of time_slice" view also shows the same 1D profile 'plugin' option, as described above. The option "Equilibrium overview" also doesn't seem to work, due to interpolation error.


The widget can be freely taken out from the main window by clicking undock (box) button.

### Use a local VIZ checkout
Steps used to use local copy of imasviz.
```
git clone ssh://git@git.iter.org/vis/viz.git
source activate <imas python env>
pip3 install the following:
PyQt5
pyqtgraph
matplotlib
Not required:
  Sphinx (pip3 install sphinx)
  Sphinx RTD theme (pip3 install sphinx_rtd_theme)

git checkout develop
# bash
export VIZ_PRODUCTION=1
export VIZ_HOME=$PWD
export UDA_DISABLED=1 #disable UDA
make

python $VIZ_HOME/imasviz/VizGUI/QtVIZ_GUI.py
```

## Troubleshooting

### SyntaxError: invalid syntax
This error could be due to missing virtual environment. Activate the venv,
$ source <path_to_venv>/bin/activate

### Python error while loading shared libraries
Load python module, or install python

### ModuleNotFoundError: No module named `imas`
Load list of modules including IMAS module.

https://wpcd-workflows.github.io/IMASviz/index.html
https://docs.psnc.pl/display/WFMS/IMASViz+v2?preview=/17146892/17147822/IMASViz.pdf

## IMAS DD4 / COCOS17 note

When writing IMAS Data Dictionary major version 4, `transp2imas` targets COCOS17.
It reads `NLJCCW` and `NLBCCW` from the TRANSP TR.DAT file to determine the
physical plasma-current and toroidal-field directions.  See
the COCOS 17 sign-policy section for the full sign table.  For old runs without those
controls, use `--ip-direction {ccw,cw}` and `--bt-direction {ccw,cw}`.

Run `python examples/verify_cocos17.py <run>.CDF` for a source-side sign sanity check.

## Selecting the IMAS Data Dictionary (v0.3)

TRANSP2IMAS can now populate either a legacy DD3.39 schema or a DD4 schema
explicitly.  Use, for example:

```bash
python convert.py -i examples/data/nstx/141711Y01/141711Y01.CDF -new -b HDF5 --dd-version 3.39.0
python convert.py -i examples/data/nstx/141711Y01/141711Y01.CDF -new -b HDF5 --dd-version 4.1.1
```

With `--dd-version auto` (default), the converter follows the IMAS environment
or the default Data Dictionary bundled with IMAS-Python.  The exact selected
version is printed at startup.  IDSs are created from the DBEntry factory for
that same version rather than creating one DD and relying on cross-major
conversion at write time.

The EC launcher mapping is explicitly schema-aware: DD3 uses
`beam.identifier` + descriptive `beam.name` + `beam.mode`; DD4 uses the short
identifier in `beam.name`, descriptive text in `beam.description`, and (DD4.1+)
`beam.polarization.o_mode_fraction`.  See the Data Dictionary compatibility section for
the complete mapping table and run
`python examples/verify_dd_compat.py 3.39.0 4.1.1` to inspect the schemas available in
your IMAS-Python installation.

## Extended CDF mapping and provenance (v0.4)

Version v0.4 expands `core_profiles`, `equilibrium`, and `summary` with
high-confidence TRANSP CDF quantities and embeds the original TRANSP/model input
files in IMAS `code.parameters` XML when the IDS has a `code` structure.

Run the source-side checks for both active regression cases before conversion:

```bash
make verify-all
make verify-waves
make verify-dd
```

The source-only verifiers do not require IMAS. `examples/verify_dd_compat.py` exercises
the target DD object structure and therefore must be run in an environment with
the requested IMAS Data Dictionary versions installed.  Detailed source
mappings are documented in the core sources and transport mapping section,
the radiation mapping section, the edge sources mapping section, and the waves mapping section.

### Validated RF `waves` mapping

`make verify-waves` exercises three source-side RF regression CDFs without an
IMAS installation:

- DIII-D `147634B13`: six-gyrotron TORBEAM EC deposition/current drive;
- DIII-D `147634G14`: six-gyrotron GENRAY EC deposition/current drive;
- C-Mod `89846`: two-frequency ICRF/TORIC absorbed-power output.

For EC, `PEECHj(TIME3,X)` is the coherent-wave absorbed power density and
thermal-electron deposition. Its finite-volume integral provides wave power.
`ECCURj(TIME3,X)` is integrated over `DAREA` to obtain cumulative toroidal
driven current.

For IC, `PIC_Fj(TIME3,X)` is the authoritative frequency/antenna-resolved
absorbed wave-power profile and `RFPWRABj(TIME)` is the corresponding absorbed
power trace. `QICHA` is used only as a single-wave fallback if a resolved
`PIC_Fj` profile is unavailable. `PEICH`, `PIICH`, and `QRFBTOT` remain
plasma-heating/source quantities in `core_sources`; they are not summed to
construct wave absorption.

The C-Mod case has two IC frequencies but only an aggregate `XGRID_NPHI` /
`PIC_NPHI` spectrum, so `n_phi` is intentionally left empty rather than
duplicating an unresolved spectrum onto both coherent waves.

Referenced namelist templates (for example `PT_TEMPLATE`, `MMM_TEMPLATE`,
`TGLF_TEMPLATE`, etc.) are resolved relative to the TRANSP input directory.  If
present, they are embedded verbatim and also represented as parsed namelist XML;
if missing, their reference and missing status remain in the provenance payload.



### Audit CDF mapping coverage

To inventory which variables in a TRANSP CDF are referenced, validation-only,
intentionally deferred, or still unclassified, run:

```bash
make audit CASE=141711Y01
```

The command writes full Markdown, CSV, and JSON reports under `build/audit/`.
For a CDF outside the bundled regression-data tree, call the tool directly:

```bash
python examples/audit_mapping_coverage.py /path/to/run.CDF \
    --markdown build/audit/run.md --csv build/audit/run.csv
```

The audit is static and intentionally conservative: it is a planning tool for
choosing the next mapping group, not a claim that every referenced variable has
a complete one-to-one IMAS mapping. See
[the mapping workflow](docs/mapping-workflow.md) for interpretation and policy.

## Passive-structure currents and geometry (`pf_passive`)

The converter can export TRANSP `PS_n` passive-loop current histories, using
**each loop's own time axis** (not just top-level `pf_passive.time`). Without
an authoritative device description, it retains a clearly logged **current-only
legacy mode** and does not invent geometry, turns, resistance or connectivity.

With an explicitly sourced geometry manifest, the converter additionally writes
all `pf_passive/loop/element/geometry` shapes supported by the DD and maps a
current to a physical loop only after its identity and sign have been specified:

```bash
python3 convert.py -i /path/to/RUN.CDF -b HDF5 -new \
    --pf-passive-geometry /path/to/machine-pf-passive.json
make verify-pf-passive
```

The only manifest included in the repository is **synthetic test data**, not a
machine file. All current bundled CDF fixtures lack `PS_n`, so actual device
current/geometry correlation is not yet physically validated. See
[the complete manifest specification, DD3/DD4 outline rules and validation
status](docs/pf-passive.md).

## Pre-0.9.0 development milestones

The following 8.x labels were internal development snapshots created before the
first public pre-1.0 baseline. They are retained here only to document the
evolution of the converter.

### v0.8.3: NBI parallel current

v0.8.3 populates `core_sources/source[nbi]/profiles_1d/j_parallel` from the
TRANSP beam-driven toroidal current `CURB`. The converter changes A/cm² to
A/m², applies the physical plasma-current sign, and converts to the IMAS
`<J.B>/B0` representation using `CONPLJB/BZ`. The same shared conversion is
used for the bootstrap-current source.

This release also appends completed `core_sources` entries instead of shrinking
a preallocated IMAS array of structures. IMAS-Python recreates elements during
`resize()`, which previously left the correct number of source slots but erased
their populated contents before writing.

For `174783Q32`, `CURB` agrees with the sum of the fourteen `BDCnn` beam-current
profiles with correlation 0.9999993, and the resulting parallel-current profile
is finite.

### v0.8.2: workflow and provenance API fix

v0.8.2 fixes the metadata/provenance workflow path for DD 4.1.1.  The
workflow writer now uses `build_code_parameters_xml()` (the actual nc2imas
API) and the new `get_ncattr_text()` helper instead of the nonexistent
`build_transp_parameters_xml()` and `_decode_attr()` names.  The latter bug
had silently suppressed TRANSP/TGLF/GLF23 component metadata because `_attr()`
caught the exception; the former aborted conversion when `workflow` was
reached.  `dataset_fair` also goes through the common code-property helper when
the target DD supplies a code node.

Run `python examples/verify_writer_api.py` as a no-IMAS regression check; it statically
verifies every `self.im.<method>()` call made by writer modules against the
methods implemented by `nc2imas` and its local base classes.

### v0.8.1: distribution_sources namelist accessor fix

v0.8.1 fixes the per-beam species lookup in `imas_distribution_sources.py`.
The v0.8.0 writer accidentally called a nonexistent `get_namelist_value()`
method. It now uses the established `get_namelist_var()` API with the correct
zero-based index (`ib-1`) and treats missing species namelist metadata as
optional rather than aborting conversion. `distribution_sources` remains
independent of the full fast-ion `distributions` IDS: NBI birth/deposition
source moments can be exported from `BDEPnn_*` even when FBM phase-space output
was not saved.

### v0.8.0: NBI distribution sources and metadata/provenance

v0.8.0 adds `distribution_sources` for NUBEAM birth/deposition moments when
`BDEPnn_*` is present, including per-beam particle, energy, toroidal-momentum,
global power/torque, and shinethrough bookkeeping.  It does not synthesize a
full `distributions` IDS if FBM phase-space arrays are absent.

The converter also writes `workflow` and `dataset_fair` where supported, and
attaches whole-IDS source provenance to `ids_properties/provenance`.  Full
TR.DAT/template text and parsed XML remain in `code.parameters`.  See
the NBI kinetic-source mapping section and the metadata and provenance section.



## Documentation

Detailed mapping audits, compatibility notes, provenance design, and migration history are maintained in [the technical reference](docs/technical-reference.md).

## Versioning and releases

The public baseline version is `0.9.0`. Runtime provenance uses
`nc2imas.__version__`; `VERSION` and `TRANSP2IMAS_VERSION` are aliases.

Release Please uses the `simple` strategy because this repository is a
script-based Python project rather than an installable package with
`pyproject.toml`/`setup.py`. The canonical release file is `version.txt`; the
annotated version in `nc2imas.py` is updated in the same release PR.

- Conventional Commits such as `fix: ...` and `feat: ...` determine the next
  semantic version.
- Pushes to `main` update a release pull request and `CHANGELOG.md`.
- Merging that release pull request updates `version.txt`, `nc2imas.py`, and the
  manifest, then creates the corresponding `vX.Y.Z` GitHub release/tag.
- `0.9.0` is the bootstrap baseline. After the one-time `v0.9.0` GitHub release
  exists, routine development should not edit version files or the manifest
  manually.

### Model input XML in workflow and physics IDSs

The complete source TR.DAT and recursively referenced external model namelists
are embedded as XML under `workflow.time_loop.component[0].parameters`.
Individual configured wave codes (TORIC, GENRAY, TORBEAM, CQL3D, FPP, etc.)
are additionally recorded in separate workflow components and in
`waves.code.library[].parameters`. Referenced transport/NUBEAM model input
subsets are associated with their corresponding IDS code libraries. Missing
files are explicitly flagged; comments do not count as active model inputs.
Input provenance is written before physics IDS conversion starts.

See [docs/input-provenance.md](docs/input-provenance.md). Test using
`make verify-inputs` (source-only) and `make verify-dd` (IMAS-Python).


### Unified anomalous theory transport: MMM07 and GLF23

`core_transport.model[]` now preserves the combined MMM07 transport coefficients
and separate Weiland (source `W19`), DRIBM/DRBM, ETG and reduced-MTM
contributions as **model outputs**, not experimental diagnostics. Diffusivities
are mapped in SI units on the original `TIME3/XB` grid, while model growth
rates, frequencies, selected `ky*rho_s`, and untyped aggregate coefficients
remain in model-specific XML, but their `TIME3` and `X` coordinates are now
stored only in the native `core_transport.time` and
`model[].profiles_1d[].grid_flux.rho_tor_norm` fields. The transport
coefficient grid remains `model[].profiles_1d[].grid_d.rho_tor_norm` (`XB`). The combined result and its
components must not be added together. Only physically nonzero model families are instantiated; most examples
contain inactive, zero-filled placeholders. ITER `38530P94` has nonzero
GLF23 `XKEGLF23`, `XKIGLF23` and `ETPHIGLF` on `TIME3/XB` and is now
converted into its own `core_transport.model[glf23]`. The aggregate
`DIFFIGLF` remains supplementary rather than assigned to an invented species.
TGLF awaits a validated output example.

Verify with `make verify-theory`, followed by `make verify-all` and (where IMAS
is installed) `make verify-dd`. See [unified theory mapping details](docs/theory-transport.md).

### Separate neoclassical theory transport and neutral-state profiles

`core_transport.model[]` preserves neoclassical particle and heat fluxes from
`GFLNC_*`, `QFLNC_*`, and `QFLNCC_*` divergences with finite-volume
reconstruction. The same generic theory interface also exports nonzero
GLF23 heat and momentum transport from ITER 38530P94, and the MMM07 models
from DIII-D T16. TGLF requires a validated example. Neoclassical presence
alone does not prove that PT_SOLVER enabled it.

Independently, `core_profiles.neutral[]` contains wall, volume, and
recombination neutral densities and eV temperatures; DD4 toroidal neutral
rotation remains supplemental XML. These are not theory transport models.

Verify separately using `make verify-theory`, `make verify-neoclassical`,
`make verify-neutrals` and `make verify-all`. See the
[theory transport](docs/theory-transport.md),
[neoclassical](docs/neoclassical-transport.md), and
[neutral profiles](docs/neutral-profiles.md) guides.

### Sawtooth history

The `sawteeth` writer exports `TLASTSAW(TIME)` into the native IMAS
`diagnostics.previous_crash_time` history. C-Mod 89846 has 38 observed
updates; valid `XIQ1(TIME)` samples are mapped to the native `q=1` radius.
Trigger mechanisms, mixing/inversion radii, and post-crash profiles are **not**
inferred. The first recorded `TLASTSAW` value is not counted as a newly
observed crash. Run `make verify-sawteeth` or convert C-Mod with
`--ids sawteeth`. See [details](docs/sawteeth.md).
