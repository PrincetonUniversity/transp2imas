## transp2imas
 
The application transp2imas converts a TRANSP NetCDF output file (\<runid\>.CDF), aided
by information and data from the namelist file, to populate an IMAS tree.
The namelist is used if available in the current directroy in the format  \<runid\>TR.DAT.

There are two backends supported for IMAS IDS storage: MDSPLUS and HDF5.
The MDSPLUS storage has 3 types of files: .tree, .datafile, and .characteristics.
The default location is ~/public/imasdb/\<machine\>/\<imas\_version\>/\<occur\>/.

For HDF5, the default location is the current directory, and each IDS is stored separately
as an .h5 file within \<machine\>/\<version\>/\<shot\>/\<run\> directory.

As of now, there is also a fortran code in this repo. For compiling that, TRANSP libraries are required.

The python tool can be used as follows.
Please note that this is tested using only Python 3.10.7 and 3.11.4. Other versions are being tested.

he environment modules valid on PPPL portal are listed here.

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

#### Convert CDF to IDS
TRANSP run output CDF file or its symlink is required to be present in the current directory. Only one CDF file is converted at a time. The input \<runid\>.CDF file and corresponding namelist \<runid\>TR.DAT (if available) are used in some data structures. An IDS is written out only if relevant input is present as valid CDF file; otherwise a message is added in a log.out file.


Python virtual environment with numpy, scipy, netcdf_class, and imass_class or system installations of these are required. Here the \<python_venv\> is the name of the virtual environment directory. The machine name (eg. ITER, D3D) and shot number are retrieved from the input CDF file, in order to avoid errors of passing in mismatched arguments.

The first step is to load modules.

source \<intel_or_gcc_imas_modules\>

If python virtual environment is already setup, skip this this section (of 6 lines).
Please note that this is tested using only Python 3.10.7 and 3.11.4. Other versions are being tested.

```
cd ~/
python3 -m venv  python_venv
(Now we can cd  to anywhere we want)
source ~/python_venv/bin/activate
pip install --upgrade pip
pip install flake8 matplotlib netCDF4 f90nml numpy  PyYAML scipy
#for imas3.39, need this too
pip install pyqt5

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


All options:
 convert.py [-h] [-i INPUT_CDF] [-d INPUT_DIR] [-b {MDS,MDSPLUS,HDF,HDF5}] [-r RUN] [-new] [-o OCCUR] [-u USER]
                  [-v] [-l {debug,info,warn,error,critical}]

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


### IMAS VIZ for visualization

#### use module in pppl portal

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

#### Use local VIZ
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

### Troubleshooting
#### SyntaxError: invalid syntax
This error could be due to missing virtual environment. Activate the venv,
$ source <path_to_venv>/bin/activate

#### python: error while loading shared libraries
Load python module, or install python

#### ModuleNotFoundError: No module named 'imas'
Load list of modules including IMAS module.

https://wpcd-workflows.github.io/IMASviz/index.html
https://docs.psnc.pl/display/WFMS/IMASViz+v2?preview=/17146892/17147822/IMASViz.pdf
