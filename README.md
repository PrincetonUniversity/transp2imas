## transp2imas
 
The application transp2imas converts a TRANSP NetCDF output file (\<runid\>.CDF), aided
by information and data from the namelist file, to data in an IMAS tree. 
The namelist is used if available in the current directroy in the format  \<runid\>TR.DAT.

There are two backends supported for IMAS IDS storage: MDSPLUS and HDF5.
The MDSPLUS storage has 3 types of files: .tree, .datafile, and .characteristics.
The default location is ~/public/imasdb/\<machine\>/\<imas\_version\>/\<occur\>/.

For HDF5, the default location is the current directory, and each IDS is stored separately
as an .h5 file within \<machine\>/\<version\>/\<shot\>/\<run\> directory.

As of now there is also a fortran code in this repo. For compiling that, TRANSP libraries are required.

The python scripts can be used as follows:

The environment modules are listed here.

```
#env_gcc.txt
module load gcc/9.3.0
module load openmpi/4.0.3
module load netcdf-c/4.7.3
module load szip/2.1.1
module load imas/3.37.1-4.11.0

#env_intel.txt
module load flexlm
module load intel/2019.u3
module load openmpi
module load netcdf-c
module load netcdf-fortran
module load imas/3.37.1-4.11.0
```

#### Convert CDF to IDS
Only one CDF file is converted at a time. Input CDF file and corresponding namelist (if available) are used.

```
source <python_venv>/bin/activate
.   env_gcc.txt  (or . env_intel.txt)
python convert.py -m <machine> -i <runid> [ -b MDS/HDF -r <ids_run_num> -new ]
# machine : ITER/D3D etc; runid is TRANSP runid; MDS/HDF for backend
#  <runid>.CDF is the CDF file required.
#  Also used if available <runid>TR.DAT


All options:
  python convert.py [-h] [-i INPUT_CDF] [-b BACKEND] [-m MACHINE] [-r RUN] [-new] [-o OCCUR] [-d] [-u USER] [-v] [-l LOG]

options:
  -h, --help            show this help message and exit
  -i INPUT_CDF, --input-cdf INPUT_CDF
                        CDF TRANSP output file, with or without .CDF
  -b BACKEND, --backend BACKEND
                        backend: MDS or HDF
  -m MACHINE, --machine MACHINE
                        machine
  -r RUN, --run RUN     run number
  -new, --new_db        new database
  -o OCCUR, --occur OCCUR
                        occurrence
  -d, --default-val     use default
  -u USER, --user USER  username or path of database
  -v, --verbose
  -l LOG, --log LOG     Provide logging level

```
The above steps will create MDS+ backend if HDF is not passed. A log file log.out is written including which variables were not found in the input CDF file, for each IDS written out.

### IMAS VIZ for visualization

#### use module in pppl portal

```
source <python_venv>/bin/activate
.   env_gcc.txt (or . env_intel.txt)
module load imas-viz/2.6.1
python  /usr/pppl/imas-viz/2.6.1/imasviz/VizGUI/QtVIZ_GUI.py

```

After opening the VIZ GUI, the database selector shows the available IDS in ~/public/imasdb,
for HDF and MDS+ separately. Opening IDS from other directories seem to be not possible.
Selecting and double clicking the desired IDS will just populate the
fields on the upper section. Click 'open' in order to open the IDS. After opening the filled IDS
will be shown in different color. Double-click those in order to load. Upon double clicking the
message window on the bottom part shows status, and errors if any in red color. This is one way for
checking if IDS have valid entries.


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
https://docs.psnc.pl/display/WFMS/IMASViz+v2?preview=/17146892/17147822/IMASViz.pdf
