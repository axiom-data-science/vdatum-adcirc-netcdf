# Tidal Constituent NetCDF files from ADCIRC

### Load conda environment

```bash
$ conda env create -f environment.yml
source activate walcc35
```

### Compile program

This will create an `adcirc2netcdf` binary you can copy into any directory with ADCIRC data files and run

```bash
$ make
```

### Create NetCDF

Copy `adcirc2netcdf` and `fort.*` files into their own directory. `adcirc2netcdf` requires `fort.14`, `fort.15` and `fort.53` files.

```bash
$ mkdir -p mydir
$ cp ./adcirc2netcdf mydir
$ cp /foo/bar/fort.* mydir
$ cd mydir
```

### Run `adcirc2netcdf` inside data folder

```bash
$ export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib
$ ./adcirc2netcdf
```

An output file called `adcirc53.nc` will be created. Rejoice.
