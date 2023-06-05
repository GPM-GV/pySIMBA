# pySIMBA
System for Integrating Multiplatform Data to Build the Atmospheric Column (SIMBA)

pySIMBA is a software that integrates multiple precipitation sensors from the 
ground and space into a user-defined Cartesian grid and outputs a netCDF data file.
This product file is the initial step to establish an easy method to 
conduct precipitation science and cross-validate sensors on a common coordinate system. 
Python based SIMBA was adapted from Interactive Data Language (IDL) code developed by 
Wingo et al. 2018. Please refer to https://gpm-gv.gsfc.nasa.gov/SIMBA/ for more information. 
Here are several peer-reviewed literature that demonstrates the framework of SIMBA:

Pabla, C. S., D. B. Wolff, D. A. Marks, S. M. Wingo, and J. L. Pippitt, 2022:
    GPM Ground Validation at NASA Wallops Precipitation Research Facility.
    J. Atmos. Oceanic Technol., 39, 1199-1215. https://doi.org/10.1175/JTECH-D-21-0122.1 

Wingo, S. M., W. A. Petersen, P. N. Gatlin, C. S. Pabla, D. A. Marks, and 
    D. B. Wolff, 2018: The System for Integrating Multiplatform Data to Build 
    the Atmospheric Column (SIMBA) precipitation observation fusion framework. 
    J. Atmos. Oceanic Technol., 35, 1353–1374. https://doi.org/10.1175/JTECH-D-17-0187.1

# Dependencies
* [PyART](https://arm-doe.github.io/pyart/)
* [NumPy](https://www.numpy.org)
* [pandas](https://pandas.pydata.org/) 
* [Xarray](https://docs.xarray.dev/en/stable/)
* [SciPy](https://www.scipy.org)
* [netCDF4](https://github.com/Unidata/netcdf4-python)
* [h5py](https://docs.h5py.org/en/stable/) 

# Installing pySIMBA
We suggest creating an environment, pysimba. Install the dependencies in this environment.

Clone pysimba repository:
    
    git clone https://github.com/GPM-GV/pySIMBA.git
    
If you would like to test with data in the notebook, download [sample_data](https://gpm-gv.gsfc.nasa.gov/SIMBA/pySIMBA/sample_data/SIMBA_WFF_2020_1029.tgz).

# Running pysimba

On command line,
    
    pysimba.py YYYY MM DD [--params_dict PARAMS_DICT]
    
    pysimba.py 2022 10 29 --params_dict wff_center_params_dict.txt
    
Required Parameters:

    YYYY - 4 digit year

    MM - 2 digit month

    DD - 2 digit day

Optional Parameters:

    --params_dict params_dict.txt - Parameters dictionary to define the column grid

# Output from pysimba:

There are two files:

    column_[main_plat_name]_[center_on]_YYYYMMDD_HHMM.nc

    column_[main_plat_name]_[center_on]_YYYYMMDD_HHMM.nc_ncHeader.txt

Example output:

    column_NPOL_WFFPad_20201029_0736.nc

    column_NPOL_WFFPad_20201029_0736.nc_ncHeader.txt
