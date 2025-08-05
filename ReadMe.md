# EFDC+

## Introduction
Environmental Fluid Dynamics Code (EFDC) is a multifunctional surface water modeling engine that includes hydrodynamic, sediment-contaminant, and eutrophication components designed to simulate aquatic systems in one, two, and three dimensions. 

Over the years, DSI has continued developing EFDC and into what is now known as EFDC+. To learn more about different versions of EFDC, please visit [A Review of EFDC Versions - EEMS Blog](https://www.eemodelingsystem.com/efdc-insider-blog/a-review-of-efdc-versions) 

## Getting Started


### Linux

After you [Prepare your environment](https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-intel-oneapi-hpc-linux/top/before-you-begin.html#before-you-begin_HPCCMAKE).

Execute the following:
```bash
git clone https://github.com/dsi-llc/EFDC_Plus.git
cd EFDC_Plus
chmod -x toolkit-setup.sh
```

Run the setup script and pass in your package manager parameter (ubuntu users would use apt).  
If you want to be able to compile the code, use option `-c`. 
If you just want to be able to run efdc, use option `-r`.

```bash
toolkit-setup.sh <-r or -c> <package manager name> 
# example: toolkit-setup -r apt
```

>Alternatively, you can download and follow the installation steps from intel:  
[Intel OneApi Base Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html)  
[Intel OneApi HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html)


Load the intel environment variables.
```bash
source /opt/intel/oneapi/setvars.sh
```
Install NetCDF library (for Ubuntu users) \
Run the command to update the package lists for upgrades for packages that need upgrading
```bash
sudo apt update
```
Install hdf5 and hdf5-devel libraries using the command: 
```bash
sudo apt install hdf5-tools hdf5-helpers libhdf5-dev libhdf5-doc libhdf5-serial-dev
```
Install the libnetcdf-dev:
```bash
sudo apt install libnetcdf-dev
```
Install m4 package:
```bash
sudo apt-get install m4
```
Download netcdf-fortran-4.5.2 from the website https://github.com/Unidata/netcdf-fortran/releases/tag/v4.5.2 \
Go to the netcdf-fortran-4.5.2 folder and create config-intel.sh as the following:
```bash
export FC=mpiifort
export F77=mpiifort
export F90=mpiifort
./configure
```
Run config.sh
```bash
./config.sh
```
Build NetCDF-Foxtran \
Go to the NetCDF-Fortran source code folder and run the following commands
```bash
make check
make install
```
Run the below command for checking whether NetCDF and NetCDF-Foxtran are installed propertly or not
```bash
nc-config --all
```
Add the DNCOUT flag to the FFLAGS section in the makefile for supporting NetCDF

#### _Build EFDC_
```bash
make -f Makefile
chmod -x efdc.x  # this allows the built file to be executed as a program.
```

#### _Run EFDC_

The run command structure is 
```bash
cd /path/to/project
mpiexec -n <number of nodes> path/to/efdc.x -NT<number of omp threads>
```
>An example of a run command for a model with 4 mpi domains, running with 6 omp threads per domain (24 cpu cores total) would look like:  
`mpiexec -n 4 ~/code/efdc/efdc.x -NT6`

<hr>

### Windows

* Install Visual Studio:  
[Visual Studio Downloads](https://visualstudio.microsoft.com/downloads/)

* Install the Intel Toolkits:  
[Intel OneApi Base Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html)  
[Intel OneApi HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html)

* Clone the EFDC_Plus repository: https://github.com/dsi-llc/EFDC_Plus.git
* Open the .sln file at the root of the repository.

## Contribute
The open source availability of this code will make it easier for scientists, researchers, and developers to contribute to the code and build more trust in their models. We welcome all the opportunities to collaborate. If you would like to contribute to the source code development, please clone the repository and submit pull requests as needed. For more active contribution and role, please email admin@ds-intl.biz
