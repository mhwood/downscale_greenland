
# Generate L0_W_Greenland Input Binaries

## Overview
To generate the input binaries for the regional model, use the following steps:
1. Generate the mitgrid file for the subdomain (manual editing required)
2. Generate the bathymetry file (manual editing required)
3. Generate the reference grid file
4. Generate the initial conditions
5. Generate the external forcing conditions
6. Generate the boundary conditions

## Description
A detailed description of each step is provided here:

### Step 1: Generate the mitgrid file
Script: create_L0_W_Greenland_mitgrid.py

An mitgrid file is a convenient way to store many fields that describe the geometry of the computational grid. Specifically, it stores the following 16 fields: XC, YC, DXF, DYF, RAC, XG, YG, DXV, DYU, RAX, DXC, DYC, RAW, RAS, DXG, DYG.

Since we are generating this regional model on a subset of the existing global model grid, we can subset the existing ECCO mitgrid files to generate the mitgrid file for our subdomain. 

First, we must choose where we want out subdomain to be - this is the only manual step in this process. As can be seen on the following plot, West Greenland is located on faces 3 and 5 of the LLC grid. We will take the top left corner of face 3 and the top left corner of face 5 in order to generate our grid. 

### Step 2: Generate the bathymetry file
Script: create_L0_W_Greenland_bathymetry.py

Next, we will use the migrid to subset the global bathymetry file onto our domain. Further, we will make some modification to the regional bathymetry to avoid numerical issues. Specifically, we will close off Hudson Bay and the Gulf of St Lawrence, and then ensure that all wet areas of the domain are connected. If any areas are not connected, they are filled in. This process ensures that there are no "lakes" or other regions which might otherwise create numerical issues.

### Step 3: Generate a reference grid
For the subsequent binaries (initial conditions, external forcings, and boundary conditions) we will need further information about our grid - namely the angles of the curvilinear grid and the vertical fraction of each cell which is wet. For this, we we will leverage the nice capabilities of the ```mnc``` package to generate a convient reference grid. Here, we will run the model for 1 timestep using the data_for_grid file inside the namelist directory. To run with this file, use the following steps:

```
cd run
ln -s ../namelist/* .
ln -s ../input/* .
ln -s ../build/* .
mpirun -np 12 ./mitgcmuv
```
This run will generate 12 tiles, each containing a piece of the grid. To put them together, use the following code:
```
python3 init_file_creation/stitch_L0_W_Greenland_nc_grid_files_for_ref.py -d ../../../../
```

### Step 4: Generate the inital conditions
Script: create_L0_W_Greenland_pickup.py
Sub-script: create_L0_ECCO_pickup.py

Next, we will create an initial conditions file by subsetting the initial conditions files for the ECCOv5 model. The main tricky part of this procedure is to correctly map the "Lat-Lon-Cap" geometry files onto the regular grid of the regional model. The initial condittions are conveniently stored in the "pickup" file. 

### Step 5: Generate the sea ice inital conditions
Script: create_L0_W_Greenland_seaice_pickup.py
Sub-script: create_L0_ECCO_seaice_pickup.py

The seaice variables are stored in a different file than the usual variables. Heere, we will create an initial conditions file for the sea ice variables by subsetting the initial conditions files for the ECCOv5 model. 

### Step 6: Generate the external forcing conditions
Script: create_L0_W_Greenland_exf.py
Sub-script: create_L0_ECCO_exf.py

### Step 7: Generate the boundary condtions
Script: create_L0_W_Greenland_BCs.py
Sub-script: create_L0_ECCO_BCs.py

In this example, the boundary conditions are generated using monthly mean output from the ECCOv5 solution. Using the files downloaded previously (see the "Pertinent Files" instructions"), the global solution is subsetted onto the regional domain and stored in annual files.

