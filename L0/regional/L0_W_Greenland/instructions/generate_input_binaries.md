
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
Next, we will leverage the 
