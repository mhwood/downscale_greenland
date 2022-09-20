# Running the ECCOv5 Darwin model with diagnostics_vec


### Step 1: Make masks and add them to a run/dv dir
The `diagnostics_vec` pkg ingests 2D masks which delineate where the model output is requested. The masks are stored in compact format, identical to the format of the bathymetry file. The masks are identically zero except where output is resquested, in which case the mask is numbered sequentially (e.g. along a boundary). 

### Step 2: Add the data.diagnostics_vec file
After the masks are constructed 

### Step 2: Add diagnostics_vec to data.pkg

