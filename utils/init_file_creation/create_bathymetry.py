
import os
import numpy as np
import netCDF4 as nc4
from pyproj import Transformer
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import ast

def create_bathymetry_file(config_dir,level_name,model_name):

    print('Generating the bathymetry for the ' + model_name + ' model from GEBCO and BedMachine')

    f = open(os.path.join(config_dir, 'domain_sizes.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)
    L_size = size_dict[model_name]
    n_rows = L_size[0]
    n_cols = L_size[1]

    file_path = os.path.join(config_dir, 'mitgrids', model_name + '.mitgrid')
    print(file_path)
    entire_grid = np.fromfile(file_path, dtype='>f8')
    entire_grid = np.reshape(entire_grid, (16, n_rows + 1, n_cols + 1))
    print(np.shape(entire_grid))

    Lon_C = entire_grid[0, :, :]
    Lat_C = entire_grid[1, :, :]

    Lon_C = Lon_C[:-1, :-1]
    Lat_C = Lat_C[:-1, :-1]

    # x2, y2 = transformer.transform(polygon_array[:, y_column], polygon_array[:, x_column])
    transformer = Transformer.from_crs('EPSG:' + str(4326), 'EPSG:' + str(3413))
    XC, YC = transformer.transform(Lat_C.ravel(), Lon_C.ravel())

    XC = XC.reshape(np.shape(Lon_C))
    YC = YC.reshape(np.shape(Lat_C))

    ###############################################################
    # First interpolate with BedMachine

    print('   - Interpolating from the BedMachine dataset')

    bathy = -10000*np.ones_like(XC)

    bedMachine_file = '/Users/michwood/Documents/Research/Data Repository/Greenland/Bathymetry/' \
                      'Bed Machine/Versions/BedMachineGreenland-2021-04-20.nc'

    ds = nc4.Dataset(bedMachine_file)
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    bed = ds.variables['bed'][:,:]
    surface = ds.variables['surface'][:,:]
    ds.close()

    min_x_index = np.argmin(np.abs(x - np.min(XC)))
    max_x_index = np.argmin(np.abs(x - np.max(XC)))
    max_y_index = np.argmin(np.abs(y - np.min(YC)))
    min_y_index = np.argmin(np.abs(y - np.max(YC)))

    x = x[min_x_index:max_x_index]
    y = y[min_y_index:max_y_index]
    bed = bed[min_y_index:max_y_index, min_x_index:max_x_index]
    surface = surface[min_y_index:max_y_index, min_x_index:max_x_index]

    # this is a manual nearest neighbor method
    bedmachine_resolution = 150
    for i in range(np.shape(XC)[0]):
        for j in range(np.shape(YC)[1]):
            x_index = np.argmin(np.abs(x - XC[i, j]))
            y_index = np.argmin(np.abs(y - YC[i, j]))
            dist = ((x[x_index]-XC[i,j])**2 + (y[y_index]-YC[i,j])**2)**0.5
            if dist<=bedmachine_resolution*np.sqrt(2)/2:
                if surface[y_index,x_index]>0:
                    bathy[i,j] = 0
                else:
                    bathy[i, j] = bed[y_index, x_index]

    # # this uses scipy
    # X,Y = np.meshgrid(x,y)
    # bed_interp = griddata(np.column_stack([X.ravel(),Y.ravel()]),bed.ravel(),(XC,YC))
    # bathy[~np.isnan(bed_interp)] = bed[~np.isnan(bed_interp)]

    # C = plt.imshow(bathy,origin='lower')
    # plt.colorbar(C)
    # plt.show()

    ###############################################################
    if np.any(bathy<-9999):
        print('   - Interpolating from the GEBCO dataset')
        gebco_file = '/Users/michwood/Documents/Research/Data Repository/Global/Bathymetry/GEBCO/gebco_2021/GEBCO_2021.nc'
        ds = nc4.Dataset(gebco_file)
        lon = ds.variables['lon'][:]
        lat = ds.variables['lat'][:]
        elev = ds.variables['elevation'][:,:]
        ds.close()

        min_lon_index = np.argmin(np.abs(lon - np.min(Lon_C)))
        max_lon_index = np.argmin(np.abs(lon - np.max(Lon_C)))
        min_lat_index = np.argmin(np.abs(lat - np.min(Lat_C)))
        max_lat_index = np.argmin(np.abs(lat - np.max(Lat_C)))

        lon = lon[min_lon_index:max_lon_index]
        lat = lat[min_lat_index:max_lat_index]
        elev = elev[min_lat_index:max_lat_index, min_lon_index:max_lon_index]

        # plt.imshow(elev,origin='lower')
        # plt.show()

        # this is a manual nearest neighbor method
        for i in range(np.shape(XC)[0]):
            for j in range(np.shape(YC)[1]):
                lon_index = np.argmin(np.abs(lon - Lon_C[i, j]))
                lat_index = np.argmin(np.abs(lat - Lat_C[i, j]))
                if elev[lat_index,lon_index]<0:
                    bathy[i,j] = elev[lat_index,lon_index]
                else:
                    bathy[i,j] = 0

    C = plt.imshow(bathy, origin='lower')
    plt.colorbar(C)
    plt.title(np.shape(bathy))
    plt.show()

    output_dir = os.path.join(config_dir,level_name,model_name)
    if 'input' not in os.listdir(output_dir):
        os.mkdir(os.path.join(output_dir,'input'))

    output_dir = os.path.join(output_dir,'input')

    output_file = os.path.join(output_dir, model_name+'_bathymetry.bin')
    bathy.ravel('C').astype('>f4').tofile(output_file)

