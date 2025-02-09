
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import shapefile as sf
from scipy.interpolate import interp2d, griddata

def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    drF = ds.variables['drF'][:]
    ds.close()
    return(XC, YC, drF)


def read_bathymetry_to_transect(config_dir,model_name,L2_XC,L2_YC,transect_lon, transect_lat):

    bathy_file = os.path.join(config_dir,'L2',model_name,'input',model_name+'_bathymetry.bin')
    grid = np.fromfile(bathy_file, '>f4')
    grid = np.reshape(grid, np.shape(L2_XC))

    points = np.column_stack([L2_XC.ravel(), L2_YC.ravel()])
    values = grid.ravel()

    bathy_transect = -1*griddata(points,values,(transect_lon,transect_lat))

    return(bathy_transect)

def read_dv_dict_from_nc(config_dir,model_name,fjord_name):
    ds = nc4.Dataset(os.path.join(config_dir,'L2',model_name,'input','L2_dv_mask_reference_dict.nc'))
    grp = ds.groups[fjord_name]
    rows = grp.variables['source_rows'][:]
    cols = grp.variables['source_cols'][:]
    ds.close()
    return(rows, cols)


def great_circle_distance(lon_ref, lat_ref, Lon, Lat):
    earth_radius = 6371000
    lon_ref_radians = np.radians(lon_ref)
    lat_ref_radians = np.radians(lat_ref)
    lons_radians = np.radians(Lon)
    lats_radians = np.radians(Lat)
    lat_diff = lats_radians - lat_ref_radians
    lon_diff = lons_radians - lon_ref_radians
    d = np.sin(lat_diff * 0.5) ** 2 + np.cos(lat_ref_radians) * np.cos(lats_radians) * np.sin(lon_diff * 0.5) ** 2
    h = 2 * earth_radius * np.arcsin(np.sqrt(d))
    return(h)


def series_to_N_points(series,N):
    #find the total length of the series
    totalDistance=0
    for s in range(len(series[:,0])-1):
        totalDistance+=((series[s,0]-series[s+1,0])**2+(series[s,1]-series[s+1,1])**2)**0.5
    intervalDistance=totalDistance/(N-1)

    #make the list of points
    newSeries=series[0,:]
    currentS = 0
    currentPoint1=series[currentS,:]
    currentPoint2=series[currentS+1,:]
    for p in range(N-2):
        distanceAccrued = 0
        while distanceAccrued<intervalDistance:
            currentLineDistance=((currentPoint1[0]-currentPoint2[0])**2+(currentPoint1[1]-currentPoint2[1])**2)**0.5
            if currentLineDistance<intervalDistance-distanceAccrued:
                distanceAccrued+=currentLineDistance
                currentS+=1
                currentPoint1 = series[currentS, :]
                currentPoint2 = series[currentS + 1, :]
            else:
                distance=intervalDistance-distanceAccrued
                newX=currentPoint1[0]+(distance/currentLineDistance)*(currentPoint2[0]-currentPoint1[0])
                newY = currentPoint1[1] + (distance / currentLineDistance) * (currentPoint2[1] - currentPoint1[1])
                distanceAccrued=intervalDistance+1
                newSeries=np.vstack([newSeries,np.array([newX,newY])])
                currentPoint1=np.array([newX,newY])
    newSeries = np.vstack([newSeries, series[-1,:]])
    return(newSeries)


def subsample_XY_to_transect(L2_XC, L2_YC, rows, cols):
    longitude = np.ones((len(rows),))
    latitude = np.ones((len(rows),))
    for i in range(len(rows)):
        longitude[i] = L2_XC[rows[i], cols[i]]
        latitude[i] = L2_YC[rows[i], cols[i]]

    # distance = np.zeros((len(rows),))
    # for i in range(len(rows)-1):
    #     dist = great_circle_distance(longitude[i],latitude[i],longitude[i+1], latitude[i+1])
    #     distance[i+1] = distance[i]+dist

    return(longitude, latitude)#, distance)


def read_dv_output(config_dir, run_dir, model_name, N, Nr, fjord_name, var_name):

    file_names = []
    for file_name in os.listdir(os.path.join(config_dir,'L2',model_name,run_dir,'dv',fjord_name)):
        if fjord_name+'_mask_'+var_name in file_name and file_name[0]!='.':
            file_names.append(file_name)
    file_names.sort()

    counter = 0
    for file_name in file_names:
        print('        - Reading from file '+file_name)
        grid = np.fromfile(os.path.join(config_dir,'L2',model_name,run_dir,'dv',fjord_name,file_name), '>f4')
        print('                - min: '+str(np.min(grid)))
        print('                - max: ' + str(np.max(grid)))

        if var_name in ['THETA','SALT','UVEL','VVEL'] or var_name[:6]=='PTRACE':
            time_steps = int(np.size(grid)/(N*Nr))
            grid = np.reshape(grid, (time_steps, Nr, N))
        else:
            time_steps = int(np.size(grid)/N)
            grid = np.reshape(grid, (time_steps, N))

        first_iteration = int(file_name.split('.')[1])-1
        iter_step = (6*60*60)/30
        iterations = np.arange(first_iteration,first_iteration+iter_step*time_steps,iter_step)

        if counter==0:
            full_grid = grid
            full_iterations = iterations
            counter+=1
        else:
            full_grid = np.concatenate([full_grid, grid], axis=0)
            full_iterations = np.concatenate([full_iterations, iterations])
            # raise ValueError('Only implemented the first dv file')

    # find indices to remove overlap
    indices = np.ones((len(full_iterations),))
    counter = 0
    for i in range(len(indices)-1):
        if full_iterations[i] in full_iterations[i+1:].tolist():
            indices[i]=0
            counter += 1
    indices = indices.astype(bool)
    print('            -Removed '+str(counter)+' duplicate points')

    full_iterations = full_iterations[indices]
    full_grid = full_grid[indices, :, :]

    full_iterations = full_iterations[::4]
    full_grid = full_grid[::4,:,:]

    print(np.shape(full_iterations))
    print(np.shape(full_grid))

    return(full_grid, full_iterations)


def read_transect_points_from_shp(config_dir, model_name, transect_name):
    shp_file = os.path.join(config_dir,'L2',model_name,'input','dv_shp',transect_name)
    r = sf.Reader(shp_file)
    shapes = r.shapes()
    points = np.array(shapes[0].points)

    total_dist = 0
    for i in range(len(points)-1):
        total_dist += great_circle_distance(points[i,0], points[i,1],
                                            points[i+1,0], points[i+1, 1])
    N = int(total_dist/100)
    print('        - Subsampling the transect to '+str(N)+' points')

    points = series_to_N_points(points,N)

    XC_boundary = points[:,0]
    YC_boundary = points[:,1]
    distance = np.zeros((len(XC_boundary),))
    for i in range(len(YC_boundary) - 1):
        dist = great_circle_distance(XC_boundary[i], YC_boundary[i], XC_boundary[i + 1], YC_boundary[i + 1])
        distance[i + 1] = distance[i] + dist

    return(XC_boundary, YC_boundary, distance)


def interpolate_dv_grid_to_transect(L2_XC_dv, L2_YC_dv, dv_grid, transect_lon, transect_lat):

    transect_grid = np.zeros((np.shape(dv_grid)[0], np.shape(dv_grid)[1], len(transect_lat)))

    # plt.plot(L2_XC_dv, L2_YC_dv, 'k.')
    # plt.plot(transect_lon,transect_lat,'g.')
    # plt.show()

    points = np.column_stack([L2_XC_dv, L2_YC_dv])

    for timestep in range(np.shape(dv_grid)[0]):
        if timestep%10==0:
            print('            - Interpolating timestep '+str(timestep)+' of '+str(np.shape(dv_grid)[0]))
        for depth in range(np.shape(dv_grid)[1]):
            transect_grid[timestep, depth, :] = griddata(points, dv_grid[timestep,depth,:].ravel(), (transect_lon, transect_lat))
            # set_int = interp2d(L2_XC_dv, L2_YC_dv, dv_grid[timestep,depth,:], fill_value=0)
            # for i in range(len(transect_lat)):
            #     transect_grid[timestep, depth, i] = set_int(transect_lon[i], transect_lat[i])

    return(transect_grid)


def write_output_to_nc(config_dir, run_dir, model_name, fjord_name, var_name, drF,
                       longitude, latitude, distance, iterations, var_grid, transect_bathymetry):

    Z_bottom = np.cumsum(drF)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2

    output_file = os.path.join(config_dir,'L2',model_name,'results'+run_dir[3:],'dv',
                               model_name+'_'+fjord_name+'_'+var_name+'.nc')

    ds = nc4.Dataset(output_file, 'w')

    if var_name in ['THETA','SALT','UVEL','VVEL'] or var_name[:6]=='PTRACE':
        ds.createDimension('depth',len(drF))
    ds.createDimension('distance',len(distance))
    ds.createDimension('iterations', len(iterations))

    lon_var = ds.createVariable('longitude','f4',('distance',))
    lon_var[:] = longitude

    lat_var = ds.createVariable('latitude', 'f4', ('distance',))
    lat_var[:] = latitude

    d_var = ds.createVariable('distance', 'f4', ('distance',))
    d_var[:] = distance

    i_var = ds.createVariable('iterations', 'f4', ('iterations',))
    i_var[:] = iterations

    b_var = ds.createVariable('bathymetry', 'f4', ('distance',))
    b_var[:] = transect_bathymetry

    if var_name in ['THETA','SALT','UVEL','VVEL'] or var_name[:6]=='PTRACE':
        z_var = ds.createVariable('depth','f4',('depth',))
        z_var[:] = Z

        var = ds.createVariable(var_name,'f4',('iterations','depth','distance'))
        var[:, :, :] = var_grid

    else:
        var = ds.createVariable(var_name, 'f4', ('iterations', 'distance'))
        var[:, :] = var_grid

    ds.close()


def store_transect_to_nc(config_dir, run_dir):

    model_name = 'L2_Kanger'

    for var_name in ['THETA']:
        for fjord_name in ['Kanger_Trough']:
            print(' - Writing the file for '+fjord_name)

            print('        - Reading in the geometry and dv files')
            rows, cols = read_dv_dict_from_nc(config_dir,model_name, fjord_name)
            N = len(rows)

            L2_XC, L2_YC, drF = read_grid_geometry_from_nc(config_dir, model_name)
            Nr = len(drF)

            L2_XC_dv, L2_YC_dv = subsample_XY_to_transect(L2_XC, L2_YC, rows, cols)

            dv_grid, dv_iterations = read_dv_output(config_dir, run_dir, model_name, N, Nr, fjord_name, var_name)

            transect_lon, transect_lat, distance = read_transect_points_from_shp(config_dir, model_name, fjord_name)

            bathymetry_transect = read_bathymetry_to_transect(config_dir, model_name, L2_XC, L2_YC, transect_lon,
                                                                transect_lat)

            # plt.plot(bathymetry_transect)
            # plt.show()

            print('        - Interpolating onto the dense transect')
            transect_grid = interpolate_dv_grid_to_transect(L2_XC_dv, L2_YC_dv, dv_grid, transect_lon, transect_lat)

            indices = ~np.isnan(transect_grid[0,0,:])

            transect_lon = transect_lon[indices]
            transect_lat = transect_lat[indices]
            distance = distance[indices]
            bathymetry_transect = bathymetry_transect[indices]
            transect_grid = transect_grid[:,:,indices]

            # C = plt.pcolormesh(transect_grid[-1, :, :], cmap='turbo')
            # plt.colorbar(C)
            # plt.gca().invert_yaxis()
            # plt.show()

            print('        - Outputting to an nc file')
            write_output_to_nc(config_dir, run_dir, model_name, fjord_name, var_name, drF,
                               transect_lon, transect_lat, distance, dv_iterations, transect_grid, bathymetry_transect)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L2, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-r", "--run_dir", action="store",
                        help="The run dir corresponding to the experiment.", dest="run_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    run_dir = args.run_dir

    store_transect_to_nc(config_dir, run_dir)






