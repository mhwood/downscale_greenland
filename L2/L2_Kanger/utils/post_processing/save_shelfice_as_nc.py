
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from pyproj import Proj, Transformer


def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1):
    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))
    # inProj = Proj(init='epsg:'+str(inputCRS))
    # outProj = Proj(init='epsg:'+str(outputCRS))
    x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon


def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    ds.close()
    return(XC, YC)


def read_dv_dict_from_nc(config_dir,model_name):
    ds = nc4.Dataset(os.path.join(config_dir,'L2',model_name,'input','L2_dv_mask_reference_dict.nc'))
    grp = ds.groups['shelfice']
    rows = grp.variables['source_rows'][:]
    cols = grp.variables['source_cols'][:]
    ds.close()
    return(rows, cols)


def read_dv_dict_from_mask(config_dir,model_name,L2_XC, L2_YC):

    grid = np.fromfile(os.path.join(config_dir,'L2',model_name,'input','dv','melange_mask.bin'), '>f4')

    grid = np.reshape(grid, np.shape(L2_XC))

    N = int(np.max(grid))
    # print(N)

    rows = np.zeros((N,)).astype(int)
    cols = np.zeros((N,)).astype(int)

    for i in range(1,N+1):
        r, c = np.where(grid==i)
        rows[i-1] = r[0]
        cols[i-1] = c[0]

    # plt.plot(rows, cols)
    # plt.show()

    return(rows,cols)


def read_dv_output_from_nc(config_dir, run_dir, model_name, N, var_name):

    fw_flux_files = []
    ht_flux_files = []

    file_names = []
    for file_name in os.listdir(os.path.join(config_dir,'L2',model_name,run_dir,'dv','melange')):
        if var_name in file_name and file_name[0]!='.':
            file_names.append(file_name)

    # index = 0
    # grid = np.fromfile(os.path.join(config_dir,'L2',model_name,'run_shelfice','dv',file_names[index]), '>f4')
    #
    # time_steps = int(np.size(grid)/N)
    # grid = np.reshape(grid, (time_steps, N))
    #
    # file_names = []
    # for file_name in os.listdir(os.path.join(config_dir, 'L2', model_name, 'run_iceplume', 'dv', fjord_name + '_Fjd')):
    #     if fjord_name + '_Fjd_mask_' + var_name in file_name and file_name[0] != '.':
    #         file_names.append(file_name)

    file_names.sort()

    counter = 0
    for file_name in file_names:
        print('            - Reading from file '+file_name)
        grid = np.fromfile(
            os.path.join(config_dir, 'L2', model_name, run_dir, 'dv', 'melange', file_name), '>f4')

        time_steps = int(np.size(grid) / N)
        grid = np.reshape(grid, (time_steps, N))

        first_iteration = int(file_name.split('.')[1]) - 1
        iter_step = (24 * 60 * 60) / 60
        iterations = np.arange(first_iteration, first_iteration + iter_step * time_steps, iter_step)

        if counter == 0:
            full_grid = grid
            full_iterations = iterations
            counter += 1
        else:
            full_grid = np.concatenate([full_grid, grid], axis=0)
            full_iterations = np.concatenate([full_iterations, iterations])

    return(full_grid, full_iterations)



def write_output_to_nc(config_dir, results_dir, model_name, var_name, L2_XC, L2_YC, L2_X, L2_Y,
                       iterations, var_grid):

    output_file = os.path.join(config_dir,'L2',model_name,results_dir,'dv',
                               model_name+'_shelfice_'+var_name+'.nc')

    ds = nc4.Dataset(output_file, 'w')

    ds.createDimension('x', np.shape(L2_XC)[1])
    ds.createDimension('y', np.shape(L2_XC)[0])
    ds.createDimension('iterations', len(iterations))

    lon_var = ds.createVariable('Longitude', 'f4', ('y', 'x'))
    lon_var[:, :] = L2_XC

    lat_var = ds.createVariable('Latitude', 'f4', ('y', 'x'))
    lat_var[:, :] = L2_YC

    lon_var = ds.createVariable('X','f4',('y','x'))
    lon_var[:,:] = L2_X

    lat_var = ds.createVariable('Y', 'f4',('y','x'))
    lat_var[:,:] = L2_Y

    lon_var = ds.createVariable('x', 'f4', ('x',))
    lon_var[:] = L2_X[0,:]

    lat_var = ds.createVariable('y', 'f4', ('y', ))
    lat_var[:] = L2_Y[:,0]

    i_var = ds.createVariable('iterations', 'f4', ('iterations',))
    i_var[:] = iterations

    var = ds.createVariable(var_name, 'f4', ('iterations','y','x'))
    var[:, :, :] = var_grid

    ds.close()



def store_shelfice_to_nc(config_dir, run_dir):

    model_name = 'L2_Kanger'
    var_name = 'SHFFWFLX'
    results_dir = 'results_'+run_dir[4:]

    print(' - Storing the shelfice dv output to nc')

    L2_XC, L2_YC = read_grid_geometry_from_nc(config_dir, model_name)

    # rows, cols = read_dv_dict_from_nc(config_dir,model_name)
    rows, cols = read_dv_dict_from_mask(config_dir, model_name, L2_XC, L2_YC)
    N = len(rows)


    print('    - Reading the files from the run dir')
    dv_grid, iterations = read_dv_output_from_nc(config_dir, run_dir, model_name, N, var_name)

    points = np.column_stack([L2_XC.ravel(), L2_YC.ravel()])
    new_points = reproject_polygon(points, inputCRS=4326, outputCRS=3413)
    L2_X = new_points[:, 0].reshape(np.shape(L2_XC))
    L2_Y = new_points[:, 1].reshape(np.shape(L2_YC))

    output_grid = np.zeros((np.shape(dv_grid)[0], np.shape(L2_XC)[0], np.shape(L2_XC)[1]))
    for timestep in range(np.shape(dv_grid)[0]):
        for i in range(len(rows)):
            output_grid[timestep,rows[i],cols[i]] = dv_grid[timestep,i]

    min_rpw = np.min(rows)-2
    max_row = np.max(rows)+3
    min_col = np.min(cols)-2
    max_col = np.max(cols)+3

    L2_XC = L2_XC[min_rpw:max_row, min_col:max_col]
    L2_YC = L2_YC[min_rpw:max_row, min_col:max_col]
    L2_X = L2_X[min_rpw:max_row, min_col:max_col]
    L2_Y = L2_Y[min_rpw:max_row, min_col:max_col]
    output_grid = output_grid[:, min_rpw:max_row, min_col:max_col]

    C = plt.imshow(output_grid[-1, :, :], origin='lower')#, vmin=0, vmax=1)
    plt.colorbar(C)
    plt.show()

    write_output_to_nc(config_dir, results_dir, model_name, var_name, L2_XC, L2_YC, L2_X, L2_Y,
                       iterations, output_grid)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L2, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-r", "--run_dir", action="store",
                        help="The run dir for the experiment.", dest="run_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    run_dir = args.run_dir

    store_shelfice_to_nc(config_dir, run_dir)






