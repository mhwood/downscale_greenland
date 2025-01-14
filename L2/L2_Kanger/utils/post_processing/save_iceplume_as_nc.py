
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import shapefile
from scipy.interpolate import interp2d, griddata


def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    Depth = ds.variables['Depth'][:, :]
    drF = ds.variables['drF'][:]
    ds.close()
    return(XC, YC,Depth, drF)



def read_glaciers_location_dict_from_shp(folder,model_name):
    sf = shapefile.Reader(os.path.join(folder,'L2',model_name,'input','domain_shp',model_name+'_iceplume_cells'))
    records = sf.records()

    loc_dict = {}
    for record in records:
        if record[4] != '':
            loc_dict[(record[0], record[1])] = record[4]

    return(loc_dict)


def read_iceplume_mask(config_dir, model_name, rows, cols):
    file_path = os.path.join(config_dir, 'L2', model_name,'input','iceplume', 'L2_iceplume_mask.bin')
    grid = np.fromfile(file_path,'>f4')
    grid = np.reshape(grid, (rows, cols))
    return(grid)

def read_melange_draft(config_dir, model_name, rows, cols):
    file_path = os.path.join(config_dir, 'L2', model_name,'input','shelfice', 'melange_draft.bin')
    grid = np.fromfile(file_path,'>f4')
    grid = np.reshape(grid, (rows, cols))
    return(grid)


def read_dv_output(config_dir, run_dir, model_name, N, Nr, var_name):

    file_names = []
    for file_name in os.listdir(os.path.join(config_dir,'L2',model_name,run_dir,'dv','iceplume')):
        if 'iceplume_mask_'+var_name in file_name and file_name[0]!='.':
            file_names.append(file_name)
    file_names.sort()

    counter = 0
    for file_name in file_names:
        print('        - Reading from file '+file_name)
        grid = np.fromfile(os.path.join(config_dir,'L2',model_name,run_dir,'dv','iceplume',file_name), '>f4')

        time_steps = int(np.size(grid)/(N*Nr))
        grid = np.reshape(grid, (time_steps, Nr, N))

        first_iteration = int(file_name.split('.')[1])-1
        iter_step = (60*60)/60
        iterations = np.arange(first_iteration,first_iteration+iter_step*time_steps,iter_step)

        if counter==0:
            full_grid = grid
            full_iterations = iterations
            counter+=1
        else:
            full_grid = np.concatenate([full_grid, grid], axis=0)
            full_iterations = np.concatenate([full_iterations, iterations])
            # raise ValueError('Only implemented the first dv file')

    full_iterations = full_iterations[::4]
    full_grid = full_grid[::4,:,:]

    print(np.shape(full_iterations))
    print(np.shape(full_grid))

    return(full_grid, full_iterations)


def write_output_to_nc(config_dir, results_dir, model_name, glacier_name, var_names, drF, iterations, output_grids,max_depth,melange_draft):

    Z_bottom = np.cumsum(drF)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2

    output_file = os.path.join(config_dir,'L2',model_name,results_dir,'dv',
                               model_name+'_iceplume_'+glacier_name+'.nc')

    ds = nc4.Dataset(output_file, 'w')

    ds.createDimension('depth',len(drF))
    ds.createDimension('iterations', len(iterations))


    i_var = ds.createVariable('iterations', 'f4', ('iterations',))
    i_var[:] = iterations

    b_var = ds.createVariable('depth', 'f4', ('depth',))
    b_var[:] = Z

    for i in range(len(var_names)):
        var = ds.createVariable(var_names[i],'f4',('iterations','depth'))
        var[:, :] = output_grids[i]

    ds.max_depth = max_depth
    ds.melange_draft = melange_draft

    ds.close()


def store_iceplume_to_nc(config_dir, run_dir):

    model_name = 'L2_Kanger'
    results_dir = 'results'+run_dir[3:]

    print('        - Reading in the geometry and dv files')
    glacier_loc_dict = read_glaciers_location_dict_from_shp(config_dir,model_name)
    glaciers = list(glacier_loc_dict.keys())
    print(glacier_loc_dict)

    L2_XC, L2_YC, Depth, drF = read_grid_geometry_from_nc(config_dir, model_name)
    Nr = len(drF)

    iceplume_mask = read_iceplume_mask(config_dir, model_name, np.shape(L2_XC)[0], np.shape(L2_XC)[1])

    if 'melange' in run_dir:
        melange_draft = read_melange_draft(config_dir, model_name, np.shape(L2_XC)[0], np.shape(L2_XC)[1])

    plume_rows, plume_cols = np.where(np.abs(iceplume_mask)==6)
    ice_indices = np.where(np.abs(iceplume_mask) != 0)
    ice_values = iceplume_mask[ice_indices]
    plume_indices = np.where(np.abs(ice_values) == 6)[0]

    N = len(ice_values)
    print(len(plume_indices))
    print(len(glaciers))

    # for i in range(len(plume_cols)):
    #     print(L2_XC[plume_rows[i],plume_cols[i]], L2_YC[plume_rows[i],plume_cols[i]])

    for i in range(len(plume_cols)):
        row = plume_rows[i]
        col = plume_cols[i]
        glacier = glacier_loc_dict[(row,col)]
        glacier_index = plume_indices[i]
        print(' - Writing the file for ' + glacier)

        max_depth = Depth[row, col]
        if 'melange' in run_dir:
            draft = melange_draft[row,col]
        else:
            draft=0

        output_grids = []

        var_names = ['ICEFRNTA','ICEFRNTM','ICEFRNTR','ICEFRNTS','ICEFRNTT','ICEFRNTW']
        var_names = ['ICEFRNTM','ICEFRNTA']
        for var_name in var_names:
            print('    - Reading grid for '+var_name)
            dv_grid, dv_iterations = read_dv_output(config_dir, run_dir, model_name, N, Nr, var_name)

            glacier_dv_grid = dv_grid[:,:,glacier_index]

            output_grids.append(glacier_dv_grid)

            # for i in range(np.shape(dv_grid)[2]):
            #     C = plt.imshow(glacier_dv_grid, cmap='turbo')
            #     plt.colorbar(C)
            #     plt.title(str(i))
            #     plt.show()

        print('        - Outputting to an nc file')
        write_output_to_nc(config_dir, results_dir, model_name, glacier, var_names, drF, dv_iterations, output_grids, max_depth,draft)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L2, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-r", "--run_dir", action="store",
                        help="The name of the run dir.", dest="run_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    run_dir = args.run_dir

    store_iceplume_to_nc(config_dir, run_dir)






