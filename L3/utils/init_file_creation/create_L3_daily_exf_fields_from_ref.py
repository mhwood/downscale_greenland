
import os
#import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import netCDF4 as nc4
import argparse
import ast
import sys

def create_src_dest_dicts_from_ref(config_dir, L2_model_name, var_name, start_year, final_year, start_month, final_month, start_day, final_day):
    dest_files = []
    start_date = datetime(start_year, start_month, start_day)
    final_date = datetime(final_year, final_month, final_day)
    for year in range(2002, 2003):
        for month in range(1, 13):
            if month in [1, 3, 5, 7, 8, 10, 12]:
                nDays = 31
            elif month in [4, 6, 9, 11]:
                nDays = 30
            else:
                if year % 4 == 0:
                    nDays = 29
                else:
                    nDays = 28
            for day in range(1, nDays + 1):
                test_date = datetime(year, month, day)
                if test_date >= start_date and test_date <= final_date:
                    dest_files.append(str(year) + '{:02d}'.format(month) + '{:02d}'.format(day))

    f = open(os.path.join(config_dir,'L2', L2_model_name, 'run','dv', 'exf_dest_ref.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)

    dest_files_out = []
    source_file_read_dict = {}
    source_file_read_index_sets = {}

    for dest_file in dest_files:
        dest_files_out.append('L3_exf_'+var_name+'.'+dest_file+'.bin')
        source_files = []
        index_set = []
        for pair in size_dict[dest_file]:
            source_files.append('surface_'+var_name+'.'+pair[0]+'.bin')
            index_set.append(pair[1])
        source_file_read_dict['L3_exf_'+var_name+'.'+dest_file+'.bin'] = source_files
        source_file_read_index_sets['L3_exf_'+var_name+'.'+dest_file+'.bin'] = index_set

    return(dest_files_out, source_file_read_dict, source_file_read_index_sets)

def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    delR = ds.variables['drF'][:]
    ds.close()
    return(XC, YC, delR)

def read_wetgrid_from_nc(config_dir, model_name, hFac):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    hFac_grid = ds.variables['HFac'+hFac][:,:]
    ds.close()

    wet_grid = np.copy(hFac_grid)
    wet_grid[wet_grid>0] = 1

    if hFac=='S':
        wet_grid = wet_grid[:,:-1,:]
    if hFac=='W':
        wet_grid = wet_grid[:,:,:-1]

    return(wet_grid)

def read_mask_reference_from_nc_dict(nc_dict_file,mask_name):
    ds = nc4.Dataset(nc_dict_file)

    grp = ds.groups[mask_name]
    source_rows = grp.variables['source_rows'][:]
    source_cols = grp.variables['source_cols'][:]

    ds.close()
    return(source_rows,source_cols)

def read_exf_variable_to_L2_points(L2_run_dir, var_name,
                                        source_files, source_file_read_indices,
                                        source_rows, source_cols,
                                        L2_XC, L2_YC, L2_wet_grid, print_level):

    n_Timesteps = 0
    for index_set in source_file_read_indices:
        n_Timesteps += int(index_set[1]-index_set[0])+1

    if print_level >= 3:
        print('            - the L2 grid for this file has '+str(n_Timesteps)+' timesteps')

    points = np.zeros((len(source_rows),2))
    values = np.zeros((n_Timesteps, 1, len(source_rows)))
    wet_grid = np.zeros((1,len(source_rows)))

    for i in range(len(source_rows)):
        x = L2_XC[source_rows[i], source_cols[i]]
        y = L2_YC[source_rows[i], source_cols[i]]
        points[i,0] = x
        points[i,1] = y
        wet_grid[0, i] = L2_wet_grid[0,source_rows[i], source_cols[i]]

    index_counter = 0
    for s in range(len(source_files)):
        source_file = source_files[s]
        file_suffix = '.'.join(source_file.split('.')[-2:])
        index_set = source_file_read_indices[s]
        if print_level >= 4:
            print('                - Adding timesteps ' + str(index_set[0]) + ' to ' + str(index_set[1]+1) + ' from file ' + source_file)

        start_file_index = int(index_set[0])
        end_file_index = int(index_set[1])

        if print_level >= 4:
            print('                - Storing at points '+str(index_counter)+' to '+
              str(index_counter+(end_file_index-start_file_index)+1)+' in the grid')

        N = len(source_rows)

        var_file = os.path.join(L2_run_dir, 'dv', 'L3_surface_mask_' + var_name + '.' + file_suffix)
        var_grid = np.fromfile(var_file, dtype='>f4')
        timesteps_in_file = int(np.size(var_grid) / (N))
        var_grid = np.reshape(var_grid, (timesteps_in_file, N))

        values[index_counter:index_counter + (end_file_index-start_file_index)+1,0, :] = var_grid[start_file_index:end_file_index+1,:]

        index_counter += (end_file_index-start_file_index)

    return(points,values,wet_grid)

def downscale_L2_exf_field_to_L3(df, surface_var_name,
                                 L2_points, L2_values, L2_wet_grid_points,
                                 L3_XC, L3_YC, L3_wet_grid, print_level):

    n_Timesteps = np.shape(L2_values)[0]
    L3_surface_var = np.zeros((n_Timesteps,np.shape(L3_XC)[0],np.shape(L3_YC)[1]))

    if print_level>=4:
        print('            -  Variable shapes entering downscale routine:')
        print('                -  L2_points: '+str(np.shape(L2_points)))
        print('                -  L2_values: ' + str(np.shape(L2_values)))
        print('                -  L2_wet_grid_points: ' + str(np.shape(L2_wet_grid_points)))
        print('                -  L3_XC: ' + str(np.shape(L3_XC)))
        print('                -  L3_YC: ' + str(np.shape(L3_YC)))
        print('                -  L3_wet_grid: ' + str(np.shape(L3_wet_grid)))

    # if surface_var_name in ['RUNOFF', 'PRECIP','SWDOWN']:
    #     spread_horizontally = False
    # else:
    #     spread_horizontally = True
    #
    # for timestep in range(n_Timesteps):
    #     if timestep%50==0:
    #         print('          Working on timestep '+str(timestep)+' of '+str(n_Timesteps))
    #     # if timestep==30:
    #     #     C = plt.imshow(L2_surface_var[timestep,:,:],origin='lower')
    #     #     plt.colorbar(C)
    #     #     plt.show()
    #     downscaled_field = df.downscale_2D_field(L2_XC_e, L2_YC_e, L2_surface_var[timestep,:,:],
    #                                              L2_wet_grid, L2_wet_grid_on_L3,
    #                                              XC_subset, YC_subset, L3_wet_grid,spread_horizontally)
    #
    #     L3_surface_var[timestep,:,:] = downscaled_field

    if surface_var_name in ['SWDOWN', 'PRECIP', 'RUNOFF']:
        remove_zeros = False
    else:
        remove_zeros = True

    for timestep in range(n_Timesteps):
        # if timestep%5==0:
        if print_level >= 5:
            print('                    -  Working on timestep '+str(timestep+1)+' of '+str(n_Timesteps))

        downscaled_field = df.downscale_3D_points(L2_points, L2_values[timestep, :, :],
                                                  L2_wet_grid_points,
                                                  L3_XC, L3_YC, L3_wet_grid,
                                                  remove_zeros=remove_zeros)

        L3_surface_var[timestep, :, :] = downscaled_field

    return(L3_surface_var)

def output_exf_variable(output_dir, dest_file, var_grid):

    output_file = os.path.join(output_dir,dest_file)
    var_grid.ravel('C').astype('>f4').tofile(output_file)

def create_exf_field(config_dir,L3_model_name, L2_model_name, var_name,
                     dest_files, source_file_read_dict, source_file_read_index_sets, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df

    if print_level >= 1:
        print('    - Creating the '+var_name+' external forcing files for the ' + L3_model_name + ' model from L2 output')

    # step 0: get the model domains
    L2_XC, L2_YC, L2_delR = read_grid_geometry_from_nc(config_dir, L2_model_name)
    L3_XC, L3_YC, L3_delR = read_grid_geometry_from_nc(config_dir, L3_model_name)

    # this is the dir where the exf output will be stored
    if 'exf' not in os.listdir(os.path.join(config_dir,'L3',L3_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L3',L3_model_name,'input','exf'))
    if var_name not in os.listdir(os.path.join(config_dir,'L3',L3_model_name,'input','exf')):
        os.mkdir(os.path.join(config_dir,'L3',L3_model_name,'input','exf',var_name))
    output_dir = os.path.join(config_dir,'L3',L3_model_name,'input','exf',var_name)

    if var_name in ['VWIND']:
        L2_wet_cells = read_wetgrid_from_nc(config_dir, L2_model_name, hFac='S')
        # L2_wet_cells_on_L3 = read_mask_from_nc(os.path.join('..', 'input', 'L2_wetgrid_on_L3.nc'), hFac='S')
        L3_wet_cells = read_wetgrid_from_nc(config_dir, L3_model_name, hFac='S')
    elif var_name in ['UWIND']:
        L2_wet_cells = read_wetgrid_from_nc(config_dir, L2_model_name, hFac='W')
        # L2_wet_cells_on_L3 = read_mask_from_nc(os.path.join('..', 'input', 'L2_wetgrid_on_L3.nc'), hFac='W')
        L3_wet_cells = read_wetgrid_from_nc(config_dir, L3_model_name, hFac='W')
    else:
        L2_wet_cells = read_wetgrid_from_nc(config_dir, L2_model_name, hFac='C')
        # L2_wet_cells_on_L3 = read_mask_from_nc(os.path.join('..', 'input', 'L2_wetgrid_on_L3.nc'), hFac='C')
        L3_wet_cells = read_wetgrid_from_nc(config_dir, L3_model_name, hFac='C')

    # L3_wet_cells = L3_wet_cells[0, :, :]

    # plt.imshow(L3_wet_cells)
    # plt.show()

    # # plt.subplot(1, 3, 1)
    # # C = plt.imshow(L2_wet_grid, origin='lower')
    # # plt.colorbar(C)
    # # plt.subplot(1, 3, 2)
    # # C = plt.imshow(L2_wet_grid_on_L3, origin='lower')
    # # plt.colorbar(C)
    # # plt.subplot(1, 3, 3)
    # # C = plt.imshow(L3_wet_grid, origin='lower')
    # # plt.colorbar(C)
    # # plt.show()
    if print_level >= 1:
        print('    - Reading the mask to reference the variable to the L2 grid')
    nc_dict_file = os.path.join(config_dir, 'L2', L2_model_name, 'namelist', 'L2_dv_mask_reference_dict.nc')
    source_rows, source_cols = read_mask_reference_from_nc_dict(nc_dict_file, 'surface')

    for dest_file in dest_files:
        if dest_file not in os.listdir(output_dir):
            if print_level >= 2:
                print('        - Downscaling the timesteps to be stored in file '+str(dest_file))
            source_files = source_file_read_dict[dest_file]
            source_file_read_indices = source_file_read_index_sets[dest_file]

            # L2_surface_var = read_surface_variable_onto_L2_grid(L2_run_dir, var_name, source_files, source_file_read_indices, n_rows_L3_e,n_cols_L3_e)

            L2_run_dir = os.path.join(config_dir, 'L2', L2_model_name, 'run')
            L2_points, L2_values, L2_wet_grid_points = read_exf_variable_to_L2_points(L2_run_dir, var_name,
                                                                                      source_files,
                                                                                      source_file_read_indices,
                                                                                      source_rows, source_cols,
                                                                                      L2_XC, L2_YC, L2_wet_cells, print_level)

            # plt.plot(L2_points[:,0],L2_points[:,1],'k.')
            # plt.plot(L3_XC[:, 0], L3_YC[:, 0], 'b-')
            # plt.plot(L3_XC[:, -1], L3_YC[:, -1], 'b-')
            # plt.plot(L3_XC[0, :], L3_YC[0, :], 'b-')
            # plt.plot(L3_XC[-1, :], L3_YC[-1, :], 'b-')
            # plt.show()

            L3_surface_var = downscale_L2_exf_field_to_L3(df, var_name,
                                                          L2_points, L2_values, L2_wet_grid_points,
                                                          L3_XC, L3_YC, L3_wet_cells, print_level)


            # plot_var = L3_surface_var[0,:,:]
            # plot_var = np.ma.masked_where(plot_var==0,plot_var)
            # C = plt.imshow(plot_var, origin='lower')
            # plt.colorbar(C)
            # plt.title(var_name)
            # plt.show()

            output_exf_variable(output_dir, dest_file, L3_surface_var)
        else:
            if print_level >= 2:
                print('        - Skipping ' + str(dest_file) + ' because it was already created')

def create_exf_fields_via_interpolation(config_dir, L3_model_name, L2_model_name, proc_id,
                                        start_year, final_year, start_month, final_month, start_day, final_day, print_level):

    var_name_list = ['UWIND','VWIND','ATEMP','AQH','PRECIP','SWDOWN','LWDOWN','RUNOFF']
    var_name = var_name_list[proc_id%10]

    if print_level>=1:
        print('    - Creating the exf field for '+var_name+' to cover days ' +
              str(start_year)+'/'+str(start_month)+'/'+str(start_day) + ' to ' + str(final_year)+'/'+str(final_month)+'/'+str(final_day))

    dest_files, source_file_read_dict, source_file_read_index_sets = \
        create_src_dest_dicts_from_ref(config_dir, L2_model_name, var_name, start_year, final_year, start_month, final_month, start_day, final_day)

    if print_level >= 1:
        print('    - Running the Downscale routine:')
    create_exf_field(config_dir, L3_model_name, L2_model_name, var_name,
                     dest_files, source_file_read_dict, source_file_read_index_sets, print_level)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L3, L3, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-i", "--proc_id", action="store",
                        help="The id of the process to run.", dest="proc_id",
                        type=int, required=True)

    parser.add_argument("-S", "--start_year", action="store",
                        help="The start year.", dest="start_year",
                        type=int, required=True)

    parser.add_argument("-s", "--start_month", action="store",
                        help="The start month.", dest="start_month",
                        type=int, required=True)

    parser.add_argument("-sd", "--start_day", action="store",
                        help="The start day.", dest="start_day",
                        type=int, required=True)

    parser.add_argument("-F", "--final_year", action="store",
                        help="The final year.", dest="final_year",
                        type=int, required=True)

    parser.add_argument("-f", "--final_month", action="store",
                        help="The final month.", dest="final_month",
                        type=int, required=True)

    parser.add_argument("-fd", "--final_day", action="store",
                        help="The final day.", dest="final_day",
                        type=int, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    proc_id = args.proc_id
    start_year = args.start_year
    final_year = args.final_year
    start_month = args.start_month
    final_month = args.final_month
    start_day = args.start_day
    final_day = args.final_day

    create_exf_fields_via_interpolation(config_dir, proc_id, start_year, final_year, start_month, final_month, start_day, final_day)
