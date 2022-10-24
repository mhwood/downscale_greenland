
import os
#import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import netCDF4 as nc4
import argparse
import ast
import sys

def create_src_dest_dicts_from_ref(config_dir, L1_model_name, mask_name,var_name, start_year, final_year, start_month, final_month, start_day, final_day):
    dest_files = []
    start_date = datetime(start_year, start_month, start_day)
    final_date = datetime(final_year, final_month, final_day)
    for year in range(1992,2022):
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

    f = open(os.path.join(config_dir,'L1_grid', L1_model_name, 'run','dv', 'L3_BC_dest_ref.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)

    dest_files_out = []
    source_file_read_dict = {}
    source_file_read_index_sets = {}

    if var_name in ['UVEL','VVEL','UICE','VICE']:
        suffix = '_rotated.bin'
    else:
        suffix = '.bin'

    for dest_file in dest_files:
        dest_files_out.append('L3_BC_'+mask_name+'_'+var_name+'.'+dest_file+suffix)
        source_files = []
        index_set = []
        for pair in size_dict[dest_file]:
            source_files.append(mask_name+'_'+var_name+'.'+pair[0]+'.bin')
            index_set.append(pair[1])
        source_file_read_dict['L3_BC_'+mask_name+'_'+var_name+'.'+dest_file+suffix] = source_files
        source_file_read_index_sets['L3_BC_'+mask_name+'_'+var_name+'.'+dest_file+suffix] = index_set

    return(dest_files_out, source_file_read_dict, source_file_read_index_sets)


def read_grid_geometry_from_nc(config_dir,model_name,var_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    AngleSN = ds.variables['AngleSN'][:, :]
    if var_name in ['UVEL','UICE']:
        hfac = ds.variables['HFacW'][:,:,:]
        hfac = hfac[:,:,:-1]
    elif var_name in ['VVEL','VICE']:
        hfac = ds.variables['HFacS'][:,:,:]
        hfac = hfac[:, :-1, :]
    else:
        hfac = ds.variables['HFacC'][:, :, :]
    delR = np.array(ds.variables['drF'][:])
    ds.close()
    mask = np.copy(hfac)
    mask[mask>0]=1
    return(XC,YC,AngleCS,AngleSN,mask,delR)


def read_mask_reference_from_nc_dict(nc_dict_file,mask_name):
    ds = nc4.Dataset(nc_dict_file)

    grp = ds.groups['L3_'+mask_name]
    source_rows = grp.variables['source_rows'][:]
    source_cols = grp.variables['source_cols'][:]

    ds.close()
    return(source_rows,source_cols)


def read_boundary_variable_to_L1_points(L1_run_dir, mask_name, var_name,
                                        source_files, source_file_read_indices,
                                        source_rows, source_cols,
                                        Nr_in, L1_XC, L1_YC,
                                        L1_AngleCS, L1_AngleSN, L1_wet_grid_3D_subset, print_level):

    n_Timesteps = 0
    for index_set in source_file_read_indices:
        n_Timesteps += int(index_set[1]-index_set[0])+1
    # n_Timesteps = 1

    if print_level>=5:
        print('                    - the L1 grid for this file has '+str(n_Timesteps)+' timesteps')

    points = np.zeros((len(source_rows),2))
    angles = np.zeros((len(source_rows),2))
    if var_name in ['ETAN','UICE','VICE','HSNOW','HEFF','AREA']:
        values = np.zeros((n_Timesteps, 1, len(source_rows)))
        wet_grid = np.zeros((1, len(source_rows)))
    else:
        values = np.zeros((n_Timesteps, Nr_in, len(source_rows)))
        wet_grid = np.zeros((Nr_in, len(source_rows)))

    for i in range(len(source_rows)):
        x = L1_XC[source_rows[i], source_cols[i]]
        y = L1_YC[source_rows[i], source_cols[i]]
        points[i,0] = x
        points[i,1] = y
        angles[i,0] = L1_AngleCS[source_rows[i], source_cols[i]]
        angles[i,1] = L1_AngleSN[source_rows[i], source_cols[i]]
        wet_column = L1_wet_grid_3D_subset[:, source_rows[i], source_cols[i]]
        wet_column[wet_column>0]=1
        if var_name in ['ETAN','UICE','VICE','HSNOW','HEFF','AREA']:
            wet_grid[0, i] = wet_column[0]
        else:
            wet_grid[:, i] = wet_column

    index_counter = 0
    for s in range(len(source_files)):
        source_file = source_files[s]
        file_suffix = '.'.join(source_file.split('.')[-2:])
        index_set = source_file_read_indices[s]
        if print_level>=3:
            print('        - Adding timesteps ' + str(index_set[0]) + ' to ' + str(index_set[1]) + ' from file ' + source_file)

        start_file_index = int(index_set[0])
        end_file_index = int(index_set[1])

        if print_level >= 4:
            print('                - Storing at points '+str(index_counter)+' to '+
              str(index_counter+(end_file_index-start_file_index))+' in the grid')

        N = len(source_rows)

        if var_name in ['UVEL','VVEL']:
            u_var_file = os.path.join(L1_run_dir,'dv', 'L3_'+mask_name, 'L3_'+mask_name + '_BC_mask_UVEL.'+file_suffix)
            u_var_grid = np.fromfile(u_var_file, dtype='>f4')
            timesteps_in_file = int(np.size(u_var_grid) / (Nr_in * N))
            u_var_grid = np.reshape(u_var_grid, (timesteps_in_file, Nr_in, N))

            v_var_file = os.path.join(L1_run_dir, 'dv', 'L3_'+mask_name, 'L3_' + mask_name + '_BC_mask_VVEL.' + file_suffix)
            v_var_grid = np.fromfile(v_var_file, dtype='>f4')
            v_var_grid = np.reshape(v_var_grid, (timesteps_in_file, Nr_in, N))

            var_grid = np.zeros_like(u_var_grid)
            for timestep in range(timesteps_in_file):
                for k in range(Nr_in):
                    if var_name=='UVEL':
                        var_grid[timestep, k, :] = angles[:,0] * u_var_grid[timestep, k, :] - angles[:,1] * v_var_grid[timestep, k, :]
                    if var_name=='VVEL':
                        var_grid[timestep, k, :] = angles[:,1] * u_var_grid[timestep, k, :] + angles[:,0] * v_var_grid[timestep, k, :]

        elif var_name in ['UICE','VICE']:
            u_var_file = os.path.join(L1_run_dir, 'dv', 'L3_'+mask_name, 'L3_' + mask_name + '_BC_mask_UICE.' + file_suffix)
            u_var_grid = np.fromfile(u_var_file, dtype='>f4')
            timesteps_in_file = int(np.size(u_var_grid) / (N))
            u_var_grid = np.reshape(u_var_grid, (timesteps_in_file, 1, N))

            v_var_file = os.path.join(L1_run_dir, 'dv', 'L3_'+mask_name, 'L3_' + mask_name + '_BC_mask_VICE.' + file_suffix)
            v_var_grid = np.fromfile(v_var_file, dtype='>f4')
            v_var_grid = np.reshape(v_var_grid, (timesteps_in_file, 1, N))

            var_grid = np.zeros_like(u_var_grid)
            for timestep in range(timesteps_in_file):
                if var_name == 'UICE':
                    var_grid[timestep, 0, :] = angles[:, 0] * u_var_grid[timestep, 0, :] - angles[:, 1] * v_var_grid[timestep, 0, :]
                if var_name == 'VICE':
                    var_grid[timestep, 0, :] = angles[:, 1] * u_var_grid[timestep, 0, :] + angles[:, 0] * v_var_grid[timestep, 0, :]
        else:
            # var_file = os.path.join(L1_run_dir, 'dv', mask_name + '_i' + str(layer), var_name,
            #                         mask_name + '_i' + str(layer) + '_mask_' + var_name + '.' + file_suffix)
            var_file = os.path.join(L1_run_dir,'dv','L3_'+mask_name, 'L3_'+mask_name +'_BC_mask_' + var_name +'.'+ file_suffix)
            var_grid = np.fromfile(var_file, dtype='>f4')
            if var_name in ['ETAN','UICE','VICE','HSNOW','HEFF','AREA']:
                timesteps_in_file = int(np.size(var_grid) / (N))
                var_grid = np.reshape(var_grid, (timesteps_in_file, 1, N))
            else:
                timesteps_in_file = int(np.size(var_grid) / (Nr_in * N))
                var_grid = np.reshape(var_grid, (timesteps_in_file, Nr_in, N))

        if var_name in ['ETAN','UICE','VICE','HSNOW','HEFF','AREA']:
            values[index_counter:index_counter + (end_file_index-start_file_index)+1,0, :] = var_grid[start_file_index:end_file_index+1,0,:]
        else:
            values[index_counter:index_counter + (end_file_index-start_file_index)+1,:, :] = var_grid[start_file_index:end_file_index+1,:,:]

        index_counter += (end_file_index-start_file_index)

    # plt.imshow(var_grid_out[0,0,:,:],origin='lower')
    # plt.show()

    return(points,values,wet_grid)


def downscale_L1_bc_field_to_L3(df, mask_name, var_name,
                                L1_points, L1_values,
                                L1_wet_grid_3D_points, L1_wet_grid_on_L3_3D_points,
                                L3_XC_subset, L3_YC_subset, L3_wet_grid_subset,print_level):

    nTimestepsOut = np.shape(L1_values)[0]

    if mask_name=='south' or mask_name=='north':
        if var_name in ['ETAN','UICE','VICE','HSNOW','HEFF','AREA']:
            L3_boundary_var = np.zeros((nTimestepsOut, 1, np.shape(L3_YC_subset)[0], np.shape(L3_YC_subset)[1]))
        else:
            L3_boundary_var = np.zeros((nTimestepsOut, np.shape(L3_wet_grid_subset)[0], np.shape(L3_YC_subset)[0], np.shape(L3_YC_subset)[1]))

    if mask_name=='west' or mask_name=='east':
        if var_name in ['ETAN','UICE','VICE','HSNOW','HEFF','AREA']:
            L3_boundary_var = np.zeros((nTimestepsOut, 1, np.shape(L3_YC_subset)[0], np.shape(L3_YC_subset)[1]))
        else:
            L3_boundary_var = np.zeros((nTimestepsOut, np.shape(L3_wet_grid_subset)[0], np.shape(L3_YC_subset)[0], np.shape(L3_YC_subset)[1]))

    if print_level>=5:
        print('                -  Variable shapes entering downscale routine:')
        print('                    -  L1_point: '+str(np.shape(L1_points)))
        print('                    -  L1_values: ' + str(np.shape(L1_values)))
        print('                    -  L1_wet_grid: ' + str(np.shape(L1_wet_grid_3D_points)))
        # print('          -  L1_wet_grid_on_L3_subset: ' + str(np.shape(L1_wet_grid_on_L3_3D_points)))
        print('                    -  L3_XC_subset: ' + str(np.shape(L3_XC_subset)))
        print('                    -  L3_YC_subset: ' + str(np.shape(L3_YC_subset)))
        print('                    -  L3_wet_grid_subset: ' + str(np.shape(L3_wet_grid_subset)))

    for timestep in range(nTimestepsOut):
        if timestep%50==0:
            if print_level >= 5:
                print('                  - Working on timestep '+str(timestep)+' of '+str(nTimestepsOut)+' for '+var_name+' on mask '+mask_name)
                print('                  - The bounds at this timestep are '+str(np.min(L1_values[timestep, :, :]))+' to '+str(np.max(L1_values[timestep, :, :])))

        if var_name in ['UVEL','VVEL','THETA','SALT','UICE','VICE']:
            downscaled_field = df.downscale_3D_points(L1_points, L1_values[timestep, :, :],
                                                       L1_wet_grid_3D_points,
                                                       L3_XC_subset, L3_YC_subset, L3_wet_grid_subset,
                                                       remove_zeros=True, printing = False)
        else:
            downscaled_field = df.downscale_3D_points_with_zeros(L1_points, L1_values[timestep, :, :],
                                                                 L1_wet_grid_3D_points,
                                                                 L3_XC_subset, L3_YC_subset, L3_wet_grid_subset,
                                                                 printing=False)
        L3_boundary_var[timestep, :, :, :] = downscaled_field

    # if var_name == 'ETAN':
    #     plt.subplot(2,1,1)
    #     C = plt.imshow(L0_boundary_var_subset[0,:,:])
    #     # plt.colorbar(C)
    #     plt.title(var_name)
    #     plt.subplot(2,1,2)
    #     plt.imshow(L1_boundary_var[0,:,:])
    #     plt.show()
    # else:
    #     plt.subplot(2, 1, 1)
    #     C = plt.imshow(L0_boundary_var_subset[0, 0, :, :])
    #     # plt.colorbar(C)
    #     plt.title(var_name)
    #     plt.subplot(2, 1, 2)
    #     plt.imshow(L1_boundary_var[0, 0, :, :])
    #     plt.show()

    return(L3_boundary_var)


########################################################################################################################


def create_bc_field(config_dir, L1_model_name, L3_model_name, mask_name,
                    boundary_var_name, dest_files, source_file_read_dict, source_file_read_index_sets, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import downscale_functions as df

    # this is the dir where the obcs output will be stored
    if 'obcs' not in os.listdir(os.path.join(config_dir,'L3',L3_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L3',L3_model_name,'input','obcs'))
    output_dir = os.path.join(config_dir,'L3',L3_model_name,'input','obcs')

    if mask_name not in os.listdir(output_dir):
        os.mkdir(os.path.join(output_dir,mask_name))
    if boundary_var_name not in os.listdir(os.path.join(output_dir,mask_name)):
        os.mkdir(os.path.join(output_dir,mask_name,boundary_var_name))

    print('    - Reading in the geometry of the L1 domain')

    L1_XC, L1_YC, L1_AngleCS, L1_AngleSN, L1_wet_grid_3D, delR_in = read_grid_geometry_from_nc(config_dir, L1_model_name, boundary_var_name)
    Nr_in = len(delR_in)

    L3_XC, L3_YC, _, _, L3_wet_grid_3D, delR_out = read_grid_geometry_from_nc(config_dir,L3_model_name,boundary_var_name)
    Nr_out = len(delR_out)

    if mask_name=='north':
        L3_XC_subset = L3_XC[-1:, :]
        L3_YC_subset = L3_YC[-1:, :]
        L3_wet_grid_3D_subset = L3_wet_grid_3D[:,-1:,:]
    if mask_name=='south':
        L3_XC_subset = L3_XC[:1, :]
        L3_YC_subset = L3_YC[:1, :]
        L3_wet_grid_3D_subset = L3_wet_grid_3D[:, :1, :]
    if mask_name=='east':
        L3_XC_subset = L3_XC[:, -1:]
        L3_YC_subset = L3_YC[:, -1:]
        L3_wet_grid_3D_subset = L3_wet_grid_3D[:, :, -1:]
    if mask_name=='west':
        L3_XC_subset = L3_XC[:, :1]
        L3_YC_subset = L3_YC[:, :1]
        L3_wet_grid_3D_subset = L3_wet_grid_3D[:, :, :1]

    if boundary_var_name in ['ETAN','UICE','VICE','HSNOW','HEFF','AREA']:
        L3_wet_grid_3D_subset = L3_wet_grid_3D_subset[:1,:,:]

    print('    - Reading the mask to reference the variable to the llc grid')
    nc_dict_file = os.path.join(config_dir,'L1_grid',L1_model_name,'input','L1_dv_mask_reference_dict.nc')
    source_rows, source_cols = read_mask_reference_from_nc_dict(nc_dict_file, mask_name)

    for dest_file in dest_files:
        if dest_file not in os.listdir(os.path.join(output_dir,mask_name,boundary_var_name)):
            # try:
            print('    - Downscaling the timesteps to be stored in file ' + str(dest_file))
            source_files = source_file_read_dict[dest_file]
            source_file_read_indices = source_file_read_index_sets[dest_file]

            print('        - Reading in the L1 diagnostics_vec output')
            L1_run_dir = os.path.join(config_dir,'L1_grid',L1_model_name,'run')
            L1_points, L1_values, L1_wet_grid_3D_points = \
                read_boundary_variable_to_L1_points(L1_run_dir, mask_name, boundary_var_name,
                                                    source_files, source_file_read_indices,
                                                    source_rows, source_cols,
                                                    Nr_in, L1_XC, L1_YC,
                                                    L1_AngleCS, L1_AngleSN, L1_wet_grid_3D, print_level)

            # plt.plot(L3_XC_subset,L3_YC_subset,'g.')
            # plt.plot(L1_points[:,0],L1_points[:,1],'k.')
            # plt.show()

            if boundary_var_name in ['THETA', 'SALT', 'UVEL', 'VVEL']:
                if Nr_out!=Nr_in:
                    L1_values, L1_wet_grid_3D_points = df.interpolate_var_points_timeseries_to_new_depth_levels(
                        L1_values, L1_wet_grid_3D_points, delR_in, delR_out)

            L1_wet_grid_on_L3_3D_points = np.copy(L1_wet_grid_3D_points)

            print('        - Downscaling the output to the new boundary')
            L3_boundary_var = downscale_L1_bc_field_to_L3(df, mask_name, boundary_var_name,
                                                          L1_points, L1_values,
                                                          L1_wet_grid_3D_points, L1_wet_grid_on_L3_3D_points,
                                                          L3_XC_subset, L3_YC_subset, L3_wet_grid_3D_subset,print_level)

            # if 'north' in mask_name or 'south' in mask_name:
            #     plt.subplot(1,2,1)
            #     plt.imshow(L3_wet_grid_3D_subset[:,0,:]==0,vmin=-0.1,vmax=0.1,cmap='seismic_r')
            #     plt.title('Mask')
            #     plt.subplot(1, 2, 2)
            #     plt.imshow(L3_boundary_var[0, :, 0, :],cmap='viridis')#,vmin=-0.1,vmax=0.1
            #     plt.title('L3 (nan values: '+str(np.sum(np.isnan(L3_boundary_var[0,:,0,:]))))
            #     plt.show()
            #
            # if 'east' in mask_name or 'west' in mask_name:
            #     plt.subplot(1, 2, 1)
            #     plt.imshow(L3_wet_grid_3D_subset[:, :, 0] == 0, vmin=-0.1, vmax=0.1, cmap='seismic_r')
            #     plt.title('Mask')
            #     plt.subplot(1, 2, 2)
            #     plt.imshow(L3_boundary_var[0, :, :, 0], cmap='viridis')  # ,vmin=-0.1,vmax=0.1
            #     plt.title('L3 (nan values: ' + str(np.sum(np.isnan(L3_boundary_var[0, :, :, 0]))))
            #     plt.show()

            if mask_name not in os.listdir(output_dir):
                os.mkdir(os.path.join(output_dir, mask_name))
            if boundary_var_name not in os.listdir(os.path.join(output_dir, mask_name)):
                os.mkdir(os.path.join(output_dir, mask_name, boundary_var_name))

            # if boundary_var_name in ['UVEL','VVEL','UICE','VICE']:
            #     output_file = os.path.join(output_dir, mask_name, boundary_var_name, dest_file[:-4]+'_rotated.bin')
            # else:
            output_file = os.path.join(output_dir, mask_name, boundary_var_name, dest_file)
            L3_boundary_var.ravel(order='C').astype('>f4').tofile(output_file)

        else:
            print('    - Skipping ' + str(dest_file) + ' because it was already created')

def create_bc_fields_via_interpolation(config_dir, L3_model_name, L1_model_name, proc_id,
                                       start_year, final_year, start_month,
                                       final_month, start_day, final_day, print_level):

    var_name_list = ['THETA','THETA','THETA',
                     'SALT','SALT','SALT',
                     'UVEL','UVEL','UVEL',
                     'VVEL','VVEL','VVEL',
                     'UICE','UICE','UICE',
                     'VICE','VICE','VICE',
                     'HSNOW','HSNOW','HSNOW',
                     'HEFF','HEFF','HEFF',
                     'AREA','AREA','AREA',
                     'ETAN','ETAN','ETAN']
    var_name = var_name_list[proc_id % len(var_name_list)]

    mask_name_list = ['north','south','east',
                      'north','south','east',
                      'north','south','east',
                      'north','south','east',
                      'north','south','east',
                      'north','south','east',
                      'north','south','east',
                      'north','south','east',
                      'north','south','east',
                      'north','south','east']
    mask_name = mask_name_list[proc_id % len(var_name_list)]

    print('Creating the bc field for ' + var_name + ' on mask ' +mask_name+ ' to cover year/month/days ' +
          str(start_year)+'/'+str(start_month)+'/'+str(start_day) + ' to ' +
          str(final_year)+'/'+str(final_month)+'/'+str(final_day))

    dest_files, source_file_read_dict, source_file_read_index_sets = \
        create_src_dest_dicts_from_ref(config_dir, L1_model_name, mask_name, var_name, start_year, final_year, start_month, final_month, start_day, final_day)

    print('  Running the Downscale routine:')
    create_bc_field(config_dir, L1_model_name, L3_model_name, mask_name,
                    var_name, dest_files,
                    source_file_read_dict, source_file_read_index_sets, print_level)

