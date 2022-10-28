
import os
from datetime import datetime
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import ast
import sys

def read_grid_geometry(config_dir,model_name):

    file_path = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    AngleSN = ds.variables['AngleSN'][:, :]
    delR = ds.variables['drF'][:]
    ds.close()

    return(XC, YC, AngleCS, AngleSN, delR)

def get_dest_file_list(model_name, boundary, var_name, n_timesteps_per_day, Nr, n_rows, n_cols, year):

    prefix = '_'.join(model_name.split('_')[:2])

    dest_files = []
    dest_file_shapes = {}
    total_timesteps = 0
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
            if var_name in ['ETAN','AREA','HEFF','HSNOW','UICE','VICE']:
                depth_levels = 1
            else:
                depth_levels = Nr
            dest_file = prefix+'_' + boundary + '_' + var_name + '.' + str(year) + '{:02d}'.format(month) + '{:02d}'.format(day) + '.bin'
            dest_files.append(dest_file)
            nTimesteps = n_timesteps_per_day
            total_timesteps += nTimesteps
            if boundary == 'west' or boundary=='east':
                dest_file_shapes[dest_file] = (nTimesteps, depth_levels, n_rows, 1)
            else:
                dest_file_shapes[dest_file] = (nTimesteps, depth_levels, 1, n_cols)

    return(dest_files,dest_file_shapes,total_timesteps)

def subset_tile_angles_to_boundary(boundary, AngleCS, AngleSN):

    # subset to the boundary
    if boundary=='south':
        AngleCS_subset = AngleCS[:1, :]
        AngleSN_subset = AngleSN[:1, :]

    if boundary=='west':
        AngleCS_subset = AngleCS[:, :1]
        AngleSN_subset = AngleSN[:, :1]

    if boundary=='north':
        AngleCS_subset = AngleCS[-1:, :]
        AngleSN_subset = AngleSN[-1:, :]

    if boundary=='east':
        AngleCS_subset = AngleCS[:, -1:]
        AngleSN_subset = AngleSN[:, -1:]

    return(AngleCS_subset, AngleSN_subset)

def rotate_velocity_vectors_to_domain(angle_cos, angle_sin, zonal_vel, meridional_vel):
    uvel = np.zeros_like(zonal_vel)
    vvel = np.zeros_like(meridional_vel)
    for k in range(np.shape(uvel)[0]):
        tmp_uvel = angle_cos.ravel() * zonal_vel[k, :, :].ravel() + angle_sin.ravel() * meridional_vel[k, :, :].ravel()
        tmp_vvel = -1 * angle_sin.ravel() * zonal_vel[k, :, :].ravel() + angle_cos.ravel() * meridional_vel[k, :, :].ravel()
        uvel[k, :, :] = tmp_uvel.reshape(np.shape(uvel[k, :, :]))
        vvel[k, :, :] = tmp_vvel.reshape(np.shape(vvel[k, :, :]))
    return (uvel, vvel)

def rotate_velocity_grid_to_domain(var_name,u_var_grid,v_var_grid,AngleCS_subset, AngleSN_subset):
    rotated_u = np.zeros_like(u_var_grid)
    rotated_v = np.zeros_like(v_var_grid)
    for timestep in range(np.shape(u_var_grid)[0]):
        u, v = rotate_velocity_vectors_to_domain(AngleCS_subset, AngleSN_subset,
                                                 u_var_grid[timestep, :, :, :], v_var_grid[timestep, :, :, :])
        rotated_u[timestep, :, :, :] = u
        rotated_v[timestep, :, :, :] = v
    if var_name[0]=='U':
        return(rotated_u)
    if var_name[0]=='V':
        return(rotated_v)

def stack_daily_bc_files_to_one(config_dir, model_name, boundary, var_name,
                                dest_files, dest_file_shapes, total_timesteps,
                                AngleCS_subset, AngleSN_subset, print_level):
    depth_levels = dest_file_shapes[dest_files[0]][1]
    rows = dest_file_shapes[dest_files[0]][2]
    cols = dest_file_shapes[dest_files[0]][3]

    output_grid = np.zeros((total_timesteps, depth_levels, rows, cols))
    timesteps_added = 0
    for dest_file in dest_files:
        if print_level>=3:
            print('            - Adding timesteps from file ' + dest_file)
        if var_name in ['UICE','VICE','UVEL','VVEL']:
            if var_name=='UVEL':
                u_dest_file = dest_file[:-4]+'_rotated.bin'
                v_dest_file = u_dest_file.replace('UVEL','VVEL')
            if var_name=='VVEL':
                v_dest_file = dest_file[:-4]+'_rotated.bin'
                u_dest_file = v_dest_file.replace('VVEL','UVEL')
            if var_name=='UICE':
                u_dest_file = dest_file[:-4]+'_rotated.bin'
                v_dest_file = u_dest_file.replace('UICE','VICE')
            if var_name=='VICE':
                v_dest_file = dest_file[:-4]+'_rotated.bin'
                u_dest_file = v_dest_file.replace('VICE','UICE')
            u_var_grid = np.fromfile(os.path.join(config_dir, 'L1_grid',model_name, 'input', 'obcs', boundary, 'U'+var_name[1:], u_dest_file),'>f4')
            v_var_grid = np.fromfile(os.path.join(config_dir, 'L1_grid',model_name, 'input', 'obcs', boundary, 'V'+var_name[1:], v_dest_file),'>f4')
            u_var_grid = np.reshape(u_var_grid, dest_file_shapes[dest_file])
            v_var_grid = np.reshape(v_var_grid, dest_file_shapes[dest_file])
            if print_level >= 4:
                print('                - Rotating the grid to the domain orientation')
            var_grid = rotate_velocity_grid_to_domain(var_name,u_var_grid,v_var_grid,AngleCS_subset, AngleSN_subset)
        else:
            var_grid = np.fromfile(os.path.join(config_dir, 'L1_grid', model_name, 'input', 'obcs', boundary, var_name, dest_file), '>f4')
            var_grid = np.reshape(var_grid, dest_file_shapes[dest_file])

        if print_level >= 4:
            print('                - Timesteps added to levels ' + str(timesteps_added) + ' through ' + str(timesteps_added + np.shape(var_grid)[0] - 1))
        if print_level >= 5:
            print('                    - The mean value on this day is ' + str(np.mean(var_grid)))
            print('                    - There are ' + str(np.sum(np.isnan(var_grid[0, :, :]))) + ' nan values')
        # rows = np.where(np.isnan(var_grid[0,:,:]))
        # print(rows)
        # print(cols)
        output_grid[timesteps_added:timesteps_added + np.shape(var_grid)[0], :, :, :] = var_grid
        timesteps_added += np.shape(var_grid)[0]

    return (output_grid)

def combine_L1_daily_BC_files_annual(config_dir, model_name, var_name,
                                     start_year, final_year, n_timesteps_per_day,
                                     print_level):

    XC, YC, AngleCS, AngleSN, delR = read_grid_geometry(config_dir,model_name)
    Nr = len(delR)
    n_rows = np.shape(XC)[0]
    n_cols = np.shape(XC)[1]

    for boundary in ['east','west','south','north']:
        if boundary in os.listdir(os.path.join(config_dir,'L1_grid',model_name,'input','obcs')):

            if print_level>=1:
                print('    - Combining the daily files for ' + var_name+' on the '+boundary+' boundary')

            for year in range(start_year,final_year+1):

                if print_level >= 2:
                    print('        - Working on year '+str(year))

                if print_level >= 2:
                    print('        - Determining the size of the input/output files')
                dest_files, dest_file_shapes, total_timesteps = get_dest_file_list(model_name, boundary, var_name, n_timesteps_per_day, Nr, n_rows, n_cols, year)

                # for d in range(len(dest_files)):
                #     print(dest_files[d],dest_file_shapes[dest_files[d]])

                if var_name in ['UVEL','VVEL','UICE','VICE']:
                    if print_level >= 2:
                        print('        - Reading the grid orientation')
                    AngleCS_subset, AngleSN_subset = \
                        subset_tile_angles_to_boundary(boundary, AngleCS, AngleSN)
                else:
                    AngleCS_subset = AngleSN_subset = []

                if print_level >= 1:
                    print('    - Stacking all of the daily files into the '+str(year)+' file')
                output_grid = stack_daily_bc_files_to_one(config_dir, model_name, boundary, var_name,
                                                          dest_files, dest_file_shapes, total_timesteps,
                                                          AngleCS_subset, AngleSN_subset, print_level)

                output_file = os.path.join(config_dir, 'L1_grid', model_name, 'input', 'obcs', 'L1_BC_'+boundary+'_' + var_name+'_' + str(year))
                if print_level >= 1:
                    print('    - Outputting to ' + output_file)
                output_grid.ravel('C').astype('>f4').tofile(output_file)


