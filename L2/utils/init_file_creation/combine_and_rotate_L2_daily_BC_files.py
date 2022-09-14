
import os
from datetime import datetime
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import ast
import sys

def get_dest_file_list(mask_name, var_name, Nr, n_rows_L2,n_cols_L2,
                       start_year, final_year, start_month, final_month, start_day, final_day):

    start_date = datetime(start_year, start_month, start_day)
    final_date = datetime(final_year, final_month, final_day)

    dest_files = []
    dest_file_shapes = {}
    total_timesteps = 0
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
                    if var_name in ['UVEL','VVEL','UICE','VICE']:
                        dest_file = 'L2_BC_'+mask_name+'_'+ var_name + '.' + str(year) + '{:02d}'.format(month)+ '{:02d}'.format(day) + '_rotated.bin'
                    else:
                        dest_file = 'L2_BC_'+mask_name+'_' + var_name + '.' + str(year) + '{:02d}'.format(month) + '{:02d}'.format(day) + '.bin'
                    dest_files.append(dest_file)
                    nTimesteps = 24
                    total_timesteps += nTimesteps
                    if mask_name in ['west','east']:
                        dest_file_shapes[dest_file] = (nTimesteps, Nr, n_rows_L2, 1)
                    if mask_name in ['north','south']:
                        dest_file_shapes[dest_file] = (nTimesteps, Nr, 1, n_cols_L2)

    return(dest_files,dest_file_shapes,total_timesteps)

def stack_daily_bc_files_to_one(config_dir, config_name, mask_name, var_name, AngleCS, AngleSN, dest_files, dest_file_shapes,total_timesteps,print_level):

    depth_levels = dest_file_shapes[dest_files[0]][1]
    rows = dest_file_shapes[dest_files[0]][2]
    cols = dest_file_shapes[dest_files[0]][3]

    # the 2 is added because we will duplicate the first and last field
    output_grid = np.zeros((total_timesteps+2,depth_levels,rows,cols))
    timesteps_added = 1
    for dest_file in dest_files:

        if var_name in ['UVEL','VVEL','UICE','VICE']:
            if 'VEL' in var_name:
                if var_name=='UVEL':
                    u_dest_file = dest_file
                else:
                    u_dest_file = dest_file.replace('VVEL','UVEL')
                if var_name=='VVEL':
                    v_dest_file = dest_file
                else:
                    v_dest_file = dest_file.replace('UVEL','VVEL')
                u_var_grid = np.fromfile(os.path.join(config_dir, 'L2', config_name, 'input', 'obcs',mask_name, 'UVEL', u_dest_file), '>f4')
                u_var_grid = np.reshape(u_var_grid, dest_file_shapes[dest_file])
                v_var_grid = np.fromfile(os.path.join(config_dir, 'L2', config_name, 'input', 'obcs',mask_name, 'VVEL', v_dest_file),'>f4')
                v_var_grid = np.reshape(v_var_grid, dest_file_shapes[dest_file])

            if 'ICE' in var_name:
                if var_name=='UICE':
                    u_dest_file = dest_file
                else:
                    u_dest_file = dest_file.replace('VICE','UICE')
                if var_name=='VICE':
                    v_dest_file = dest_file
                else:
                    v_dest_file = dest_file.replace('UICE','VICE')
                u_var_grid = np.fromfile(os.path.join(config_dir, 'L2', config_name, 'input', 'obcs',mask_name, 'UICE', u_dest_file), '>f4')
                u_var_grid = np.reshape(u_var_grid, dest_file_shapes[dest_file])
                v_var_grid = np.fromfile(os.path.join(config_dir, 'L2', config_name, 'input', 'obcs',mask_name, 'VICE', v_dest_file),'>f4')
                v_var_grid = np.reshape(v_var_grid, dest_file_shapes[dest_file])

            var_grid = np.zeros_like(u_var_grid)
            for timestep in range(np.shape(var_grid)[0]):
                for k in range(np.shape(var_grid)[1]):
                    if var_name=='UICE' or var_name=='UVEL':
                        var_grid[timestep, k, :, :] = AngleCS*u_var_grid[timestep, k, :, :] + AngleSN*v_var_grid[timestep, k, :, :]
                    if var_name=='VICE' or var_name=='VVEL':
                        var_grid[timestep, k, :, :] = -1*AngleSN*u_var_grid[timestep, k, :, :] + AngleCS*v_var_grid[timestep, k, :, :]

        else:
            var_grid = np.fromfile(os.path.join(config_dir,'L2',config_name, 'input','obcs',mask_name,var_name,dest_file),'>f4')
            var_grid = np.reshape(var_grid,dest_file_shapes[dest_file])

        if print_level>=2:
            print('        - Adding timesteps from file '+dest_file+' to levels '+str(timesteps_added)+' to '+str(timesteps_added+np.shape(var_grid)[0]))
        if print_level >= 4:
            print('                - The mean value on this day is ' + str(np.mean(var_grid)))
            print('                - There are '+str(np.sum(np.isnan(var_grid[0,:,:])))+' nan values')
        # rows = np.where(np.isnan(var_grid[0,:,:]))
        # print(rows)
        # print(cols)
        output_grid[timesteps_added:timesteps_added+np.shape(var_grid)[0],:,:,:] = var_grid
        timesteps_added += np.shape(var_grid)[0]

    # here we duplicate the first and last field
    # this is done so that all timesteps can be interpolated by the model
    output_grid[0,:,:,:] = output_grid[1,:,:,:]
    output_grid[-1,:,:,:] = output_grid[-2,:,:,:]
    return(output_grid)


def combine_and_rotate_L2_daily_bcs(config_dir, config_name, proc_id,
                                         start_year, final_year, start_month, final_month, start_day, final_day, print_level):

    var_name_list = ['THETA', 'THETA', 'THETA',
                     'SALT', 'SALT', 'SALT',
                     'UVEL', 'UVEL', 'UVEL',
                     'VVEL', 'VVEL', 'VVEL',
                     'UICE', 'UICE', 'UICE',
                     'VICE', 'VICE', 'VICE',
                     'HSNOW', 'HSNOW', 'HSNOW',
                     'HEFF', 'HEFF', 'HEFF',
                     'AREA', 'AREA', 'AREA',
                     'ETAN', 'ETAN', 'ETAN']
    var_name = var_name_list[proc_id % len(var_name_list)]

    mask_name_list = ['north', 'south', 'east',
                      'north', 'south', 'east',
                      'north', 'south', 'east',
                      'north', 'south', 'east',
                      'north', 'south', 'east',
                      'north', 'south', 'east',
                      'north', 'south', 'east',
                      'north', 'south', 'east',
                      'north', 'south', 'east',
                      'north', 'south', 'east']
    mask_name = mask_name_list[proc_id % len(var_name_list)]

    grid_file = os.path.join(config_dir,'nc_grids',config_name+'_grid.nc')
    ds = nc4.Dataset(grid_file)
    AngleCS = ds.variables['AngleCS'][:,:]
    AngleSN = ds.variables['AngleSN'][:, :]
    n_rows_L2 = np.shape(AngleCS)[0]
    n_cols_L2 = np.shape(AngleCS)[1]
    drF = ds.variables['drF'][:]
    Nr = len(drF)
    ds.close()

    if mask_name=='north':
        AngleCS = AngleCS[-1:,:]
        AngleSN = AngleSN[-1:,:]
    if mask_name=='south':
        AngleCS = AngleCS[:1,:]
        AngleSN = AngleSN[:1,:]
    if mask_name=='west':
        AngleCS = AngleCS[:,:1]
        AngleSN = AngleSN[:,:1]
    if mask_name=='east':
        AngleCS = AngleCS[:,-1:]
        AngleSN = AngleSN[:,-1:]

    if var_name in ['ETAN','UICE','VICE','AREA','HSNOW','HEFF']:Nr = 1

    if print_level >=1:
        print('    - Combining the monthly files for ' + var_name)
    dest_files, dest_file_shapes, total_timesteps = get_dest_file_list(mask_name, var_name, Nr, n_rows_L2, n_cols_L2,
                                                                       start_year, final_year, start_month, final_month, start_day, final_day)

    if print_level >= 1:
        print('    - Stacking all of the daily files into a big global file')
    output_grid = stack_daily_bc_files_to_one(config_dir, config_name, mask_name, var_name,
                                              AngleCS, AngleSN, dest_files, dest_file_shapes, total_timesteps, print_level)


    output_file = os.path.join(config_dir, 'L2', config_name, 'input', 'obcs', 'L2_BC_'+mask_name+'_' + var_name + '.bin')
    if print_level >= 1:
        print('    - Outputting to ' + output_file)
    output_grid.ravel('C').astype('>f4').tofile(output_file)