
import os
from datetime import datetime
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import ast
import sys


def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    delR = ds.variables['drF'][:]
    ds.close()
    return(XC, YC, delR)

def get_dest_file_list(config_dir, config_name, var_name, n_rows_L3, n_cols_L3,
                       start_year, final_year, start_month, final_month, start_day, final_day, print_level):

    start_date = datetime(start_year, start_month, start_day)
    final_date = datetime(final_year, final_month, final_day)

    exf_dir = os.path.join(config_dir, 'L3', config_name,'input','exf',var_name)

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
                    dest_file = 'L3_exf_' + var_name + '.' + str(year) + '{:02d}'.format(month)+ '{:02d}'.format(day) + '.bin'
                    dest_files.append(dest_file)
                    n_Timesteps = int(np.size(np.fromfile(os.path.join(exf_dir, dest_file),'>f4') )/(n_rows_L3*n_cols_L3))
                    total_timesteps += n_Timesteps
                    dest_file_shapes[dest_file] = (n_Timesteps, n_rows_L3, n_cols_L3)

    return(dest_files,dest_file_shapes,total_timesteps)

def stack_daily_exf_files_to_one(config_dir, config_name, var_name, dest_files, dest_file_shapes, total_timesteps, print_level):

    rows = dest_file_shapes[dest_files[0]][1]
    cols = dest_file_shapes[dest_files[0]][2]

    # the 2 is added because we will duplicate the first and last field
    output_grid = np.zeros((total_timesteps+2,rows,cols))
    timesteps_added = 1
    for dest_file in dest_files:
        var_grid = np.fromfile(os.path.join(config_dir,'L3',config_name, 'input','exf',var_name,dest_file),'>f4')

        var_grid = np.reshape(var_grid,dest_file_shapes[dest_file])
        if print_level>=3:
            print('            - Adding timesteps from file '+dest_file+' to levels '+str(timesteps_added)+' to '+str(timesteps_added+np.shape(var_grid)[0]))
        if print_level >= 4:
            print('                - The mean value on this day is ' + str(np.mean(var_grid)))
            print('                - There are '+str(np.sum(np.isnan(var_grid[0,:,:])))+' nan values')
        # rows = np.where(np.isnan(var_grid[0,:,:]))
        # print(rows)
        # print(cols)
        output_grid[timesteps_added:timesteps_added+np.shape(var_grid)[0],:,:] = var_grid
        timesteps_added += np.shape(var_grid)[0]

    # here we duplicate the first and last field
    # this is done so that all timesteps can be interpolated by the model
    output_grid[0,:,:] = output_grid[1,:,:]
    output_grid[-1,:,:] = output_grid[-2,:,:]
    return(output_grid)


def combine_L3_daily_exf_files(config_dir, config_name, proc_id,
                               start_year, final_year, start_month, final_month, start_day, final_day, print_level):

    var_name_list = ['UWIND','VWIND','ATEMP','AQH','PRECIP','SWDOWN','LWDOWN','RUNOFF']
    var_name = var_name_list[proc_id % len(var_name_list)]

    XC, _, _ = read_grid_geometry_from_nc(config_dir, config_name)
    n_rows_L3 = np.shape(XC)[0]
    n_cols_L3 = np.shape(XC)[1]

    if print_level>=1:
        print('    - Combining the daily external forcing files for ' + var_name)

    dest_files, dest_file_shapes, total_timesteps = get_dest_file_list(config_dir, config_name, var_name,
                                                                       n_rows_L3, n_cols_L3,
                                                                       start_year, final_year, start_month,
                                                                       final_month, start_day, final_day,
                                                                       print_level)

    if print_level >= 1:
        print('    - Stacking all of the daily files into a big global file')
    output_grid = stack_daily_exf_files_to_one(config_dir, config_name, var_name, dest_files, dest_file_shapes, total_timesteps, print_level)

    output_file = os.path.join(config_dir, 'L3', config_name, 'input', 'exf', 'L3_exf_' + var_name + '.bin')
    if print_level >= 1:
        print('    - Outputting to ' + output_file)
    output_grid.ravel('C').astype('>f4').tofile(output_file)
