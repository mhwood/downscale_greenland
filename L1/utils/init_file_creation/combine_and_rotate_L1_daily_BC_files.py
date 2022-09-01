
import os
from datetime import datetime
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import ast
import sys

def get_dest_file_list(model_name, boundary, var_name, Nr, sNx, sNy, tile_numbers,
                       start_year, final_year, start_month, final_month, start_day, final_day):

    start_date = datetime(start_year, start_month, start_day)
    final_date = datetime(final_year, final_month, final_day)

    prefix = '_'.join(model_name.split('_')[:2])

    dest_files = []
    dest_file_shapes = {}
    total_timesteps = 0
    for year in range(1992, 2023):
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
                test_date = datetime(year, month, day)
                if test_date >= start_date and test_date <= final_date:
                    dest_file = prefix+'_' + boundary + '_' + var_name + '.' + str(year) + '{:02d}'.format(
                        month) + '{:02d}'.format(day) + '.bin'
                    dest_files.append(dest_file)
                    nTimesteps = 24
                    total_timesteps += nTimesteps
                    if boundary == 'west':
                        dest_file_shapes[dest_file] = (nTimesteps, depth_levels, sNy*len(tile_numbers), 1)
                    else:
                        dest_file_shapes[dest_file] = (nTimesteps, depth_levels, 1, sNx*len(tile_numbers))

    return(dest_files,dest_file_shapes,total_timesteps)

def read_grid_tile_angles(config_dir,model_name,ordered_nonblank_tiles):

    ordered_AngleCS_tiles = []
    ordered_AngleSN_tiles = []

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L1', model_name, 'run_for_grid')


    for r in range(len(ordered_nonblank_tiles)):
        row_AngleCSs = []
        row_AngleSNs = []
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    AngleCS = ds.variables['AngleCS'][:, :]
                    AngleSN = ds.variables['AngleSN'][:, :]
                    ds.close()
                    row_AngleCSs.append(AngleCS)
                    row_AngleSNs.append(AngleSN)
        ordered_AngleCS_tiles.append(row_AngleCSs)
        ordered_AngleSN_tiles.append(row_AngleSNs)

    return(ordered_AngleCS_tiles, ordered_AngleSN_tiles)

def subset_tile_angles_to_boundary(boundary, sNx, sNy, tile_numbers, ordered_nonblank_tiles,
                                     ordered_AngleCS_tiles, ordered_AngleSN_tiles):

    angles_started = False

    for tile_number in tile_numbers:

        if tile_number>0:

            # get the geometry for this tile
            for r in range(len(ordered_nonblank_tiles)):
                for c in range(len(ordered_nonblank_tiles[0])):
                    if ordered_nonblank_tiles[r][c] == tile_number:
                        AngleCS = ordered_AngleCS_tiles[r][c]
                        AngleSN = ordered_AngleSN_tiles[r][c]

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

        else:
            if boundary in ['north','south']:
                AngleCS_subset = np.zeros((1, sNx))
                AngleSN_subset = np.zeros((1, sNx))
            if boundary in ['east','west']:
                AngleCS_subset = np.zeros((sNy, 1))
                AngleSN_subset = np.zeros((sNy, 1))

        if not angles_started:
            angles_started = True
            full_AngleCS_subset = AngleCS_subset
            full_AngleSN_subset = AngleSN_subset
        else:
            full_AngleCS_subset = np.vstack([full_AngleCS_subset,AngleCS_subset])
            full_AngleSN_subset = np.vstack([full_AngleSN_subset,AngleSN_subset])

    full_AngleCS_subset = full_AngleCS_subset.ravel()
    full_AngleSN_subset = full_AngleSN_subset.ravel()

    return(full_AngleCS_subset, full_AngleSN_subset)

def rotate_velocity_vectors_to_domain(angle_cos, angle_sin, zonal_vel, meridional_vel):
    uvel = np.zeros_like(zonal_vel)
    vvel = np.zeros_like(meridional_vel)
    for k in range(np.shape(uvel)[0]):
        tmp_uvel = angle_cos * zonal_vel[k, :, :].ravel() + angle_sin * meridional_vel[k, :, :].ravel()
        tmp_vvel = -1 * angle_sin * zonal_vel[k, :, :].ravel() + angle_cos * meridional_vel[k, :, :].ravel()
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

    # the 2 is added because we will duplicate the first and last field
    output_grid = np.zeros((total_timesteps + 2, depth_levels, rows, cols))
    timesteps_added = 1
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
            u_var_grid = np.fromfile(os.path.join(config_dir, 'L1',model_name, 'input', 'obcs', boundary, 'U'+var_name[1:], u_dest_file),'>f4')
            v_var_grid = np.fromfile(os.path.join(config_dir, 'L1',model_name, 'input', 'obcs', boundary, 'V'+var_name[1:], v_dest_file),'>f4')
            u_var_grid = np.reshape(u_var_grid, dest_file_shapes[dest_file])
            v_var_grid = np.reshape(v_var_grid, dest_file_shapes[dest_file])
            if print_level >= 3:
                print('            - Rotating the grid to the domain orientation')
            var_grid = rotate_velocity_grid_to_domain(var_name,u_var_grid,v_var_grid,AngleCS_subset, AngleSN_subset)
        else:
            var_grid = np.fromfile(os.path.join(config_dir, 'L1', model_name, 'input', 'obcs', boundary, var_name, dest_file), '>f4')
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

    # here we duplicate the first and last field
    # this is done so that all timesteps can be interpolated by the model
    output_grid[0, :, :, :] = output_grid[1, :, :, :]
    output_grid[-1, :, :, :] = output_grid[-2, :, :, :]
    return (output_grid)



def combine_L1_daily_BC_files(config_dir, model_name, var_name,
                              Nr, sNx, sNy, ordered_nonblank_tiles,
                              northern_tiles, southern_tiles, eastern_tiles, western_tiles,
                              start_year, final_year, start_month, final_month, start_day, final_day,
                              print_level):

    for boundary in ['south','east','west','north']:
        if boundary in os.listdir(os.path.join(config_dir,'L1',model_name,'input','obcs')):

            if print_level>=1:
                print('    - Combining the daily files for ' + var_name+' on the '+boundary+' boundary')

            if boundary == 'south':
                tile_numbers = southern_tiles
            if boundary == 'west':
                tile_numbers = western_tiles
            if boundary == 'north':
                tile_numbers = northern_tiles
            if boundary == 'east':
                tile_numbers = eastern_tiles

            if print_level >= 2:
                print('        - Determining the size of the input/output files')
            dest_files, dest_file_shapes, total_timesteps = get_dest_file_list(model_name, boundary, var_name, Nr, sNx, sNy, tile_numbers,
                                                                               start_year, final_year, start_month, final_month, start_day, final_day)

            if var_name in ['UVEL','VVEL','UICE','VICE']:
                if print_level >= 2:
                    print('        - Reading the grid orientationn')
                ordered_AngleCS_tiles, ordered_AngleSN_tiles = read_grid_tile_angles(config_dir, model_name, ordered_nonblank_tiles)
                AngleCS_subset, AngleSN_subset = \
                    subset_tile_angles_to_boundary(boundary, sNx, sNy, tile_numbers, ordered_nonblank_tiles,
                                                   ordered_AngleCS_tiles, ordered_AngleSN_tiles)
            else:
                AngleCS_subset = AngleSN_subset = []

            if print_level >= 1:
                print('    - Stacking all of the daily files into a big global file')
            output_grid = stack_daily_bc_files_to_one(config_dir, model_name, boundary, var_name,
                                                      dest_files, dest_file_shapes, total_timesteps,
                                                      AngleCS_subset, AngleSN_subset, print_level)


            output_file = os.path.join(config_dir, 'L1', model_name, 'input', 'obcs', 'L1_BC_'+boundary+'_' + var_name + '.bin')
            if print_level >= 1:
                print('    - Outputting to ' + output_file)
            output_grid.ravel('C').astype('>f4').tofile(output_file)


