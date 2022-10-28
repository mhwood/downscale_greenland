
import os
from datetime import datetime
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import ast
import sys

def get_dest_file_list(model_name, var_name, sNx, sNy, ordered_nonblank_tiles, year):

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
            dest_file = prefix+'_surface_' + var_name + '.' + str(year) + '{:02d}'.format(month) + '{:02d}'.format(day) + '.bin'
            dest_files.append(dest_file)
            nTimesteps = 4
            total_timesteps += nTimesteps
            dest_file_shapes[dest_file] = (nTimesteps, sNy * len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0]), sNx)

    return(dest_files,dest_file_shapes,total_timesteps)

def read_stitched_grid_tile_angles(config_dir,model_name,sNx,sNy,ordered_nonblank_tiles,ordered_nonblank_rotations):

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    stitched_AngleCS = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    stitched_AngleSN = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    grid_dir = os.path.join(config_dir, 'L1', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    AngleCS = ds.variables['AngleCS'][:, :]
                    AngleSN = ds.variables['AngleSN'][:, :]
                    ds.close()

                    for i in range(ordered_nonblank_rotations[r][c]):
                        AngleCS = np.rot90(AngleCS)
                        AngleSN = np.rot90(AngleSN)

                    stitched_AngleCS[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = AngleCS
                    stitched_AngleSN[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = AngleSN

    return(stitched_AngleCS, stitched_AngleSN)

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
    # uvel = np.zeros_like(zonal_vel)
    # vvel = np.zeros_like(meridional_vel)
    uvel = angle_cos * zonal_vel + angle_sin * meridional_vel
    vvel = -1 * angle_sin * zonal_vel + angle_cos * meridional_vel
    return (uvel, vvel)

def rotate_velocity_grid_to_domain(var_name,u_var_grid,v_var_grid,AngleCS_subset, AngleSN_subset):
    rotated_u = np.zeros_like(u_var_grid)
    rotated_v = np.zeros_like(v_var_grid)
    for timestep in range(np.shape(u_var_grid)[0]):
        u, v = rotate_velocity_vectors_to_domain(AngleCS_subset, AngleSN_subset,
                                                 u_var_grid[timestep, :, :], v_var_grid[timestep, :, :])
        rotated_u[timestep, :, :] = u
        rotated_v[timestep, :, :] = v
    if var_name[0]=='U':
        return(rotated_u)
    if var_name[0]=='V':
        return(rotated_v)

def stack_daily_exf_files_to_one(config_dir, model_name, var_name,
                                dest_files, dest_file_shapes, total_timesteps,
                                AngleCS_subset, AngleSN_subset, print_level):

    rows = dest_file_shapes[dest_files[0]][1]
    cols = dest_file_shapes[dest_files[0]][2]

    output_grid = np.zeros((total_timesteps, rows, cols))
    timesteps_added = 0
    for dest_file in dest_files:
        if print_level>=3:
            print('            - Adding timesteps from file ' + dest_file)
        if var_name in ['UWIND','VWIND','USTRESS','VSTRESS']:
            if var_name=='USTRESS':
                u_dest_file = dest_file[:-4]+'_rotated.bin'
                v_dest_file = u_dest_file.replace('USTRESS','VSTRESS')
            if var_name=='VSTRESS':
                v_dest_file = dest_file[:-4]+'_rotated.bin'
                u_dest_file = v_dest_file.replace('VSTRESS','USTRESS')
            if var_name=='UWIND':
                u_dest_file = dest_file[:-4]+'_rotated.bin'
                v_dest_file = u_dest_file.replace('UWIND','VWIND')
            if var_name=='VWIND':
                v_dest_file = dest_file[:-4]+'_rotated.bin'
                u_dest_file = v_dest_file.replace('VWIND','UWIND')
            u_var_grid = np.fromfile(os.path.join(config_dir, 'L1',model_name, 'input', 'exf', 'U'+var_name[1:], u_dest_file),'>f4')
            v_var_grid = np.fromfile(os.path.join(config_dir, 'L1',model_name, 'input', 'exf', 'V'+var_name[1:], v_dest_file),'>f4')
            u_var_grid = np.reshape(u_var_grid, dest_file_shapes[dest_file])
            v_var_grid = np.reshape(v_var_grid, dest_file_shapes[dest_file])
            if print_level >= 4:
                print('                - Rotating the grid to the domain orientation')
            var_grid = rotate_velocity_grid_to_domain(var_name,u_var_grid,v_var_grid,AngleCS_subset, AngleSN_subset)
        else:
            var_grid = np.fromfile(os.path.join(config_dir, 'L1', model_name, 'input', 'exf', var_name, dest_file), '>f4')
            var_grid = np.reshape(var_grid, dest_file_shapes[dest_file])

        if print_level >= 4:
            print('                - Timesteps added to levels ' + str(timesteps_added) + ' through ' + str(timesteps_added + np.shape(var_grid)[0] - 1))
        if print_level >= 5:
            print('                    - The mean value on this day is ' + str(np.mean(var_grid)))
            print('                    - There are ' + str(np.sum(np.isnan(var_grid[0, :, :]))) + ' nan values')
        # rows = np.where(np.isnan(var_grid[0,:,:]))
        # print(rows)
        # print(cols)
        output_grid[timesteps_added:timesteps_added + np.shape(var_grid)[0], :, :] = var_grid
        timesteps_added += np.shape(var_grid)[0]

    return (output_grid)



def combine_L1_daily_exf_files_annual(Lf, config_dir, model_name, var_name,
                                     sNx, sNy, ordered_nonblank_tiles,ordered_nonblank_rotations,
                                     start_year, final_year,
                                     print_level):

    for year in range(start_year,final_year+1):

        if print_level >= 2:
            print('        - Determining the size of the input/output files')
        dest_files, dest_file_shapes, total_timesteps = get_dest_file_list(model_name, var_name, sNx, sNy, ordered_nonblank_tiles, year)

        if var_name in ['UWIND','VWIND','USTRESS','VSTRESS']:
            if print_level >= 2:
                print('        - Reading the grid orientationn')
            AngleCS, AngleSN = read_stitched_grid_tile_angles(config_dir,model_name,sNx,sNy,ordered_nonblank_tiles,ordered_nonblank_rotations)

            AngleCS_faces = Lf.read_stitched_grid_to_faces(AngleCS, sNx, sNy, dim=2)
            AngleSN_faces = Lf.read_stitched_grid_to_faces(AngleSN, sNx, sNy, dim=2)

            AngleCS_subset = Lf.read_faces_to_compact(AngleCS_faces, sNx, sNy, dim=2)
            AngleSN_subset = Lf.read_faces_to_compact(AngleSN_faces, sNx, sNy, dim=2)

        else:
            AngleCS_subset = AngleSN_subset = []

        if print_level >= 1:
            print('    - Stacking all of the daily files into the '+str(year)+') file')
        output_grid = stack_daily_exf_files_to_one(config_dir, model_name, var_name,
                                                  dest_files, dest_file_shapes, total_timesteps,
                                                  AngleCS_subset, AngleSN_subset, print_level)

        output_file = os.path.join(config_dir, 'L1', model_name, 'input', 'exf', 'L1_exf_' + var_name+ '_' + str(year))
        if print_level >= 1:
            print('    - Outputting to ' + output_file)
        if print_level >= 2:
            print('        - Output shape: '+str(np.shape(output_grid)))
        output_grid.ravel('C').astype('>f4').tofile(output_file)


