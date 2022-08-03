
import os
import numpy as np
import netCDF4 as nc4
import ast
import matplotlib.pyplot as plt
import sys


def read_grid_tile_angles(config_dir,model_name,ordered_nonblank_tiles):

    ordered_AngleCS_tiles = []
    ordered_AngleSN_tiles = []

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L05', model_name, 'run_for_grid')


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

def subset_tile_angles_to_boundary(boundary, tile_numbers, ordered_nonblank_tiles,
                                     ordered_AngleCS_tiles, ordered_AngleSN_tiles):

    angles_started = False

    for tile_number in tile_numbers:

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
        uvel[k, :] = angle_cos * zonal_vel[k, :] + angle_sin * meridional_vel[k, :]
        vvel[k, :] = -1 * angle_sin * zonal_vel[k, :] + angle_cos * meridional_vel[k, :]
    return (uvel, vvel)


def rotate_BCs(config_dir,model_name,sNx,sNy,ordered_nonblank_tiles,tile_face_index_dict,
               northern_tiles, southern_tiles, eastern_tiles, western_tiles):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import aste_functions as af

    print('    - Rotating the UVEL, VVEL, SIuice, and SIvice BCs for the '+model_name+' model from ASTE data')

    # step 0: get the model domain
    print('    - Reading in the model geometry')
    ordered_AngleCS_tiles, ordered_AngleSN_tiles = read_grid_tile_angles(config_dir,model_name,ordered_nonblank_tiles)

    pairs = [['UVEL','VVEL'],['SIuice', 'SIvice']]

    for boundary in ['west','south','east','north']:
        if boundary == 'south':
            tile_numbers = southern_tiles
        if boundary == 'west':
            tile_numbers = western_tiles
        if boundary == 'north':
            tile_numbers = northern_tiles
        if boundary == 'east':
            tile_numbers = eastern_tiles

        if len(tile_numbers) > 0:

            for pair in pairs:

                # subset the geometry to the boundary
                AngleCS_subset, AngleSN_subset = \
                    subset_tile_angles_to_boundary(boundary, tile_numbers, ordered_nonblank_tiles,
                                                     ordered_AngleCS_tiles, ordered_AngleSN_tiles)

                u_file = os.path.join(config_dir, 'L05', model_name, 'input', 'obcs',
                                           'L05_BC_' + boundary + '_'+pair[0]+'_rotated.bin')
                u_grid = np.fromfile(u_file,'>f4')

                v_file = os.path.join(config_dir, 'L05', model_name, 'input', 'obcs',
                                         'L05_BC_' + boundary + '_'+pair[1]+'_rotated.bin')
                v_grid = np.fromfile(v_file, '>f4')

                if pair[0]=='UVEL':
                    Nr = 50
                    print('     - Rotating VEL fields on the ' + boundary + ' boundary')
                else:
                    Nr = 1
                    print('     - Rotating SI fields on the ' + boundary + ' boundary')

                if boundary in ['north','south']:
                    n_timesteps = int(np.size(u_grid)/(Nr*len(tile_numbers)*sNx))
                    u_grid = u_grid.reshape((n_timesteps, Nr, len(tile_numbers) * sNx))
                    v_grid = v_grid.reshape((n_timesteps, Nr, len(tile_numbers) * sNy))
                else:
                    n_timesteps = int(np.size(u_grid) / (Nr * len(tile_numbers) * sNy))
                    u_grid = u_grid.reshape((n_timesteps, Nr, len(tile_numbers) * sNx))
                    v_grid = v_grid.reshape((n_timesteps, Nr, len(tile_numbers) * sNy))

                rotated_u = np.zeros_like(u_grid)
                rotated_v = np.zeros_like(v_grid)

                for timestep in range(np.shape(u_grid)[0]):
                    u, v = rotate_velocity_vectors_to_domain(AngleCS_subset, AngleSN_subset,
                                                                   u_grid[timestep,:,:], v_grid[timestep,:,:])
                    rotated_u[timestep,:,:] = u
                    rotated_v[timestep, :, :] = v

                u_output_file = os.path.join(config_dir, 'L05', model_name, 'input', 'obcs',
                                                'L05_BC_' + boundary + '_'+pair[0]+'.bin')
                rotated_u.ravel(order='C').astype('>f4').tofile(u_output_file)

                print('       - Output file shape: '+str(np.shape(rotated_u)))

                v_output_file = os.path.join(config_dir, 'L05', model_name, 'input', 'obcs',
                                                'L05_BC_' + boundary + '_'+pair[1]+'.bin')
                rotated_v.ravel(order='C').astype('>f4').tofile(v_output_file)



