
import os
import numpy as np
import netCDF4 as nc4
import ast
import matplotlib.pyplot as plt
import sys

def read_compact_to_stitched_grid(file_path,sNx,sNy,Nr,compact_tile_size,
                                  ordered_nonblank_tiles, tile_face_index_dict):

    compact_grid = np.fromfile(file_path,'>f4')
    n_rows = int(np.size(compact_grid)/(Nr*compact_tile_size))
    compact_grid = np.reshape(compact_grid,(Nr,n_rows,compact_tile_size))

    grid_faces = {}

    face_1 = compact_grid[:, :3 * sNx, :]
    face_1 = np.reshape(face_1, (Nr, sNx, 3 * sNx))
    grid_faces[1] = face_1
    # print('face_1', np.shape(face_1))

    # plt.imshow(face_1[0,:,:],origin='lower')
    # plt.show()

    # face_3 = np.rot90(compact_grid[:, 3 * sNx:, :], axes=(2, 1))
    face_3 = compact_grid[:, 3 * sNx:, :]
    grid_faces[3] = face_3
    # print('face_3',np.shape(face_3))

    # plt.imshow(face_3[0,:,:], origin='lower')
    # plt.show()


    stitched_grid = np.zeros((Nr,sNy*len(ordered_nonblank_tiles), sNy*len(ordered_nonblank_tiles[0])))
    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            face = tile_face_index_dict[tile_number][0]
            row = tile_face_index_dict[tile_number][1]
            col = tile_face_index_dict[tile_number][2]
            rotations = tile_face_index_dict[tile_number][3]

            tile_grid = grid_faces[face][:,row:row+sNy,col:col+sNx]
            for i in range(rotations):
                tile_grid = np.rot90(tile_grid,axes=(2,1))


            stitched_grid[:,r*sNy:(r+1)*sNy,c*sNx:(c+1)*sNx] =tile_grid

    return(stitched_grid)

def read_grid_tile_angles(config_dir,model_name,sNx,sNy,ordered_nonblank_tiles,ordered_nonblank_rotations):

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L05', model_name, 'run_for_grid')

    AngleCS_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    AngleSN_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            rotations = ordered_nonblank_rotations[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    AngleCS = ds.variables['AngleCS'][:, :]
                    AngleSN = ds.variables['AngleSN'][:, :]
                    ds.close()

                    for i in range(rotations):
                        AngleCS = np.rot90(AngleCS)
                        AngleSN = np.rot90(AngleSN)

                    AngleCS_grid[r*sNy:(r+1)*sNy, c*sNx:(c+1)*sNx] = AngleCS
                    AngleSN_grid[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = AngleSN

    return(AngleCS_grid, AngleSN_grid)

def rotate_velocity_vectors_to_domain(angle_cos, angle_sin, zonal_vel, meridional_vel):
    uvel = np.zeros_like(zonal_vel)
    vvel = np.zeros_like(meridional_vel)
    for k in range(np.shape(uvel)[0]):
        uvel[k, :, :] = angle_cos * zonal_vel[k, :, :] + angle_sin * meridional_vel[k, :, :]
        vvel[k, :, :] = -1 * angle_sin * zonal_vel[k, :, :] + angle_cos * meridional_vel[k, :, :]
    return (uvel, vvel)

def store_grid_as_compact(output_file,output_grid,compact_tile_size,
                          sNx, sNy, Nr, faces,face_shapes,
                          ordered_nonblank_tiles,tile_face_index_dict):

    face_grids = {}
    for f in range(len(faces)):
        face = faces[f]
        face_grids[face] = np.zeros((Nr,face_shapes[f][0], face_shapes[f][1]))

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            tile_face = tile_face_index_dict[tile_number][0]
            min_row = tile_face_index_dict[tile_number][1]
            min_col = tile_face_index_dict[tile_number][2]
            tile_subset = output_grid[:,r*sNy:(r+1)*sNy,c*sNx:(c+1)*sNx]
            face_grids[tile_face][:,min_row:min_row+sNy,min_col:min_col+sNx] = tile_subset

    for f in range(len(faces)):
        face = faces[f]
        grid = face_grids[face]

        # plt.imshow(grid[0,:,:])
        # plt.show()

        n_rows = int((np.shape(grid)[1] * np.shape(grid)[2]) / (compact_tile_size))
        grid = np.reshape(grid, (Nr, n_rows, compact_tile_size))
        if face == 1:
            compact_stack = grid
        else:
            compact_stack = np.concatenate([compact_stack, grid], axis=1)

        print(np.shape(compact_stack))

    compact_stack.ravel(order='C').astype('>f4').tofile(output_file)


def rotate_ICs(config_dir,model_name,sNx,sNy,faces, face_shapes,ordered_nonblank_tiles,ordered_nonblank_rotations,tile_face_index_dict):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import aste_functions as af

    print('    - Rotating the UVEL, VVEL, SIuice, and SIvice BCs for the '+model_name+' model from ASTE data')

    # step 0: get the model domain
    print('    - Reading in the model geometry')
    AngleCS, AngleSN = read_grid_tile_angles(config_dir,model_name,sNx,sNy,
                                             ordered_nonblank_tiles,ordered_nonblank_rotations)

    pairs = [['UVEL','VVEL'],['SIuice', 'SIvice']]
    pairs = [['UVEL', 'VVEL']]

    compact_tile_size = 90

    for pair in pairs:

        if pair[0]=='UVEL':
            Nr = 90
            print('     - Rotating VEL fields')
        else:
            Nr = 1
            print('     - Rotating SI fields')

        u_file = os.path.join(config_dir, 'L05', model_name, 'input', 'L05_IC_'+pair[0]+'_rotated.bin')
        u_grid = read_compact_to_stitched_grid(u_file,sNx,sNy,Nr,compact_tile_size,
                                  ordered_nonblank_tiles, tile_face_index_dict)

        v_file = os.path.join(config_dir, 'L05', model_name, 'input', 'L05_IC_'+pair[1]+'_rotated.bin')
        v_grid = read_compact_to_stitched_grid(v_file, sNx, sNy, Nr, compact_tile_size,
                                               ordered_nonblank_tiles, tile_face_index_dict)

        # plt.imshow(v_grid[0,:,:],origin='lower',vmin=-0.5,vmax=0.5,cmap='seismic')
        # plt.show()

        rotated_u, rotated_v = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                       u_grid, v_grid)

        # plt.imshow(rotated_v[0, :, :], origin='lower', vmin=-0.5, vmax=0.5, cmap='seismic')
        # plt.show()

        u_output_file = os.path.join(config_dir, 'L05', model_name, 'input', 'L05_IC_' + pair[0] + '.bin')
        store_grid_as_compact(u_output_file, rotated_u, compact_tile_size,
                              sNx, sNy, Nr, faces, face_shapes,
                              ordered_nonblank_tiles, tile_face_index_dict)

        v_output_file = os.path.join(config_dir, 'L05', model_name, 'input', 'L05_IC_' + pair[1] + '.bin')
        store_grid_as_compact(v_output_file, rotated_v, compact_tile_size,
                              sNx, sNy, Nr, faces, face_shapes,
                              ordered_nonblank_tiles, tile_face_index_dict)

        # u_output_file = os.path.join(config_dir, 'L05', model_name, 'input', 'L05_IC_'+pair[0]+'.bin')
        # rotated_u.ravel(order='C').astype('>f4').tofile(u_output_file)
        #
        # print('       - Output file shape: '+str(np.shape(rotated_u)))
        #
        # v_output_file = os.path.join(config_dir, 'L05', model_name, 'input', 'L05_IC_'+pair[1]+'.bin')
        # rotated_v.ravel(order='C').astype('>f4').tofile(v_output_file)



