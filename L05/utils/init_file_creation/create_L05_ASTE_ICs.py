
import os
import numpy as np
import netCDF4 as nc4
import ast
import matplotlib.pyplot as plt
from MITgcmutils import mds
from scipy.interpolate import griddata
import sys


def read_grid_tile_geometry(config_dir,model_name,var_name,ordered_nonblank_tiles,tile_face_index_dict):
    ordered_XC_tiles = []
    ordered_YC_tiles = []
    ordered_AngleCS_tiles = []
    ordered_AngleSN_tiles = []
    ordered_hfac_tiles = []

    delR = np.array([1.00, 1.14, 1.30, 1.49, 1.70,
                     1.93, 2.20, 2.50, 2.84, 3.21,
                     3.63, 4.10, 4.61, 5.18, 5.79,
                     6.47, 7.20, 7.98, 8.83, 9.73,
                     10.69, 11.70, 12.76, 13.87, 15.03,
                     16.22, 17.45, 18.70, 19.97, 21.27,
                     22.56, 23.87, 25.17, 26.46, 27.74,
                     29.00, 30.24, 31.45, 32.65, 33.82,
                     34.97, 36.09, 37.20, 38.29, 39.37,
                     40.45, 41.53, 42.62, 43.73, 44.87,
                     46.05, 47.28, 48.56, 49.93, 51.38,
                     52.93, 54.61, 56.42, 58.38, 60.53,
                     62.87, 65.43, 68.24, 71.33, 74.73,
                     78.47, 82.61, 87.17, 92.21, 97.79,
                     103.96, 110.79, 118.35, 126.73, 136.01,
                     146.30, 157.71, 170.35, 184.37, 199.89,
                     217.09, 236.13, 257.21, 280.50, 306.24,
                     334.64, 365.93, 400.38, 438.23, 479.74, ])

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L05', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        row_XCs = []
        row_YCs = []
        row_AngleCSs = []
        row_AngleSNs = []
        row_HFacs = []
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            tile_face = tile_face_index_dict[tile_number][0]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
                    AngleCS = ds.variables['AngleCS'][:, :]
                    AngleSN = ds.variables['AngleSN'][:, :]
                    if tile_face < 3:
                        if var_name == 'VVEL':
                            hFac = 'S'
                        elif var_name == 'UVEL':
                            hFac = 'W'
                        else:
                            hFac = 'C'
                    else:
                        if var_name == 'VVEL':
                            hFac = 'W'
                        elif var_name == 'UVEL':
                            hFac = 'S'
                        else:
                            hFac = 'C'
                    HFac = ds.variables['HFac' + hFac][:, :, :]
                    if hFac == 'W':
                        HFac = HFac[:,:,1:]
                    if hFac == 'S':
                        HFac = HFac[:,1:,:]
                    ds.close()
                    row_XCs.append(XC)
                    row_YCs.append(YC)
                    row_AngleCSs.append(AngleCS)
                    row_AngleSNs.append(AngleSN)
                    row_HFacs.append(HFac)
        ordered_XC_tiles.append(row_XCs)
        ordered_YC_tiles.append(row_YCs)
        ordered_AngleCS_tiles.append(row_AngleCSs)
        ordered_AngleSN_tiles.append(row_AngleSNs)
        ordered_hfac_tiles.append(row_HFacs)

    return(ordered_XC_tiles, ordered_YC_tiles, ordered_AngleCS_tiles, ordered_AngleSN_tiles, ordered_hfac_tiles, delR)

def subset_aste_to_tile_region(XC_subset, YC_subset,
                                   aste_XC, aste_YC, aste_hfacC, aste_grid):

    ll_dist = ((aste_XC-np.min(XC_subset))**2 + (aste_YC-np.min(YC_subset))**2)**0.5
    ll_index_row, ll_index_col = np.where(ll_dist == np.min(ll_dist))

    ul_dist = ((aste_XC - np.min(XC_subset)) ** 2 + (aste_YC - np.max(YC_subset)) ** 2) ** 0.5
    ul_index_row, ul_index_col = np.where(ul_dist == np.min(ul_dist))

    lr_dist = ((aste_XC - np.max(XC_subset)) ** 2 + (aste_YC - np.min(YC_subset)) ** 2) ** 0.5
    lr_index_row, lr_index_col = np.where(lr_dist == np.min(lr_dist))

    ur_dist = ((aste_XC - np.max(XC_subset)) ** 2 + (aste_YC - np.max(YC_subset)) ** 2) ** 0.5
    ur_index_row, ur_index_col = np.where(ur_dist == np.min(ur_dist))


    min_y_index = np.min([ll_index_row[0],ul_index_row[0],lr_index_row[0],ur_index_row[0]])
    max_y_index = np.max([ll_index_row[0], ul_index_row[0], lr_index_row[0], ur_index_row[0]])
    min_x_index = np.min([ll_index_col[0], ul_index_col[0], lr_index_col[0], ur_index_col[0]])
    max_x_index = np.max([ll_index_col[0], ul_index_col[0], lr_index_col[0], ur_index_col[0]])

    buffer = 1
    if min_x_index - buffer>=0:
        min_x_index -= buffer
    if max_x_index + buffer<= np.shape(aste_XC)[1]-1:
        max_x_index += buffer
    if min_y_index - buffer>=0:
        min_y_index -= buffer
    if max_y_index + buffer<= np.shape(aste_XC)[0]-1:
        max_y_index += buffer

    buffer = 1
    if min_x_index - buffer >= 0:
        min_x_index -= buffer
    if max_x_index + buffer <= np.shape(aste_XC)[1] - 1:
        max_x_index += buffer
    if min_y_index - buffer >= 0:
        min_y_index -= buffer
    if max_y_index + buffer <= np.shape(aste_XC)[0] - 1:
        max_y_index += buffer

    aste_XC_subset = aste_XC[min_y_index:max_y_index,min_x_index:max_x_index]
    aste_YC_subset = aste_YC[min_y_index:max_y_index,min_x_index:max_x_index]
    aste_hfacC_subset = aste_hfacC[:, min_y_index:max_y_index, min_x_index:max_x_index]
    aste_grid_subset = aste_grid[:, min_y_index:max_y_index, min_x_index:max_x_index]

    # plt.plot(aste_XC_subset,aste_YC_subset,'k.')
    # plt.plot(XC_subset, YC_subset, 'b.')
    # plt.show()

    return(aste_XC_subset, aste_YC_subset, aste_hfacC_subset, aste_grid_subset)

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

def create_L05_ICs(config_dir,model_name,var_name,sNx,sNy,faces,face_shapes,
                   ordered_nonblank_tiles,tile_face_index_dict,
                   ordered_aste_tiles, ordered_aste_tile_rotations):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import aste_functions as af

    compact_tile_size = 90

    print('    - Creating the '+var_name+' IC file for the '+model_name+' model from ASTE data')

    # step 0: get the model domain
    print('    - Reading in the model geometry')
    ordered_XC_tiles, ordered_YC_tiles, ordered_AngleCS_tiles, ordered_AngleSN_tiles, ordered_hfac_tiles, delR = \
        read_grid_tile_geometry(config_dir,model_name,var_name,ordered_nonblank_tiles,tile_face_index_dict)

    # step 1: get the aste faces geometry
    aste_dir = '/Users/michwood/Documents/Research/Data Repository/Greenland/Ocean Properties/Models/ASTE'
    print('    - Reading in the aste geometry')
    aste_XC, aste_YC, aste_AngleCS, aste_AngleSN, aste_hfacC, aste_delR = \
        af.read_aste_grid_geometry(aste_dir, ordered_aste_tiles, ordered_aste_tile_rotations)

    print('    - Reading in the aste grid for variable ' + str(var_name))
    if var_name in ['UVEL','VVEL','SIuice','SIvice']:
        aste_grid = af.read_aste_field_to_stiched_grid(aste_dir, var_name, ordered_aste_tiles,
                                                       ordered_aste_tile_rotations,
                                                       aste_AngleCS, aste_AngleSN, rotate_velocity = True)
    else:
        aste_grid = af.read_aste_field_to_stiched_grid(aste_dir, var_name, ordered_aste_tiles,
                                                       ordered_aste_tile_rotations)

    # get the initial conditions
    aste_grid = aste_grid[0,:,:,:]

    # plt.imshow(aste_grid[0,0,:,:],origin='lower',cmap='seismic',vmin=-1,vmax=1)
    # plt.show()

    if np.shape(aste_grid)[0]>1:
        Nr = np.size(delR)
    else:
        Nr = 1

    output_grid = np.zeros((Nr, sNy * len(ordered_aste_tiles), sNx*len(ordered_aste_tiles[0])))

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            print('      - Downscaling tile '+str(tile_number))

            tile_XC = ordered_XC_tiles[r][c]
            tile_YC = ordered_YC_tiles[r][c]
            tile_hfac = ordered_hfac_tiles[r][c]

            print('          - Subsetting the ASTE data to the tile domain')
            # subset aste to the region around the tile to run a quicker script
            aste_XC_subset, aste_YC_subset, aste_hfacC_subset, aste_grid_subset = \
                subset_aste_to_tile_region(tile_XC, tile_YC,
                                               aste_XC, aste_YC, aste_hfacC, aste_grid)

            # use the hfacs to make the mask
            print('          - Making the masks')
            aste_wet_cells_subset = np.copy(aste_hfacC_subset)
            aste_wet_cells_subset[aste_wet_cells_subset>0]=1
            wetgrid = np.copy(tile_hfac)
            wetgrid[wetgrid > 0] = 1
            aste_wet_cells_on_domain = np.copy(wetgrid)

            # interpolate to the new depths of this model
            print('          - Interpolating the depths to the new depths of the model')
            if Nr>1:
                aste_grid_subset, aste_wet_cells_subset = df.interpolate_var_grid_faces_to_new_depth_levels(
                    aste_grid_subset, aste_wet_cells_subset, aste_delR, delR)

            print('          - Running the downscale routine')
            interp_field = df.downscale_3D_field(aste_XC_subset, aste_YC_subset,
                                                 aste_grid_subset, aste_wet_cells_subset,
                                                 aste_wet_cells_on_domain,
                                                 tile_XC, tile_YC, wetgrid,
                                                 mean_vertical_difference=0, fill_downward=True, remove_zeros=True,
                                                 printing=False)

            # plt.subplot(2, 2, 1)
            # plt.imshow(aste_grid_subset[0,:,:])
            # plt.subplot(2, 2, 2)
            # plt.imshow(aste_wet_cells_subset[0, :, :])
            # plt.subplot(2, 2, 3)
            # plt.imshow(interp_field[0,:,:])
            # plt.subplot(2, 2, 4)
            # plt.imshow(wetgrid[0, :, :])
            # plt.show()

            output_grid[:,r*sNy:(r+1)*sNy, c*sNx:(c+1)*sNx] = interp_field

    if var_name in ['UVEL','VVEL','SIuice','SIvice']:
        output_file = os.path.join(config_dir,'L05',model_name,'input',
                                   'L05_IC_'+var_name+'_rotated.bin')
    else:
        output_file = os.path.join(config_dir, 'L05', model_name, 'input',
                                   'L05_IC_' + var_name + '.bin')

    store_grid_as_compact(output_file, output_grid, compact_tile_size,
                          sNx, sNy, Nr, faces, face_shapes,
                          ordered_nonblank_tiles, tile_face_index_dict)




