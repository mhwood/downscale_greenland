
import os
import numpy as np
import netCDF4 as nc4
import ast
import matplotlib.pyplot as plt
import sys

def read_grid_geometry(config_dir,model_name, n_rows, n_cols, boundary):
    file_path = os.path.join(config_dir, 'mitgrids', model_name + '.mitgrid')
    entire_grid = np.fromfile(file_path, dtype='>f8')
    entire_grid = np.reshape(entire_grid, (16, n_rows + 1, n_cols + 1))
    XC = entire_grid[0, :, :]
    YC = entire_grid[1, :, :]
    XC = XC[:-1, :-1]
    YC = YC[:-1, :-1]

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

    return(XC, YC, delR)

def read_aste_data_to_stiched_grid(aste_dir,var_name,ordered_aste_tiles,ordered_aste_tile_rotations, L_XC, L_YC):

    aste_Nr = 50
    aste_sNx = 90
    aste_sNy = 90
    aste_timesteps = 192
    aste_XC = np.zeros((aste_sNy*len(ordered_aste_tiles),aste_sNx*len(ordered_aste_tiles[0])))
    aste_YC = np.zeros((aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))
    aste_grid = np.zeros((aste_timesteps,aste_Nr,aste_sNy*len(ordered_aste_tiles),aste_sNx*len(ordered_aste_tiles[0])))
    aste_hfacC = np.zeros((aste_Nr, aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))

    for r in range(len(ordered_aste_tiles)):
        aste_tile_row = ordered_aste_tiles[r]
        aste_rotation_row = ordered_aste_tile_rotations[r]
        for c in range(len(ordered_aste_tiles[r])):

            # get the variable grid
            file_name = os.path.join(aste_dir,var_name,var_name+'.'+'{:04d}'.format(aste_tile_row[c]))+'.nc'
            ds = nc4.Dataset(file_name)
            var_grid = ds.variables[var_name][:,:,:,:]
            XC = ds.variables['lon'][:,:]
            YC = ds.variables['lat'][:, :]
            ds.close()

            # get the hfac grid
            file_name = os.path.join(aste_dir, 'nctiles_grid', 'GRID.' + '{:04d}'.format(aste_tile_row[c])) + '.nc'
            ds = nc4.Dataset(file_name)
            hfac_grid = ds.variables['hFacC'][:, :, :]
            ds.close()

            # rotate things as necessary
            for n in range(aste_rotation_row[c]):
                var_grid = np.rot90(var_grid, axes=(2, 3))
                hfac_grid = np.rot90(hfac_grid, axes=(1, 2))
                XC = np.rot90(XC)
                YC = np.rot90(YC)

            # put it into the big grid
            aste_grid[:,:, r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = var_grid
            aste_hfacC[:, r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = hfac_grid
            aste_XC[r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = XC
            aste_YC[r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = YC

    # C = plt.imshow(aste_grid[0,:,:],origin='lower')
    # C = plt.imshow(aste_XC, origin='lower')
    # C = plt.imshow(aste_YC, origin='lower')
    # plt.colorbar(C)
    # plt.show()

    ll_dist = ((L_XC[0, 0] - aste_XC) ** 2 + (L_YC[0, 0] - aste_YC) ** 2) ** 0.5
    ll_row, ll_col = np.where(ll_dist == np.min(ll_dist))
    ul_dist = ((L_XC[-1, 0] - aste_XC) ** 2 + (L_YC[-1, 0] - aste_YC) ** 2) ** 0.5
    ul_row, ul_col = np.where(ul_dist == np.min(ul_dist))

    lr_dist = ((L_XC[0, -1] - aste_XC) ** 2 + (L_YC[0, -1] - aste_YC) ** 2) ** 0.5
    lr_row, lr_col = np.where(lr_dist == np.min(lr_dist))
    ur_dist = ((L_XC[-1, -1] - aste_XC) ** 2 + (L_YC[-1, -1] - aste_YC) ** 2) ** 0.5
    ur_row, ur_col = np.where(ur_dist == np.min(ur_dist))

    min_row = np.min([ll_row[0], ul_row[0], lr_row[0], ur_row[0]])
    max_row = np.max([ll_row[0], ul_row[0], lr_row[0], ur_row[0]])
    min_col = np.min([ll_col[0], ul_col[0], lr_col[0], ur_col[0]])
    max_col = np.max([ll_col[0], ul_col[0], lr_col[0], ur_col[0]])

    dist_buffer = 3
    aste_XC = aste_XC[min_row - dist_buffer:max_row + dist_buffer,
              min_col - dist_buffer:max_col + dist_buffer]
    aste_YC = aste_YC[min_row - dist_buffer:max_row + dist_buffer,
              min_col - dist_buffer:max_col + dist_buffer]
    aste_grid = aste_grid[:,:, min_row - dist_buffer:max_row + dist_buffer,
                min_col - dist_buffer:max_col + dist_buffer]
    aste_hfacC = aste_hfacC[:, min_row - dist_buffer:max_row + dist_buffer,
                min_col - dist_buffer:max_col + dist_buffer]

    aste_delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                        10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                        31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                        93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                        139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                        341.50, 364.50, 387.50, 410.50, 433.50, 456.50])

    return(aste_XC,aste_YC,aste_grid,aste_hfacC,aste_delR)

def read_mask_from_nc(nc_file, hFac='C'):
    ds = nc4.Dataset(nc_file)
    mask = ds.variables['wet_grid_'+hFac][:,:,:]
    ds.close()
    return(mask)

def read_mask_from_grid_nc(nc_file, hFac='C'):
    ds = nc4.Dataset(nc_file)
    mask = ds.variables['HFac'+hFac][:,:,:]
    mask[mask>0]=1
    mask[mask<=0]=0
    ds.close()
    return(mask)

def downscale_ASTE_boundary_field(df,boundary, var_name,
                                  aste_XC, aste_YC, aste_grid, ASTE_wet_cells, ASTE_wet_cells_on_domain_3D,
                                  XC, YC, domain_wet_cells_3D):

    nTimestepsOut = np.shape(aste_grid)[0]

    if boundary=='south':
        ASTE_wet_cells_on_domain_3D_subset = ASTE_wet_cells_on_domain_3D[:,:1, :]
        domain_wet_cells_3D = domain_wet_cells_3D[:,:1 , :]
        if var_name=='ETAN':
            L1_boundary_var = np.zeros((nTimestepsOut, np.shape(YC)[0], np.shape(YC)[1]))
        else:
            L1_boundary_var = np.zeros((nTimestepsOut, np.shape(domain_wet_cells_3D)[0], np.shape(YC)[0], np.shape(YC)[1]))

    if boundary=='north':
        ASTE_wet_cells_on_domain_3D_subset = ASTE_wet_cells_on_domain_3D[:, -1:, :]
        domain_wet_cells_3D = domain_wet_cells_3D[:, -1:, :]
        if var_name=='ETAN':
            L1_boundary_var = np.zeros((nTimestepsOut, np.shape(YC)[0], np.shape(YC)[1]))
        else:
            L1_boundary_var = np.zeros((nTimestepsOut, np.shape(domain_wet_cells_3D)[0], np.shape(YC)[0], np.shape(YC)[1]))

    if boundary=='west':
        ASTE_wet_cells_on_domain_3D_subset = ASTE_wet_cells_on_domain_3D[:, :, :1]
        domain_wet_cells_3D = domain_wet_cells_3D[:, :, :1]
        if var_name=='ETAN':
            L1_boundary_var = np.zeros((nTimestepsOut, np.shape(YC)[0], np.shape(YC)[1]))
        else:
            L1_boundary_var = np.zeros((nTimestepsOut, np.shape(domain_wet_cells_3D)[0], np.shape(YC)[0], np.shape(YC)[1]))

    if boundary=='east':
        ASTE_wet_cells_on_domain_3D_subset = ASTE_wet_cells_on_domain_3D[:, :, -1:]
        domain_wet_cells_3D = domain_wet_cells_3D[:, :, -1:]
        if var_name=='ETAN':
            L1_boundary_var = np.zeros((nTimestepsOut, np.shape(YC)[0], np.shape(YC)[1]))
        else:
            L1_boundary_var = np.zeros((nTimestepsOut, np.shape(domain_wet_cells_3D)[0], np.shape(YC)[0], np.shape(YC)[1]))

    print('        -  Variable shapes entering downscale routine:')
    print('          -  aste_XC: '+str(np.shape(aste_XC)))
    print('          -  aste_YC: ' + str(np.shape(aste_YC)))
    print('          -  aste_grid: ' + str(np.shape(aste_grid)))
    print('          -  ASTE_wet_cells: ' + str(np.shape(ASTE_wet_cells)))
    print('          -  ASTE_wet_cells_on_domain_3D_subset: ' + str(np.shape(ASTE_wet_cells_on_domain_3D_subset)))
    print('          -  XC: ' + str(np.shape(XC)))
    print('          -  YC: ' + str(np.shape(YC)))
    print('          -  domain_wet_cells_3D: ' + str(np.shape(domain_wet_cells_3D)))

    # if np.shape(domain_wet_cells_3D_subset)[-1]==np.shape(L1_XC)[-1]+1:
    #     domain_wet_cells_3D_subset = domain_wet_cells_3D_subset[:,:,:-1]
    #
    # if np.shape(domain_wet_cells_3D_subset)[-2]==np.shape(L1_XC)[-2]+1:
    #     domain_wet_cells_3D_subset = domain_wet_cells_3D_subset[:,:-1,:]
    #
    # for timestep in range(nTimestepsOut):
    #     if timestep%50==0:
    #         print('          Working on timestep '+str(timestep)+' of '+str(nTimestepsOut)+' for '+var_name+' on mask '+boundary)
    #     if var_name=='ETAN':
    #         downscaled_field = df.downscale_2D_field(aste_XC, aste_YC, aste_grid[timestep,:,:],
    #                                                  ASTE_wet_cells[0,:,:], ASTE_wet_cells_on_domain_3D_subset[0,:,:],
    #                                                  L1_XC_subset, YC, domain_wet_cells_3D_subset[0,:,:])
    #         L1_boundary_var[timestep, :, :] = downscaled_field
    #     else:
    #         # if var_name in ['UVEL','VVEL']:
    #         #     remove_zeros=False
    #         # else:
    #         remove_zeros = True
    #         downscaled_field = df.downscale_3D_field(aste_XC, aste_YC, aste_grid[timestep, :, :, :],
    #                                                  ASTE_wet_cells, ASTE_wet_cells_on_domain_3D_subset,
    #                                                  L1_XC_subset, YC, domain_wet_cells_3D_subset,remove_zeros=remove_zeros)
    #         # if var_name in ['UVEL','VVEL'] and relax_velocity:
    #         #     if timestep<440:
    #         #         relaxation_factor = timestep/440
    #         #         downscaled_field *= relaxation_factor
    #         L1_boundary_var[timestep, :, :, :] = downscaled_field


    # if var_name == 'ETAN':
    #     plt.subplot(2,1,1)
    #     C = plt.imshow(aste_grid[0,:,:])
    #     # plt.colorbar(C)
    #     plt.title(var_name)
    #     plt.subplot(2,1,2)
    #     plt.imshow(L1_boundary_var[0,:,:])
    #     plt.show()
    # else:
    #     plt.subplot(2, 1, 1)
    #     C = plt.imshow(aste_grid[0, 0, :, :])
    #     # plt.colorbar(C)
    #     plt.title(var_name)
    #     plt.subplot(2, 1, 2)
    #     plt.imshow(L1_boundary_var[0, 0, :, :])
    #     plt.show()

    return(L1_boundary_var)


def create_L3_ASTE_boundary_condition(config_dir,model_name,var_name,boundary,ordered_aste_tiles,ordered_aste_tile_rotations):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df

    print('    - Creating the '+var_name+' boundary conditions on the '+boundary+' boundary for the '+model_name+' model from ASTE data')

    f = open(os.path.join(config_dir, 'domain_sizes.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)
    L_size = size_dict[model_name]
    n_rows = L_size[0]
    n_cols = L_size[1]

    # step 0: get the model domain
    XC, YC, delR = read_grid_geometry(config_dir,model_name, n_rows, n_cols, boundary)

    # step 1: create a stitched grid around the domain using ASTE data
    aste_dir = '/Users/michwood/Documents/Research/Data Repository/Greenland/Ocean Properties/Models/ASTE'
    aste_XC, aste_YC, aste_grid, aste_hfacC, aste_delR = read_aste_data_to_stiched_grid(aste_dir, var_name,
                                                                                       ordered_aste_tiles, ordered_aste_tile_rotations, XC, YC)

    ASTE_wet_cells = np.copy(aste_hfacC)
    ASTE_wet_cells[ASTE_wet_cells>0]=1

    if boundary=='north':
        XC = XC[-1:,:]
        YC = YC[-1:, :]
    elif boundary=='south':
        XC = XC[:1,:]
        YC = YC[:1, :]
    elif boundary=='west':
        XC = XC[:,:1]
        YC = YC[:, :1]
    elif boundary=='east':
        XC = XC[:, -1:]
        YC = YC[:, -1:]
    else:
        raise ValueError('Boundary not recognized')

    ASTE_grid_on_domain_file = os.path.join(config_dir, 'L3', model_name, 'input',
                                            'ASTE_wetgrid_on_' + model_name + '.nc')
    domain_grid_file = os.path.join(config_dir, 'L3', model_name, 'input', model_name + '_grid.nc')

    if var_name in ['VVEL']:
        ASTE_wet_cells_on_domain_3D = read_mask_from_nc(ASTE_grid_on_domain_file, hFac='S')
        domain_wet_cells_3D = read_mask_from_grid_nc(domain_grid_file, hFac='S')
        domain_wet_cells_3D = domain_wet_cells_3D[:, 1:, :]
    elif var_name in ['UVEL']:
        ASTE_wet_cells_on_domain_3D = read_mask_from_nc(ASTE_grid_on_domain_file, hFac='W')
        domain_wet_cells_3D = read_mask_from_grid_nc(domain_grid_file, hFac='W')
        domain_wet_cells_3D = domain_wet_cells_3D[:, :, 1:]
    else:
        ASTE_wet_cells_on_domain_3D = read_mask_from_nc(ASTE_grid_on_domain_file, hFac='C')
        domain_wet_cells_3D = read_mask_from_grid_nc(domain_grid_file, hFac='C')


    subset_copy = np.copy(aste_grid)

    print('          - Shape boundary variable subset: ' + str(np.shape(aste_grid)))

    print('    - Interpolating vertically')
    aste_grid, ASTE_wet_cells = df.interpolate_var_grid_faces_to_new_depth_levels(
        aste_grid, ASTE_wet_cells, aste_delR, delR)

    print('          - Shape boundary variable subset: ' + str(np.shape(aste_grid)))
    # print('          - Shape L0 grid subset: ' + str(np.shape(aste_XC)))

    # ASTE_wet_cells = np.copy(L0_wet_grid_3D_vertically_interpolated)
    # min_row, min_col = np.where(np.logical_and(L0_grid_XC == aste_XC[0, 0], L0_grid_YC == aste_YC[0, 0]))
    # max_row, max_col = np.where(np.logical_and(L0_grid_XC == aste_XC[-1, -1], L0_grid_YC == aste_YC[-1, -1]))
    # ASTE_wet_cells = ASTE_wet_cells[:, min_row[0]:max_row[0] + 1, min_col[0]:max_col[0] + 1]
    #
    # print('          - Shape: ' + str(np.shape(ASTE_wet_cells)))

    print('    - Downscaling the output to the new boundary')
    L1_boundary_var = downscale_ASTE_boundary_field(df, boundary, var_name,
                                                   aste_XC, aste_YC, aste_grid,
                                                   ASTE_wet_cells, ASTE_wet_cells_on_domain_3D,
                                                   XC, YC, domain_wet_cells_3D)



    # # if 'north' in mask_name or 'south' in mask_name:
    # #     for time_step in range(5):#np.shape(L1_boundary_var)):
    # #         plt.subplot(1,2,1)
    # #         plt.imshow(subset_copy[time_step,:,1,:],vmin=-0.1,vmax=0.1,cmap='seismic_r')
    # #         plt.title('L0')
    # #         plt.subplot(1, 2, 2)
    # #         plt.imshow(L1_boundary_var[time_step, :, 0, :],vmin=-0.1,vmax=0.1,cmap='seismic_r')
    # #         plt.title('L1 (nan values: '+str(np.sum(np.isnan(L1_boundary_var[time_step,:,0,:]))))
    # #         plt.show()
    # # else:
    # #     plt.subplot(1, 2, 1)
    # #     plt.imshow(subset_copy[0, :, :, 1])
    # #     plt.subplot(1, 2, 2)
    # #     plt.imshow(L1_boundary_var[0, :, :, 0])
    # #     plt.show()
    #
    # output_obcs_field(output_dir, dest_file, mask_name, boundary_var_name, L1_boundary_var)






