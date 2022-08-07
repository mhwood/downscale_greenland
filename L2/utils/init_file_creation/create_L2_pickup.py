
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
    for year in range(2002,2003):
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

    f = open(os.path.join(config_dir,'L1', L1_model_name, 'run','dv', 'daily_bc_dest_ref.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)

    dest_files_out = []
    source_file_read_dict = {}
    source_file_read_index_sets = {}

    for dest_file in dest_files:
        dest_files_out.append('L2_BC_'+mask_name+'_'+var_name+'.'+dest_file+'.bin')
        source_files = []
        index_set = []
        for pair in size_dict[dest_file]:
            source_files.append(mask_name+'_'+var_name+'.'+pair[0]+'.bin')
            index_set.append(pair[1])
        source_file_read_dict['L2_BC_'+mask_name+'_'+var_name+'.'+dest_file+'.bin'] = source_files
        source_file_read_index_sets['L2_BC_'+mask_name+'_'+var_name+'.'+dest_file+'.bin'] = index_set

    return(dest_files_out, source_file_read_dict, source_file_read_index_sets)


def read_L1_grid_tile_geometry(config_dir,model_name, sNx, sNy, Nr,
                               ordered_nonblank_tiles, ordered_nonblank_rotations,
                               faces, ordered_tiles_faces_dict):

    XC_faces = {}
    YC_faces = {}
    HFac_faces = {}
    for face in faces:
        face_XC_grid = np.zeros((sNy*len(ordered_tiles_faces_dict[face]),
                                 sNx*len(ordered_tiles_faces_dict[face][0])))
        face_YC_grid = np.zeros((sNy * len(ordered_tiles_faces_dict[face]),
                                 sNx * len(ordered_tiles_faces_dict[face][0])))
        face_HFac_grid = np.zeros((Nr,
                                   sNy * len(ordered_tiles_faces_dict[face]),
                                   sNx * len(ordered_tiles_faces_dict[face][0])))
        XC_faces[face] = face_XC_grid
        YC_faces[face] = face_YC_grid
        HFac_faces[face] = face_HFac_grid

    # stitched_XC_grid = np.zeros((sNy*len(ordered_tiles_faces_dict), sNx*len(ordered_nonblank_tiles[1])))
    # stitched_YC_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[1])))

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L1', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
                    HFac = ds.variables['HFacC'][:,:,:]
                    ds.close()

                    for face in faces:
                        for fr in range(len(ordered_tiles_faces_dict[face])):
                            for fc in range(len(ordered_tiles_faces_dict[face][0])):
                                if ordered_tiles_faces_dict[face][fr][fc] == tile_number:
                                    XC_faces[face][fr * sNy:(fr + 1) * sNy, fc * sNx: (fc + 1) * sNx] = XC
                                    YC_faces[face][fr * sNy:(fr + 1) * sNy, fc * sNx: (fc + 1) * sNx] = YC
                                    HFac_faces[face][:, fr * sNy:(fr + 1) * sNy, fc * sNx: (fc + 1) * sNx] = HFac

                    # for i in range(ordered_nonblank_rotations[r][c]):
                    #     XC = np.rot90(XC)
                    #     YC = np.rot90(YC)
                    #
                    # stitched_XC_grid[r * sNy:(r + 1) * sNy,c * sNx: (c + 1) * sNx] = XC
                    # stitched_YC_grid[r * sNy:(r + 1) * sNy,c * sNx: (c + 1) * sNx] = YC

    return(XC_faces, YC_faces, HFac_faces)


def read_grid_geometry_from_nc(config_dir,model_name,var_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    if var_name=='UVEL':
        hfac = ds.variables['HFacW'][:,:,:]
        hfac = hfac[:,:,:-1]
    elif var_name=='VVEL':
        hfac = ds.variables['HFacS'][:,:,:]
        hfac = hfac[:, 1:, :]
    else:
        hfac = ds.variables['HFacC'][:, :, :]
    ds.close()
    mask = np.copy(hfac)
    mask[mask>0]=1
    return(XC,YC,mask)


def read_mask_reference_from_nc_dict(nc_dict_file,mask_name):
    ds = nc4.Dataset(nc_dict_file)

    grp = ds.groups[mask_name]
    source_faces = grp.variables['source_faces'][:]
    source_rows = grp.variables['source_rows'][:]
    source_cols = grp.variables['source_cols'][:]

    ds.close()
    return(source_faces,source_rows,source_cols)


def read_boundary_variable_to_L1_points(L1_run_dir, mask_name, var_name,
                                        source_files, source_file_read_indices,
                                        source_faces,source_rows,source_cols, Nr_in,
                                        faces, L1_XC_faces, L1_YC_faces, L1_Hfac_faces):

    # n_Timesteps = 0
    # for index_set in source_file_read_indices:
    #     n_Timesteps += int(index_set[1]-index_set[0])
    n_Timesteps = 41

    print('       + the L1 grid for this file will have '+str(n_Timesteps)+' timesteps')

    points = np.zeros((len(source_faces),2))
    if var_name == 'ETAN':
        values = np.zeros((n_Timesteps, len(source_faces)))
        wet_grid = np.zeros((len(source_faces),))
    else:
        values = np.zeros((n_Timesteps, Nr_in, len(source_faces)))
        wet_grid = np.zeros((Nr_in, len(source_faces)))

    for i in range(len(source_faces)):
        x = L1_XC_faces[source_faces[i]][source_rows[i], source_cols[i]]
        y = L1_YC_faces[source_faces[i]][source_rows[i], source_cols[i]]
        points[i,0] = x
        points[i,1] = y
        wet_column = L1_Hfac_faces[source_faces[i]][:, source_rows[i], source_cols[i]]
        wet_column[wet_column>0]=1
        if var_name == 'ETAN':
            wet_grid[i] = wet_column[0]
        else:
            wet_grid[:, i] = wet_column

    index_counter = 0
    for s in range(len(source_files)):
        source_file = source_files[s]
        file_suffix = '.'.join(source_file.split('.')[-2:])
        index_set = source_file_read_indices[s]
        print('         - Adding timesteps ' + str(index_set[0]) + ' to ' + str(index_set[1]) + ' from file ' + source_file)

        start_file_index = int(index_set[0])
        end_file_index = int(index_set[1])

        print('           - Storing at points '+str(index_counter)+' to '+
              str(index_counter+(end_file_index-start_file_index))+' in the grid')

        N = len(source_faces)

        if var_name=='UVEL':
            var_file = os.path.join(L1_run_dir,'dv', 'L2_'+mask_name + '_BC_mask_VVEL.'+file_suffix)
            var_grid = np.fromfile(var_file, dtype='>f4')
            timesteps_in_file = int(np.size(var_grid) / (Nr_in*N))
            var_grid = np.reshape(var_grid, (timesteps_in_file, Nr_in, N))
        elif var_name=='VVEL':
            # var_file = os.path.join(L1_run_dir, 'dv', mask_name + '_i' + str(layer), 'UVEL',
            #                         mask_name + '_i' + str(layer) + '_mask_UVEL.' + file_suffix)
            var_file = os.path.join(L1_run_dir,'dv', 'L2_'+mask_name +'_BC_mask_UVEL.'+file_suffix)
            var_grid = np.fromfile(var_file, dtype='>f4')
            timesteps_in_file = int(np.size(var_grid) / (Nr_in * N))
            var_grid = -1*np.reshape(var_grid, (timesteps_in_file, Nr_in, N))
        else:
            # var_file = os.path.join(L1_run_dir, 'dv', mask_name + '_i' + str(layer), var_name,
            #                         mask_name + '_i' + str(layer) + '_mask_' + var_name + '.' + file_suffix)
            var_file = os.path.join(L1_run_dir,'dv', 'L2_'+mask_name +'_BC_mask_' + var_name +'.'+ file_suffix)
            var_grid = np.fromfile(var_file, dtype='>f4')
            if var_name=='ETAN':
                timesteps_in_file = int(np.size(var_grid) / (N))
                var_grid = np.reshape(var_grid, (timesteps_in_file, N))
            else:
                timesteps_in_file = int(np.size(var_grid) / (Nr_in * N))
                var_grid = np.reshape(var_grid, (timesteps_in_file, Nr_in, N))

        if var_name== 'ETAN':
            values[index_counter:index_counter + (end_file_index-start_file_index),:] = var_grid
        else:
            values[index_counter:index_counter + (end_file_index-start_file_index),:, :] = var_grid


        index_counter += (end_file_index-start_file_index)

    # plt.imshow(var_grid_out[0,0,:,:],origin='lower')
    # plt.show()

    return(points,values,wet_grid)


def downscale_L1_bc_field_to_L2(df, mask_name, var_name,
                                L1_points, L1_values,
                                L1_wet_grid_3D_points, L1_wet_grid_on_L2_3D_points,
                                L2_XC_subset, L2_YC_subset, L2_wet_grid_subset):

    nTimestepsOut = np.shape(L1_values)[0]

    if mask_name=='south' or mask_name=='north':
        if var_name=='ETAN':
            L2_boundary_var = np.zeros((nTimestepsOut, np.shape(L2_YC_subset)[0], np.shape(L2_YC_subset)[1]))
        else:
            L2_boundary_var = np.zeros((nTimestepsOut, np.shape(L2_wet_grid_subset)[0], np.shape(L2_YC_subset)[0], np.shape(L2_YC_subset)[1]))

    if mask_name=='west' or mask_name=='east':
        if var_name=='ETAN':
            L2_boundary_var = np.zeros((nTimestepsOut, np.shape(L2_YC_subset)[0], np.shape(L2_YC_subset)[1]))
        else:
            L2_boundary_var = np.zeros((nTimestepsOut, np.shape(L2_wet_grid_subset)[0], np.shape(L2_YC_subset)[0], np.shape(L2_YC_subset)[1]))

    print('        -  Variable shapes entering downscale routine:')
    print('          -  L1_point: '+str(np.shape(L1_points)))
    print('          -  L1_values: ' + str(np.shape(L1_values)))
    print('          -  L1_wet_grid: ' + str(np.shape(L1_wet_grid_3D_points)))
    # print('          -  L1_wet_grid_on_L2_subset: ' + str(np.shape(L1_wet_grid_on_L2_3D_points)))
    print('          -  L2_XC_subset: ' + str(np.shape(L2_XC_subset)))
    print('          -  L2_YC_subset: ' + str(np.shape(L2_YC_subset)))
    print('          -  L2_wet_grid_subset: ' + str(np.shape(L2_wet_grid_subset)))

    for timestep in range(1):#nTimestepsOut):
        if timestep%50==0:
            print('          Working on timestep '+str(timestep)+' of '+str(nTimestepsOut)+' for '+var_name+' on mask '+mask_name)
        remove_zeros = True
        downscaled_field = df.downscale_3D_points(L1_points, L1_values[timestep, :, :],
                                                 L1_wet_grid_3D_points,
                                                 L2_XC_subset, L2_YC_subset, L2_wet_grid_subset,
                                                  remove_zeros=remove_zeros, printing = False)
        L2_boundary_var[timestep, :, :, :] = downscaled_field


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

    return(L2_boundary_var)


########################################################################################################################


def create_pickup(config_dir, L1_model_name, L2_model_name,
                     sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                     faces, ordered_tiles_faces_dict):

    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import downscale_functions as df

    # this is the dir where the output will be stored
    output_dir = os.path.join(config_dir,'L2',L2_model_name,'input')

    Nr = 90

    delR_in = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                        10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                        31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                        93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                        139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                        341.50, 364.50, 387.50, 410.50, 433.50, 456.50])

    delR_out = np.array([1.00, 1.14, 1.30, 1.49, 1.70,
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


    print('    - Reading in the geometry of the L1 domain')
    L1_XC_faces, L1_YC_faces, L1_Hfac_faces = read_L1_grid_tile_geometry(config_dir,L1_model_name,sNx, sNy, Nr,
                                                                        ordered_nonblank_tiles, ordered_nonblank_rotations,
                                                                        faces, ordered_tiles_faces_dict)

    # for face in faces:
    #     plt.subplot(1,2,1)
    #     C = plt.imshow(L1_XC_faces[face],origin='lower')
    #     plt.colorbar(C)
    #     plt.subplot(1, 2, 2)
    #     C = plt.imshow(L1_YC_faces[face], origin='lower')
    #     plt.colorbar(C)
    #     plt.show()

    L2_XC, L2_YC, L2_wet_grid_3D = read_grid_geometry_from_nc(config_dir,L2_model_name,boundary_var_name)


#     print('    - Reading the mask to reference the variable to the llc grid')
#     nc_dict_file = os.path.join(config_dir,'L1',L1_model_name,'namelist','L1_dv_mask_reference_dict.nc')
#     source_faces,source_rows,source_cols = read_mask_reference_from_nc_dict(nc_dict_file, mask_name)
#
#     for dest_file in dest_files:
#         if dest_file not in os.listdir(os.path.join(output_dir,mask_name,boundary_var_name)):
#             # try:
#             print('    - Downscaling the timesteps to be stored in file ' + str(dest_file))
#             source_files = source_file_read_dict[dest_file]
#             source_file_read_indices = source_file_read_index_sets[dest_file]
#
#             print('    - Reading in the L1 diagnostics_vec output')
#             L1_run_dir = os.path.join(config_dir,'L1',L1_model_name,'run')
#             L1_points, L1_values, L1_wet_grid_3D_points = read_boundary_variable_to_L1_points(L1_run_dir, mask_name, boundary_var_name,
#                                                                        source_files, source_file_read_indices,
#                                                                        source_faces,source_rows, source_cols,
#                                                                        Nr, faces, L1_XC_faces, L1_YC_faces, L1_Hfac_faces)
#
#             print(np.shape(L1_points),np.shape(L1_values))
#
#             L1_wet_grid_on_L2_3D_points = np.copy(L1_wet_grid_3D_points)
#
#             print('    - Downscaling the output to the new boundary')
#             L2_boundary_var = downscale_L1_bc_field_to_L2(df, mask_name, boundary_var_name,
#                                                           L1_points, L1_values,
#                                                           L1_wet_grid_3D_points, L1_wet_grid_on_L2_3D_points,
#                                                           L2_XC_subset, L2_YC_subset, L2_wet_grid_3D_subset)
#
#             if 'north' in mask_name or 'south' in mask_name:
#                 plt.subplot(1,2,1)
#                 plt.imshow(L2_wet_grid_3D_subset[:,0,:]==0,vmin=-0.1,vmax=0.1,cmap='seismic_r')
#                 plt.title('Mask')
#                 plt.subplot(1, 2, 2)
#                 plt.imshow(L2_boundary_var[0, :, 0, :],cmap='viridis')#,vmin=-0.1,vmax=0.1
#                 plt.title('L2 (nan values: '+str(np.sum(np.isnan(L2_boundary_var[0,:,0,:]))))
#                 plt.show()
#
#             if 'east' in mask_name or 'west' in mask_name:
#                 plt.subplot(1, 2, 1)
#                 plt.imshow(L2_wet_grid_3D_subset[:, :, 0] == 0, vmin=-0.1, vmax=0.1, cmap='seismic_r')
#                 plt.title('Mask')
#                 plt.subplot(1, 2, 2)
#                 plt.imshow(L2_boundary_var[0, :, :, 0], cmap='viridis')  # ,vmin=-0.1,vmax=0.1
#                 plt.title('L2 (nan values: ' + str(np.sum(np.isnan(L2_boundary_var[0, :, :, 0]))))
#                 plt.show()
#
#             if mask_name not in os.listdir(output_dir):
#                 os.mkdir(os.path.join(output_dir, mask_name))
#             if boundary_var_name not in os.listdir(os.path.join(output_dir, mask_name)):
#                 os.mkdir(os.path.join(output_dir, mask_name, boundary_var_name))
#
#             output_file = os.path.join(output_dir, mask_name, boundary_var_name, dest_file)
#             L2_boundary_var.ravel(order='C').astype('>f4').tofile(output_file)
#
#         else:
#             print('    - Skipping ' + str(dest_file) + ' because it was already created')


