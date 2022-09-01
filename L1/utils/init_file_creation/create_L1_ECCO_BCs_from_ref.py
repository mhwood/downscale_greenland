
import os
import numpy as np
import netCDF4 as nc4
import ast
from datetime import datetime
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

    delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                     10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                     31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                     93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                     139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                     341.50, 364.50, 387.50, 410.50, 433.50, 456.50])

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L1', model_name, 'run_for_grid')

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
                    # if tile_face < 3:
                    if var_name in ['VVEL','VICE']:
                        hFac = 'S'
                    elif var_name in ['UVEL','UICE']:
                        hFac = 'W'
                    else:
                        hFac = 'C'
                    # else:
                    #     if var_name == 'VVEL':
                    #         hFac = 'W'
                    #     elif var_name == 'UVEL':
                    #         hFac = 'S'
                    #     else:
                    #         hFac = 'C'
                    HFac = ds.variables['HFac' + hFac][:, :, :]
                    if hFac == 'W':
                        HFac = HFac[:,:,:-1]
                        if tile_number==1:
                            HFac[:,:,0] = HFac[:,:,1]
                    if hFac == 'S':
                        HFac = HFac[:,:-1,:]
                        if tile_number in [1,2,3,4]:
                            HFac[:,0,:] = HFac[:,1,:]
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

def create_src_dest_dicts_from_ref(config_dir, model_name, boundary, var_name,
                                   start_year, final_year, start_month, final_month, start_day, final_day):
    prefix = '_'.join(model_name.split('_')[:2])

    dest_files = []
    start_date = datetime(start_year, start_month, start_day)
    final_date = datetime(final_year, final_month, final_day)
    for year in range(1992, 1993):
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

    f = open(os.path.join(config_dir,'L0', 'run','dv',prefix, model_name+'_BC_dest_ref.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)

    dest_files_out = []
    source_file_read_dict = {}
    source_file_read_index_sets = {}

    if var_name in ['UVEL','VVEL','UICE','VICE']:
        suffix = '_rotated.bin'
    else:
        suffix = '.bin'


    for dest_file in dest_files:
        dest_files_out.append(prefix+'_'+boundary+'_'+var_name+'.'+dest_file+suffix)
        source_files = []
        index_set = []
        for pair in size_dict[dest_file]:
            source_files.append(prefix+'_'+boundary+'_'+var_name+'.'+pair[0]+'.bin')
            index_set.append(pair[1])
        source_file_read_dict[prefix+'_'+boundary+'_'+var_name+'.'+dest_file+suffix] = source_files
        source_file_read_index_sets[prefix+'_'+boundary+'_'+var_name+'.'+dest_file+suffix] = index_set

    return(dest_files_out, source_file_read_dict, source_file_read_index_sets)

def read_L0_boundary_variable_points(L0_run_dir, model_name, boundary, var_name,
                                     source_files,source_file_read_indices,
                                     llc, Nr, faces, rows, cols,
                                     ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces,
                                     print_level):
    n_Timesteps = 0
    for index_set in source_file_read_indices:
        n_Timesteps += index_set[1] - index_set[0] + 1

    # print('           + the L0 grid for this file will have ' + str(n_Timesteps) + ' timesteps')

    # make a blank grid of zeros
    points = np.zeros((np.size(faces), 2))
    hfac_points = np.zeros((Nr,np.size(faces)))
    if var_name in ['ETAN','AREA','HEFF','HSNOW','UICE','VICE']:
        values = np.zeros((n_Timesteps, np.size(faces)))
    else:
        values = np.zeros((n_Timesteps, Nr, np.size(faces)))

    index_counter = 0
    for s in range(len(source_files)):
        source_file = source_files[s]
        file_suffix = '.'.join(source_file.split('.')[-2:])
        index_set = source_file_read_indices[s]
        if print_level >= 4:
            print('                - Adding timesteps ' + str(index_set[0]) + ' to ' + str(index_set[1]) + ' from file ' + source_file)

        start_file_index = index_set[0]
        end_file_index = index_set[1]

        if print_level >= 4:
            print('                - Storing at points ' + str(index_counter) + ' to ' + str(
            index_counter + (end_file_index - start_file_index)) + ' in the grid')

        N = len(faces)

        prefix = '_'.join(model_name.split('_')[:2])
        if var_name in ['UVEL','VVEL']:
            u_var_file = os.path.join(L0_run_dir, 'dv',prefix, prefix+'_' + boundary + '_mask_UVEL.' + file_suffix)
            u_var_grid = np.fromfile(u_var_file, dtype='>f4')
            v_var_file = os.path.join(L0_run_dir, 'dv',prefix, prefix+'_' + boundary + '_mask_VVEL.' + file_suffix)
            v_var_grid = np.fromfile(v_var_file, dtype='>f4')
        elif var_name in ['UICE','VICE']:
            u_var_file = os.path.join(L0_run_dir, 'dv',prefix, prefix+'_' + boundary + '_mask_UICE.' + file_suffix)
            u_var_grid = np.fromfile(u_var_file, dtype='>f4')
            v_var_file = os.path.join(L0_run_dir, 'dv',prefix, prefix+'_' + boundary + '_mask_VICE.' + file_suffix)
            v_var_grid = np.fromfile(v_var_file, dtype='>f4')
        else:
            var_file = os.path.join(L0_run_dir, 'dv',prefix,prefix+'_'+ boundary + '_mask_' + var_name +'.' + file_suffix)
            var_grid = np.fromfile(var_file, dtype='>f4')

        if var_name in ['ETAN','AREA','HEFF','HSNOW','UICE','VICE']:
            if var_name in ['UICE','VICE']:
                timesteps_in_file = int(np.size(u_var_grid) / (N))
                u_var_grid = np.reshape(u_var_grid, (timesteps_in_file, N))
                u_var_grid = u_var_grid[start_file_index:end_file_index+1, :]
                v_var_grid = np.reshape(v_var_grid, (timesteps_in_file, N))
                v_var_grid = v_var_grid[start_file_index:end_file_index+1, :]
            else:
                timesteps_in_file = int(np.size(var_grid) / (N))
                var_grid = np.reshape(var_grid, (timesteps_in_file, N))
                var_grid = var_grid[start_file_index:end_file_index+1, :]

            for n in range(N):
                if faces[n]!=0:
                    points[n, 0] = ecco_XC_faces[faces[n]][rows[n], cols[n]]
                    points[n, 1] = ecco_YC_faces[faces[n]][rows[n], cols[n]]
                    hfac_points[:, n] = ecco_hFacC_faces[faces[n]][:, rows[n], cols[n]]
                    if var_name in ['UICE', 'VICE']:
                        angle_cos = ecco_AngleCS_faces[faces[n]][rows[n], cols[n]]
                        angle_sin = ecco_AngleSN_faces[faces[n]][rows[n], cols[n]]
                        if var_name=='UICE':
                            zonal_velocity = angle_cos * u_var_grid[:, n] - angle_sin * v_var_grid[:, n]
                            values[index_counter:index_counter + (end_file_index - start_file_index)+1, n] = zonal_velocity
                        if var_name=='VICE':
                            meridional_velocity = angle_sin * u_var_grid[:, n] + angle_cos * v_var_grid[:, n]
                            values[index_counter:index_counter + (end_file_index - start_file_index)+1,n] = meridional_velocity
                    else:
                        values[index_counter:index_counter + (end_file_index - start_file_index)+1, n] = var_grid[:, n]
        else:
            if var_name in ['UVEL','VVEL']:
                timesteps_in_file = int(np.size(u_var_grid) / (Nr * N))
                u_var_grid = np.reshape(u_var_grid, (timesteps_in_file, Nr, N))
                u_var_grid = u_var_grid[start_file_index:end_file_index+1, :, :]
                v_var_grid = np.reshape(v_var_grid, (timesteps_in_file, Nr, N))
                v_var_grid = v_var_grid[start_file_index:end_file_index+1, :, :]
            else:
                timesteps_in_file = int(np.size(var_grid) / (Nr * N))
                var_grid = np.reshape(var_grid, (timesteps_in_file, Nr, N))
                var_grid = var_grid[start_file_index:end_file_index+1, :, :]
            for n in range(N):
                if faces[n] != 0:
                    points[n, 0] = ecco_XC_faces[faces[n]][rows[n], cols[n]]
                    points[n, 1] = ecco_YC_faces[faces[n]][rows[n], cols[n]]
                    hfac_points[:, n] = ecco_hFacC_faces[faces[n]][:, rows[n], cols[n]]
                    if var_name in ['UVEL', 'VVEL']:
                        angle_cos = ecco_AngleCS_faces[faces[n]][rows[n], cols[n]]
                        angle_sin = ecco_AngleSN_faces[faces[n]][rows[n], cols[n]]
                        if var_name=='UVEL':
                            zonal_velocity = np.zeros((np.shape(u_var_grid)[0],np.shape(u_var_grid)[1]))
                            for k in range(np.shape(zonal_velocity)[1]):
                                zonal_velocity[:,k] = angle_cos * u_var_grid[:, k, n] - angle_sin * v_var_grid[:, k, n]
                            values[index_counter:index_counter + (end_file_index - start_file_index)+1, :,n] = zonal_velocity
                        if var_name=='VVEL':
                            meridional_velocity = np.zeros((np.shape(u_var_grid)[0],np.shape(u_var_grid)[1]))
                            for k in range(np.shape(meridional_velocity)[1]):
                                meridional_velocity[:,k] = angle_sin * u_var_grid[:, k, n] + angle_cos * v_var_grid[:, k, n]
                            values[index_counter:index_counter + (end_file_index - start_file_index)+1,:,n] = meridional_velocity
                    else:
                        values[index_counter:index_counter + (end_file_index - start_file_index)+1, :, n] = var_grid[:, :, n]

        index_counter += (end_file_index - start_file_index)+1

    return(points,values,hfac_points)

def subset_tile_geometry_to_boundary(boundary, tile_number,ordered_nonblank_tiles,
                                     ordered_XC_tiles, ordered_YC_tiles,
                                     ordered_AngleCS_tiles, ordered_AngleSN_tiles, ordered_hfac_tiles):

    # get the geometry for this tile
    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            if ordered_nonblank_tiles[r][c] == tile_number:
                XC = ordered_XC_tiles[r][c]
                YC = ordered_YC_tiles[r][c]
                AngleCS = ordered_AngleCS_tiles[r][c]
                AngleSN = ordered_AngleSN_tiles[r][c]
                hFac = ordered_hfac_tiles[r][c]

    # subset to the boundary
    if boundary=='south':
        XC_subset = XC[:1, :]
        YC_subset = YC[:1, :]
        AngleCS_subset = AngleCS[:1, :]
        AngleSN_subset = AngleSN[:1, :]
        hFac_subset = hFac[:,:1, :]

    if boundary=='west':
        XC_subset = XC[:, :1]
        YC_subset = YC[:, :1]
        AngleCS_subset = AngleCS[:, :1]
        AngleSN_subset = AngleSN[:, :1]
        hFac_subset = hFac[:,:, :1]

    if boundary=='north':
        XC_subset = XC[-1:, :]
        YC_subset = YC[-1:, :]
        AngleCS_subset = AngleCS[-1:, :]
        AngleSN_subset = AngleSN[-1:, :]
        hFac_subset = hFac[:,-1:, :]

    if boundary=='east':
        XC_subset = XC[:, -1:]
        YC_subset = YC[:, -1:]
        AngleCS_subset = AngleCS[:, -1:]
        AngleSN_subset = AngleSN[:, -1:]
        hFac_subset = hFac[:,:, -1:]

    return(XC_subset, YC_subset, AngleCS_subset, AngleSN_subset, hFac_subset)

def create_L1_BCs(config_dir,model_name,var_name,
                  Nr, sNx, sNy, ordered_nonblank_tiles, tile_face_index_dict,
                  ecco_dir, llc,
                  northern_tiles, southern_tiles, eastern_tiles, western_tiles,
                  start_year, final_year, start_month,final_month, start_day, final_day, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import ecco_functions as ef

    if 'obcs' not in os.listdir(os.path.join(config_dir, 'L1', model_name, 'input')):
        os.mkdir(os.path.join(config_dir, 'L1', model_name, 'input', 'obcs'))

    if print_level>=1:
        print('    - Creating the '+var_name+' BC files for the '+model_name+' model from ECCOv5 data')

    # llc = 270
    n_timesteps = 24

    # step 0: get the model domain
    if print_level >= 1:
        print('    - Reading in the model geometry')
    ordered_XC_tiles, ordered_YC_tiles, ordered_AngleCS_tiles, ordered_AngleSN_tiles, ordered_hfac_tiles, delR = \
        read_grid_tile_geometry(config_dir,model_name,var_name,ordered_nonblank_tiles,tile_face_index_dict)

    # step 1: get the ecco faces geometry
    if print_level >= 1:
        print('    - Reading in the ECCO geometry')
    # ecco_XC, ecco_YC, ecco_AngleCS, ecco_AngleSN, ecco_hfacC, ecco_delR = \
    #     ef.read_ecco_grid_geometry(ecco_dir, llc, ordered_ecco_tiles, ordered_ecco_tile_rotations)
    ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces = \
        ef.read_ecco_geometry_to_faces(ecco_dir, llc)

    L0_run_dir = os.path.join(config_dir, 'L0', 'run')

    ####################################################################################################################
    # Loop through the boundaries

    for boundary in ['south','west','east','north']:

        if boundary == 'south':
            tile_numbers = southern_tiles
        if boundary == 'west':
            tile_numbers = western_tiles
        if boundary == 'north':
            tile_numbers = northern_tiles
        if boundary == 'east':
            tile_numbers = eastern_tiles

        if len(tile_numbers) > 0:

            if print_level >= 2:
                print('        - Running the downscale routine for the '+boundary+' boundary')


            if boundary not in os.listdir(os.path.join(config_dir,'L1',model_name,'input','obcs')):
                os.mkdir(os.path.join(config_dir,'L1',model_name,'input','obcs',boundary))
            if var_name not in os.listdir(os.path.join(config_dir, 'L1', model_name, 'input', 'obcs', boundary)):
                os.mkdir(os.path.join(config_dir, 'L1', model_name, 'input', 'obcs', boundary, var_name))

            dest_files, source_file_read_dict, source_file_read_index_sets = \
                create_src_dest_dicts_from_ref(config_dir, model_name, boundary, var_name,
                                               start_year, final_year, start_month,
                                               final_month, start_day, final_day)

            if print_level >= 3:
                print('            - Reading the mask to reference the variable to the llc grid')
            nc_dict_file = os.path.join(config_dir, 'L0', 'input', 'L0_dv_mask_reference_dict.nc')
            ds = nc4.Dataset(nc_dict_file)
            boundary_group = '_'.join(model_name.split('_')[:2]) + '_'+boundary
            grp = ds.groups[boundary_group]
            faces = grp.variables['source_faces'][:]
            rows = grp.variables['source_rows'][:]
            cols = grp.variables['source_cols'][:]
            ds.close()

            ############################################################################################################
            # Loop through the destination files

            for dest_file in dest_files:
                if dest_file not in os.listdir(os.path.join(config_dir, 'L1', model_name, 'input', 'obcs', boundary, var_name)):
                    if print_level >= 3:
                        print('            - Downscaling the timesteps to be stored in file ' + str(dest_file))
                    source_files = source_file_read_dict[dest_file]
                    source_file_read_indices = source_file_read_index_sets[dest_file]

                    if print_level >= 3:
                        print('            - Reading in the L0 diagnostics_vec output')
                    L0_boundary_points, L0_boundary_values, L0_boundary_point_hFacC = \
                        read_L0_boundary_variable_points(L0_run_dir, model_name, boundary, var_name,
                                                         source_files, source_file_read_indices,
                                                         llc, Nr, faces, rows, cols,
                                                         ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces,
                                                         print_level)

                    if var_name in ['AREA', 'HEFF', 'HSNOW', 'UICE', 'VICE', 'ETAN']:
                        if boundary in ['north', 'south']:
                            output_grid = np.zeros((n_timesteps, sNx * len(tile_numbers)))
                        else:
                            output_grid = np.zeros((n_timesteps, sNy * len(tile_numbers)))
                    else:
                        if boundary in ['north', 'south']:
                            output_grid = np.zeros((n_timesteps, Nr, sNx * len(tile_numbers)))
                        else:
                            output_grid = np.zeros((n_timesteps, Nr, sNy * len(tile_numbers)))

                    if var_name in ['AREA', 'HEFF', 'HSNOW', 'UICE', 'VICE', 'ETAN']:
                        L0_boundary_values = L0_boundary_values[:, L0_boundary_points[:, 0] != 0]
                        L0_boundary_values = np.reshape(L0_boundary_values, (
                        np.shape(L0_boundary_values)[0], 1, np.shape(L0_boundary_values)[1]))
                        L0_boundary_point_hFacC = L0_boundary_point_hFacC[:1, :]
                    else:
                        L0_boundary_values = L0_boundary_values[:, :, L0_boundary_points[:, 0] != 0]
                    L0_boundary_point_hFacC = L0_boundary_point_hFacC[:, L0_boundary_points[:, 0] != 0]
                    L0_boundary_points = L0_boundary_points[L0_boundary_points[:, 0] != 0, :]

                    L0_boundary_point_mask = np.copy(L0_boundary_point_hFacC)
                    L0_boundary_point_mask[L0_boundary_point_mask > 0] = 1

                    for tn in range(len(tile_numbers)):
                        tile_number = tile_numbers[tn]
                        if tile_number>0:
                            if print_level >= 4:
                                print('                - Downscaling data for tile number '+str(tile_number))

                            XC_subset, YC_subset, AngleCS_subset, AngleSN_subset, hFac_subset = \
                                subset_tile_geometry_to_boundary(boundary, tile_number, ordered_nonblank_tiles,
                                                                 ordered_XC_tiles, ordered_YC_tiles,
                                                                 ordered_AngleCS_tiles, ordered_AngleSN_tiles, ordered_hfac_tiles)

                            mask_subset = np.copy(hFac_subset)
                            mask_subset[mask_subset>0]=1
                            if var_name in ['AREA', 'HEFF', 'HSNOW', 'UICE', 'VICE', 'ETAN']:
                                mask_subset = mask_subset[:1, :, :]

                            # plt.plot(L0_boundary_points[:, 0], L0_boundary_points[:, 1], 'g.')
                            # plt.plot(XC_subset.ravel(),YC_subset.ravel(),'k-')
                            # plt.show()

                            for timestep in range(np.shape(L0_boundary_values)[0]):
                                # if timestep%50 == 0:
                                #     print('        - Downscaling timestep '+str(timestep))

                                if print_level >= 5:
                                    if timestep == 0:
                                        print('                - L0_boundary_points shape: '+str(np.shape(L0_boundary_points)))
                                        print('                - L0_boundary_values shape: ' + str(np.shape(L0_boundary_values)))
                                        print('                - L0_boundary_point_mask shape: ' + str(np.shape(L0_boundary_point_mask)))
                                        print('                - XC_subset shape: ' + str(np.shape(XC_subset)))
                                        print('                - YC_subset shape: ' + str(np.shape(YC_subset)))
                                        print('                - mask_subset shape: ' + str(np.shape(mask_subset)))

                                # interp_field = df.downscale_3D_boundary_points(L0_boundary_points, L0_boundary_values[timestep,:,:], L0_boundary_point_mask,
                                #                              XC_subset, YC_subset, mask_subset,
                                #                              mean_vertical_difference=0, fill_downward=True,
                                #                              printing=False, remove_zeros=remove_zeros)

                                interp_field = df.downscale_3D_points_with_zeros(L0_boundary_points,
                                                                                  L0_boundary_values[timestep, :, :],
                                                                                  L0_boundary_point_mask,
                                                                                  XC_subset, YC_subset, mask_subset,
                                                                                  mean_vertical_difference=0,
                                                                                  fill_downward=True,
                                                                                  printing=False)

                                if var_name in ['AREA', 'HEFF', 'HSNOW', 'UICE', 'VICE', 'ETAN']:
                                    if boundary in ['north', 'south']:
                                        output_grid[timestep, tn * sNx:(tn + 1) * sNx] = interp_field[0, 0, :]
                                    else:
                                        output_grid[timestep, tn * sNy:(tn + 1) * sNy] = interp_field[0, :, 0]
                                else:
                                    if boundary in ['north', 'south']:
                                        output_grid[timestep, :, tn * sNx:(tn + 1) * sNx] = interp_field[:,0,:]
                                    else:
                                        output_grid[timestep, :, tn * sNy:(tn + 1) * sNy] = interp_field[:,:,0]

                            # if timestep==0 and tile_number == 4:
                            #     if boundary in ['north', 'south']:
                            #         plt.subplot(1,2,1)
                            #         plt.imshow(mask_subset[:,0,:])
                            #         plt.subplot(1, 2, 2)
                            #         plt.imshow(interp_field[:,0,:])
                            #     else:
                            #         plt.subplot(1, 2, 1)
                            #         plt.imshow(mask_subset[:, :, 0])
                            #         plt.subplot(1, 2, 2)
                            #         plt.imshow(interp_field[:, :, 0])
                            #     plt.show()

                    # if var_name in ['AREA', 'HEFF', 'HSNOW', 'UICE', 'VICE', 'ETAN']:
                    #     plt.plot(output_grid[0, :])
                    #     plt.show()
                    # else:
                    #     plt.imshow(output_grid[0, :, :],vmin=-0.1,vmax=0.1,cmap='seismic')
                    #     plt.show()

                    # if var_name in ['UVEL','VVEL','UICE','VICE']:
                    #     output_file = os.path.join(config_dir, 'L1', model_name, 'input', 'obcs', boundary, var_name,dest_file[:-4]+'_rotated.bin')
                    # else:
                    output_file = os.path.join(config_dir, 'L1', model_name, 'input', 'obcs', boundary, var_name, dest_file)
                    output_grid.ravel(order='C').astype('>f4').tofile(output_file)





