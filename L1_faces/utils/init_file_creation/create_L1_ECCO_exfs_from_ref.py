
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
#from MITgcmutils import mds
#from scipy.interpolate import griddata
import sys
from datetime import datetime
import ast

def read_stitched_grid_tile_geometry(config_dir,model_name,var_name,sNx,sNy,ordered_nonblank_tiles,ordered_nonblank_rotations):

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])
    Nr = 50

    stitched_XC = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    stitched_YC = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    stitched_HFac = np.zeros((Nr, sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    grid_dir = os.path.join(config_dir, 'L1', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
                    if var_name in ['VSTRESS','VWIND']:
                        hFac = 'S'
                    elif var_name in ['USTRESS','UWIND']:
                        hFac = 'W'
                    else:
                        hFac = 'C'
                    HFac = ds.variables['HFac' + hFac][:, :, :]
                    if hFac == 'W':
                        HFac = HFac[:, :, :-1]
                        if tile_number == 1:
                            HFac[:, :, 0] = HFac[:, :, 1]
                    if hFac == 'S':
                        HFac = HFac[:, :-1, :]
                        if tile_number in [1, 2, 3, 4]:
                            HFac[:, 0, :] = HFac[:, 1, :]
                    ds.close()

                    for i in range(ordered_nonblank_rotations[r][c]):
                        XC = np.rot90(XC)
                        YC = np.rot90(YC)
                        HFac = np.rot90(HFac,axes=(1,2))

                    stitched_XC[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = XC
                    stitched_YC[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = YC
                    stitched_HFac[:, r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = HFac

    # plt.subplot(1,3,1)
    # C = plt.imshow(stitched_XC,origin='lower')
    # plt.colorbar(C)
    #
    # plt.subplot(1, 3, 2)
    # C = plt.imshow(stitched_YC, origin='lower')
    # plt.colorbar(C)
    #
    # plt.show()

    return(stitched_XC, stitched_YC, stitched_HFac)

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
                    if var_name in ['VSTRESS','VWIND']:
                        hFac = 'S'
                    elif var_name in ['USTRESS','UWIND']:
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

def create_src_dest_dicts_from_ref(config_dir, model_name, var_name,
                                   start_year, final_year, start_month, final_month, start_day, final_day):
    prefix = '_'.join(model_name.split('_')[:2])

    dest_files = []
    start_date = datetime(start_year, start_month, start_day)
    final_date = datetime(final_year, final_month, final_day)
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
                test_date = datetime(year, month, day)
                if test_date >= start_date and test_date <= final_date:
                    dest_files.append(str(year) + '{:02d}'.format(month) + '{:02d}'.format(day))

    f = open(os.path.join(config_dir,'L0', 'run','dv',prefix, model_name+'_exf_dest_ref.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)

    dest_files_out = []
    source_file_read_dict = {}
    source_file_read_index_sets = {}

    if var_name in ['UWIND','VWIND','USTRESS','VSTRESS']:
        suffix = '_rotated.bin'
    else:
        suffix = '.bin'

    for dest_file in dest_files:
        dest_files_out.append(prefix+'_surface_'+var_name+'.'+dest_file+suffix)
        source_files = []
        index_set = []
        for pair in size_dict[dest_file]:
            source_files.append(prefix+'_surface_'+var_name+'.'+pair[0]+'.bin')
            index_set.append(pair[1])
        source_file_read_dict[prefix+'_surface_'+var_name+'.'+dest_file+suffix] = source_files
        source_file_read_index_sets[prefix+'_surface_'+var_name+'.'+dest_file+suffix] = index_set

    return(dest_files_out, source_file_read_dict, source_file_read_index_sets)

def read_L0_surface_variable_points(L0_run_dir, model_name, var_name,
                                     source_files, source_file_read_indices,
                                     faces, rows, cols,
                                     ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces,
                                     print_level):
    n_Timesteps = 0
    for index_set in source_file_read_indices:
        n_Timesteps += index_set[1] - index_set[0]

    # make a blank grid of zeros
    points = np.zeros((np.size(faces), 2))
    hfac_points = np.zeros((1,np.size(faces)))
    values = np.zeros((n_Timesteps, np.size(faces)))

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
        boundary = 'surface'
        if var_name in ['UWIND','VWIND']:
            u_var_file = os.path.join(L0_run_dir, 'dv',prefix, prefix+'_' + boundary + '_mask_UWIND.' + file_suffix)
            u_var_grid = np.fromfile(u_var_file, dtype='>f4')
            v_var_file = os.path.join(L0_run_dir, 'dv',prefix, prefix+'_' + boundary + '_mask_VWIND.' + file_suffix)
            v_var_grid = np.fromfile(v_var_file, dtype='>f4')
        elif var_name in ['USTRESS','VSTRESS']:
            u_var_file = os.path.join(L0_run_dir, 'dv',prefix, prefix+'_' + boundary + '_mask_USTRESS.' + file_suffix)
            u_var_grid = np.fromfile(u_var_file, dtype='>f4')
            v_var_file = os.path.join(L0_run_dir, 'dv',prefix, prefix+'_' + boundary + '_mask_VSTRESS.' + file_suffix)
            v_var_grid = np.fromfile(v_var_file, dtype='>f4')
        else:
            var_file = os.path.join(L0_run_dir, 'dv',prefix,prefix+'_'+ boundary + '_mask_' + var_name +'.' + file_suffix)
            var_grid = np.fromfile(var_file, dtype='>f4')

        if var_name in ['UWIND','VWIND','USTRESS','VSTRESS']:
            timesteps_in_file = int(np.size(u_var_grid) / N)
            u_var_grid = np.reshape(u_var_grid, (timesteps_in_file, N))
            u_var_grid = u_var_grid[start_file_index:end_file_index, :]
            v_var_grid = np.reshape(v_var_grid, (timesteps_in_file, N))
            v_var_grid = v_var_grid[start_file_index:end_file_index, :]
        else:
            timesteps_in_file = int(np.size(var_grid) / N)
            var_grid = np.reshape(var_grid, (timesteps_in_file, N))
            var_grid = var_grid[start_file_index:end_file_index, :]

        for n in range(N):
            if faces[n] != 0:
                points[n, 0] = ecco_XC_faces[faces[n]][rows[n], cols[n]]
                points[n, 1] = ecco_YC_faces[faces[n]][rows[n], cols[n]]
                hfac_points[0, n] = ecco_hFacC_faces[faces[n]][0, rows[n], cols[n]]
                if var_name in ['UWIND','VWIND','USTRESS', 'VSTRESS']:
                    angle_cos = ecco_AngleCS_faces[faces[n]][rows[n], cols[n]]
                    angle_sin = ecco_AngleSN_faces[faces[n]][rows[n], cols[n]]
                    if var_name=='UWIND' or var_name=='USTRESS':
                        zonal_velocity = angle_cos * u_var_grid[:, n] - angle_sin * v_var_grid[:, n]
                        values[index_counter:index_counter + (end_file_index - start_file_index), n] = zonal_velocity
                    if var_name=='VWIND' or var_name=='VSTRESS':
                        meridional_velocity = angle_sin * u_var_grid[:, n] + angle_cos * v_var_grid[:, n]
                        values[index_counter:index_counter + (end_file_index - start_file_index), n] = meridional_velocity
                else:
                    values[index_counter:index_counter + (end_file_index - start_file_index), n] = var_grid[:, n]

        index_counter += (end_file_index - start_file_index)

    return(points,values,hfac_points)

def create_L1_exfs(Lf, config_dir, model_name, var_name, ecco_dir, llc,
                   ordered_ecco_tiles, ordered_ecco_tile_rotations,
                   start_year, final_year, start_month, final_month, start_day, final_day,
                   sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations, tile_face_index_dict, face_size_dict, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import ecco_functions as ef

    if 'exf' not in os.listdir(os.path.join(config_dir, 'L1', model_name, 'input')):
        os.mkdir(os.path.join(config_dir, 'L1', model_name, 'input', 'exf'))
    if var_name not in os.listdir(os.path.join(config_dir, 'L1', model_name, 'input', 'exf')):
        os.mkdir(os.path.join(config_dir, 'L1', model_name, 'input', 'exf', var_name))

    if print_level >= 1:
        print('    - Creating the ' + var_name + ' exf files for the ' + model_name + ' model from ECCOv5 output data')

    # llc = 270
    n_timesteps = 4

    # step 0: get the model domain
    if print_level >= 1:
        print('    - Reading in the model geometry')
    # ordered_XC_tiles, ordered_YC_tiles, ordered_AngleCS_tiles, ordered_AngleSN_tiles, ordered_hfac_tiles, delR = \
    #     read_grid_tile_geometry(config_dir, model_name, var_name, ordered_nonblank_tiles, tile_face_index_dict)
    L1_XC, L1_YC, L1_mask = read_stitched_grid_tile_geometry(config_dir, model_name, var_name, sNx, sNy, ordered_nonblank_tiles,
                                     ordered_nonblank_rotations)
    L1_mask[L1_mask>0]=1
    L1_mask = L1_mask[:1,:,:]

    # step 1: get the ecco faces geometry
    if print_level >= 1:
        print('    - Reading in the ECCO geometry')
    ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces = \
        ef.read_ecco_geometry_to_faces(ecco_dir, llc)

    L0_run_dir = os.path.join(config_dir, 'L0', 'run')

    dest_files, source_file_read_dict, source_file_read_index_sets = \
        create_src_dest_dicts_from_ref(config_dir, model_name, var_name,
                                       start_year, final_year, start_month,
                                       final_month, start_day, final_day)

    if print_level >= 3:
        print('            - Reading the mask to reference the variable to the llc grid')
    nc_dict_file = os.path.join(config_dir, 'L0', 'input', 'L0_dv_mask_reference_dict.nc')
    ds = nc4.Dataset(nc_dict_file)
    grp = ds.groups['L1_CE_surface']
    faces = grp.variables['source_faces'][:]
    rows = grp.variables['source_rows'][:]
    cols = grp.variables['source_cols'][:]
    ds.close()

    llc = 270

    ####################################################################################################################
    # Loop through the daily files

    for dest_file in dest_files:
        if dest_file not in os.listdir(os.path.join(config_dir, 'L1', model_name, 'input', 'exf', var_name)):

            if print_level >= 3:
                print('            - Downscaling the timesteps to be stored in file ' + str(dest_file))

            source_files = source_file_read_dict[dest_file]
            source_file_read_indices = source_file_read_index_sets[dest_file]

            output_grid = np.zeros((n_timesteps, sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

            L0_surface_points, L0_surface_values, L0_surface_hFacC = \
                read_L0_surface_variable_points(L0_run_dir, model_name, var_name,
                                            source_files, source_file_read_indices,
                                            faces, rows, cols,
                                            ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces,
                                            ecco_hFacC_faces,
                                            print_level)

            L0_surface_values = np.reshape(L0_surface_values, (
                np.shape(L0_surface_values)[0], 1, np.shape(L0_surface_values)[1]))
            L0_surface_mask = np.copy(L0_surface_hFacC)
            L0_surface_mask[L0_surface_mask > 0] = 1

            for timestep in range(np.shape(L0_surface_values)[0]):
                # if timestep%50 == 0:
                #     print('        - Downscaling timestep '+str(timestep))

                if print_level >= 5:
                    if timestep == 0:
                        print('                - L0_surface_points shape: ' + str(np.shape(L0_surface_points)))
                        print('                - L0_surface_values shape: ' + str(np.shape(L0_surface_values)))
                        print('                - L0_surface_mask: ' + str(np.shape(L0_surface_mask)))
                        print('                - L1_XC shape: ' + str(np.shape(L1_XC)))
                        print('                - L1_YC shape: ' + str(np.shape(L1_YC)))
                        print('                - L1_mask shape: ' + str(np.shape(L1_mask)))

                interp_field = df.downscale_3D_points_with_zeros(L0_surface_points,
                                                                 L0_surface_values[timestep, :, :],
                                                                 L0_surface_mask,
                                                                 L1_XC, L1_YC, L1_mask,
                                                                 mean_vertical_difference=0,
                                                                 fill_downward=True,
                                                                 printing=False)

                output_grid[timestep, :, :] = interp_field[0, :, :]

                # if timestep==0:
                #     plt.subplot(1,2,1)
                #     plt.imshow(L1_mask[0,:,:])
                #     plt.subplot(1, 2, 2)
                #     plt.imshow(interp_field[0,:,:])
                #     plt.show()

            output_faces = Lf.read_stitched_grid_to_faces(output_grid, sNx, sNy, dim=3)

            # plt.imshow(runoff_faces[1][0, :, :])
            # plt.show()
            #
            # plt.imshow(runoff_faces[5][0, :, :])
            # plt.show()

            output_compact = Lf.read_faces_to_compact(output_faces, sNx, sNy)

            if print_level >= 4:
                print('                - Compact output shape: ' + str(np.shape(output_compact)))

            # plt.imshow(interp_runoff[0, :, :],origin='lower')
            # plt.show()

            output_file = os.path.join(config_dir, 'L1', model_name, 'input', 'exf', var_name, dest_file)
            output_compact.ravel(order='C').astype('>f4').tofile(output_file)





