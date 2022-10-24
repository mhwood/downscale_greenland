
import os
import numpy as np
import netCDF4 as nc4
import ast
from datetime import datetime
import matplotlib.pyplot as plt
from MITgcmutils import mds
from scipy.interpolate import griddata
import sys


def read_grid_geometry(config_dir,model_name,var_name):

    file_path = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    AngleSN = ds.variables['AngleSN'][:, :]

    if var_name in ['VVEL', 'VICE']:
        hFac = 'S'
    elif var_name in ['UVEL', 'UICE']:
        hFac = 'W'
    else:
        hFac = 'C'

    HFac = ds.variables['HFac' + hFac][:, :, :]
    if hFac == 'W':
        HFac = HFac[:, :, :-1]
    if hFac == 'S':
        HFac = HFac[:, :-1, :]

    delR = ds.variables['drF'][:]

    ds.close()

    return(XC, YC, AngleCS, AngleSN, HFac, delR)

def create_src_dest_dicts_from_ref(config_dir, model_name, boundary, var_name,
                                   start_year, final_year, start_month, final_month, start_day, final_day, read_darwin):
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

    if read_darwin:
        f = open(os.path.join(config_dir,'L0', 'run_darwin','dv',prefix, model_name+'_BC_dest_ref.txt'))
    else:
        f = open(os.path.join(config_dir, 'L0', 'run', 'dv', prefix, model_name + '_BC_dest_ref.txt'))
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

def read_L0_boundary_variable_points(config_dir, L0_run_dir, model_name, boundary, var_name,
                                     source_files,source_file_read_indices,
                                     llc, Nr, boundary_read_groups, n_dv_points,
                                     ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces,
                                     print_level):
    # if print_level >= 3:
    #     print('            - Reading the mask to reference the variable to the llc grid')
    nc_dict_file = os.path.join(config_dir, 'L0', 'input', 'L0_dv_mask_reference_dict.nc')

    n_Timesteps = 0
    for index_set in source_file_read_indices:
        n_Timesteps += index_set[1] - index_set[0]

    # print('           + the L0 grid for this file will have ' + str(n_Timesteps) + ' timesteps')

    # make a blank grid of zeros
    points = np.zeros((n_dv_points, 2))
    hfac_points = np.zeros((Nr,n_dv_points))
    if var_name in ['ETAN','AREA','HEFF','HSNOW','UICE','VICE']:
        values = np.zeros((n_Timesteps, n_dv_points))
    else:
        values = np.zeros((n_Timesteps, Nr, n_dv_points))

    index_counter = 0
    for s in range(len(source_files)):
        source_file = source_files[s]
        file_suffix = '.'.join(source_file.split('.')[-2:])
        index_set = source_file_read_indices[s]
        if print_level >= 4:
            print('                - Adding timesteps ' + str(index_set[0]) + ' to ' + str(index_set[1]))

        start_file_index = index_set[0]
        end_file_index = index_set[1]

        if print_level >= 4:
            print('                - Storing at points ' + str(index_counter) + ' to ' + str(
            index_counter + (end_file_index - start_file_index)) + ' in the grid')

        boundary_points_counted = 0
        for read_boundary in boundary_read_groups:

            read_file = source_file.replace(boundary,read_boundary)

            if print_level >= 4:
                print('                - Adding timesteps ' + str(index_set[0]) + ' to ' + str(
                    index_set[1]) + ' from source files ' + read_file)

            ds = nc4.Dataset(nc_dict_file)
            boundary_group = '_'.join(model_name.split('_')[:2]) + '_'+read_boundary
            grp = ds.groups[boundary_group]
            faces = grp.variables['source_faces'][:]
            rows = grp.variables['source_rows'][:]
            cols = grp.variables['source_cols'][:]
            ds.close()

            N = len(faces)

            prefix = '_'.join(model_name.split('_')[:2])
            if var_name in ['UVEL','VVEL']:
                u_var_file = os.path.join(L0_run_dir, 'dv',prefix, prefix+'_' + read_boundary + '_mask_UVEL.' + file_suffix)
                u_var_grid = np.fromfile(u_var_file, dtype='>f4')
                v_var_file = os.path.join(L0_run_dir, 'dv',prefix, prefix+'_' + read_boundary + '_mask_VVEL.' + file_suffix)
                v_var_grid = np.fromfile(v_var_file, dtype='>f4')
            elif var_name in ['UICE','VICE']:
                u_var_file = os.path.join(L0_run_dir, 'dv',prefix, prefix+'_' + read_boundary + '_mask_UICE.' + file_suffix)
                u_var_grid = np.fromfile(u_var_file, dtype='>f4')
                v_var_file = os.path.join(L0_run_dir, 'dv',prefix, prefix+'_' + read_boundary + '_mask_VICE.' + file_suffix)
                v_var_grid = np.fromfile(v_var_file, dtype='>f4')
            else:
                var_file = os.path.join(L0_run_dir, 'dv',prefix,prefix+'_'+ read_boundary + '_mask_' + var_name +'.' + file_suffix)
                var_grid = np.fromfile(var_file, dtype='>f4')

            if var_name in ['ETAN','AREA','HEFF','HSNOW','UICE','VICE']:
                if var_name in ['UICE','VICE']:
                    timesteps_in_file = int(np.size(u_var_grid) / (N))
                    u_var_grid = np.reshape(u_var_grid, (timesteps_in_file, N))
                    u_var_grid = u_var_grid[start_file_index:end_file_index, :]
                    v_var_grid = np.reshape(v_var_grid, (timesteps_in_file, N))
                    v_var_grid = v_var_grid[start_file_index:end_file_index, :]
                else:
                    timesteps_in_file = int(np.size(var_grid) / (N))
                    var_grid = np.reshape(var_grid, (timesteps_in_file, N))
                    var_grid = var_grid[start_file_index:end_file_index, :]

                for n in range(N):
                    if faces[n]!=0:
                        points[boundary_points_counted+n, 0] = ecco_XC_faces[faces[n]][rows[n], cols[n]]
                        points[boundary_points_counted+n, 1] = ecco_YC_faces[faces[n]][rows[n], cols[n]]
                        hfac_points[:, boundary_points_counted+n] = ecco_hFacC_faces[faces[n]][:, rows[n], cols[n]]
                        if var_name in ['UICE', 'VICE']:
                            angle_cos = ecco_AngleCS_faces[faces[n]][rows[n], cols[n]]
                            angle_sin = ecco_AngleSN_faces[faces[n]][rows[n], cols[n]]
                            if var_name=='UICE':
                                zonal_velocity = angle_cos * u_var_grid[:, n] - angle_sin * v_var_grid[:, n]
                                values[index_counter:index_counter + (end_file_index - start_file_index),boundary_points_counted+n] = zonal_velocity
                            if var_name=='VICE':
                                meridional_velocity = angle_sin * u_var_grid[:, n] + angle_cos * v_var_grid[:, n]
                                values[index_counter:index_counter + (end_file_index - start_file_index),boundary_points_counted+n] = meridional_velocity
                        else:
                            values[index_counter:index_counter + (end_file_index - start_file_index), boundary_points_counted+n] = var_grid[:,n]

                boundary_points_counted += N
            else:
                if var_name in ['UVEL','VVEL']:
                    timesteps_in_file = int(np.size(u_var_grid) / (Nr * N))
                    u_var_grid = np.reshape(u_var_grid, (timesteps_in_file, Nr, N))
                    u_var_grid = u_var_grid[start_file_index:end_file_index, :, :]
                    v_var_grid = np.reshape(v_var_grid, (timesteps_in_file, Nr, N))
                    v_var_grid = v_var_grid[start_file_index:end_file_index, :, :]
                else:
                    timesteps_in_file = int(np.size(var_grid) / (Nr * N))
                    var_grid = np.reshape(var_grid, (timesteps_in_file, Nr, N))
                    var_grid = var_grid[start_file_index:end_file_index, :, :]

                for n in range(N):
                    if faces[n] != 0:
                        points[boundary_points_counted+n, 0] = ecco_XC_faces[faces[n]][rows[n], cols[n]]
                        points[boundary_points_counted+n, 1] = ecco_YC_faces[faces[n]][rows[n], cols[n]]
                        hfac_points[:, boundary_points_counted+n] = ecco_hFacC_faces[faces[n]][:, rows[n], cols[n]]
                        if var_name in ['UVEL', 'VVEL']:
                            angle_cos = ecco_AngleCS_faces[faces[n]][rows[n], cols[n]]
                            angle_sin = ecco_AngleSN_faces[faces[n]][rows[n], cols[n]]
                            if var_name=='UVEL':
                                zonal_velocity = np.zeros((np.shape(u_var_grid)[0],np.shape(u_var_grid)[1]))
                                for k in range(np.shape(zonal_velocity)[1]):
                                    zonal_velocity[:,k] = angle_cos * u_var_grid[:, k, n] - angle_sin * v_var_grid[:, k, n]
                                values[index_counter:index_counter + (end_file_index - start_file_index), :,boundary_points_counted+n] = zonal_velocity
                            if var_name=='VVEL':
                                meridional_velocity = np.zeros((np.shape(u_var_grid)[0],np.shape(u_var_grid)[1]))
                                for k in range(np.shape(meridional_velocity)[1]):
                                    meridional_velocity[:,k] = angle_sin * u_var_grid[:, k, n] + angle_cos * v_var_grid[:, k, n]
                                values[index_counter:index_counter + (end_file_index - start_file_index), :,boundary_points_counted+n] = meridional_velocity
                        else:
                            values[index_counter:index_counter + (end_file_index - start_file_index), :, boundary_points_counted+n] = var_grid[:, :, n]

                boundary_points_counted += N

        index_counter += (end_file_index - start_file_index)

    return(points,values,hfac_points)

def subset_tile_geometry_to_boundary(boundary,  XC, YC, AngleCS, AngleSN, hFac):

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
                  ecco_dir, llc, n_timesteps_per_day, boundary_read_group_dict,
                  start_year, final_year, start_month,final_month, start_day, final_day, print_level, read_darwin = False):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import ecco_functions as ef

    if 'obcs' not in os.listdir(os.path.join(config_dir, 'L1_grid', model_name, 'input')):
        os.mkdir(os.path.join(config_dir, 'L1_grid', model_name, 'input', 'obcs'))

    if print_level>=1:
        print('    - Creating the '+var_name+' BC files for the '+model_name+' model from ECCOv5 data')

    # llc = 270
    n_timesteps = n_timesteps_per_day

    # step 0: get the model domain
    if print_level >= 1:
        print('    - Reading in the model geometry')
    XC, YC, AngleCS, AngleSN, hfac, delR = read_grid_geometry(config_dir,model_name,var_name)
    Nr = len(delR)

    # step 1: get the ecco faces geometry
    if print_level >= 1:
        print('    - Reading in the ECCO geometry')
    ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces = \
        ef.read_ecco_geometry_to_faces(ecco_dir, llc)

    if read_darwin:
        L0_run_dir = os.path.join(config_dir, 'L0', 'run_darwin')
    else:
        L0_run_dir = os.path.join(config_dir, 'L0', 'run')

    ####################################################################################################################
    # Loop through the boundaries

    for boundary in ['south']:#,'west','east','north']:

        if print_level >= 2:
            print('        - Running the downscale routine for the '+boundary+' boundary')

        if boundary not in os.listdir(os.path.join(config_dir,'L1_grid',model_name,'input','obcs')):
            os.mkdir(os.path.join(config_dir,'L1_grid',model_name,'input','obcs',boundary))
        if var_name not in os.listdir(os.path.join(config_dir, 'L1_grid', model_name, 'input', 'obcs', boundary)):
            os.mkdir(os.path.join(config_dir, 'L1_grid', model_name, 'input', 'obcs', boundary, var_name))

        dest_files, source_file_read_dict, source_file_read_index_sets = \
            create_src_dest_dicts_from_ref(config_dir, model_name, boundary, var_name,
                                           start_year, final_year, start_month,
                                           final_month, start_day, final_day, read_darwin)

        # get the points just on the boundary
        XC_subset, YC_subset, AngleCS_subset, AngleSN_subset, hFac_subset = subset_tile_geometry_to_boundary(boundary,  XC, YC, AngleCS, AngleSN, hfac)

        # these convert to the old reference, which was stored as faces
        boundary_read_groups = boundary_read_group_dict[boundary]

        if print_level >= 3:
            print('            - Reading the mask to determine the number of dv points')

        n_dv_points = 0
        nc_dict_file = os.path.join(config_dir, 'L0', 'input', 'L0_dv_mask_reference_dict.nc')
        ds = nc4.Dataset(nc_dict_file)

        for b in boundary_read_groups:
            boundary_group = '_'.join(model_name.split('_')[:2]) + '_' + b
            grp = ds.groups[boundary_group]
            faces = grp.variables['source_faces'][:]
            n_dv_points += len(faces)
        ds.close()

        if print_level >= 6:
            print('                   - n_dv_points: '+str(n_dv_points))

        ############################################################################################################
        # Loop through the destination files

        for dest_file in dest_files:
            if dest_file not in os.listdir(os.path.join(config_dir, 'L1_grid', model_name, 'input', 'obcs', boundary, var_name)):
                if print_level >= 3:
                    print('            - Downscaling the timesteps to be stored in file ' + str(dest_file))
                source_files = source_file_read_dict[dest_file]
                source_file_read_indices = source_file_read_index_sets[dest_file]

                if print_level >= 4:
                    print('                - Reading in the L0 diagnostics_vec output')
                L0_boundary_points, L0_boundary_values, L0_boundary_point_hFacC = \
                    read_L0_boundary_variable_points(config_dir, L0_run_dir, model_name, boundary, var_name,
                                                     source_files, source_file_read_indices,
                                                     llc, Nr, boundary_read_groups,n_dv_points,
                                                     ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces,
                                                     print_level)

                # plt.plot(L0_boundary_points[:,0],L0_boundary_points[:,1], 'k.')
                # plt.plot(XC_subset,YC_subset,'g.')
                # plt.title(dest_file)
                # plt.show()

                if var_name in ['AREA', 'HEFF', 'HSNOW', 'UICE', 'VICE', 'ETAN']:
                    if boundary in ['north', 'south']:
                        output_grid = np.zeros((n_timesteps, np.shape(XC)[1]))
                    else:
                        output_grid = np.zeros((n_timesteps, np.shape(XC)[0]))
                else:
                    if boundary in ['north', 'south']:
                        output_grid = np.zeros((n_timesteps, Nr, np.shape(XC)[1]))
                    else:
                        output_grid = np.zeros((n_timesteps, Nr, np.shape(XC)[0]))

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

                ####################################################################################################################
                if print_level >= 4:
                    print('                - Downscaling data')

                mask_subset = np.copy(hFac_subset)
                mask_subset[mask_subset>0]=1
                if var_name in ['AREA', 'HEFF', 'HSNOW', 'UICE', 'VICE', 'ETAN']:
                    mask_subset = mask_subset[:1, :, :]

                # plt.plot(L0_boundary_points[:, 0], L0_boundary_points[:, 1], 'g.')
                # plt.plot(XC_subset.ravel(),YC_subset.ravel(),'k-')
                # plt.show()

                plt.imshow(mask_subset[:,0,:])
                plt.show()

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
                            output_grid[timestep, :] = interp_field[0, 0, :]
                        else:
                            output_grid[timestep, :] = interp_field[0, :, 0]
                    else:
                        if boundary in ['north', 'south']:
                            output_grid[timestep, :, :] = interp_field[:,0,:]
                        else:
                            output_grid[timestep, :, :] = interp_field[:,:,0]

                # if var_name in ['AREA', 'HEFF', 'HSNOW', 'UICE', 'VICE', 'ETAN']:
                #     plt.plot(output_grid[0, :])
                #     plt.show()
                # else:
                #     plt.imshow(output_grid[0, :, :],cmap='plasma')
                #     plt.show()

                # if var_name in ['UVEL','VVEL','UICE','VICE']:
                #     output_file = os.path.join(config_dir, 'L1_grid', model_name, 'input', 'obcs', boundary, var_name,dest_file[:-4]+'_rotated.bin')
                # else:
                output_file = os.path.join(config_dir, 'L1_grid', model_name, 'input', 'obcs', boundary, var_name, dest_file)
                output_grid.ravel(order='C').astype('>f4').tofile(output_file)





