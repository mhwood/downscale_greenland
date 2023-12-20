
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
#from MITgcmutils import mds
#from scipy.interpolate import griddata
import sys
from datetime import datetime
import ast

def read_grid_geometry(config_dir,model_name,var_name):

    file_path = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    AngleSN = ds.variables['AngleSN'][:, :]

    if var_name in ['VWIND']:
        hFac = 'S'
    elif var_name in ['UWIND']:
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

def create_src_dest_dicts_from_ref(config_dir, model_name, var_name,
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
        f = open(os.path.join(config_dir,'L0', 'run_darwin','dv',prefix, model_name+'_exf_dest_ref.txt'))
    else:
        f = open(os.path.join(config_dir, 'L0', 'run', 'dv', prefix, model_name + '_exf_dest_ref.txt'))
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

def create_L1_exfs(config_dir, model_name, var_name, ecco_dir, llc,
                   start_year, final_year, start_month, final_month, start_day, final_day, print_level, read_darwin = False):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import ecco_functions as ef

    if 'exf' not in os.listdir(os.path.join(config_dir, 'L1_grid', model_name, 'input')):
        os.mkdir(os.path.join(config_dir, 'L1_grid', model_name, 'input', 'exf'))
    if var_name not in os.listdir(os.path.join(config_dir, 'L1_grid', model_name, 'input', 'exf')):
        os.mkdir(os.path.join(config_dir, 'L1_grid', model_name, 'input', 'exf', var_name))

    if print_level >= 1:
        print('    - Creating the ' + var_name + ' exf files for the ' + model_name + ' model from ECCOv5 output data')

    if print_level >= 1:
        print('    - Reading in the model geometry')
    L1_XC, L1_YC, L1_AngleCS, L1_AngleSN, L1_hfac, delR = read_grid_geometry(config_dir, model_name, var_name)
    L1_mask = np.copy(L1_hfac)
    L1_mask[L1_mask>0]=1
    L1_mask = L1_mask[:1,:,:]

    # step 1: get the ecco faces geometry
    if print_level >= 1:
        print('    - Reading in the ECCO geometry')
    ecco_XC_faces, ecco_YC_faces, ecco_AngleCS_faces, ecco_AngleSN_faces, ecco_hFacC_faces = \
        ef.read_ecco_geometry_to_faces(ecco_dir, llc)

    if read_darwin:
        L0_run_dir = os.path.join(config_dir, 'L0', 'run_darwin')
    else:
        L0_run_dir = os.path.join(config_dir, 'L0', 'run')

    dest_files, source_file_read_dict, source_file_read_index_sets = \
        create_src_dest_dicts_from_ref(config_dir, model_name, var_name,
                                       start_year, final_year, start_month,
                                       final_month, start_day, final_day, read_darwin)

    if print_level >= 3:
        print('            - Reading the mask to reference the variable to the llc grid')
    nc_dict_file = os.path.join(config_dir, 'L0', 'input', 'L0_dv_mask_reference_dict.nc')
    ds = nc4.Dataset(nc_dict_file)
    grp = ds.groups['_'.join(model_name.split('_')[:2])+'_surface']
    faces = grp.variables['source_faces'][:]
    rows = grp.variables['source_rows'][:]
    cols = grp.variables['source_cols'][:]
    ds.close()

    n_timesteps = 4

    ####################################################################################################################
    # Loop through the daily files

    for dest_file in dest_files:
        if dest_file not in os.listdir(os.path.join(config_dir, 'L1_grid', model_name, 'input', 'exf', var_name)):

            if print_level >= 3:
                print('            - Downscaling the timesteps to be stored in file ' + str(dest_file))

            source_files = source_file_read_dict[dest_file]
            source_file_read_indices = source_file_read_index_sets[dest_file]

            output_grid = np.zeros((n_timesteps, np.shape(L1_XC)[0], np.shape(L1_XC)[1]))

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

            # plt.imshow(output_grid[0, :, :],origin='lower')
            # plt.show()

            output_file = os.path.join(config_dir, 'L1_grid', model_name, 'input', 'exf', var_name, dest_file)
            output_grid.ravel(order='C').astype('>f4').tofile(output_file)





