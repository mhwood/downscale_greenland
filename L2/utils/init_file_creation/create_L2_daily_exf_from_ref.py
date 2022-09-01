
import os
#import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import netCDF4 as nc4
import argparse
import ast
import sys

def create_src_dest_dicts_from_ref(config_dir, L1_model_name, var_name, start_year, final_year, start_month, final_month, start_day, final_day):
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

    f = open(os.path.join(config_dir,'L1', L1_model_name, 'run','dv', 'exf_dest_ref.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)

    dest_files_out = []
    source_file_read_dict = {}
    source_file_read_index_sets = {}

    if var_name in ['UWIND','VWIND']:
        suffix = '_rotated.bin'
    else:
        suffix = '.bin'

    for dest_file in dest_files:
        dest_files_out.append('L2_exf_'+var_name+'.'+dest_file+suffix)
        source_files = []
        index_set = []
        for pair in size_dict[dest_file]:
            source_files.append('L2_surface_mask_'+var_name+'.'+pair[0]+'.bin')
            index_set.append(pair[1])
        source_file_read_dict['L2_exf_'+var_name+'.'+dest_file+suffix] = source_files
        source_file_read_index_sets['L2_exf_'+var_name+'.'+dest_file+suffix] = index_set

    return(dest_files_out, source_file_read_dict, source_file_read_index_sets)


def read_L1_grid_tile_geometry(config_dir,model_name, sNx, sNy,
                               ordered_nonblank_tiles, ordered_nonblank_rotations,
                               faces, ordered_tiles_faces_dict):

    XC_faces = {}
    YC_faces = {}
    AngleCS_faces = {}
    AngleSN_faces = {}
    Mask_faces = {}
    for face in faces:
        face_XC_grid = np.zeros((sNy*len(ordered_tiles_faces_dict[face]),
                                 sNx*len(ordered_tiles_faces_dict[face][0])))
        face_YC_grid = np.zeros((sNy * len(ordered_tiles_faces_dict[face]),
                                 sNx * len(ordered_tiles_faces_dict[face][0])))
        face_AngleCS_grid = np.zeros((sNy * len(ordered_tiles_faces_dict[face]),
                                 sNx * len(ordered_tiles_faces_dict[face][0])))
        face_AngleSN_grid = np.zeros((sNy * len(ordered_tiles_faces_dict[face]),
                                 sNx * len(ordered_tiles_faces_dict[face][0])))
        face_Mask_grid = np.zeros((sNy * len(ordered_tiles_faces_dict[face]),
                                   sNx * len(ordered_tiles_faces_dict[face][0])))
        XC_faces[face] = face_XC_grid
        YC_faces[face] = face_YC_grid
        AngleCS_faces[face] = face_AngleCS_grid
        AngleSN_faces[face] = face_AngleSN_grid
        Mask_faces[face] = face_Mask_grid

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
                    AngleCS = ds.variables['AngleCS'][:, :]
                    AngleSN = ds.variables['AngleSN'][:, :]
                    Mask = ds.variables['HFacC'][:,:,:]
                    Mask = Mask[0,:,:]
                    Mask[Mask>0]=1
                    ds.close()

                    for face in faces:
                        for fr in range(len(ordered_tiles_faces_dict[face])):
                            for fc in range(len(ordered_tiles_faces_dict[face][0])):
                                if ordered_tiles_faces_dict[face][fr][fc] == tile_number:
                                    XC_faces[face][fr * sNy:(fr + 1) * sNy, fc * sNx: (fc + 1) * sNx] = XC
                                    YC_faces[face][fr * sNy:(fr + 1) * sNy, fc * sNx: (fc + 1) * sNx] = YC
                                    AngleCS_faces[face][fr * sNy:(fr + 1) * sNy, fc * sNx: (fc + 1) * sNx] = AngleCS
                                    AngleSN_faces[face][fr * sNy:(fr + 1) * sNy, fc * sNx: (fc + 1) * sNx] = AngleSN
                                    Mask_faces[face][fr * sNy:(fr + 1) * sNy, fc * sNx: (fc + 1) * sNx] = Mask

                    # for i in range(ordered_nonblank_rotations[r][c]):
                    #     XC = np.rot90(XC)
                    #     YC = np.rot90(YC)
                    #
                    # stitched_XC_grid[r * sNy:(r + 1) * sNy,c * sNx: (c + 1) * sNx] = XC
                    # stitched_YC_grid[r * sNy:(r + 1) * sNy,c * sNx: (c + 1) * sNx] = YC

    return(XC_faces, YC_faces, AngleCS_faces, AngleSN_faces, Mask_faces)


def read_grid_geometry_from_nc(config_dir,model_name,var_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    AngleSN = ds.variables['AngleSN'][:, :]
    if var_name=='UWIND':
        hfac = ds.variables['HFacW'][:,:,:]
        hfac = hfac[:,:,:-1]
    elif var_name=='VWIND':
        hfac = ds.variables['HFacS'][:,:,:]
        hfac = hfac[:, :-1, :]
    else:
        hfac = ds.variables['HFacC'][:, :, :]
    ds.close()
    hfac = hfac[:1,:,:]
    mask = np.copy(hfac)
    mask[mask>0]=1
    return(XC,YC,AngleCS,AngleSN,mask)


def read_mask_reference_from_nc_dict(nc_dict_file,mask_name):
    ds = nc4.Dataset(nc_dict_file)

    grp = ds.groups[mask_name]
    source_faces = grp.variables['source_faces'][:]
    source_rows = grp.variables['source_rows'][:]
    source_cols = grp.variables['source_cols'][:]

    ds.close()
    return(source_faces,source_rows,source_cols)


def read_exf_variable_to_L1_points(L1_run_dir, var_name,
                                        source_files, source_file_read_indices,
                                        source_faces,source_rows,source_cols,
                                        L1_XC_faces, L1_YC_faces, L1_AngleCS_faces, L1_AngleSN_faces, L1_Mask_faces,print_level):

    n_Timesteps = 0
    for index_set in source_file_read_indices:
        n_Timesteps += int(index_set[1]-index_set[0])+1

    if print_level>=4:
        print('                - the L1 grid for this file will have '+str(n_Timesteps)+' timesteps')

    points = np.zeros((len(source_faces),2))
    angles = np.zeros((len(source_faces),2))
    values = np.zeros((n_Timesteps, 1, len(source_faces)))
    wet_grid = np.zeros((1,len(source_faces)))

    for i in range(len(source_faces)):
        x = L1_XC_faces[source_faces[i]][source_rows[i], source_cols[i]]
        y = L1_YC_faces[source_faces[i]][source_rows[i], source_cols[i]]
        points[i,0] = x
        points[i,1] = y
        angles[i, 0] = L1_AngleCS_faces[source_faces[i]][source_rows[i], source_cols[i]]
        angles[i, 1] = L1_AngleSN_faces[source_faces[i]][source_rows[i], source_cols[i]]
        wet_grid[0,i] = L1_Mask_faces[source_faces[i]][source_rows[i], source_cols[i]]

    index_counter = 0
    for s in range(len(source_files)):
        source_file = source_files[s]
        file_suffix = '.'.join(source_file.split('.')[-2:])
        index_set = source_file_read_indices[s]
        if print_level >= 3:
            print('            - Adding timesteps ' + str(index_set[0]) + ' to ' + str(index_set[1]) + ' from file ' + source_file)

        start_file_index = int(index_set[0])
        end_file_index = int(index_set[1])

        if print_level >= 3:
            print('            - Storing at points '+str(index_counter)+' to '+
              str(index_counter+(end_file_index-start_file_index))+' in the grid')

        N = len(source_faces)

        if var_name in ['UWIND','VWIND']:
            u_var_file = os.path.join(L1_run_dir, 'dv', 'L2_surface_mask_UWIND.' + file_suffix)
            u_var_grid = np.fromfile(u_var_file, dtype='>f4')
            timesteps_in_file = int(np.size(u_var_grid) / (N))
            # u_var_grid = np.reshape(u_var_grid, (timesteps_in_file, N))
            u_var_grid = np.reshape(u_var_grid, (timesteps_in_file, N + 109))
            u_var_grid = u_var_grid[:, :-109]

            v_var_file = os.path.join(L1_run_dir, 'dv', 'L2_surface_mask_VWIND.' + file_suffix)
            v_var_grid = np.fromfile(v_var_file, dtype='>f4')
            timesteps_in_file = int(np.size(v_var_grid) / (N))
            # v_var_grid = np.reshape(v_var_grid, (timesteps_in_file, N))
            v_var_grid = np.reshape(v_var_grid, (timesteps_in_file, N + 109))
            v_var_grid = v_var_grid[:, :-109]

            if print_level >= 3:
                print('             - Rotating vectors to \"natural\" coords')
            var_grid = np.zeros_like(u_var_grid)
            for timestep in range(np.shape(var_grid)[0]):
                if var_name=='UWIND':
                    var_grid[timestep,:] = angles[:, 0] * u_var_grid[timestep, :] - angles[:, 1] * v_var_grid[timestep, :]
                if var_name=='VWIND':
                    var_grid[timestep,:] = angles[:, 1] * u_var_grid[timestep, :] + angles[:, 0] * v_var_grid[timestep, :]

        else:
            var_file = os.path.join(L1_run_dir, 'dv', 'L2_surface_mask_' + var_name + '.' + file_suffix)
            var_grid = np.fromfile(var_file, dtype='>f4')
            timesteps_in_file = int(np.size(var_grid) / (N))
            # var_grid = np.reshape(var_grid, (timesteps_in_file, N))
            print(timesteps_in_file, N)
            var_grid = np.reshape(var_grid, (timesteps_in_file, N+109))
            var_grid = var_grid[:,:-109]

            # plt.plot(var_grid[-1,:])
            # plt.show()
            # raise ValueError('Stop')

        values[index_counter:index_counter + (end_file_index-start_file_index)+1,0,:] = var_grid[start_file_index:end_file_index+1,:]

        index_counter += (end_file_index-start_file_index)

    return(points,values,wet_grid)


def downscale_L1_exf_field_to_L2(df, var_name,
                                 L1_points, L1_values,
                                 L1_wet_grid_points, L1_wet_grid_on_L2_points,
                                 L2_XC, L2_YC, L2_wet_grid, print_level):

    nTimestepsOut = np.shape(L1_values)[0]

    L2_exf_var = np.zeros((nTimestepsOut, np.shape(L2_YC)[0], np.shape(L2_YC)[1]))

    if print_level >= 4:
        print('            -  Variable shapes entering downscale routine:')
        print('                -  L1_point: '+str(np.shape(L1_points)))
        print('                -  L1_values: ' + str(np.shape(L1_values)))
        print('                -  L1_wet_grid: ' + str(np.shape(L1_wet_grid_points)))
        # print('          -  L1_wet_grid_on_L2_subset: ' + str(np.shape(L1_wet_grid_on_L2_3D_points)))
        print('                -  L2_XC: ' + str(np.shape(L2_XC)))
        print('                -  L2_YC: ' + str(np.shape(L2_YC)))
        print('                -  L2_wet_grid: ' + str(np.shape(L2_wet_grid)))

    for timestep in range(nTimestepsOut):
        # if timestep%5==0:
        #     print('          Working on timestep '+str(timestep)+' of '+str(nTimestepsOut)+' for '+var_name+' on mask '+mask_name)

        remove_zeros = True
        downscaled_field = df.downscale_3D_points_with_zeros(L1_points, L1_values[timestep, :],
                                                  L1_wet_grid_points,
                                                  L2_XC, L2_YC, L2_wet_grid)
        L2_exf_var[timestep, :, :] = downscaled_field

    # plt.imshow(L2_exf_var[0,:,:],origin='lower')
    # plt.show()

    return(L2_exf_var)


########################################################################################################################


def create_exf_field(config_dir, L1_model_name, L2_model_name,
                    sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                    faces, ordered_tiles_faces_dict,
                    var_name, dest_files, source_file_read_dict, source_file_read_index_sets, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import downscale_functions as df

    # this is the dir where the exf output will be stored
    if 'exf' not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,'input','exf'))
    output_dir = os.path.join(config_dir,'L2',L2_model_name,'input','exf')

    if var_name not in os.listdir(output_dir):
        os.mkdir(os.path.join(output_dir,var_name))

    if print_level>=2:
        print('        - Reading in the geometry of the L1 domain')
    L1_XC_faces, L1_YC_faces, L1_AngleCS_faces, L1_AngleSN_faces, L1_Mask_faces = read_L1_grid_tile_geometry(config_dir,L1_model_name,sNx, sNy,
                                                                        ordered_nonblank_tiles, ordered_nonblank_rotations,
                                                                        faces, ordered_tiles_faces_dict)

    # # for face in faces:
    # #     plt.subplot(1,2,1)
    # #     C = plt.imshow(L1_XC_faces[face],origin='lower')
    # #     plt.colorbar(C)
    # #     plt.subplot(1, 2, 2)
    # #     C = plt.imshow(L1_YC_faces[face], origin='lower')
    # #     plt.colorbar(C)
    # #     plt.show()

    L2_XC, L2_YC, _, _, L2_wet_grid = read_grid_geometry_from_nc(config_dir,L2_model_name,var_name)

    if print_level >= 2:
        print('        - Reading the mask to reference the variable to the llc grid')
    nc_dict_file = os.path.join(config_dir,'L1',L1_model_name,'input','L1_dv_mask_reference_dict.nc')
    source_faces,source_rows,source_cols = read_mask_reference_from_nc_dict(nc_dict_file, 'L2_surface')
    print(np.min(source_faces),np.min(source_rows),np.min(source_cols))
    print(np.max(source_faces), np.max(source_rows), np.max(source_cols))

    for dest_file in dest_files:
        if dest_file not in []:#os.listdir(os.path.join(output_dir,var_name)):
            # try:
            if print_level >= 3:
                print('            - Downscaling the timesteps to be stored in file ' + str(dest_file))
            source_files = source_file_read_dict[dest_file]
            source_file_read_indices = source_file_read_index_sets[dest_file]

            if print_level >= 3:
                print('            - Reading in the L1 diagnostics_vec output')
            L1_run_dir = os.path.join(config_dir,'L1',L1_model_name,'run')
            L1_points, L1_values, L1_wet_grid_points = read_exf_variable_to_L1_points(L1_run_dir, var_name,
                                                                       source_files, source_file_read_indices,
                                                                       source_faces, source_rows, source_cols,
                                                                       L1_XC_faces, L1_YC_faces, L1_AngleCS_faces, L1_AngleSN_faces,
                                                                       L1_Mask_faces, print_level)

            L1_wet_grid_on_L2_points = np.copy(L1_wet_grid_points)

            # print('    - Downscaling the output to the new boundary')
            L2_exf_var = downscale_L1_exf_field_to_L2(df, var_name,
                                                           L1_points, L1_values,
                                                           L1_wet_grid_points, L1_wet_grid_on_L2_points,
                                                           L2_XC, L2_YC, L2_wet_grid, print_level)

            output_file = os.path.join(output_dir, var_name, dest_file)
            L2_exf_var.ravel(order='C').astype('>f4').tofile(output_file)

        else:
            if print_level >= 3:
                print('            - Skipping ' + str(dest_file) + ' because it was already created')

def create_exf_fields_via_interpolation(config_dir, L2_model_name, L1_model_name, proc_id,
                                        sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                                        faces, ordered_tiles_faces_dict,
                                        start_year, final_year, start_month,
                                        final_month, start_day, final_day, print_level):

    var_name_list = ['ATEMP','AQH','LWDOWN','SWDOWN','UWIND','VWIND','PRECIP','RUNOFF']
    var_name = var_name_list[proc_id % len(var_name_list)]

    if print_level>=1:
        print('    - Creating the exf fields for ' + var_name + ' cover year/month/days ' +
          str(start_year)+'/'+str(start_month)+'/'+str(start_day) + ' to ' +
          str(final_year)+'/'+str(final_month)+'/'+str(final_day))

    dest_files, source_file_read_dict, source_file_read_index_sets = \
        create_src_dest_dicts_from_ref(config_dir, L1_model_name, var_name, start_year, final_year, start_month, final_month, start_day, final_day)

    # print(dest_files, source_file_read_dict, source_file_read_index_sets)

    if print_level >= 1:
        print('    - Running the Downscale routine:')
    create_exf_field(config_dir, L1_model_name, L2_model_name,
                    sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                    faces, ordered_tiles_faces_dict,
                    var_name, dest_files,
                    source_file_read_dict, source_file_read_index_sets, print_level)

