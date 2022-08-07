
import os
#import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import netCDF4 as nc4
import argparse
import ast
import sys

def create_src_dest_dicts_from_ref(config_dir, L05_model_name, var_name, start_year, final_year, start_month, final_month, start_day, final_day):
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

    f = open(os.path.join(config_dir,'L05', L05_model_name, 'run','dv', 'daily_exf_dest_ref.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)

    dest_files_out = []
    source_file_read_dict = {}
    source_file_read_index_sets = {}

    for dest_file in dest_files:
        dest_files_out.append('L2_exf_'+var_name+'.'+dest_file+'.bin')
        source_files = []
        index_set = []
        for pair in size_dict[dest_file]:
            source_files.append('L2_surface_mask_'+var_name+'.'+pair[0]+'.bin')
            index_set.append(pair[1])
        source_file_read_dict['L2_exf_'+var_name+'.'+dest_file+'.bin'] = source_files
        source_file_read_index_sets['L2_exf_'+var_name+'.'+dest_file+'.bin'] = index_set

    return(dest_files_out, source_file_read_dict, source_file_read_index_sets)


def read_L05_grid_tile_geometry(config_dir,model_name, sNx, sNy,
                               ordered_nonblank_tiles, ordered_nonblank_rotations,
                               faces, ordered_tiles_faces_dict):

    XC_faces = {}
    YC_faces = {}
    Mask_faces = {}
    for face in faces:
        face_XC_grid = np.zeros((sNy*len(ordered_tiles_faces_dict[face]),
                                 sNx*len(ordered_tiles_faces_dict[face][0])))
        face_YC_grid = np.zeros((sNy * len(ordered_tiles_faces_dict[face]),
                                 sNx * len(ordered_tiles_faces_dict[face][0])))
        face_Mask_grid = np.zeros((sNy * len(ordered_tiles_faces_dict[face]),
                                   sNx * len(ordered_tiles_faces_dict[face][0])))
        XC_faces[face] = face_XC_grid
        YC_faces[face] = face_YC_grid
        Mask_faces[face] = face_Mask_grid

    # stitched_XC_grid = np.zeros((sNy*len(ordered_tiles_faces_dict), sNx*len(ordered_nonblank_tiles[1])))
    # stitched_YC_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[1])))

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L05', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
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
                                    Mask_faces[face][fr * sNy:(fr + 1) * sNy, fc * sNx: (fc + 1) * sNx] = Mask

                    # for i in range(ordered_nonblank_rotations[r][c]):
                    #     XC = np.rot90(XC)
                    #     YC = np.rot90(YC)
                    #
                    # stitched_XC_grid[r * sNy:(r + 1) * sNy,c * sNx: (c + 1) * sNx] = XC
                    # stitched_YC_grid[r * sNy:(r + 1) * sNy,c * sNx: (c + 1) * sNx] = YC

    return(XC_faces, YC_faces, Mask_faces)


def read_grid_geometry_from_nc(config_dir,model_name,var_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    if var_name=='UWIND':
        hfac = ds.variables['HFacW'][:,:,:]
        hfac = hfac[:,:,:-1]
    elif var_name=='VWIND':
        hfac = ds.variables['HFacS'][:,:,:]
        hfac = hfac[:, :-1, :]
    else:
        hfac = ds.variables['HFacC'][:, :, :]
    ds.close()
    hfac = hfac[0,:,:]
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


def read_exf_variable_to_L05_points(L05_run_dir, var_name,
                                        source_files, source_file_read_indices,
                                        source_faces,source_rows,source_cols,
                                        L05_XC_faces, L05_YC_faces, L05_Mask_faces):

    n_Timesteps = 0
    for index_set in source_file_read_indices:
        n_Timesteps += int(index_set[1]-index_set[0])+1

    print('       + the L05 grid for this file will have '+str(n_Timesteps)+' timesteps')

    points = np.zeros((len(source_faces),2))
    values = np.zeros((n_Timesteps, len(source_faces)))
    wet_grid = np.zeros((len(source_faces),))

    for i in range(len(source_faces)):
        x = L05_XC_faces[source_faces[i]][source_rows[i], source_cols[i]]
        y = L05_YC_faces[source_faces[i]][source_rows[i], source_cols[i]]
        points[i,0] = x
        points[i,1] = y
        wet_grid[i] = L05_Mask_faces[source_faces[i]][source_rows[i], source_cols[i]]

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

        var_file = os.path.join(L05_run_dir, 'dv', 'L2_surface_mask_' + var_name + '.' + file_suffix)
        var_grid = np.fromfile(var_file, dtype='>f4')
        timesteps_in_file = int(np.size(var_grid) / (N))
        var_grid = np.reshape(var_grid, (timesteps_in_file, N))

        values[index_counter:index_counter + (end_file_index-start_file_index)+1,:] = var_grid[start_file_index:end_file_index+1,:]

        index_counter += (end_file_index-start_file_index)

    return(points,values,wet_grid)


def downscale_L05_exf_field_to_L2(df, var_name,
                                 L05_points, L05_values,
                                 L05_wet_grid_points, L05_wet_grid_on_L2_points,
                                 L2_XC, L2_YC, L2_wet_grid):

    nTimestepsOut = np.shape(L05_values)[0]

    L2_exf_var = np.zeros((nTimestepsOut, np.shape(L2_YC)[0], np.shape(L2_YC)[1]))

    # print('        -  Variable shapes entering downscale routine:')
    # print('          -  L05_point: '+str(np.shape(L05_points)))
    # print('          -  L05_values: ' + str(np.shape(L05_values)))
    # print('          -  L05_wet_grid: ' + str(np.shape(L05_wet_grid_points)))
    # # print('          -  L05_wet_grid_on_L2_subset: ' + str(np.shape(L05_wet_grid_on_L2_3D_points)))
    # print('          -  L2_XC: ' + str(np.shape(L2_XC)))
    # print('          -  L2_YC: ' + str(np.shape(L2_YC)))
    # print('          -  L2_wet_grid: ' + str(np.shape(L2_wet_grid)))

    if var_name in ['SWDOWN','PRECIP','RUNOFF']:
        remove_zeros = False
    else:
        remove_zeros = True

    for timestep in range(nTimestepsOut):
        # if timestep%5==0:
        #     print('          Working on timestep '+str(timestep)+' of '+str(nTimestepsOut)+' for '+var_name+' on mask '+mask_name)

        remove_zeros = True
        downscaled_field = df.downscale_exf_field(L05_points, L05_values[timestep, :],
                                                  L05_wet_grid_points,
                                                  L2_XC, L2_YC, L2_wet_grid,
                                                  remove_zeros=remove_zeros)
        L2_exf_var[timestep, :, :] = downscaled_field

    # plt.imshow(L2_exf_var[0,:,:],origin='lower')
    # plt.show()

    return(L2_exf_var)


########################################################################################################################


def create_exf_field(config_dir, L05_model_name, L2_model_name,
                    sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                    faces, ordered_tiles_faces_dict,
                    var_name, dest_files, source_file_read_dict, source_file_read_index_sets):

    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import downscale_functions as df

    # this is the dir where the exf output will be stored
    if 'exf' not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,'input','exf'))
    output_dir = os.path.join(config_dir,'L2',L2_model_name,'input','exf')

    if var_name not in os.listdir(output_dir):
        os.mkdir(os.path.join(output_dir,var_name))

    print('    - Reading in the geometry of the L05 domain')
    L05_XC_faces, L05_YC_faces, L05_Mask_faces = read_L05_grid_tile_geometry(config_dir,L05_model_name,sNx, sNy,
                                                                        ordered_nonblank_tiles, ordered_nonblank_rotations,
                                                                        faces, ordered_tiles_faces_dict)

    # # for face in faces:
    # #     plt.subplot(1,2,1)
    # #     C = plt.imshow(L05_XC_faces[face],origin='lower')
    # #     plt.colorbar(C)
    # #     plt.subplot(1, 2, 2)
    # #     C = plt.imshow(L05_YC_faces[face], origin='lower')
    # #     plt.colorbar(C)
    # #     plt.show()

    L2_XC, L2_YC, L2_wet_grid = read_grid_geometry_from_nc(config_dir,L2_model_name,var_name)

    print('    - Reading the mask to reference the variable to the llc grid')
    nc_dict_file = os.path.join(config_dir,'L05',L05_model_name,'namelist','L05_dv_mask_reference_dict.nc')
    source_faces,source_rows,source_cols = read_mask_reference_from_nc_dict(nc_dict_file, 'surface')

    for dest_file in dest_files:
        if dest_file not in []:#os.listdir(os.path.join(output_dir,var_name)):
            # try:
            print('    - Downscaling the timesteps to be stored in file ' + str(dest_file))
            source_files = source_file_read_dict[dest_file]
            source_file_read_indices = source_file_read_index_sets[dest_file]

            print('    - Reading in the L05 diagnostics_vec output')
            L05_run_dir = os.path.join(config_dir,'L05',L05_model_name,'run')
            L05_points, L05_values, L05_wet_grid_points = read_exf_variable_to_L05_points(L05_run_dir, var_name,
                                                                       source_files, source_file_read_indices,
                                                                       source_faces, source_rows, source_cols,
                                                                       L05_XC_faces, L05_YC_faces, L05_Mask_faces)

            L05_wet_grid_on_L2_points = np.copy(L05_wet_grid_points)

            # print('    - Downscaling the output to the new boundary')
            L2_exf_var = downscale_L05_exf_field_to_L2(df, var_name,
                                                           L05_points, L05_values,
                                                           L05_wet_grid_points, L05_wet_grid_on_L2_points,
                                                           L2_XC, L2_YC, L2_wet_grid)

            output_file = os.path.join(output_dir, var_name, dest_file)
            L2_exf_var.ravel(order='C').astype('>f4').tofile(output_file)

        else:
            print('    - Skipping ' + str(dest_file) + ' because it was already created')

def create_exf_fields_via_interpolation(config_dir, L05_model_name, L2_model_name, proc_id,
                                       sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                                       faces, ordered_tiles_faces_dict,
                                       start_year, final_year, start_month,
                                       final_month, start_day, final_day):

    var_name_list = ['ATEMP','AQH','LWDOWN','SWDOWN','UWIND','VWIND','PRECIP','RUNOFF']
    var_name = var_name_list[proc_id % len(var_name_list)]

    print('Creating the exf fields for ' + var_name + ' cover year/month/days ' +
          str(start_year)+'/'+str(start_month)+'/'+str(start_day) + ' to ' +
          str(final_year)+'/'+str(final_month)+'/'+str(final_day))

    dest_files, source_file_read_dict, source_file_read_index_sets = \
        create_src_dest_dicts_from_ref(config_dir, L05_model_name, var_name, start_year, final_year, start_month, final_month, start_day, final_day)

    # print(dest_files, source_file_read_dict, source_file_read_index_sets)

    print('  Running the Downscale routine:')
    create_exf_field(config_dir, L05_model_name, L2_model_name,
                    sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                    faces, ordered_tiles_faces_dict,
                    var_name, dest_files,
                    source_file_read_dict, source_file_read_index_sets)

