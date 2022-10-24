
import os
#import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import netCDF4 as nc4
import argparse
import ast
import sys


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

def create_annual_list_of_files(year,var_name):
    if year % 4 == 0:
        n_days = 366
    else:
        n_days = 365

    day = 0
    month = 1

    output_files = []
    start_indices = []
    end_indices = []

    for d in range(n_days):
        day+=1
        if var_name not in ['UWIND','VWIND']:
            output_file = 'L3_exf_'+var_name+'_'+str(year)+'{:02d}'.format(month)+'{:02d}'.format(day)+'.bin'
        else:
            output_file = 'L3_exf_' + var_name + '_' + str(year) + '{:02d}'.format(month) + '{:02d}'.format(day) + '_rotated.bin'
        output_files.append(output_file)
        start_indices.append(d*4)
        end_indices.append((d+1)*4)

        if month in [1,3,5,7,8,10,12]:
            if day==31:
                day = 0
                month += 1
        elif month in [4,6,9,11]:
            if day == 30:
                day = 0
                month += 1
        else:
            if year%4==0:
                if day==29:
                    day = 0
                    month += 1
            else:
                if day==28:
                    day = 0
                    month += 1

    return(output_files,start_indices,end_indices)

def read_exf_variable_from_L1(L1f, L1_exf_dir, var_name, year,
                              sNx, sNy, faces, face_size_dict,
                              L1_XC_faces, L1_YC_faces, L1_AngleCS_faces, L1_AngleSN_faces,
                              L1_Mask_faces, L3_XC, L3_YC, print_level):

    L1_file_name = 'L1_exf_'+var_name+'_'+str(year)

    N = 0
    for face in faces:
        N += face_size_dict[face][0] * face_size_dict[face][1]

    XC = L1f.read_grid_faces_to_stitched_grid(L1_XC_faces, dim=2)
    YC = L1f.read_grid_faces_to_stitched_grid(L1_YC_faces, dim=2)
    AngleCS = L1f.read_grid_faces_to_stitched_grid(L1_AngleCS_faces, dim=2)
    AngleSN = L1f.read_grid_faces_to_stitched_grid(L1_AngleSN_faces, dim=2)
    Mask = L1f.read_grid_faces_to_stitched_grid(L1_Mask_faces, dim=2)

    if var_name not in ['UWIND','VWIND']:
        if print_level>=3:
            print('            - Reading '+L1_file_name)
        compact = np.fromfile(os.path.join(L1_exf_dir,L1_file_name),'>f4')
        n_rows = int(N / sNx)
        n_timesteps = int(np.size(compact) / (n_rows * sNx))
        compact = np.reshape(compact, (n_timesteps, n_rows, sNx))
        exf_field = L1f.read_compact_grid_to_stitched_grid(compact, sNx, sNy, faces, face_size_dict, dim=3)
    else:
        if var_name=='UWIND':
            u_file = L1_file_name
            v_file = L1_file_name.replace('UWIND','VWIND')
        if var_name=='VWIND':
            v_file = L1_file_name
            u_file = L1_file_name.replace('VWIND', 'UWIND')

        if print_level>=3:
            print('            - Reading '+u_file)
        u_compact = np.fromfile(os.path.join(L1_exf_dir, u_file), '>f4')
        n_rows = int(N / sNx)
        n_timesteps = int(np.size(u_compact) / (n_rows * sNx))
        u_compact = np.reshape(u_compact, (n_timesteps, n_rows, sNx))
        u_exf_field = L1f.read_compact_grid_to_stitched_grid(u_compact, sNx, sNy, faces, face_size_dict, dim=3)

        if print_level>=3:
            print('            - Reading '+v_file)
        v_compact = np.fromfile(os.path.join(L1_exf_dir, v_file), '>f4')
        v_compact = np.reshape(v_compact, (n_timesteps, n_rows, sNx))
        v_exf_field = L1f.read_compact_grid_to_stitched_grid(v_compact, sNx, sNy, faces, face_size_dict, dim=3)

        exf_field = np.zeros_like(u_exf_field)
        for i in range(np.shape(u_exf_field)[0]):
            if var_name=='UWIND':
                exf_field[i,:,:] = AngleCS * u_exf_field[i,:,:] - AngleSN * v_exf_field[i,:,:]
            if var_name == 'VWIND':
                exf_field[i,:,:] = AngleSN * u_exf_field[i,:,:] + AngleCS * v_exf_field[i,:,:]

        del u_compact
        del v_compact
        del u_exf_field
        del v_exf_field

    if print_level>=4:
        print('                - the L1 grid for this file will have '+str(n_timesteps)+' timesteps')

    ##############################################################
    # these lines will chop down the L1 data to hurry up the interpolation

    ll_dist = ((XC - L3_XC[0, 0]) ** 2 + (YC - L3_YC[0, 0]) ** 2) ** 0.5
    ll_row, ll_col = np.where(ll_dist == np.min(ll_dist))

    ur_dist = ((XC - L3_XC[-1, -1]) ** 2 + (YC - L3_YC[-1, -1]) ** 2) ** 0.5
    ur_row, ur_col = np.where(ur_dist == np.min(ur_dist))

    lr_dist = ((XC - L3_XC[0, -1]) ** 2 + (YC - L3_YC[0, -1]) ** 2) ** 0.5
    lr_row, lr_col = np.where(lr_dist == np.min(lr_dist))

    ul_dist = ((XC - L3_XC[-1, 0]) ** 2 + (YC - L3_YC[-1, 0]) ** 2) ** 0.5
    ul_row, ul_col = np.where(ul_dist == np.min(ul_dist))

    min_row = np.min([ll_row[0], ul_row[0], ur_row[0], lr_row[0]])
    max_row = np.max([ll_row[0], ul_row[0], ur_row[0], lr_row[0]])
    min_col = np.min([ll_col[0], ul_col[0], ur_col[0], lr_col[0]])
    max_col = np.max([ll_col[0], ul_col[0], ur_col[0], lr_col[0]])

    XC = XC[min_row:max_row, min_col:max_col]
    YC = YC[min_row:max_row, min_col:max_col]
    Mask = Mask[min_row:max_row, min_col:max_col]
    exf_field = exf_field[:, min_row:max_row, min_col:max_col]

    ##############################################################

    points = np.column_stack([XC.ravel(),YC.ravel()])
    values = np.zeros((n_timesteps,np.size(XC)))
    for timestep in range(n_timesteps):
        values[timestep,:] = exf_field[timestep,:,:].ravel()
    mask = Mask.ravel()

    return(points,values,mask)

def downscale_L1_exf_field_to_L3(df, var_name,
                                 L1_points, L1_values,
                                 L1_wet_grid_points, L1_wet_grid_on_L3_points,
                                 L3_XC, L3_YC, L3_wet_grid, print_level):

    nTimestepsOut = np.shape(L1_values)[0]

    L3_exf_var = np.zeros((nTimestepsOut, np.shape(L3_YC)[0], np.shape(L3_YC)[1]))

    if print_level >= 5:
        print('                -  Variable shapes entering downscale routine:')
        print('                    -  L1_point: '+str(np.shape(L1_points)))
        print('                    -  L1_values: ' + str(np.shape(L1_values)))
        print('                    -  L1_wet_grid: ' + str(np.shape(L1_wet_grid_points)))
        # print('          -  L1_wet_grid_on_L3_subset: ' + str(np.shape(L1_wet_grid_on_L3_3D_points)))
        print('                    -  L3_XC: ' + str(np.shape(L3_XC)))
        print('                    -  L3_YC: ' + str(np.shape(L3_YC)))
        print('                    -  L3_wet_grid: ' + str(np.shape(L3_wet_grid)))

    for timestep in range(nTimestepsOut):
        # if timestep%50==0:
        if print_level >= 5:
            print('                    - Working on timestep '+str(timestep)+' of '+str(nTimestepsOut))

        downscaled_field = df.downscale_2D_points_with_zeros(L1_points, L1_values[timestep, :],
                                                  L1_wet_grid_points,
                                                  L3_XC, L3_YC, L3_wet_grid[0,:,:])
        L3_exf_var[timestep, :, :] = downscaled_field

    # plt.imshow(L3_exf_var[0,:,:],origin='lower')
    # plt.show()

    return(L3_exf_var)


########################################################################################################################


def create_exf_field(L1f, config_dir, L1_model_name, L3_model_name,
                     sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                     faces, face_size_dict, ordered_tiles_faces_dict,
                     var_name, years, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import downscale_functions as df

    # this is the dir where the exf output will be stored
    if 'exf' not in os.listdir(os.path.join(config_dir,'L3',L3_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L3',L3_model_name,'input','exf'))
    output_dir = os.path.join(config_dir,'L3',L3_model_name,'input','exf')

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

    L3_XC, L3_YC, _, _, L3_wet_grid = read_grid_geometry_from_nc(config_dir,L3_model_name,var_name)

    for year in years:
        if print_level >= 2:
            print('        - Downscaling the timesteps for year ' + str(year))

        output_files, start_indices, end_indices = create_annual_list_of_files(year, var_name)

        all_files_completed = True
        for output_file_name in output_files:
            if output_file_name not in os.listdir(os.path.join(output_dir,var_name)):
                all_files_completed = False

        if not all_files_completed:

            if print_level >= 3:
                print('            - Reading in the L1 file')

            L1_exf_dir = os.path.join(config_dir,'L1',L1_model_name,'input','exf')
            L1_points, L1_values, L1_wet_grid_points = read_exf_variable_from_L1(L1f, L1_exf_dir, var_name, year,
                                                                                 sNx, sNy, faces, face_size_dict,
                                                                                 L1_XC_faces, L1_YC_faces, L1_AngleCS_faces, L1_AngleSN_faces,
                                                                                 L1_Mask_faces, L3_XC, L3_YC, print_level)


            for ff in range(90):#len(output_files)):
                output_file_name = output_files[ff]
                print(output_file_name)
                if output_file_name not in os.listdir(os.path.join(output_dir, var_name)):
                    start_index = start_indices[ff]
                    end_index = end_indices[ff]

                    if print_level >= 4:
                        print('                - Creating file '+output_file_name+' (indices '+str(start_index)+' to '+str(end_index)+')')

                    L1_values_subset = L1_values[start_index:end_index,:]

                    L1_wet_grid_on_L3_points = np.copy(L1_wet_grid_points)

                    if print_level >= 5:
                        print('                    - Downscaling the output to the new boundary')
                    L3_exf_var = downscale_L1_exf_field_to_L3(df, var_name,
                                                              L1_points, L1_values_subset,
                                                              L1_wet_grid_points, L1_wet_grid_on_L3_points,
                                                              L3_XC, L3_YC, L3_wet_grid, print_level)

                    output_file = os.path.join(output_dir, var_name, output_file_name)
                    L3_exf_var.ravel(order='C').astype('>f4').tofile(output_file)
                else:
                    if print_level >= 4:
                        print('                - Skipping file '+output_file_name+' (already created)')

        else:
            if print_level >= 3:
                print('            - Skipping ' + str(year) + ' - all daily files completed for this year')

def create_exf_fields_via_interpolation(L1f, config_dir, L3_model_name, L1_model_name, proc_id,
                                        sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                                        faces, face_size_dict, ordered_tiles_faces_dict,
                                        start_year, final_year, print_level):

    var_name_list = ['ATEMP','AQH','LWDOWN','SWDOWN','UWIND','VWIND','PRECIP','RUNOFF']
    var_name = var_name_list[proc_id % len(var_name_list)]

    if print_level>=1:
        print('    - Creating the exf fields for ' + var_name + ' to cover years ' +
          str(start_year) + ' to ' + str(final_year))
    years = np.arange(start_year,final_year+1).tolist()

    if print_level >= 1:
        print('    - Running the Downscale routine:')
    create_exf_field(L1f, config_dir, L1_model_name, L3_model_name,
                     sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                     faces, face_size_dict, ordered_tiles_faces_dict,
                     var_name, years, print_level)

