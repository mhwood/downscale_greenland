
import os
import simplegrid as sg
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import argparse
import ast

def read_grid_geometry_from_nc(config_dir,model_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:]
    ds.close()
    return(XC,YC,Depth)

def great_circle_distance(lon_ref, lat_ref, Lon, Lat):
    earth_radius = 6371000
    lon_ref_radians = np.radians(lon_ref)
    lat_ref_radians = np.radians(lat_ref)
    lons_radians = np.radians(Lon)
    lats_radians = np.radians(Lat)
    lat_diff = lats_radians - lat_ref_radians
    lon_diff = lons_radians - lon_ref_radians
    d = np.sin(lat_diff * 0.5) ** 2 + np.cos(lat_ref_radians) * np.cos(lats_radians) * np.sin(lon_diff * 0.5) ** 2
    h = 2 * earth_radius * np.arcsin(np.sqrt(d))
    return(h)

def create_mask(resolution, parent_XC, parent_YC, parent_Depth, XC_boundary, YC_boundary):
    mask_indices = []
    counter = 1
    mask = np.zeros_like(parent_XC)
    for i in range(len(XC_boundary)):
        # dist = ((parent_XC-XC_boundary[i])**2 + (parent_YC-YC_boundary[i])**2)**0.5
        dist = great_circle_distance(XC_boundary[i], YC_boundary[i], parent_XC, parent_YC)
        rows, cols = np.where(dist<resolution*np.sqrt(2))
        for i in range(len(rows)):
            if mask[rows[i],cols[i]]==0:
                if parent_Depth[rows[i],cols[i]]>0:
                    mask[rows[i],cols[i]] = counter
                    mask_indices.append([rows[i],cols[i]])
                    counter +=1
    mask_indices = np.array(mask_indices)
    return(mask, mask_indices, counter)

def create_surface_mask(resolution, parent_XC, parent_YC, parent_Depth, XC_boundary, YC_boundary):
    mask_indices = []
    counter = 1
    mask = np.zeros_like(parent_XC)

    ll_dist = great_circle_distance(np.min(XC_boundary),np.min(YC_boundary),parent_XC,parent_YC)
    ll_dist = np.reshape(ll_dist,np.shape(parent_XC))
    ll_row, ll_col = np.where(ll_dist==np.min(ll_dist))
    ll_row = ll_row[0]
    ll_col = ll_col[0]

    lr_dist = great_circle_distance(np.max(XC_boundary), np.min(YC_boundary), parent_XC, parent_YC)
    lr_dist = np.reshape(lr_dist, np.shape(parent_XC))
    lr_row, lr_col = np.where(lr_dist == np.min(lr_dist))
    lr_row = lr_row[0]
    lr_col = lr_col[0]

    ur_dist = great_circle_distance(np.max(XC_boundary), np.max(YC_boundary), parent_XC, parent_YC)
    ur_dist = np.reshape(ur_dist, np.shape(parent_XC))
    ur_row, ur_col = np.where(ur_dist == np.min(ur_dist))
    ur_row = ur_row[0]
    ur_col = ur_col[0]

    ul_dist = great_circle_distance(np.min(XC_boundary), np.max(YC_boundary), parent_XC, parent_YC)
    ul_dist = np.reshape(ul_dist, np.shape(parent_XC))
    ul_row, ul_col = np.where(ul_dist == np.min(ul_dist))
    ul_row = ul_row[0]
    ul_col = ul_col[0]

    min_row = np.min([ll_row, lr_row, ur_row, ul_row])
    min_col = np.min([ll_col, lr_col, ur_col, ul_col])
    max_row = np.max([ll_row, lr_row, ur_row, ul_row])
    max_col = np.max([ll_col, lr_col, ur_col, ul_col])

    if min_row>0:
        min_row-=1
    if min_col>0:
        min_col-=1
    if max_row<np.shape(parent_XC)[0]-1:
        max_row+=1
    if max_col<np.shape(parent_YC)[1]-1:
        max_col+=1

    for i in range(min_row,max_row+1):
        for j in range(min_col,max_col+1):
            if parent_Depth[i,j]>0:
                mask[i,j] = counter
                mask_indices.append([i,j])
                counter +=1
    mask_indices = np.array(mask_indices)
    return(mask, mask_indices, counter)

def output_mask_dictionary_to_nc(output_dir,output_file_name,all_mask_dicts,mask_names_list):
    if output_file_name in os.listdir(output_dir):
        os.remove(os.path.join(output_dir,output_file_name))

    ds = nc4.Dataset(os.path.join(output_dir,output_file_name),'w')

    for m in range(len(mask_names_list)):
        grp = ds.createGroup(mask_names_list[m])
        grp.createDimension('n_points', np.shape(all_mask_dicts[m])[0])
        var = grp.createVariable('source_rows', 'i4', ('n_points',))
        var[:] = all_mask_dicts[m][:, 0].astype(int)
        var = grp.createVariable('source_cols', 'i4', ('n_points',))
        var[:] = all_mask_dicts[m][:, 1].astype(int)

    ds.close()

########################################################################################################################

def create_L2_diagnostic_vec_masks(config_dir, L2_model_name, L3_model_name, print_level):

    print_status = True

    # this is the dir where the dv masks will be stored
    if 'dv' not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,'input','dv'))
    output_dir = os.path.join(config_dir,'L2',L2_model_name,'input','dv')

    L2_XC, L2_YC, L2_Depth = read_grid_geometry_from_nc(config_dir, L2_model_name)

    L3_XC, L3_YC, _ = read_grid_geometry_from_nc(config_dir, L3_model_name)


    ###############################################################################################
    # Create the masks

    resolution = 5000

    mask_names_list = []
    all_mask_dicts = []

    for boundary in ['surface','north', 'south', 'east']:
        if print_level >= 1:
            print('    - Creating the ' + boundary + ' mask')

        if boundary == 'south':
            XC_boundary = L3_XC[0, :]
            YC_boundary = L3_YC[0, :]

        if boundary == 'north':
            XC_boundary = L3_XC[-1, :]
            YC_boundary = L3_YC[-1, :]

        if boundary == 'east':
            XC_boundary = L3_XC[:, -1]
            YC_boundary = L3_YC[:, -1]

        if boundary == 'west':
            XC_boundary = L3_XC[:, 0]
            YC_boundary = L3_YC[:, 0]

        if boundary == 'surface':
            XC_boundary = L3_XC.ravel()
            YC_boundary = L3_YC.ravel()

        if boundary!='surface':
            mask, mask_indices, counter = create_mask(resolution, L2_XC, L2_YC, L2_Depth,
                                                      XC_boundary,YC_boundary)
        else:
            mask, mask_indices, counter = create_surface_mask(resolution, L2_XC, L2_YC, L2_Depth,
                                                              XC_boundary, YC_boundary)

        if print_level >= 2:
            print('        - This mask has ' + str(counter) + ' points')
        all_mask_dicts.append(mask_indices)
        mask_names_list.append(boundary)

        mask = np.array(mask)

        if boundary == 'surface':
            output_file = os.path.join(config_dir, 'L2', L2_model_name, 'input', 'dv', 'L3_surface_mask.bin')
        else:
            output_file = os.path.join(config_dir, 'L2', L2_model_name, 'input', 'dv',
                                       'L3_' + boundary + '_BC_mask.bin')
        mask.ravel(order='C').astype('>f4').tofile(output_file)

    output_dir = os.path.join(config_dir, 'L2', L2_model_name, 'input')
    output_file_name = 'L2_dv_mask_reference_dict.nc'
    output_mask_dictionary_to_nc(output_dir, output_file_name, all_mask_dicts, mask_names_list)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_diagnostic_vec_masks(config_dir)