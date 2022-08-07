
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
    ds.close()
    return(XC,YC)

def create_mask(resolution, parent_XC, parent_YC, XC_boundary, YC_boundary):
    mask_indices = []
    counter = 1
    mask = np.zeros_like(parent_XC)
    for i in range(len(XC_boundary)):
        dist = ((parent_XC-XC_boundary[i])**2 + (parent_YC-YC_boundary[i])**2)**0.5
        rows, cols = np.where(dist<resolution*np.sqrt(2))
        for i in range(len(rows)):
            if mask[rows[i],cols[i]]==0:
                mask[rows[i],cols[i]] = counter
                mask_indices.append([rows[i],cols[i]])
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

def create_L2_diagnostic_vec_masks(config_dir, L2_config_name, L3_config_name):

    print_status = True

    # this is the dir where the dv masks will be stored
    if 'dv' not in os.listdir(os.path.join(config_dir,'L2',L2_config_name,'input')):
        os.mkdir(os.path.join(config_dir,'L2',L2_config_name,'input','dv'))
    output_dir = os.path.join(config_dir,'L2',L2_config_name,'input','dv')

    L2_XC, L2_YC = read_grid_geometry_from_nc(config_dir, L2_config_name)

    L3_XC, L3_YC = read_grid_geometry_from_nc(config_dir, L3_config_name)

    print('    - Generating the masks')

    ###############################################################################################
    # Create the masks

    resolution = 1 / 12
    sNx = 90

    mask_names_list = []
    all_mask_dicts = []

    for boundary in ['north', 'south', 'east', 'surface']:
        if print_status:
            print('    Creating the ' + boundary + ' mask')

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

        mask, mask_indices, counter = create_mask(resolution, L2_XC, L2_YC,
                                                        XC_boundary,YC_boundary)

        print('        - This mask has ' + str(counter) + ' points')
        all_mask_dicts.append(mask_indices)
        mask_names_list.append(boundary)

        C = plt.imshow(mask, origin='lower')
        plt.colorbar(C)
        plt.show()

        mask = np.array(mask)

        if boundary == 'surface':
            output_file = os.path.join(config_dir, 'L2', L2_config_name, 'input', 'dv', 'L3_surface_mask.bin')
        else:
            output_file = os.path.join(config_dir, 'L2', L2_config_name, 'input', 'dv',
                                       'L3_' + boundary + '_BC_mask.bin')
        mask.ravel(order='C').astype('>f4').tofile(output_file)

    output_dir = os.path.join(config_dir, 'L2', L2_config_name, 'namelist')
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