import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import simplegrid as sg
from scipy.interpolate import griddata
from ecco_v4_py.llc_array_conversion import llc_faces_to_compact, llc_compact_to_faces
import argparse
import ast



def read_parent_XC_YC(config_dir,model_name,faces,size_dict):

    XC_faces = {}
    YC_faces = {}
    Bathy_faces = {}

    for face in faces:
        grid_file = os.path.join(config_dir, 'L05', model_name, 'input', 'tile' + '{:03d}'.format(face) + '.mitgrid')
        grid = np.fromfile(grid_file, '>f8')
        grid = grid.reshape(16, size_dict[face][0] + 1, size_dict[face][1] + 1)
        XC_faces[face] = grid[0, :-1, :-1]
        YC_faces[face] = grid[1, :-1, :-1]


    bathy_file = os.path.join(config_dir, 'L05', model_name, 'input', 'bathymetry.bin')
    bathy_grid = np.fromfile(bathy_file,'>f4')

    points_counted = 0
    for face in faces:
        n_points = size_dict[face][0] * size_dict[face][1]
        bathy_face = bathy_grid[points_counted:points_counted + n_points]
        if face==3:
            bathy_face = bathy_face.reshape((size_dict[face][0], size_dict[face][1]))
            Bathy_faces[face] = bathy_face
        else:
            bathy_face = bathy_face.reshape((size_dict[face][0], size_dict[face][1]))
            Bathy_faces[face] = bathy_face
        points_counted += n_points

    # plt.subplot(2,3,1)
    # C = plt.imshow(XC_faces[1],origin='lower')
    # plt.colorbar(C)
    # plt.title('XC')
    #
    # plt.subplot(2, 3, 2)
    # C = plt.imshow(YC_faces[1], origin='lower')
    # plt.colorbar(C)
    # plt.title('YC')
    #
    # plt.subplot(2, 3, 3)
    # C = plt.imshow(Bathy_faces[1], origin='lower')
    # plt.colorbar(C)
    # plt.title('Bathy')
    #
    # plt.subplot(2, 3, 4)
    # C = plt.imshow(XC_faces[3], origin='lower')
    # plt.colorbar(C)
    #
    # plt.subplot(2, 3, 5)
    # C = plt.imshow(YC_faces[3], origin='lower')
    # plt.colorbar(C)
    #
    # plt.subplot(2, 3, 6)
    # C = plt.imshow(Bathy_faces[3], origin='lower')
    # plt.colorbar(C)
    # plt.title('Bathy')
    #
    # plt.show()

    return(XC_faces, YC_faces, Bathy_faces)


def read_child_XC_YC(config_dir, model_name, n_cols, n_rows):
    mitgrid_file = os.path.join(config_dir,'mitgrids', model_name+'.mitgrid')
    grid_dict = sg.gridio.read_mitgridfile(mitgrid_file, n_cols, n_rows)
    XC = grid_dict['XC'].T
    YC = grid_dict['YC'].T
    return(XC,YC)


def create_mask_faces(faces, resolution, XC_faces, YC_faces, XC_boundary, YC_boundary):
    mask_faces = {}
    mask_indices = []
    counter = 1
    for face in faces:
        mask_face = np.zeros_like(XC_faces[face])
        for i in range(len(XC_boundary)):
            dist = ((XC_faces[face]-XC_boundary[i])**2 + (YC_faces[face]-YC_boundary[i])**2)**0.5
            rows, cols = np.where(dist<resolution*np.sqrt(2))
            for i in range(len(rows)):
                if mask_face[rows[i],cols[i]]==0:
                    mask_face[rows[i],cols[i]] = counter
                    mask_indices.append([face,rows[i],cols[i]])
                    counter +=1
        mask_faces[face] = mask_face
    mask_indices = np.array(mask_indices)
    return(mask_faces, mask_indices, counter)


def output_mask_dictionary_to_nc(output_dir,output_file_name,all_mask_dicts,mask_names_list):
    if output_file_name in os.listdir(output_dir):
        os.remove(os.path.join(output_dir,output_file_name))

    ds = nc4.Dataset(os.path.join(output_dir,output_file_name),'w')

    for m in range(len(mask_names_list)):
        grp = ds.createGroup(mask_names_list[m])
        grp.createDimension('n_points', np.shape(all_mask_dicts[m])[0])
        var = grp.createVariable('source_faces', 'i4', ('n_points',))
        var[:] = all_mask_dicts[m][:, 0].astype(int)
        var = grp.createVariable('source_rows', 'i4', ('n_points',))
        var[:] = all_mask_dicts[m][:, 1].astype(int)
        var = grp.createVariable('source_cols', 'i4', ('n_points',))
        var[:] = all_mask_dicts[m][:, 2].astype(int)

    ds.close()


########################################################################################################################

def create_dv_masks(config_dir):

    print_status = True

    parent_model_name = 'L05_CE_Greenland'
    child_model_name = 'L2_CE_Greenland'

    if print_status:
        print('Creating the diagnostics_vec masks to use in the L05 domain')
    #
    # if 'input' not in os.listdir('..'):
    #     os.mkdir(os.path.join('..','input'))
    if 'dv' not in os.listdir(os.path.join(config_dir,'L05',parent_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L05',parent_model_name,'input','dv'))

    ###############################################################################################
    # Read in the grid geometries from dictionaries

    faces_dim_file = os.path.join(config_dir, 'L05', parent_model_name, 'namelist', 'face_dimensions.txt')
    f = open(faces_dim_file)
    dict_str = f.read()
    f.close()
    L05_size_dict = ast.literal_eval(dict_str)
    L05_faces = list(L05_size_dict.keys())

    f = open(os.path.join(config_dir, 'domain_sizes.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)
    L2_size = size_dict['L2_CE_Greenland']
    L2_n_rows = L2_size[0]
    L2_n_cols = L2_size[1]
    dtype = '>f4'

    ###############################################################################################
    # Read in the grids

    if print_status:
        print('    Reading in the L05 domain files')

    # read the main ECCO grid XC and YC as compacts

    XC_faces, YC_faces, Bathy_faces = read_parent_XC_YC(config_dir,parent_model_name,L05_faces,L05_size_dict)

    # read the extended subset of the XC and YC grids
    child_XC, child_YC = read_child_XC_YC(config_dir, child_model_name, L2_n_cols, L2_n_rows)
    if print_status:
        print('    The subdomain has shape '+str(np.shape(child_XC)))

    ###############################################################################################
    # Create the masks

    resolution = 1/6
    sNx = 90

    mask_names_list = []
    all_mask_dicts = []

    for boundary in ['north','south','east','surface']:
        if print_status:
            print('    Creating the '+boundary+' mask')
        if boundary == 'south':
            XC_boundary = child_XC[0,:]
            YC_boundary = child_YC[0,:]

        if boundary == 'north':
            XC_boundary = child_XC[-1,:]
            YC_boundary = child_YC[-1,:]

        if boundary == 'east':
            XC_boundary = child_XC[:,-1]
            YC_boundary = child_YC[:,-1]

        if boundary == 'west':
            XC_boundary = child_XC[:,0]
            YC_boundary = child_YC[:,0]

        if boundary == 'surface':
            XC_boundary = child_XC.ravel()
            YC_boundary = child_YC.ravel()

        mask_faces, mask_indices, counter = create_mask_faces(L05_faces, resolution, XC_faces, YC_faces, XC_boundary, YC_boundary)
        print('        - This mask has '+str(counter)+' points')
        all_mask_dicts.append(mask_indices)
        mask_names_list.append(boundary)

        plt.subplot(1, 2, 1)
        C = plt.imshow(mask_faces[1],origin='lower')
        plt.colorbar(C)
        plt.subplot(1, 2, 2)
        C = plt.imshow(mask_faces[3], origin='lower')
        plt.colorbar(C)
        plt.show()

        for f in range(len(L05_faces)):
            face = L05_faces[f]
            grid = mask_faces[face]
            n_rows = int((np.shape(grid)[0] * np.shape(grid)[1]) / sNx)
            grid = np.reshape(grid, (n_rows, sNx))
            if face == 1:
                compact_stack = grid
            else:
                compact_stack = np.concatenate([compact_stack, grid], axis=0)

        if boundary=='surface':
            output_file = os.path.join(config_dir, 'L05', parent_model_name, 'input', 'dv', 'L2_surface_mask.bin')
        else:
            output_file = os.path.join(config_dir, 'L05', parent_model_name, 'input', 'dv', 'L2_' + boundary + '_BC_mask.bin')
        compact_stack.ravel(order='C').astype('>f4').tofile(output_file)

    output_dir = os.path.join(config_dir,'L05',parent_model_name,'namelist')
    output_file_name = 'L05_dv_mask_reference_dict.nc'
    output_mask_dictionary_to_nc(output_dir,output_file_name,all_mask_dicts,mask_names_list)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L05, L05, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_dv_masks(config_dir)
