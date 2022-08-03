
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
import argparse
from pyproj import Transformer
from MITgcmutils import mds
import sys

def read_pickup_file_to_compact(pickup_file_path):

    Nr = 90
    print('      Reading from '+pickup_file_path)
    global_data, _, global_metadata = mds.rdmds(pickup_file_path, returnmeta=True)

    has_Nr = {'uvel': True, 'vvel': True, 'theta': True,
              'salt': True, 'gunm1': True, 'gvnm1': True,
              'gunm2': True, 'gvnm2': True, 'etan': False,
              'detahdt': False, 'etah': False}

    var_names = []
    row_bounds = []
    var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        if has_Nr[var_name.strip().lower()]:
            end_row = start_row + Nr
        else:
            end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]
        var_grids.append(var_grid)
        row_bounds.append([start_row,end_row])
        start_row=end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,var_grids,global_metadata)


def read_pickup_file_to_faces(config_dir,config_name,field_name,faces, face_dimensions):

    sNx = 180
    sNy = 180

    pickup_file = 'pickup.0000000001'
    pickup_file_path = os.path.join(config_dir,'L1',config_name, 'input', pickup_file)
    var_names, row_bounds, compact_var_grids, global_metadata = read_pickup_file_to_compact(pickup_file_path)
    var_name_index = var_names.index(field_name)

    plt.imshow(compact_var_grids[2][0,:,:])
    plt.show()

    print('    - Reading '+var_names[var_name_index]+' to faces')
    compact_var_grid = compact_var_grids[var_name_index]
    grid_faces = {}
    lower_index = 0
    for face in faces:
        upper_index = lower_index + int(face_dimensions[face][0]*face_dimensions[face][1]/sNx)
        face_grid = compact_var_grid[:, lower_index:upper_index, :]
        dimensions = face_dimensions[face]
        face_grid = np.reshape(face_grid, (np.shape(face_grid)[0], dimensions[0], dimensions[1]))
        grid_faces[face] = face_grid
        lower_index = upper_index

    return (var_names, grid_faces, global_metadata)


def plot_pickup_field(config_dir, model_name, field_name, faces, size_dict):

    llc = 1080
    sNx = 180

    var_names, field_faces, global_metadata = read_pickup_file_to_faces(config_dir,model_name, field_name, faces, size_dict)

    fig = plt.figure(figsize=(8,8))
    plt.style.use('dark_background')

    gs = GridSpec(4,4,left = 0.05, right=0.95)

    for face in range(1,7):
        if face in faces:

            if face==1:
                ax1 = fig.add_subplot(gs[2:,0])
                ax1.imshow(field_faces[1][0,:,:],origin='lower')#,

            if face==3:
                ax3 = fig.add_subplot(gs[1, 1])
                ax3.imshow(field_faces[3][0,:,:],origin='lower')#,

            if face == 4:
                ax4 = fig.add_subplot(gs[1, 2])
                ax4.imshow(field_faces[4][0,:,:],origin='lower')#,

            if face == 5:
                ax5 = fig.add_subplot(gs[0, 2:])
                ax5.imshow(field_faces[5][0,:,:],origin='lower')#,


    output_path = os.path.join(config_dir,'L1',model_name,'plots',model_name+'_pickup_'+field_name+'.png')
    plt.savefig(output_path)
    plt.close(fig)




