
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys


def read_bathy_from_binary(config_dir,model_name,faces,size_dict):

    bathy_file = os.path.join(config_dir,'L05',model_name,'input','bathymetry.bin')
    bathy_grid = np.fromfile(bathy_file,'>f4')

    points_counted = 0
    bathy_faces = {}
    for face in faces:
        n_points = size_dict[face][0]*size_dict[face][1]
        bathy_face = bathy_grid[points_counted:points_counted+n_points]
        bathy_face = bathy_face.reshape((size_dict[face][0], size_dict[face][1]))
        bathy_faces[face] = bathy_face
        points_counted += n_points

    return(bathy_faces)

def create_L05_CE_Greenland_masked_bathy(config_dir):

    model_name = 'L05_CE_Greenland'

    llc = 540

    sNx = 90
    sNy = 90

    faces = [1,3]
    size_dict = {1:(90,270),3:(270,90)}

    bathy_faces = read_bathy_from_binary(config_dir, model_name, faces, size_dict)

    # plt.imshow(bathy_faces[1]<0,origin='lower')
    # plt.show()

    # zero out sermlik fjord
    bathy_faces[1][38:60,:3] = 0
    bathy_faces[1][44:60,:5] = 0
    bathy_faces[1][47:60, :11] = 0

    # plt.imshow(bathy_faces[1]<0,origin='lower')
    # plt.show()

    compact_tile_size = 90

    for face in faces:
        grid = bathy_faces[face]
        n_rows = int((np.shape(grid)[0] * np.shape(grid)[1]) / compact_tile_size)
        grid = np.reshape(grid, (n_rows, compact_tile_size))
        if face == 1:
            compact_stack = grid
        else:
            compact_stack = np.concatenate([compact_stack, grid], axis=0)

    output_file = os.path.join(config_dir, 'L05', model_name, 'input', 'bathymetry_masked.bin')
    compact_stack.ravel(order='C').astype('>f4').tofile(output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L05, L05, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_L05_CE_Greenland_masked_bathy(config_dir)
   

