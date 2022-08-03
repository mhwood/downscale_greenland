
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys


def create_pickup_file(config_dir):

    sys.path.insert(1, os.path.join(config_dir,'L2', 'utils','init_file_creation'))
    import create_L2_seaice_pickup_from_L05 as csp

    L05_model_name = 'L05_CE_Greenland'
    L05_iteration = 1008
    L2_model_name = 'L2_CE_Greenland'

    sNx = 90
    sNy = 90
    ordered_nonblank_tiles = [[1, 2, 3], [6, 5, 4]]
    ordered_nonblank_rotations = [[0, 0, 0], [3, 3, 3]]

    faces = [1,3]
    ordered_tiles_faces_dict = {1:[[1,2,3]],
                                3:[[4],[5],[6]]}

    print('Creating the pickup for the '+L2_model_name+' model')

    # pass to general function to generate mitgrid
    csp.create_seaice_pickup_from_L05(config_dir, L05_model_name, L05_iteration, L2_model_name,
                      sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                      faces, ordered_tiles_faces_dict)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L05, L05, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)


    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

