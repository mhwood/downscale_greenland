
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys


def create_seaice_pickup_file(config_dir, parent_model, parent_model_pickup_iteration, L2_model_name,
                       sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                       faces, ordered_tiles_faces_dict, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L2', 'utils','init_file_creation'))
    import create_L2_seaice_pickup_from_L1 as csp

    print('Creating the pickup for the '+L2_model_name+' model')

    # pass to general function to generate mitgrid
    csp.create_seaice_pickup_from_L1(config_dir, parent_model, parent_model_pickup_iteration, L2_model_name,
                                     sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                                     faces, ordered_tiles_faces_dict, print_level)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L05, L05, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)


    args = parser.parse_args()
    config_dir = args.config_dir

    create_seaice_pickup_file(config_dir)
   

