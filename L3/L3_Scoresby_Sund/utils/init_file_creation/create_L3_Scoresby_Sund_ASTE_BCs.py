
import os
import argparse
import sys


def create_boundary_condition_files(config_dir):

    sys.path.insert(1, os.path.join(config_dir,'L3','utils','init_file_creation'))
    import create_L3_ASTE_BCs as cbc

    model_name = 'L3_Scoresby_Sund'
    ordered_aste_tiles = [[14,5]]
    ordered_aste_tile_rotations = [[2,3]] # rotations are counter-clockwise

    var_names = ['THETA']
    boundaries = ['north']

    for var_name in var_names:
        for boundary in boundaries:
            cbc.create_L3_ASTE_boundary_condition(config_dir,model_name,
                                                  var_name,boundary,
                                                  ordered_aste_tiles,ordered_aste_tile_rotations)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_boundary_condition_files(config_dir)
   

