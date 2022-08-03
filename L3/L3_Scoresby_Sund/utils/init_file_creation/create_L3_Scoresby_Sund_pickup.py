
import os
import argparse
import sys


def create_pickup_file(config_dir):

    sys.path.insert(1, os.path.join(config_dir,'L3','utils','init_file_creation'))

    L3_model_name = 'L3_Scoresby_Sund'
    source = 'L2'


    if source == 'ASTE':
        import create_L3_ASTE_pickup as cap

        ordered_aste_tiles = [[14,5]]
        ordered_aste_tile_rotations = [[2,3]] # rotations are counter-clockwise

        cap.create_L3_ASTE_pickup_file(config_dir,model_name,
                                       ordered_aste_tiles,ordered_aste_tile_rotations)

    if source == 'L2':
        import create_L3_pickup_from_L2 as c32

        L2_model_name = 'L2_CE_Greenland'
        c32.create_L3_pickup_file(config_dir, L2_model_name, L3_model_name)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

