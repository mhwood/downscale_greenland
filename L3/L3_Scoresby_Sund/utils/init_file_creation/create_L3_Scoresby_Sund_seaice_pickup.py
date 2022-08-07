
import os
import argparse
import sys


def create_seaice_pickup_file(config_dir, L3_model_name, parent_model, parent_model_pickup_iteration, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L3','utils','init_file_creation'))

    if parent_model == 'ASTE':
        import create_L3_ASTE_seaice_pickup as cap

        ordered_aste_tiles = [[14,5]]
        ordered_aste_tile_rotations = [[2,3]] # rotations are counter-clockwise

        cap.create_L3_ASTE_seaice_pickup_file(config_dir,model_name,
                                       ordered_aste_tiles,ordered_aste_tile_rotations)

    if parent_model == 'L2_CE_Greenland':
        import create_L3_seaice_pickup_from_L2 as c32
        c32.create_L3_seaice_pickup_file(config_dir, L3_model_name,
                                         parent_model, parent_model_pickup_iteration, print_level)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_seaice_pickup_file(config_dir)
   

