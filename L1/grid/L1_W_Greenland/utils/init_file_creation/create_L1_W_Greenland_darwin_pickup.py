
import os
import argparse
import sys


def create_darwin_pickup_file(config_dir, model_name,
                       ecco_dir, llc,ordered_ecco_tiles, ordered_ecco_tile_rotations,
                       parent_model_pickup_iteration, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L1_grid','utils','init_file_creation'))

    source = 'ECCO'
    import create_L1_ECCO_darwin_pickup as cep

    ###################################################################
    # these are the ecco functions

    cep.create_L1_ECCO_darwin_pickup_file(config_dir, model_name,
                                   ecco_dir, llc,ordered_ecco_tiles, ordered_ecco_tile_rotations,
                                   parent_model_pickup_iteration, print_level)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

