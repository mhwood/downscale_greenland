
import os
import argparse
import sys


def create_bathy_file(config_dir):

    sys.path.insert(1, os.path.join(config_dir,'utils','init_file_creation'))
    import create_bathymetry as cb

    model_name = 'L3_Scoresby_Sund'
    level_name = 'L3'

    cb.create_bathymetry_file(config_dir, level_name, model_name)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_bathy_file(config_dir)
   

