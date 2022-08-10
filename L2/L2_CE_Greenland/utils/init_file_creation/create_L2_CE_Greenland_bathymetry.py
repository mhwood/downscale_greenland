
import os
import argparse
import sys


def create_bathy_file(config_dir, model_name, central_wet_row, central_wet_col, hFacMinDr, hFacMin, delR, print_level):

    sys.path.insert(1, os.path.join(config_dir,'utils','init_file_creation'))
    import create_bathymetry as cb

    level_name = 'L2'

    cb.create_bathymetry_file(config_dir, level_name, model_name,
                              central_wet_row, central_wet_col, hFacMinDr, hFacMin, delR, print_level)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_bathy_file(config_dir)
   

