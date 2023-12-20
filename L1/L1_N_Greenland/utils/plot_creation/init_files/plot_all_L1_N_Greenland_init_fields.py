
import os
import netCDF4 as nc4
import numpy as np
import argparse
import sys

def plot_all_L1_W_Greenland_fields(config_dir, print_level=4):

    L1_model_name = 'L1_N_Greenland'

    sys.path.insert(1, os.path.join(config_dir, 'L1', L1_model_name, 'utils', 'plot_creation','init_files'))
    print(os.path.join(config_dir, 'L1', L1_model_name, 'utils', 'plot_creation','init_files'))

    # import plot_L1_N_Greenland_pickup_fields as ppf
    # ppf.plot_pickup(config_dir, print_level)

    # import plot_L1_N_Greenland_seaice_pickup_fields as ppf
    # ppf.plot_seaice_pickup(config_dir, print_level)

    # import plot_L1_N_Greenland_exfs as pex
    # pex.plot_exfs(config_dir, print_level)

    import plot_L1_N_Greenland_BCs as pbc
    pbc.plot_BCs(config_dir, print_level)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_all_L1_W_Greenland_fields(config_dir, print_level=4)
   

