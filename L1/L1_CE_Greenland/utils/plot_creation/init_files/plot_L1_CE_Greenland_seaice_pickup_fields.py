
import os
import argparse
import sys
import numpy as np

def plot_L1_CE_Greenland_seaice_pickup(config_dir, L1_model_name, pickup_iteration,
                                sNx, sNy, faces, face_size_dict):

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils', 'plot_creation','init_files'))

    import plot_L1_seaice_pickup_fields as ppf
    ppf.create_seaice_pickup_plot(config_dir, L1_model_name, pickup_iteration,
                           sNx, sNy, faces, face_size_dict)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_L1_CE_Greenland_init_fields(config_dir)
   

