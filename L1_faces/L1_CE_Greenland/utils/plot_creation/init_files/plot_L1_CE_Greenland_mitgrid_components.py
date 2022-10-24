
import os
import argparse
import sys
import numpy as np

def plot_L1_CE_Greenland_mitgrid_components(config_dir, L1_model_name, faces, face_size_dict):

    L1_model_name = 'L1_CE_Greenland'

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils', 'plot_creation','init_files'))

    import plot_L1_mitgrid_components as pmc
    pmc.create_mitgrid_plot(config_dir, L1_model_name, faces, face_size_dict)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_L1_CE_Greenland_init_fields(config_dir)
   

