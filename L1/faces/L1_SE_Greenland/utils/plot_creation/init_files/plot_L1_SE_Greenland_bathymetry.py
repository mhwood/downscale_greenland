
import os
import argparse
import sys
import numpy as np

def plot_L1_SE_Greenland_bathymetry(config_dir, L1_model_name, sNx, sNy, faces, face_size_dict):

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils', 'plot_creation','init_files'))

    sys.path.insert(1, os.path.join(config_dir, 'L1', L1_model_name, 'utils'))
    import L1_SE_Greenland_functions as Lf

    import plot_L1_bathymetry as pmc
    pmc.create_bathymetry_plot(Lf, config_dir, L1_model_name, sNx, sNy, faces, face_size_dict)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_L1_SE_Greenland_init_fields(config_dir)
   

