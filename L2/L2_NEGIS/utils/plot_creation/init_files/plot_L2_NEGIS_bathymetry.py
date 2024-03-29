
import os
import argparse
import sys
import numpy as np

def plot_L2_NEGIS_bathymetry(config_dir):

    L2_model_name = 'L2_NEGIS'

    sys.path.insert(1, os.path.join(config_dir, 'L2', 'utils', 'plot_creation', 'init_files'))

    import plot_L2_bathymetry as pb
    pb.create_bathymetry_plot(config_dir, L2_model_name)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_L2_NEGIS_bathymetry(config_dir)
   

