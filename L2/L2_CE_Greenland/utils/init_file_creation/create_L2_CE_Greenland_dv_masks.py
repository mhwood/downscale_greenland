import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys


def create_dv_masks(config_dir):
    sys.path.insert(1, os.path.join(config_dir, 'L2', 'utils', 'init_file_creation'))
    import create_L2_dv_masks as cdv

    L2_model_name = 'L2_CE_Greenland'
    L3_model_name = 'L3_Scoresby_Sund'

    print('Creating the diagnostics_vec masks for the ' + L2_model_name + ' model')

    # pass to general function
    cdv.create_L2_diagnostic_vec_masks(config_dir, L2_model_name, L3_model_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L05, L05, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_dv_masks(config_dir)


