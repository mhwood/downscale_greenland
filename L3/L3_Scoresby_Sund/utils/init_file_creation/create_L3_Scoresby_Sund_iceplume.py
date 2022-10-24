
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys


def create_L3_iceplume_files(config_dir, mankoff_dir, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L3', 'utils','init_file_creation'))
    import create_L3_iceplume_files as cip

    model_name = 'L3_Scoresby_Sund'

    years = [1992]

    cip.create_L3_iceplume_files(config_dir, model_name, mankoff_dir, years, print_level)








if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_L3_iceplume_files(config_dir)
   

