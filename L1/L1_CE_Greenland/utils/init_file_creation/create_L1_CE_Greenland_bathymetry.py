
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys

def create_L1_CE_Greenland_bathymetry(config_dir, L1_model_name, ecco_dir, sNx, sNy, print_level):
    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils', 'init_file_creation'))
    import create_L1_bathymetry as cb

    llc = 1080

    cb.create_bathymetry(config_dir, L1_model_name, ecco_dir, sNx, sNy, llc, print_level)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_dir", action="store",
                        help="The directory where the ECCO mitgrid files are stored.", dest="ecco_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    ecco_dir = args.ecco_dir

    create_L1_CE_Greenland_bathymetry(config_dir,ecco_dir)
   

