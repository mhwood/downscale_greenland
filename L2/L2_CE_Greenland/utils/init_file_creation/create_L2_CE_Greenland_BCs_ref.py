
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys


def create_BCs(config_dir):

    sys.path.insert(1, os.path.join(config_dir,'L2', 'utils','init_file_creation'))
    import create_L2_daily_obcs_ref as cbcr

    L05_model_name = 'L05_CE_Greenland'
    L2_model_name = 'L2_CE_Greenland'

    print('Creating the BCs ref file for the '+L2_model_name+' model')

    Nr = 50
    averaging_period = 3600
    seconds_per_iter = 600

    # pass to general function to generate mitgrid
    cbcr.create_bc_fields_reference_dict(config_dir, L05_model_name, Nr, averaging_period, seconds_per_iter)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L05, L05, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_BCs(config_dir)
   

