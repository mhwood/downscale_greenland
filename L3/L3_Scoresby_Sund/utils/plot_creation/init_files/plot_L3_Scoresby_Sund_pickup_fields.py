
import os
import argparse
import sys
import numpy as np

def plot_L3_Scoresby_Sund_pickup(config_dir, pickup_iteration):

    L3_model_name = 'L3_Scoresby_Sund'

    sys.path.insert(1, os.path.join(config_dir, 'L3', 'utils', 'plot_creation', 'init_files'))

    import plot_L3_pickup_fields as pp
    pp.create_pickup_plot(config_dir, L3_model_name, pickup_iteration)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-i", "--pickup_iteration", action="store",
                        help="The iteration number of the pickup file.", dest="pickup_iteration",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    pickup_iteration = args.pickup_iteration

    plot_L3_Scoresby_Sund_init_fields(config_dir,pickup_iteration)
   

