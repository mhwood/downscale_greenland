
import os
import netCDF4 as nc4
import numpy as np
import argparse
import sys

def plot_exfs(config_dir, print_level):

    L2_model_name = 'L2_Kanger'

    sys.path.insert(1, os.path.join(config_dir, 'L2', 'utils','plot_creation','init_files'))

    year = 1992
    month = 1
    day = 1
    hour = 3

    ds = nc4.Dataset(os.path.join(config_dir,'nc_grids',L2_model_name+'_grid.nc'))
    hFac = ds.variables['HFacC'][:,:,:]
    ds.close()
    n_rows = np.shape(hFac)[1]
    n_cols = np.shape(hFac)[2]

    # # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    import plot_L2_exf_fields as pbc
    pbc.create_exf_plot(config_dir, L2_model_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L2, L2, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_exfs(config_dir, print_level=4)
   

