
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import cmocean.cm as cm
import ast
import sys


def create_bathymetry_plot(config_dir, L1_model_name, n_rows, n_cols):

    bathy_file = os.path.join(config_dir,'L1_grid',L1_model_name,'input',L1_model_name+'_bathymetry.bin')
    bathy_grid = np.fromfile(bathy_file,'>f4')
    bathy_grid = np.reshape(bathy_grid,(n_rows,n_cols))
    depth = -1*bathy_grid

    fig = plt.figure(figsize=(8, 6))
    plt.style.use('dark_background')

    # plt.subplot(1, 2, 1)
    plt.contour(depth,levels=[0,500,1000,1500,2000],colors='k',linewidths=0.25)
    C = plt.imshow(depth,origin='lower',cmap = cm.deep)#,vmin=vmin,vmax=vmax)
    plt.colorbar(C)
    plt.title('Bathymetry (500m Contours)')
    plt.gca().set_xticks(np.arange(0, np.shape(depth)[1],90))
    plt.gca().set_yticks(np.arange(0, np.shape(depth)[0], 90))
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])

    # plt.subplot(1, 2, 2)
    # C = plt.imshow(wet_grid, origin='lower')  # ,vmin=vmin,vmax=vmax)
    # plt.colorbar(C)
    # plt.title('Wet Grid')
    # plt.gca().set_xticklabels([])
    # plt.gca().set_yticklabels([])

    plt.grid(linestyle='--',alpha=0.5,color='silver')

    output_file = os.path.join(config_dir, 'L1_grid', L1_model_name, 'plots', 'init_files', L1_model_name+'_bathymetry.png')
    plt.savefig(output_file,bbox_inches = 'tight')
    plt.close(fig)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_plot(config_dir)


