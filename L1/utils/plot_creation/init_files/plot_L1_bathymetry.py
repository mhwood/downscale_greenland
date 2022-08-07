
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import cmocean.cm as cm
import ast
import sys


def create_bathymetry_plot(config_dir, L1_model_name, sNx, sNy, faces, face_size_dict):
    sys.path.insert(1, os.path.join(config_dir, 'L1', L1_model_name, 'utils'))
    import L1_CE_Greenland_functions as Lf

    bathy_file = os.path.join(config_dir,'L1',L1_model_name,'input','bathymetry.bin')
    bathy_compact = np.fromfile(bathy_file,'>f4')
    n_rows = int(np.size(bathy_compact)/sNx)
    bathy_compact = np.reshape(bathy_compact,(n_rows,sNx))

    print(np.shape(bathy_compact))

    bathy_grid = Lf.read_compact_grid_to_stitched_grid(bathy_compact, sNx, sNy, faces, face_size_dict)
    depth = -1*bathy_grid

    fig = plt.figure(figsize=(8, 6))
    plt.style.use('dark_background')

    # plt.subplot(1, 2, 1)
    plt.contour(depth,levels=[0,500,1000,1500,2000],colors='k',linewidths=0.25)
    C = plt.imshow(depth,origin='lower',cmap = cm.deep)#,vmin=vmin,vmax=vmax)
    plt.colorbar(C)
    plt.title('Bathymetry (500m Contours)')
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])

    # plt.subplot(1, 2, 2)
    # C = plt.imshow(wet_grid, origin='lower')  # ,vmin=vmin,vmax=vmax)
    # plt.colorbar(C)
    # plt.title('Wet Grid')
    # plt.gca().set_xticklabels([])
    # plt.gca().set_yticklabels([])

    output_file = os.path.join(config_dir, 'L1', L1_model_name, 'plots', 'init_files', L1_model_name+'_bathymetry.png')
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


