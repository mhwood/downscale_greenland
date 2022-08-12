
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import cmocean.cm as cm
import ast

def read_L2_bathy_and_wetgrid_from_grid(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    depth = ds.variables['Depth'][:, :]
    wet_grid = ds.variables['HFacC'][:, :, :]
    ds.close()
    wet_grid = wet_grid[0, :, :]
    wet_grid[wet_grid>0] = 1
    return(depth, wet_grid)

def create_bathymetry_plot(config_dir, L2_model_name):

    print('    - Reading the bathymetry and wet mask from the stiched nc grid')
    depth, wet_grid = read_L2_bathy_and_wetgrid_from_grid(config_dir, L2_model_name)

    print('    - Plotting the bathymetry and wet mask')
    fig = plt.figure(figsize=(12, 6))
    plt.style.use('dark_background')

    plt.subplot(1, 2, 1)
    plt.contour(depth,levels=[0,250,500,750,1000],colors='k',linewidths=0.25)
    C = plt.imshow(depth,origin='lower',cmap = cm.deep)#,vmin=vmin,vmax=vmax)
    plt.colorbar(C)
    plt.title('Bathymetry (250m Contours)')
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])

    plt.subplot(1, 2, 2)
    C = plt.imshow(wet_grid, origin='lower')  # ,vmin=vmin,vmax=vmax)
    plt.colorbar(C)
    plt.title('Wet Grid')
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])

    output_file = os.path.join(config_dir, 'L2', L2_model_name, 'plots', 'init_files', L2_model_name+'_bathymetry.png')
    plt.savefig(output_file,bbox_inches = 'tight')
    plt.close(fig)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_plot(config_dir)


