
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import cmocean.cm as cm
import ast
import sys


def read_mask_to_domain(config_dir, model_name, mask_name, n_rows, n_cols):
    if mask_name!='surface':
        mask_file = os.path.join(config_dir, 'L2', model_name, 'input', 'dv', 'L3_'+mask_name+'_BC_mask.bin')
    else:
        mask_file = os.path.join(config_dir, 'L2', model_name, 'input', 'dv', 'L3_' + mask_name + '_mask.bin')
    mask_grid = np.fromfile(mask_file, '>f4')
    mask_grid = np.reshape(mask_grid, (n_rows, n_cols))

    mask_grid = np.ma.masked_where(mask_grid==0,mask_grid)

    return(mask_grid)

def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    Depth = ds.variables['Depth'][:,:]
    ds.close()
    return(XC, YC, Depth)

def plot_dv_masks(config_dir, L2_model_name, L3_model_name, print_level):

    L2_XC, L2_YC, L2_Depth = read_grid_geometry_from_nc(config_dir, L2_model_name)
    L3_XC, L3_YC, L3_Depth = read_grid_geometry_from_nc(config_dir, L3_model_name)

    n_rows = np.shape(L2_XC)[0]
    n_cols = np.shape(L2_XC)[1]

    for boundary in ['south','east','north','surface']:
        if print_level>=1:
            print('    - Plotting the '+boundary+' mask')
        mask_grid = read_mask_to_domain(config_dir, L2_model_name, boundary, n_rows, n_cols)
        fig = plt.figure(figsize=(15, 7))
        plt.style.use('dark_background')

        plt.subplot(1, 2, 1)
        plt.contour(L2_Depth,levels=[1],colors='k',linewidths=0.25)
        plt.imshow(L2_Depth,origin='lower',cmap = cm.deep, alpha=0.5)#,vmin=vmin,vmax=vmax)
        C = plt.imshow(mask_grid,origin='lower',cmap='jet')
        plt.colorbar(C)
        plt.title(L2_model_name+' Domain\n'+boundary+' mask')
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

        plt.subplot(1, 2, 2)
        C = plt.pcolormesh(L3_XC, L3_YC, L3_Depth, cmap=cm.deep, shading='nearest', alpha=0.7)  # ,vmin=vmin,vmax=vmax)
        plt.contour(L3_XC, L3_YC, L3_Depth, levels=[1], colors='k', linewidths=0.25)
        plt.colorbar(C)

        rows,cols = np.where(mask_grid>0)
        x = L2_XC[rows, cols]
        y = L2_YC[rows, cols]
        c = mask_grid[rows,cols]
        plt.scatter(x,y,s=3,c=c,vmin=1,vmax=np.max(mask_grid),cmap='jet')

        plt.title(L3_model_name + ' Domain')
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

        output_file = os.path.join(config_dir, 'L2', L2_model_name, 'plots', 'init_files', L2_model_name+'_dv_mask_'+boundary+'.png')
        plt.savefig(output_file,bbox_inches = 'tight')
        plt.close(fig)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_plot(config_dir)


