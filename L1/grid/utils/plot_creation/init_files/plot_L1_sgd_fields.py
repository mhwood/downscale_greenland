
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from MITgcmutils import mds
from random import randint
import cmocean.cm as cm
import sys
import argparse
from matplotlib.gridspec import GridSpec

def read_L1_geometry_from_grid(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:, :]
    drF = ds.variables['drF'][:]
    Depth = ds.variables['Depth'][:, :]
    ds.close()

    return(XC, YC, Depth, drF)


def read_sgd_mask(Lf, sgd_dir, depths, sNx, sNy, faces, face_size_dict):

    file_path = os.path.join(sgd_dir, 'L1_sgd_mask.bin')
    mask_compact = np.fromfile(file_path,'>f4')
    n_rows = int(np.size(mask_compact)/(depths*sNx))
    mask_compact = np.reshape(mask_compact,(depths,n_rows,sNx))
    mask_field = Lf.read_compact_grid_to_stitched_grid(mask_compact, sNx, sNy, faces, face_size_dict, dim=3)

    count_grid = np.zeros((np.shape(mask_field)[1], np.shape(mask_field)[2]))
    depth_grid = np.zeros((np.shape(mask_field)[1], np.shape(mask_field)[2]))

    for i in range(np.shape(mask_field)[0]):
        count_grid += (mask_field[i,:,:]!=0).astype(int)
        depths = (i+1)*(mask_field[i, :, :] != 0).astype(int)
        depth_grid[depths!=0] = depths[depths!=0]

    return(count_grid, depth_grid)


def create_sgd_plot(Lf, config_dir, L1_model_name, sNx, sNy, faces, face_size_dict):

    XC, YC, Depth, drF = read_L1_geometry_from_grid(config_dir, L1_model_name)
    rows = np.shape(XC)[0]
    cols = np.shape(XC)[1]
    depths = np.size(drF)

    sgd_dir = os.path.join(config_dir,'L1',L1_model_name,'input','sgd')
    count_grid, depth_grid = read_sgd_mask(Lf, sgd_dir, depths, sNx, sNy, faces, face_size_dict)

    indices = count_grid!=0

    discharge_points = np.column_stack([XC[indices],YC[indices]])

    fig = plt.figure(figsize=(8, 12))
    plt.style.use("dark_background")

    gs = GridSpec(8, 4, left=0.05, right=0.95, wspace=0.05, hspace=0.05)

    ax0 = fig.add_subplot(gs[:3, :2])
    plt.pcolormesh(count_grid, cmap='jet', shading='nearest')

    ax1 = fig.add_subplot(gs[:3, 2:])
    plt.pcolormesh(depth_grid, cmap='jet', shading='nearest')

    ax2 = fig.add_subplot(gs[3:6, :])
    plt.pcolormesh(XC,YC,Depth,cmap=cm.deep,shading='nearest')
    plt.plot(discharge_points[:,0],discharge_points[:,1],'r.')


    # year = 2002
    # timestep = 7
    #
    # counter = 1
    # exf_dir = os.path.join(config_dir,'L1',L1_model_name,'input','exf')
    #
    # var_grids = []
    # for ff in range(len(var_names)):
    #     var_grid_subset = read_exf_fields(Lf, exf_dir, var_names[ff], exf_shape, year, timestep,
    #                                       sNx, sNy, faces, face_size_dict)
    #     print('   - The '+var_names[ff]+' grid has values in the range '+str(np.min(var_grid_subset))+' to '+str(np.max(var_grid_subset)))
    #     var_grids.append(var_grid_subset)
    #
    # for ff in range(len(var_names)):
    #
    #     plt.subplot(2, 4, counter)
    #     cmap = 'viridis'
    #
    #     var_grid_subset = var_grids[ff]
    #
    #     if var_names[ff] in ['UWIND','VWIND']:
    #         cmap = cm.balance
    #         val = np.max(np.abs(var_grid_subset[var_grid_subset != 0]))
    #         vmin = -val
    #         vmax = val
    #     else:
    #         if np.any(var_grids[ff][ :, :])>0:
    #             vmin = np.min(var_grid_subset[var_grid_subset != 0])
    #             vmax = np.max(var_grid_subset[var_grid_subset != 0])
    #         else:
    #             vmin=-0.1
    #             vmax=0.1
    #
    #     if var_names[ff] in ['AQH']:
    #         cmap = cm.haline
    #     if var_names[ff] in ['RUNOFF']:
    #         cmap = cm.rain
    #     if var_names[ff] in ['PRECIP']:
    #         cmap = cm.rain
    #     if var_names[ff] in ['ATEMP']:
    #         cmap = cm.thermal
    #     if var_names[ff] in ['SWDOWN','LWDOWN']:
    #         cmap = cm.thermal
    #
    #     C = plt.imshow(var_grid_subset, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
    #
    #     plt.gca().set_xticklabels([])
    #     plt.gca().set_yticklabels([])
    #
    #     plt.colorbar(C)
    #     plt.title(var_names[ff])
    #     counter += 1

    output_dir = os.path.join(config_dir, 'L1', L1_model_name, 'plots', 'init_files')
    plt.savefig(os.path.join(output_dir, L1_model_name+'_sgd_fields.png'))#, bbox_inches='tight')
    plt.close(fig)


########################################################################################################################
# User inputs


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_plot(config_dir)

