
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
import argparse
from pyproj import Transformer
from MITgcmutils import mds
import sys

def read_grid_tile_angles_to_faces(config_dir,model_name,sNx,sNy,ordered_nonblank_tiles):
    ordered_AngleCS_tiles = []
    ordered_AngleSN_tiles = []

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    angleCS_grid = np.zeros((sNy*len(ordered_nonblank_tiles),
                             sNx*len(len(ordered_nonblank_tiles[0]))))
    angleSN_grid = np.zeros((sNy * len(ordered_nonblank_tiles),
                             sNx * len(len(ordered_nonblank_tiles[0]))))

    grid_dir = os.path.join(config_dir, 'L1', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        row_XCs = []
        row_YCs = []
        row_AngleCSs = []
        row_AngleSNs = []
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    AngleCS = ds.variables['AngleCS'][:, :]
                    AngleSN = ds.variables['AngleSN'][:, :]
                    ds.close()
                    row_XCs.append(XC)
                    row_YCs.append(YC)
                    row_AngleCSs.append(AngleCS)
                    row_AngleSNs.append(AngleSN)
        ordered_AngleCS_tiles.append(row_AngleCSs)
        ordered_AngleSN_tiles.append(row_AngleSNs)

    return(ordered_XC_tiles, ordered_YC_tiles, ordered_AngleCS_tiles, ordered_AngleSN_tiles, delR)


def read_pickup_file_to_compact(pickup_file_path):

    Nr = 7
    print('      Reading from '+pickup_file_path)
    global_data, _, global_metadata = mds.rdmds(pickup_file_path, returnmeta=True)

    var_names = []
    row_bounds = []
    var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        if var_name.lower() == 'sitices':
            end_row = start_row + Nr
        else:
            end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]
        var_grids.append(var_grid)
        row_bounds.append([start_row,end_row])
        start_row=end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,var_grids,global_metadata)


def read_pickup_file_to_faces(config_dir,config_name):

    sNx = 180
    sNy = 180

    pickup_file = 'pickup_seaice.0000000001'
    pickup_file_path = os.path.join(config_dir,'L1',config_name, 'input', pickup_file)
    var_names, row_bounds, compact_var_grids, global_metadata = read_pickup_file_to_compact(pickup_file_path)

    var_grids = []

    for vn in range(len(var_names)):
        compact_var_grid = compact_var_grids[vn]

        face_1 = compact_var_grid[:, :3 * sNx, :]
        face_1 = np.reshape(face_1, (np.shape(face_1)[0], sNx, 3 * sNx))

        face_3 = np.rot90(compact_var_grid[:, 3 * sNx:, :], axes=(2, 1))

        full_grid = np.concatenate([face_1, face_3], axis=1)
        var_grids.append(full_grid)

    return (var_names, var_grids)


def read_angles_to_faces():
    a=1

def plot_pickup(config_dir, model_name):

    llc = 1080
    sNx = 180

    var_names, var_grids = read_pickup_file_to_faces(config_dir,model_name)

    cmap_bounds = {'siTICES':[237,273], 'siAREA':[0, 1], 'siHEFF':[0,2],
                   'siHSNOW':[0,0.5], 'siUICE':[-1,1], 'siVICE':[-1,1]}

    fig = plt.figure(figsize=(16,10))
    plt.style.use('dark_background')

    for i in range(len(var_names)):
        if var_names[i] in ['siUICE', 'siVICE']:
            cmap = 'seismic'
        else:
            cmap = 'viridis'

        print(np.min(var_grids[i][var_grids[i]!=0]),np.max(var_grids[i][var_grids[i]!=0]))

        plt.subplot(2,3,i+1)
        C = plt.imshow(var_grids[i][0,:,:],origin='lower',cmap = cmap,
                       vmin = cmap_bounds[var_names[i]][0], vmax = cmap_bounds[var_names[i]][1])
        plt.colorbar(C,orientation='horizontal', fraction=0.05, pad=0.15)
        plt.title(var_names[i])

    plt.suptitle('Seaice Pickup File Fields')

    output_path = os.path.join(config_dir,'L1',model_name,'plots',model_name+'_seaice_pickup.png')
    plt.savefig(output_path)
    plt.close(fig)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    config_name = 'L1_CE_Greenland'

    plot_pickup(config_dir,config_name)
