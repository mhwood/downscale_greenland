
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import simplegrid as sg
from pyproj import Transformer
import sys

def create_L1_CE_Greenland_mitgrid_file(config_dir, L1_model_name, ecco_dir,
                                         sNx,sNy,global_ordered_nonblank_tiles,global_tile_face_index_dict,
                                         print_level):

    llc = 1080

    #####################################################################################
    # doing this manually for now because crunched on time
    # could easily be generalized

    XC_grid = np.zeros((len(global_ordered_nonblank_tiles) * sNy,
                        len(global_ordered_nonblank_tiles[0]) * sNx))
    YC_grid = np.zeros((len(global_ordered_nonblank_tiles) * sNy,
                        len(global_ordered_nonblank_tiles[0]) * sNx))
    XG_grid = np.zeros((len(global_ordered_nonblank_tiles) * sNy+1,
                        len(global_ordered_nonblank_tiles[0]) * sNx+1))
    YG_grid = np.zeros((len(global_ordered_nonblank_tiles) * sNy+1,
                        len(global_ordered_nonblank_tiles[0]) * sNx+1))

    #####################################################################################
    # add face 1

    grid_file = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'mitgrid_tiles',
                             'tile' + '{:03d}'.format(1) + '.mitgrid')
    grid = np.fromfile(grid_file, '>f8')
    grid = np.reshape(grid, (16, 3 * llc + 1, llc + 1))
    grid_subset = grid[:,3060:,:541]
    XC_grid[:sNy,:] = grid_subset[0,:-1,:-1]
    YC_grid[:sNy, :] = grid_subset[1, :-1, :-1]
    XG_grid[:sNy, :] = grid_subset[5, :-1, :]
    YG_grid[:sNy, :] = grid_subset[6, :-1, :]

    #####################################################################################
    # add face 3

    grid_file = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'mitgrid_tiles',
                             'tile' + '{:03d}'.format(3) + '.mitgrid')
    grid = np.fromfile(grid_file, '>f8')
    grid = np.reshape(grid, (16, llc + 1, llc + 1))
    grid_subset = grid[:, 540:, :181]
    XC_grid[sNy:, :] = np.rot90(grid_subset[0, :-1, :-1], k=3)
    YC_grid[sNy:, :] = np.rot90(grid_subset[1, :-1, :-1], k=3)
    XG_grid[sNy:, :] = np.rot90(grid_subset[5, :, :], k=3)
    YG_grid[sNy:, :] = np.rot90(grid_subset[6, :, :], k=3)

    #####################################################################################
    # plot it for sanity

    # plt.subplot(2,2,1)
    # C = plt.imshow(XC_grid,origin='lower')
    # plt.colorbar(C)
    #
    # plt.subplot(2, 2, 2)
    # C = plt.imshow(YC_grid, origin='lower')
    # plt.colorbar(C)
    #
    # plt.subplot(2, 2, 3)
    # C = plt.imshow(XG_grid, origin='lower')
    # plt.colorbar(C)
    #
    # plt.subplot(2, 2, 4)
    # C = plt.imshow(YG_grid, origin='lower')
    # plt.colorbar(C)
    #
    # plt.show()

    #####################################################################################

    XC_grid = np.flipud(XC_grid)
    YC_grid = np.flipud(YC_grid)
    XG_grid = np.flipud(XG_grid)
    YG_grid = np.flipud(YG_grid)

    mitgrid_matrices = dict()
    mitgrid_matrices['XG'] = XG_grid
    mitgrid_matrices['YG'] = YG_grid
    mitgrid_matrices['XC'] = XC_grid
    mitgrid_matrices['YC'] = YC_grid

    output_file = os.path.join(config_dir,'L1_grid',L1_model_name,'input',L1_model_name+'.mitgrid')

    mg_new_L0_grid_L1_domain, n_rows, n_cols = \
        sg.regrid.regrid(mitgrid_matrices=mitgrid_matrices,
                         lon_subscale=1, lat_subscale=1,
                         lon1=XG_grid[0, 0], lat1=YG_grid[0, 0],
                         lon2=XG_grid[-1, -1], lat2=YG_grid[-1, -1],
                         verbose=False, outfile=output_file)

    sg.gridio.write_mitgridfile(output_file, mg_new_L0_grid_L1_domain, n_rows, n_cols)





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

    create_L1_CE_Greenland_mitgrid_files(config_dir,ecco_dir)
   

