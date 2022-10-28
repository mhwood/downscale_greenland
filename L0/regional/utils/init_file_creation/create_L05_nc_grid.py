
import os
import numpy as np
import netCDF4 as nc4
import argparse
import matplotlib.pyplot as plt

def read_stitched_grid_geometry(config_dir,model_name,sNx, sNy,
                            ordered_nonblank_tiles,ordered_nonblack_rotations):

    Depth_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    XC_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    YC_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L1', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    Depth = ds.variables['Depth'][:, :]
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
                    ds.close()

                    for i in range(ordered_nonblack_rotations[r][c]):
                        Depth = np.rot90(Depth)
                        XC = np.rot90(XC)
                        YC = np.rot90(YC)

                    Depth_grid[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = Depth
                    XC_grid[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = XC
                    YC_grid[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = YC

    # C = plt.imshow(Depth_grid, origin='lower')
    # plt.colorbar(C)
    # plt.show()

    output_grids = [Depth_grid,XC_grid,YC_grid]
    output_names = ['Depth','XC','YC']

    return(output_grids,output_names)

def output_grid_to_nc(output_path, output_grids,output_names):

    ds = nc4.Dataset(output_path,'w')

    ds.createDimension('X',np.shape(output_grids[output_names.index('XC')])[1])
    ds.createDimension('Y', np.shape(output_grids[output_names.index('YC')])[0])

    Depth_var = ds.createVariable('Depth','f4',('Y','X'))
    Depth_var[:, :] = output_grids[output_names.index('Depth')]

    XC_var = ds.createVariable('XC', 'f4', ('Y', 'X'))
    XC_var[:, :] = output_grids[output_names.index('XC')]

    YC_var = ds.createVariable('YC', 'f4', ('Y', 'X'))
    YC_var[:, :] = output_grids[output_names.index('YC')]

    ds.close()

    a=1

def create_nc_grid(config_dir,model_name):

    ordered_nonblank_tiles = [[55, 56, 57], [91, 85, 79]]
    ordered_nonblack_rotations = [[0, 0, 0],[3, 3, 3]]

    sNx = 180
    sNy = 180

    output_grids,output_names = read_stitched_grid_geometry(config_dir,model_name, sNx, sNy,
                                                            ordered_nonblank_tiles,ordered_nonblack_rotations)

    output_path = os.path.join(config_dir,'L1',model_name,'input',model_name+'_grid.nc')
    output_grid_to_nc(output_path, output_grids, output_names)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-m", "--model_name", action="store",
                        help="The name of the L1 model.", dest="model_name",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    model_name = args.model_name

    create_nc_grid(config_dir, model_name)
   

