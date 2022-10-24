
import os
import numpy as np
import netCDF4 as nc4
import argparse
import ast

########################################################################################################################

def read_L1_grid_tile_geometry(config_dir, model_name, Nr, sNx, sNy, ordered_nonblank_tiles, ordered_tile_rotations):

    stitched_XC = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    stitched_YC = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    stitched_AngleCS = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    stitched_AngleSN = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    stitched_DXC = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    stitched_DYC = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    stitched_rAz = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    stitched_Depth = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    # stitched_hFacC = np.zeros((Nr, sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    # stitched_hFacS = np.zeros((Nr, sNy * len(ordered_nonblank_tiles) + 1, sNx * len(ordered_nonblank_tiles[0])))
    # stitched_hFacW = np.zeros((Nr, sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0]) + 1))

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L1', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            rotations = ordered_tile_rotations[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
                    AngleCS = ds.variables['AngleCS'][:, :]
                    AngleSN = ds.variables['AngleSN'][:, :]
                    DXC = ds.variables['dxC'][:, :]
                    DYC = ds.variables['dyC'][:, :]
                    rAz = ds.variables['rAz'][:, :]
                    # hFacC = ds.variables['HFacC'][:, :, :]
                    # hFacS = ds.variables['HFacS'][:, :, :]
                    # hFacW = ds.variables['HFacW'][:, :, :]
                    DRF = ds.variables['drF'][:]
                    Depth = ds.variables['Depth'][:, :]
                    ds.close()

                    DXC = DXC[:,:-1]
                    DYC = DYC[:-1, :]
                    rAz = rAz[:sNy, :sNx]

                    for i in range(rotations):
                        XC = np.rot90(XC)
                        YC = np.rot90(YC)
                        AngleCS = np.rot90(AngleCS)
                        AngleSN = np.rot90(AngleSN)
                        Depth = np.rot90(Depth)

                        DXC = np.rot90(DXC)
                        DYC = np.rot90(DYC)

                        rAz = np.rot90(rAz)

                        # hFacC = np.rot90(hFacC, axes=(1, 2))
                        # hFacS = np.rot90(hFacS, axes=(1, 2))
                        # hFacW = np.rot90(hFacW, axes=(1, 2))

                    stitched_XC[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = XC
                    stitched_YC[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = YC

                    stitched_AngleCS[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = AngleCS
                    stitched_AngleSN[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = AngleSN

                    stitched_DXC[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = DXC
                    stitched_DYC[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = DYC

                    stitched_rAz[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = rAz

                    # stitched_hFacC[:, r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = hFacC
                    # stitched_hFacS[:, r * sNy:(r + 1) * sNy + 1, c * sNx:(c + 1) * sNx] = hFacS
                    # stitched_hFacW[:, r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx + 1] = hFacW

                    stitched_Depth[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = Depth

    return(stitched_XC, stitched_YC, stitched_AngleCS, stitched_AngleSN, stitched_Depth, stitched_DXC, stitched_DYC, stitched_rAz, DRF)
            #stitched_DXC, stitched_DYC,
            #stitched_hFacC, stitched_hFacS, stitched_hFacW,


def write_grid_to_nc(config_dir, model_name,
                     XC, YC, AngleCS, AngleSN, DXC, DYC, RAZ, DRF, Depth):

    output_path = os.path.join(config_dir, 'nc_grids', model_name+'_grid.nc')

    ds = nc4.Dataset(output_path,'w')

    ds.createDimension('X', np.shape(XC)[1])
    ds.createDimension('Y', np.shape(XC)[0])
    # ds.createDimension('Xp1', np.shape(DXC)[1])
    # ds.createDimension('Yp1', np.shape(DYC)[0])
    ds.createDimension('Z',np.shape(DRF)[0])

    var = ds.createVariable('XC','f4',('Y','X'))
    var[:,:] = XC

    var = ds.createVariable('YC', 'f4', ('Y', 'X'))
    var[:, :] = YC

    var = ds.createVariable('AngleCS', 'f4', ('Y', 'X'))
    var[:, :] = AngleCS

    var = ds.createVariable('AngleSN', 'f4', ('Y', 'X'))
    var[:, :] = AngleSN

    var = ds.createVariable('dxC', 'f4', ('Y', 'X'))
    var[:, :] = DXC

    var = ds.createVariable('dyC', 'f4', ('Y', 'X'))
    var[:, :] = DYC

    var = ds.createVariable('rAz', 'f4', ('Y', 'X'))
    var[:, :] = RAZ

    # var = ds.createVariable('HFacC', 'f4', ('Z', 'Y', 'X'))
    # var[:, :, :] = hFacC
    #
    # var = ds.createVariable('HFacW', 'f4', ('Z', 'Y', 'Xp1'))
    # var[:, :, :] = hFacW
    #
    # var = ds.createVariable('HFacS', 'f4', ('Z', 'Yp1', 'X'))
    # var[:, :, :] = hFacS

    var = ds.createVariable('drF', 'f4', ('Z',))
    var[:] = DRF

    var = ds.createVariable('Depth', 'f4', ('Y', 'X'))
    var[:, :] = Depth

    ds.close()

def stitch_grid_files(config_dir):

    print('Stitching the nc grid files')

    model_name = 'L1_W_Greenland'
    ordered_tiles = [[6,9,12],[5,8,11],[4,7,10],[3,2,1]]
    ordered_tile_rotations = [[1,1,1],[1,1,1],[1,1,1],[2,2,2]]
    Nr = 50
    sNx = 180
    sNy = 180

    XC, YC, AngleCS, AngleSN, Depth, DXC, DYC, RAZ, DRF = read_L1_grid_tile_geometry(config_dir, model_name, Nr, sNx, sNy, ordered_tiles, ordered_tile_rotations)

    write_grid_to_nc(config_dir, model_name,
                     XC, YC, AngleCS, AngleSN,
                     DXC, DYC, RAZ,
                     DRF, Depth)

    # zero_rows = 0
    # for i in range(np.shape(hFacC)[0]):
    #     if np.all(hFacC[i,:,:]==0):
    #         zero_rows+=1
    #
    # if zero_rows>1:
    #     print('    - The grid has '+str(zero_rows)+' zero rows - consider chopping off some of them to reduce computational time')




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    stitch_grid_files(config_dir)