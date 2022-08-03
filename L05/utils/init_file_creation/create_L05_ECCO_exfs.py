
import os
import numpy as np
import netCDF4 as nc4
#import ast
import matplotlib.pyplot as plt
#from MITgcmutils import mds
#from scipy.interpolate import griddata
import sys


def read_grid_tile_geometry(config_dir,model_name,sNx,sNy,ordered_nonblank_tiles,ordered_nonblank_rotations):

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    stitched_XC = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    stitched_YC = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    # stitched_mask = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    grid_dir = os.path.join(config_dir, 'L05', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
                    # if tile_face < 3:
                    #     if var_name == 'VVEL':
                    #         hFac = 'S'
                    #     elif var_name == 'UVEL':
                    #         hFac = 'W'
                    #     else:
                    #         hFac = 'C'
                    # HFac = ds.variables['HFac' + hFac][:, :, :]
                    # if hFac == 'W':
                    #     HFac = HFac[:,:,1:]
                    # if hFac == 'S':
                    #     HFac = HFac[:,1:,:]
                    ds.close()

                    for i in range(ordered_nonblank_rotations[r][c]):
                        XC = np.rot90(XC)
                        YC = np.rot90(YC)

                    stitched_XC[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = XC
                    stitched_YC[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = YC

    # plt.subplot(1,3,1)
    # C = plt.imshow(stitched_XC,origin='lower')
    # plt.colorbar(C)
    #
    # plt.subplot(1, 3, 2)
    # C = plt.imshow(stitched_YC, origin='lower')
    # plt.colorbar(C)
    #
    # plt.show()

    return(stitched_XC, stitched_YC)

def subset_ecco_exf_to_tile_domain(XC, YC, ecco_lon, ecco_lat, ecco_exf_Grid):

    ecco_Lon, ecco_Lat = np.meshgrid(ecco_lon, ecco_lat)

    # plt.subplot(2,2,1)
    # C = plt.imshow(XC,origin='lower')
    # plt.colorbar(C)
    # plt.subplot(2, 2, 2)
    # C = plt.imshow(YC, origin='lower')
    # plt.colorbar(C)
    # plt.subplot(2, 2, 3)
    # C = plt.imshow(ecco_Lon, origin='lower')
    # plt.colorbar(C)
    # plt.subplot(2, 2, 4)
    # C = plt.imshow(ecco_Lat, origin='lower')
    # plt.colorbar(C)
    # plt.show()

    ll_dist = ((ecco_Lon-np.min(XC))**2 + (ecco_Lat-np.min(YC))**2)**0.5
    ll_index_row, ll_index_col = np.where(ll_dist == np.min(ll_dist))

    ul_dist = ((ecco_Lon - np.min(XC)) ** 2 + (ecco_Lat - np.max(YC)) ** 2) ** 0.5
    ul_index_row, ul_index_col = np.where(ul_dist == np.min(ul_dist))

    lr_dist = ((ecco_Lon - np.max(XC)) ** 2 + (ecco_Lat - np.min(YC)) ** 2) ** 0.5
    lr_index_row, lr_index_col = np.where(lr_dist == np.min(lr_dist))

    ur_dist = ((ecco_Lon - np.max(XC)) ** 2 + (ecco_Lat - np.max(YC)) ** 2) ** 0.5
    ur_index_row, ur_index_col = np.where(ur_dist == np.min(ur_dist))

    min_y_index = np.min([ll_index_row[0],ul_index_row[0],lr_index_row[0],ur_index_row[0]])
    max_y_index = np.max([ll_index_row[0], ul_index_row[0], lr_index_row[0], ur_index_row[0]])
    min_x_index = np.min([ll_index_col[0], ul_index_col[0], lr_index_col[0], ur_index_col[0]])
    max_x_index = np.max([ll_index_col[0], ul_index_col[0], lr_index_col[0], ur_index_col[0]])

    buffer = 3
    if min_x_index - buffer>=0:
        min_x_index -= buffer
    if max_x_index + buffer<= np.shape(ecco_Lon)[1]-1:
        max_x_index += buffer
    if min_y_index - buffer>=0:
        min_y_index -= buffer
    if max_y_index + buffer<= np.shape(ecco_Lon)[0]-1:
        max_y_index += buffer

    buffer = 1
    if min_x_index - buffer >= 0:
        min_x_index -= buffer
    if max_x_index + buffer <= np.shape(ecco_Lon)[1] - 1:
        max_x_index += buffer
    if min_y_index - buffer >= 0:
        min_y_index -= buffer
    if max_y_index + buffer <= np.shape(ecco_Lon)[0] - 1:
        max_y_index += buffer

    ecco_Lon_subset = ecco_Lon[min_y_index:max_y_index,min_x_index:max_x_index]
    ecco_Lat_subset = ecco_Lat[min_y_index:max_y_index,min_x_index:max_x_index]
    ecco_grid_subset = ecco_exf_Grid[:, min_y_index:max_y_index, min_x_index:max_x_index]

    # plt.pcolormesh(ecco_Lon_subset,ecco_Lat_subset,ecco_grid_subset[0,:,:])
    # plt.plot(XC[:,0],YC[:,0],'r-')
    # plt.plot(XC[:, -1], YC[:, -1], 'r-')
    # plt.plot(XC[0,:], YC[0,:], 'r-')
    # plt.plot(XC[-1, :], YC[-1, :], 'r-')
    # plt.show()

    lon_subset = ecco_lon[min_x_index:max_x_index]
    lat_subset = ecco_lat[min_y_index:max_y_index]
    print(' lon0 = '+str(lon_subset[0]))
    print(' lon_inc = '+str(np.diff(lon_subset)[0]))
    print(' nlon = '+str(len(lon_subset)))
    print(' lat0 = ' + str(lat_subset[0]))
    print(' lat_inc = ' + str(np.diff(lat_subset)[0]))
    print(' nlat = ' + str(len(lat_subset)))
    if min_y_index<6 or max_y_index>len(ecco_lat)-6:
        print('Need to be careful with the lat incs!')

    return(ecco_grid_subset)

def create_L05_exfs(config_dir,model_name,var_name,file_prefix,year,
                    sNx,sNy,ordered_nonblank_tiles,ordered_nonblank_rotations, is_runoff=False):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import ecco_functions as ec

    print('    - Creating the '+var_name+' exf files in year '+str(year)+' for the '+model_name+' model from ECCOv5 data')

    # step 0: get the model domain
    print('    - Reading in the model geometry')
    XC, YC = read_grid_tile_geometry(config_dir,model_name,sNx,sNy,ordered_nonblank_tiles,ordered_nonblank_rotations)

    print('    - Reading in the ECCO exf grid')
    ecco_dir = '/Users/michwood/Documents/Research/Projects/Ocean_Modeling/ECCO'
    ecco_lon, ecco_lat, ecco_exf_Grid = ec.read_ecco_exf_file(ecco_dir, file_prefix, year, llc=270)

    ecco_grid_subset = subset_ecco_exf_to_tile_domain(XC, YC, ecco_lon, ecco_lat, ecco_exf_Grid)

    output_file = os.path.join(config_dir,'L05',model_name,'input','exf','L05_exf_'+var_name+'_'+str(year))
    ecco_grid_subset.ravel(order='C').astype('>f4').tofile(output_file)




