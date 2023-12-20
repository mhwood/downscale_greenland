
import os
import numpy as np
import netCDF4 as nc4
import ast
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import sys

def read_grid_geometry(config_dir,model_name, n_rows, n_cols):
    file_path = os.path.join(config_dir, 'mitgrids', model_name + '.mitgrid')
    entire_grid = np.fromfile(file_path, dtype='>f8')
    entire_grid = np.reshape(entire_grid, (16, n_rows + 1, n_cols + 1))
    XC = entire_grid[0, :, :]
    YC = entire_grid[1, :, :]
    XC = XC[:-1, :-1]
    YC = YC[:-1, :-1]

    delR = np.array([1.00, 1.14, 1.30, 1.49, 1.70,
                         1.93, 2.20, 2.50, 2.84, 3.21,
                         3.63, 4.10, 4.61, 5.18, 5.79,
                         6.47, 7.20, 7.98, 8.83, 9.73,
                         10.69, 11.70, 12.76, 13.87, 15.03,
                         16.22, 17.45, 18.70, 19.97, 21.27,
                         22.56, 23.87, 25.17, 26.46, 27.74,
                         29.00, 30.24, 31.45, 32.65, 33.82,
                         34.97, 36.09, 37.20, 38.29, 39.37,
                         40.45, 41.53, 42.62, 43.73, 44.87,
                         46.05, 47.28, 48.56, 49.93, 51.38,
                         52.93, 54.61, 56.42, 58.38, 60.53,
                         62.87, 65.43, 68.24, 71.33, 74.73,
                         78.47, 82.61, 87.17, 92.21, 97.79,
                         103.96, 110.79, 118.35, 126.73, 136.01,
                         146.30, 157.71, 170.35, 184.37, 199.89,
                         217.09, 236.13, 257.21, 280.50, 306.24,
                         334.64, 365.93, 400.38, 438.23, 479.74, ])

    return(XC, YC, delR)

def read_aste_grid_geometry(aste_dir,ordered_aste_tiles,ordered_aste_tile_rotations, L_XC, L_YC):

    aste_Nr = 50
    aste_sNx = 90
    aste_sNy = 90
    aste_XC = np.zeros((aste_sNy*len(ordered_aste_tiles),aste_sNx*len(ordered_aste_tiles[0])))
    aste_YC = np.zeros((aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))
    aste_depth = np.zeros((aste_sNy*len(ordered_aste_tiles),aste_sNx*len(ordered_aste_tiles[0])))

    for r in range(len(ordered_aste_tiles)):
        aste_tile_row = ordered_aste_tiles[r]
        aste_rotation_row = ordered_aste_tile_rotations[r]
        for c in range(len(ordered_aste_tiles[r])):
            file_name = os.path.join(aste_dir,'nctiles_grid','GRID.'+'{:04d}'.format(aste_tile_row[c]))+'.nc'
            ds = nc4.Dataset(file_name)
            depth = ds.variables['Depth'][:,:]
            XC = ds.variables['XC'][:,:]
            YC = ds.variables['YC'][:, :]
            ds.close()
            for n in range(aste_rotation_row[c]):
                depth = np.rot90(depth)
                XC = np.rot90(XC)
                YC = np.rot90(YC)
            aste_depth[r*aste_sNy:(r+1)*aste_sNy,c*aste_sNx:(c+1)*aste_sNx] = depth
            aste_XC[r*aste_sNy:(r+1)*aste_sNy,c*aste_sNx:(c+1)*aste_sNx] = XC
            aste_YC[r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = YC

    # C = plt.imshow(aste_depth,origin='lower')
    # # C = plt.imshow(aste_XC, origin='lower')
    # # C = plt.imshow(aste_YC, origin='lower')
    # plt.colorbar(C)
    # plt.show()

    ll_dist = ((L_XC[0, 0] - aste_XC) ** 2 + (L_YC[0, 0] - aste_YC) ** 2) ** 0.5
    ll_row, ll_col = np.where(ll_dist == np.min(ll_dist))
    ul_dist = ((L_XC[-1, 0] - aste_XC) ** 2 + (L_YC[-1, 0] - aste_YC) ** 2) ** 0.5
    ul_row, ul_col = np.where(ul_dist == np.min(ul_dist))

    lr_dist = ((L_XC[0, -1] - aste_XC) ** 2 + (L_YC[0, -1] - aste_YC) ** 2) ** 0.5
    lr_row, lr_col = np.where(lr_dist == np.min(lr_dist))
    ur_dist = ((L_XC[-1, -1] - aste_XC) ** 2 + (L_YC[-1, -1] - aste_YC) ** 2) ** 0.5
    ur_row, ur_col = np.where(ur_dist == np.min(ur_dist))

    min_row = np.min([ll_row[0], ul_row[0], lr_row[0], ur_row[0]])
    max_row = np.max([ll_row[0], ul_row[0], lr_row[0], ur_row[0]])
    min_col = np.min([ll_col[0], ul_col[0], lr_col[0], ur_col[0]])
    max_col = np.max([ll_col[0], ul_col[0], lr_col[0], ur_col[0]])

    dist_buffer = 3
    aste_XC = aste_XC[min_row - dist_buffer:max_row + dist_buffer,
              min_col - dist_buffer:max_col + dist_buffer]
    aste_YC = aste_YC[min_row - dist_buffer:max_row + dist_buffer,
              min_col - dist_buffer:max_col + dist_buffer]
    aste_depth = aste_depth[min_row - dist_buffer:max_row + dist_buffer,
                min_col - dist_buffer:max_col + dist_buffer]

    aste_delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                        10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                        31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                        93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                        139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                        341.50, 364.50, 387.50, 410.50, 433.50, 456.50])

    return(aste_XC,aste_YC,aste_depth,aste_delR)

def get_ASTE_bathymetry_on_domain(XC,YC,aste_XC, aste_YC, aste_depth):

    points = np.hstack([np.reshape(aste_XC, (np.size(aste_XC), 1)),
                        np.reshape(aste_YC, (np.size(aste_YC), 1))])
    values = np.reshape(aste_depth, (np.size(aste_depth), 1))
    grid = griddata(points, values, (XC, YC), method='nearest')[:,:,0]
    return(grid)


def create_L3_ASTE_wetgrid_file(config_dir,model_name,ordered_aste_tiles,ordered_aste_tile_rotations):

    print('    - Creating the ASTE wetgrid on the '+model_name+' model domain')

    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import downscale_functions as df

    f = open(os.path.join(config_dir, 'domain_sizes.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)
    L_size = size_dict[model_name]
    n_rows = L_size[0]
    n_cols = L_size[1]

    # step 0: get the model domain
    XC, YC, delR = read_grid_geometry(config_dir,model_name, n_rows, n_cols)

    Z_bottom = np.cumsum(delR)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2

    # step 1: stitch the ASTE tiles together around the domain
    aste_dir = '/Users/michwood/Documents/Research/Data Repository/Greenland/Ocean Properties/Models/ASTE'
    aste_XC, aste_YC, aste_depth, _ = read_aste_grid_geometry(aste_dir, ordered_aste_tiles, ordered_aste_tile_rotations, XC, YC)
    aste_depth *= -1

    # step 2: interpolate the ASTE bathymetry on to the domain
    aste_bathy_on_domain = get_ASTE_bathymetry_on_domain(XC, YC, aste_XC, aste_YC, aste_depth)

    # plt.subplot(1,2,1)
    # plt.pcolormesh(aste_XC,aste_YC,aste_depth)
    # plt.plot(XC[:, 0], YC[:, 0], 'k-')
    # plt.plot(XC[:, -1], YC[:, -1], 'k-')
    # plt.plot(XC[0, :], YC[0, :], 'k-')
    # plt.plot(XC[-1, :], YC[-1, :], 'k-')
    #
    # plt.subplot(1,2,2)
    # plt.pcolormesh(XC,YC,aste_bathy_on_domain)
    # plt.show()

    # step 3: create the wet grids
    print('   - Calculating the hFacC wetgrid')
    wet_grid_C = df.create_3D_wet_grid(aste_bathy_on_domain, delR, hFac='C')
    print('   - Calculating the hFacS wetgrid')
    wet_grid_S = df.create_3D_wet_grid(aste_bathy_on_domain, delR, hFac='S')
    print('   - Calculating the hFacW wetgrid')
    wet_grid_W = df.create_3D_wet_grid(aste_bathy_on_domain, delR, hFac='W')


    output_file = os.path.join(config_dir,'L3', model_name, 'input', 'ASTE_wetgrid_on_'+model_name+'.nc')
    ds = nc4.Dataset(output_file, 'w')

    ds.createDimension('x', np.shape(XC)[1])
    ds.createDimension('y', np.shape(XC)[0])
    ds.createDimension('z', len(delR))

    var = ds.createVariable('XC', 'f4', ('y', 'x'))
    var[:] = XC
    var = ds.createVariable('YC', 'f4', ('y', 'x'))
    var[:] = YC
    var = ds.createVariable('z', 'f4', ('z',))
    var[:] = Z

    # var = ds.createVariable('bathymetry','f4',('y','x'))
    # var[:,:] = L0_bathy_on_L1
    var = ds.createVariable('wet_grid_C', 'f4', ('z', 'y', 'x'))
    var[:, :, :] = wet_grid_C
    var = ds.createVariable('wet_grid_W', 'f4', ('z', 'y', 'x'))
    var[:, :, :] = wet_grid_W
    var = ds.createVariable('wet_grid_S', 'f4', ('z', 'y', 'x'))
    var[:, :, :] = wet_grid_S

    ds.close()


   

