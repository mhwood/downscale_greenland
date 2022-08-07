
import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from PIL import Image, ImageDraw
import cmocean.cm as cm
import argparse


def read_geometry_from_grid_nc(config_dir,level_name,model_name):

    grid_path = os.path.join(config_dir,level_name,model_name,'input',model_name+'_grid.nc')

    ds = nc4.Dataset(grid_path)
    Depth = ds.variables['Depth'][:, :]
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    ds.close()

    bottom = np.column_stack([XC[0, :], YC[0, :]])
    top = np.column_stack([XC[-1, :], YC[-1, :]])
    left = np.column_stack([XC[:, 0], YC[:, 0]])
    right = np.column_stack([XC[:, -1], YC[:, -1]])
    polygon = np.vstack([bottom, right, np.flipud(top), np.flipud(left)])

    return(XC, YC, Depth, polygon)

def generate_global_plot(output_path,
                  center_lon, rotation_lon, center_lat, rotation_lat,
                  L1_polygon):

    plot_mode = 'light'

    with_land = True

    fig = plt.figure(figsize=(12,12))
    if plot_mode=='dark':
        plt.style.use("dark_background")

    m = Basemap(projection='ortho', resolution=None,
                lat_0=center_lat + rotation_lat, lon_0=center_lon + rotation_lon)
    if with_land:
        m.bluemarble(scale=0.5)

    # calval_lon, calval_lat = m(center_lon, center_lat)
    # m.plot(calval_lon,calval_lat,'kD',markersize=4)
    # m.plot(calval_lon,calval_lat,'rD',markersize=3)

    # for face in [1,2,3,4,5]:
    #     xc_grid = L0_xc_grid_faces[face-1]
    #     yc_grid = L0_yc_grid_faces[face-1]
    #     var_grid = L0_var_grid_faces[face-1]
    #     var_grid[var_grid>0] = 0
    #
    #     if with_land:
    #         var_grid = np.ma.masked_where(var_grid==0, var_grid)
    #
    #     if 'XY_corrected_face_'+str(face)+'.nc' in os.listdir(os.path.join(config_dir, 'plots','Basemap Plot Reference')):
    #         ds = nc4.Dataset(os.path.join(config_dir, 'plots', 'Basemap Plot Reference', 'XY_corrected_face_' + str(face) + '.nc'))
    #         X = ds.variables['X'][:,:]
    #         Y = ds.variables['Y'][:,:]
    #         ds.close()
    #     else:
    #         X, Y = m(xc_grid, yc_grid)
    #         if np.any(np.abs(X) > 1e10):
    #             rows = np.arange(np.shape(X)[0])
    #             cols = np.arange(np.shape(Y)[1])
    #             Cols, Rows = np.meshgrid(cols, rows)
    #             Cols = Cols.ravel()
    #             Rows = Rows.ravel()
    #             X_ravel = X.ravel()
    #             Y_ravel = Y.ravel()
    #             Cols = Cols[np.abs(X_ravel) < 1e10]
    #             Rows = Rows[np.abs(X_ravel) < 1e10]
    #             X_ravel = X_ravel[np.abs(X_ravel) < 1e10]
    #             Y_ravel = Y_ravel[np.abs(Y_ravel) < 1e10]
    #             nan_rows, nan_cols = np.where(np.abs(X) > 1e10)
    #             for i in range(len(nan_rows)):
    #                 if i%100==0:
    #                     print(i,len(nan_rows),round(100*i/len(nan_rows),2))
    #                 row = nan_rows[i]
    #                 col = nan_cols[i]
    #                 closest_index = np.argmin((Cols - col) ** 2 + (Rows - row) ** 2)
    #                 X[row, col] = X_ravel[closest_index]
    #                 Y[row, col] = Y_ravel[closest_index]
    #         ds = nc4.Dataset(os.path.join(config_dir, 'plots','Basemap Plot Reference','XY_corrected_face_'+str(face)+'.nc'),'w')
    #         ds.createDimension('rows',np.shape(X)[0])
    #         ds.createDimension('cols', np.shape(X)[1])
    #         Xvar = ds.createVariable('X','f4',('rows','cols'))
    #         Yvar = ds.createVariable('Y', 'f4', ('rows', 'cols'))
    #         Xvar[:] = X
    #         Yvar[:] = Y
    #         ds.close()
    #     m.pcolormesh(X, Y, var_grid, vmin=-5000, vmax = 5000, cmap = cm.topo)

    polygon_lon, polygon_lat = m(L1_polygon[:, 0], L1_polygon[:, 1])
    m.plot(polygon_lon, polygon_lat, 'r-', linewidth=2)

    axicon = fig.add_axes([0.22, 0.9, 0.5, 0.1])
    plt.text(0, 0, 'Global Source Model (1/3$^{\circ}$)', fontsize=20)
    axicon.axis('off')
    axicon.set_xticks([])
    axicon.set_yticks([])

    plt.savefig(output_path)
    plt.close(fig)

def create_globe_domain_plot(config_dir, L1_model_name):

    L1_XC, L1_YC, L1_Depth, L1_polygon = read_geometry_from_grid_nc(config_dir, 'L1', L1_model_name)

    ##################################
    # some metadata

    center_lon = np.mean(L1_XC)
    center_lat = np.mean(L1_YC)

    rotation_lon = 0
    rotation_lat = -20

    # ################################################################################################
    # # read in the global grid
    #
    # file_path = os.path.join(config_dir,'L0_540','input','bathy_llc540')
    # L0_var_grid_faces = read_L0_var(file_path,'EtaN')
    #
    ##################################
    # make the plot

    output_path = os.path.join(config_dir,'plots','bathymetry','L0_Global_bathymetry.png')

    generate_global_plot(output_path,
                  center_lon, rotation_lon, center_lat, rotation_lat,
                  L1_polygon)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-L1", "--L1_model_name", action="store",
                        help="The name of the L1 model.", dest="L1_model_name",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    L1_model_name = args.L1_model_name

    create_globe_domain_plot(config_dir, L1_model_name)