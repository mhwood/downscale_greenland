
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys


def create_mitgrid(config_dir, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L2', 'utils','init_file_creation'))
    import create_L2_mitgrid as cm

    model_name = 'L2_CE_Greenland'

    if print_level>=1:
        print('    - Generating the grid in polar coordinates')

    # make the point grid in polar coordinates
    # these were the original L3 coords
    min_x = 543736
    max_x = 897457
    max_y = -1873700
    min_y = -2092319

    # expand L3 coords to L2
    x_range = max_x - min_x
    y_range = max_y - min_y

    # min_x -= x_range
    max_x += x_range
    min_y -= y_range
    max_y += y_range

    resolution = 3000

    min_x = int(min_x/resolution)*resolution
    max_x = (int(max_x / resolution)+1) * resolution
    min_y = int(min_y / resolution) * resolution
    max_y = (int(max_y / resolution)+1) * resolution

    extra_top_rows = 10
    extra_bottom_rows = 9
    extra_right_cols = 2

    x = np.arange(min_x,max_x+2*resolution+extra_right_cols*resolution,resolution)
    y = np.arange(min_y-extra_bottom_rows*resolution,
                  max_y+2*resolution+extra_top_rows*resolution,
                  resolution)
    XC, YC = np.meshgrid(x,y)

    x = np.arange(min_x, max_x + 2*resolution + extra_right_cols * resolution, resolution)
    y = np.arange(min_y - extra_bottom_rows * resolution,
                  max_y + 2*resolution + extra_top_rows * resolution,
                  resolution)
    XG, YG = np.meshgrid(x-resolution/2,y-resolution/2)

    if print_level >= 2:
        print('    - The C grid has '+str(np.shape(XC)[0])+' rows and '+str(np.shape(XC)[1])+' cols')
        print('    - The G grid has ' + str(np.shape(XG)[0]) + ' rows and ' + str(np.shape(XG)[1]) + ' cols')

    # reproject the grid to lon, lat
    if print_level >= 2:
        print('        - Reprojecting the grid to lat/lon')
    transformer = Transformer.from_crs('EPSG:' + str(3413), 'EPSG:' + str(4326))

    Lat_C, Lon_C = transformer.transform(XC.ravel(), YC.ravel())
    Lat_C = np.reshape(Lat_C,np.shape(XC))
    Lon_C = np.reshape(Lon_C, np.shape(XC))

    Lat_G, Lon_G = transformer.transform(XG.ravel(), YG.ravel())
    Lat_G = np.reshape(Lat_G, np.shape(XG))
    Lon_G = np.reshape(Lon_G, np.shape(XG))

    # fig = plt.figure(figsize=(10,5))
    # plt.subplot(1,2,1)
    # C = plt.pcolormesh(Lon_G,Lat_G,Lon_G)
    # plt.colorbar(C)
    # plt.title('Longitude')
    # plt.subplot(1, 2, 2)
    # C = plt.pcolormesh(Lon_G, Lat_G, Lat_G)
    # plt.colorbar(C)
    # plt.title('Latitude')
    # plt.savefig(os.path.join(config_dir,'L3','L3_Scoresby_Sund','plots','L3_Scoresby_Sund_Lat_Lon.png'))
    # plt.close(fig)
    # plt.show()

    # pass to general function to generate mitgrid
    cm.create_L2_mitgrid_file(config_dir, model_name, Lat_C, Lon_C, Lat_G, Lon_G, print_level)








if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_mitgrid(config_dir)
   

