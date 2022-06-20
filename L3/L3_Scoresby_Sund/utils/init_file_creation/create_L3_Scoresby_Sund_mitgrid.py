
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys


def create_mitgrid(config_dir):

    sys.path.insert(1, os.path.join(config_dir,'L3', 'utils','init_file_creation'))
    import create_L3_mitgrid as cm

    model_name = 'L3_Scoresby_Sund'

    # pass this a grid on the subdomain in EPSG 3413 coordinates
    print('Creating the mitgrid file for the '+model_name+' model')

    print('    - Generating the grid in polar coordinates')

    # make the point grid in polar coordinates
    min_x = 543736
    max_x = 897457
    max_y = -1873700
    min_y = -2092319

    resolution = 2000

    min_x = int(min_x/resolution)*resolution
    max_x = (int(max_x / resolution)+1) * resolution
    min_y = int(min_y / resolution) * resolution
    max_y = (int(max_y / resolution)+1) * resolution

    x = np.arange(min_x,max_x+resolution,resolution)
    y = np.arange(min_y,max_y+resolution,resolution)
    XC, YC = np.meshgrid(x,y)

    x = np.arange(min_x, max_x + 2 * resolution, resolution)
    y = np.arange(min_y, max_y + 2 * resolution, resolution)
    XG, YG = np.meshgrid(x-resolution/2,y-resolution/2)

    print('    - The C grid has '+str(np.shape(XC)[0])+' rows and '+str(np.shape(XC)[1])+' cols')
    print('    - The G grid has ' + str(np.shape(XG)[0]) + ' rows and ' + str(np.shape(XG)[1]) + ' cols')

    # reproject the grid to lon, lat
    print('    - Reprojecting the grid to lat/lon')
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
    cm.create_L3_mitgrid_file(config_dir,model_name,Lat_C, Lon_C,Lat_G, Lon_G)








if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_mitgrid(config_dir)
   

