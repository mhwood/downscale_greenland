
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import netCDF4 as nc4
from pyproj import Transformer
import sys


def create_mitgrid(config_dir):

    sys.path.insert(1, os.path.join(config_dir,'L3', 'utils','init_file_creation'))

    model_name = 'L3_Scoresby_Sund'
    read_from_file = True

    # pass this a grid on the subdomain in EPSG 3413 coordinates
    print('Creating the stretched mitgrid file for the '+model_name+' model')

    print('    - Generating the grid in polar coordinates')

    if not read_from_file:

        # make the point grid in polar coordinates
        min_x = 543236
        max_x = 897457
        max_y = -1870200
        min_y = -2094319

        x_resolution_right = 2000
        x_resolution_left = 500
        y_resolution = 500

        min_x = int(min_x/x_resolution_right)*x_resolution_right
        # max_x = (int(max_x / x_resolution_right)+1) * x_resolution_right

        min_y = int(min_y / y_resolution) * y_resolution
        max_y = (int(max_y / y_resolution)+1) * y_resolution

        x_hat = np.arange(421)

        k = 0.05
        dx = (x_resolution_right-x_resolution_left)/(1+np.exp(-k*(x_hat-300))) + x_resolution_left
        x = min_x + np.cumsum(dx)

        # plt.subplot(1,2,1)
        # plt.title('dx')
        # plt.plot(x_hat,dx)
        # plt.subplot(1,2,2)
        # plt.plot(x_hat,x)
        # plt.title('x')
        # plt.plot([np.min(x_hat),np.max(x_hat)],[898000,898000],'k--')
        # plt.title(str(np.max(x))+', '+str(898000))
        # plt.show()

        y = np.arange(min_y,max_y+2*y_resolution,y_resolution)
        XC, YC = np.meshgrid(x,y)

        xg = x - dx/2
        yg = y - y_resolution/2
        XG, YG = np.meshgrid(xg, yg)

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

        ds = nc4.Dataset(os.path.join(config_dir, 'L3', 'L3_Scoresby_Sund',
                                      'input_for_grid', 'L3_Scoresby_Sund_Lat_Lon.nc'), 'w')

        ds.createDimension('rows',np.shape(Lon_C)[0])
        ds.createDimension('cols',np.shape(Lon_C)[1])

        xcvar = ds.createVariable('XC', 'f4', ('rows', 'cols'))
        xcvar[:, :] = XC
        ycvar = ds.createVariable('YC', 'f4', ('rows', 'cols'))
        ycvar[:, :] = YC

        xgvar = ds.createVariable('XG', 'f4', ('rows', 'cols'))
        xgvar[:, :] = XG
        ygvar = ds.createVariable('YG', 'f4', ('rows', 'cols'))
        ygvar[:, :] = YG

        xcvar = ds.createVariable('Lon_C','f4',('rows','cols'))
        xcvar[:,:] = Lon_C
        ycvar = ds.createVariable('Lat_C', 'f4', ('rows', 'cols'))
        ycvar[:, :] = Lat_C

        xgvar = ds.createVariable('Lon_G', 'f4', ('rows', 'cols'))
        xgvar[:, :] = Lon_G
        ygvar = ds.createVariable('Lat_G', 'f4', ('rows', 'cols'))
        ygvar[:, :] = Lat_G

        ds.close()

    else:

        import create_L3_mitgrid as cm

        ds = nc4.Dataset(os.path.join(config_dir, 'L3', 'L3_Scoresby_Sund',
                                      'input_for_grid', 'L3_Scoresby_Sund_Lat_Lon.nc'))

        Lon_C = ds.variables['Lon_C'][:, :]
        Lat_C = ds.variables['Lat_C'][:, :]

        Lon_G = ds.variables['Lon_G'][:, :]
        Lat_G = ds.variables['Lat_G'][:, :]

        ds.close()

        # pass to general function to generate mitgrid
        cm.create_L3_mitgrid_file(config_dir,model_name, Lat_C, Lon_C, Lat_G, Lon_G)








if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_mitgrid(config_dir)
   

