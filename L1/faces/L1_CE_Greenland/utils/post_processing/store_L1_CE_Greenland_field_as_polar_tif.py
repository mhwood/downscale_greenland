import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from datetime import datetime, timedelta
from scipy.interpolate import griddata
import argparse
import ast
import sys
from osgeo import gdal
from osgeo import osr
from pyproj import Transformer

def read_field_from_results_nc(file_path, field_name, timestep):

    ds= nc4.Dataset(file_path)
    field = ds.variables[field_name][:, :, :]
    ds.close()
    field = field[timestep,:,:]

    return(field)

def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    Depth = ds.variables['Depth'][:,:]
    ds.close()
    return(XC, YC, Depth)

def reproject_points(points,inputCRS,outputCRS,x_column=0,y_column=1):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(points[:, x_column], points[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
        run_test = False
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon = np.copy(points)
    output_polygon[:, x_column] = x2
    output_polygon[:, y_column] = y2
    return output_polygon

def stack_files_to_nc(Lf, output_dir, output_file, subset_folder, subset,
                      iteration_subset, seconds_per_iter, sNx, sNy, XC, YC, faces, face_size_dict):

    if subset == 'surfDiag':
        var_names = ['EtaN']
    elif subset == 'seaiceDiag':
        var_names = ['Area','Heff','Hsnow','Uice','Vice']
    elif subset == 'dynDiag':
        var_names = ['Uvel','Vvel','Theta','Salt']
    elif subset == 'awDiag':
        var_names = ['Theta','Salt']
    else:
        raise ValueError('Variables names not defined for this subset')

    iterations = np.array(iteration_subset)
    time = np.zeros((len(iteration_subset),))

    output_array = np.zeros((len(var_names),len(iteration_subset),np.shape(XC)[0],np.shape(XC)[1]))

    counter = 0
    for iter_number in iteration_subset:
        file_path = os.path.join(subset_folder,subset+'.'+'{:010d}'.format(iter_number)+'.data')

        grid_compact = np.fromfile(file_path, '>f4')

        if subset=='surfDiag':
            N = int(np.size(grid_compact) / sNx)
            grid_compact = np.reshape(grid_compact, (N, sNx))
            grid = Lf.read_compact_grid_to_stitched_grid(grid_compact, sNx, sNy, faces, face_size_dict, dim=2)

            # plt.imshow(grid, origin='lower')
            # plt.show()

            output_array[0,counter,:,:] = grid
        else:
            N = int(np.size(grid_compact) / (sNx*len(var_names)))
            grid_compact = np.reshape(grid_compact, (len(var_names), N, sNx))
            grid = Lf.read_compact_grid_to_stitched_grid(grid_compact, sNx, sNy, faces, face_size_dict, dim=3)

            # plt.imshow(grid, origin='lower')
            # plt.show()

            output_array[:, counter, :, :] = grid

        time[counter] = iter_number*seconds_per_iter
        counter+=1

    ds = nc4.Dataset(os.path.join(output_dir,output_file),'w')
    ds.createDimension('iterations',len(iterations))
    ds.createDimension('rows',np.shape(XC)[0])
    ds.createDimension('cols',np.shape(XC)[1])

    tvar = ds.createVariable('time','f4',('iterations',))
    ivar = ds.createVariable('iterations', 'f4', ('iterations',))

    for vn in range(len(var_names)):
        evar = ds.createVariable(var_names[vn], 'f4', ('iterations','rows','cols'))
        evar[:, :, :] = output_array[vn,:, :, :]

    tvar[:] = time
    ivar[:] = iterations

    xvar = ds.createVariable('longitude', 'f4', ('rows', 'cols'))
    xvar[:, :] = XC

    yvar = ds.createVariable('latitude', 'f4', ('rows', 'cols'))
    yvar[:, :] = YC

    # xvar.long_name = 'Cartesian x-coordinate'
    # xvar.standard_name = 'projection_x_coordinate'
    # xvar.units = 'degrees'
    # xvar.axis = 'X'
    # xvar.coverage_content_type = 'coordinate'
    # xvar.valid_min = np.min(XC)
    # xvar.valid_max = np.max(XC)
    # xvar.comment = 'Projected horizontal coordinates of the grid'
    #
    # yvar.long_name = 'Cartesian y-coordinate'
    # yvar.standard_name = 'projection_y_coordinate'
    # yvar.units = 'degrees'
    # yvar.axis = 'Y'
    # yvar.coverage_content_type = 'coordinate'
    # yvar.valid_min = np.min(YC)
    # yvar.valid_max = np.max(YC)
    # yvar.comment = 'Projected vertical coordinates of the grid'
    #
    # pvar = ds.createVariable('projection','S1')
    # pvar.spatial_ref = 'GEOGCS["WGS 84",' \
    #                     'DATUM["WGS_1984",' \
    #                     'SPHEROID["WGS 84",6378137,298.257223563]],' \
    #                     'PRIMEM["Greenwich",0],' \
    #                     'UNIT["degree",0.01745329251994328]]'
    # pvar.GeoTransform = str(np.min(XC)) + ' ' + str(3000) + ' 0 ' + str(
    #     np.max(YC)) + ' 0 ' + str(-3000)
    ds.close()

def write_data_to_tif(output_file, epsg, x,y, interp_grid):

    geotransform = (np.min(x), x[1]-x[0], 0, np.max(y), 0, y[0]-y[1])

    output_raster = gdal.GetDriverByName('GTiff').Create(output_file, len(x), len(y), 1,
                                                         gdal.GDT_Float32)  # Open the file
    output_raster.SetGeoTransform(geotransform)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)

    output_raster.SetProjection(srs.ExportToWkt())
    output_raster.GetRasterBand(1).WriteArray(np.flipud(interp_grid))  # Writes my array to the raster

    output_raster.FlushCache()



########################################################################################################################

def stack_data_to_nc(config_dir, subset, field_name):

    L1_model_name = 'L1_CE_Greenland'
    timestep = 17

    XC, YC, Depth = read_grid_geometry_from_nc(config_dir,L1_model_name)

    resolution = 1000

    points = np.column_stack([XC.ravel(), YC.ravel()])
    points = reproject_points(points, inputCRS=4326, outputCRS=3413)
    x = np.arange(np.min(points[:, 0]) - resolution, np.max(points[:, 0]) + 2 * resolution, resolution)
    y = np.arange(np.min(points[:, 1]) - resolution, np.max(points[:, 1]) + 2 * resolution, resolution)
    X, Y = np.meshgrid(x, y)

    file_path = os.path.join(config_dir, 'L1', L1_model_name, 'results',
                             subset, 'L1_'+subset+'.200201.1052070_1057854.nc')

    field = read_field_from_results_nc(file_path, field_name, timestep)

    interp_grid = griddata(points, Depth.ravel(), (X, Y), fill_value=0)

    # plt.imshow(interp_grid,origin='lower')
    # plt.show()

    output_file = os.path.join(config_dir,'L1',L1_model_name,'plots','tif_files','L1_Depth.tif')#'L1_'+subset+'_'+field_name+'_'+str(timestep)+'.tif')
    epsg = 3413
    write_data_to_tif(output_file, epsg, x, y, interp_grid)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-s", "--subset", action="store",
                        help="The subset to stack (e.g. surfDiag, awDiag, seaiceDiag, dynDiag).", dest="subset",
                        type=str, required=True)

    parser.add_argument("-f", "--field_name", action="store",
                        help="The field to plot (e.g. Theta).", dest="field_name",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    subset = args.subset
    field_name = args.field_name

    stack_data_to_nc(config_dir, subset, field_name)
