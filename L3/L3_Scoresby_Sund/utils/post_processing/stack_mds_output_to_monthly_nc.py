import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from datetime import datetime, timedelta
import argparse
import ast
import sys

def iter_number_to_year_month(iter_number,seconds_per_iter):
    total_seconds = iter_number*seconds_per_iter
    date = datetime(1992,1,1) + timedelta(seconds=total_seconds)
    year = date.year
    month=date.month
    year_month = str(year)+'{:02d}'.format(month)
    return(year_month)

def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    ds.close()
    return(XC, YC)

def stack_files_to_nc(output_dir, output_file, subset_folder, subset,
                      iteration_subset, seconds_per_iter, XC, YC):

    if subset == 'EtaN_day_snap':
        var_names = ['EtaN']
    elif subset == 'EXF_day_snap':
        var_names = ['EXFaqh','EXFatemp','EXFpreci','EXFroff','EXFlwdn','EXFswdn','EXFuwind','EXFvwind']
    elif subset == 'SI_daily_snap':
        var_names = ['SIarea','SIheff','SIhsnow','SIuice','SIvice']
    elif subset == 'TS_surf_daily_snap':
        var_names = ['Uvel','Vvel']
    elif subset == 'vel_surf_daily_snap':
        var_names = ['Uvel','Vvel']
    elif subset == 'TS_AW_daily_snap':
        var_names = ['Theta','Salt']
    elif subset == 'state_3D_mon_snap':
        var_names = ['Theta','Salt']
    elif subset == 'state_3D_mon_mean':
        var_names = ['Theta','Salt']
    elif subset == 'vel_3D_mon_snap':
        var_names = ['Uvel','Vvel']
    elif subset == 'vel_3D_mon_mean':
        var_names = ['Uvel','Vvel']
    else:
        raise ValueError('Variables names not defined for this subset')

    iterations = np.array(iteration_subset)
    time = np.zeros((len(iteration_subset),))

    output_array = np.zeros((len(var_names),len(iteration_subset),np.shape(XC)[0],np.shape(XC)[1]))
    # first col is for the different variables
    # second col is for the iter number
    # third and fourth are the sizes

    counter = 0
    for iter_number in iteration_subset:
        file_path = os.path.join(subset_folder,subset+'.'+'{:010d}'.format(iter_number)+'.data')

        grid = np.fromfile(file_path, '>f4')

        if subset in ['EtaN_day_snap','EXF_day_snap','SI_daily_snap','TS_surf_daily_snap','TS_AW_daily_snap','vel_surf_daily_snap']:
            grid = np.reshape(grid,(len(var_names),np.shape(XC)[0],np.shape(XC)[1]))
            output_array[:,counter,:,:] = grid
        else:
            grid = np.reshape(grid,(len(var_names)*Nr,np.shape(XC)[0],np.shape(XC)[1]))
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
        if var_names[vn] not in ['Theta_1','Salt_1']:
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




########################################################################################################################

def stack_data_to_nc(config_dir, subset):

    L1_model_name = 'L3_Scoresby_Sund'

    if subset=='All':
        subsets = ['EtaN_day_snap','vel_surf_daily_snap','SI_daily_snap','TS_surf_daily_snap','TS_AW_daily_snap']
    else:
        subsets = [subset]

    for subset in subsets:

        print('    - Creating datasets for the '+subset+' subset')

        # create a list of iteration bounds (in seconds)
        subset_folder = os.path.join(config_dir, 'L3', L1_model_name, 'run', 'diags', subset)
        iter_numbers = []
        for file_name in os.listdir(subset_folder):
            if file_name[-4:] == 'data':
                iter_number = int(file_name.split('.')[-2])
                iter_numbers.append(iter_number)

        XC, YC = read_grid_geometry_from_nc(config_dir,L1_model_name)

        iter_numbers = sorted(iter_numbers)

        seconds_per_iter = 60
        unique_year_months = []
        year_months = []
        for iter_number in iter_numbers:
            year_month = iter_number_to_year_month(iter_number, seconds_per_iter)
            year_months.append(year_month)
            if year_month not in unique_year_months:
                unique_year_months.append(year_month)

        print('        - Found data for the following year-months: '+str(unique_year_months))

        if 'results' not in os.listdir(os.path.join(config_dir,'L3',L1_model_name)):
            os.mkdir(os.path.join(config_dir,'L3',L1_model_name,'results'))
        if subset not in os.listdir(os.path.join(config_dir, 'L3',L1_model_name,'results')):
            os.mkdir(os.path.join(config_dir, 'L3', L1_model_name, 'results', subset))
        output_dir = os.path.join(config_dir, 'L3', L1_model_name, 'results', subset)

        for year_month in unique_year_months:
            print('            - Creating the file for '+str(year_month))
            # get the iter bounds
            start_index = year_months.index(year_month)
            end_index = len(year_months) - 1 - year_months[::-1].index(year_month)
            iteration_subset = iter_numbers[start_index:end_index+1]
            min_iter = np.min(np.array(iteration_subset))
            max_iter = np.max(np.array(iteration_subset))
            output_file = 'L1_'+subset+'.' + year_month+'.' +str(int(min_iter))+'_'+str(int(max_iter))+'.nc'

            # if output_file not in os.listdir(output_dir):
            stack_files_to_nc(output_dir,output_file,subset_folder,subset,
                              iteration_subset,seconds_per_iter, XC, YC)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-s", "--subset", action="store",
                        help="The subset to stack (e.g. surfDiag, awDiag, seaiceDiag, dynDiag).", dest="subset",
                        type=str, required=False, default='All')

    args = parser.parse_args()
    config_dir = args.config_dir
    subset = args.subset

    stack_data_to_nc(config_dir, subset)
