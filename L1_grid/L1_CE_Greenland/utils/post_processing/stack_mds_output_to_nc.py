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
    return(year_month,year,month)

def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    AngleCS = ds.variables['AngleCS'][:, :]
    AngleSN = ds.variables['AngleSN'][:, :]
    ds.close()
    return(XC, YC, AngleCS, AngleSN)

def stack_files_to_monthly_nc(output_dir, output_file, subset_folder, subset,
                              iteration_subset, seconds_per_iter, XC, YC):

    if subset == 'EtaN_day_snap':
        var_names = ['EtaN']
    elif subset == 'SI_daily_snap':
        var_names = ['SIarea','SIheff','SIhsnow','SIuice','SIvice']
    elif subset == 'TS_surf_daily_snap':
        var_names = ['Theta','Theta_1','Salt','Salt_1']
    elif subset == 'TS_AW_daily_snap':
        var_names = ['Theta','Theta_1','Salt','Salt_1']
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
        grid = np.reshape(grid,(len(var_names),np.shape(XC)[0],np.shape(XC)[1]))
        output_array[:,counter,:,:] = grid

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

    ds.close()

def stack_files_to_annual_nc(output_dir, output_file, subset_folder, subset,
                              iteration_subset, seconds_per_iter, XC, YC, Nr):

    if subset == 'EtaN_mon_mean':
        var_names = ['EtaN']
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

    output_array = np.zeros((len(var_names)*Nr,len(iteration_subset),np.shape(XC)[0],np.shape(XC)[1]))
    # first col is for the different variables
    # second col is for the iter number
    # third and fourth are the sizes

    counter = 0
    for iter_number in iteration_subset:
        file_path = os.path.join(subset_folder,subset+'.'+'{:010d}'.format(iter_number)+'.data')

        grid = np.fromfile(file_path, '>f4')

        if subset in ['EtaN_mon_mean']:
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

    ds.close()



########################################################################################################################

def stack_data_to_nc(config_dir, subset):

    L1_model_name = 'L1_CE_Greenland'

    if subset=='All':
        subsets = ['EtaN_day_snap','SI_daily_snap','TS_surf_daily_snap','TS_AW_daily_snap',
                   'EtaN_mon_mean','state_3D_mon_mean','state_3D_mon_snap','vel_3D_mon_mean','vel_3D_mon_snap']
    else:
        subsets = [subset]

    for subset in subsets:

        print('    - Creating datasets for the '+subset+' subset')

        # get all of the iteration numbers
        subset_folder = os.path.join(config_dir, 'L1_grid', L1_model_name, 'run', 'diags', subset)
        iter_numbers = []
        for file_name in os.listdir(subset_folder):
            if file_name[-4:] == 'data':
                iter_number = int(file_name.split('.')[-2])
                iter_numbers.append(iter_number)
        iter_numbers = sorted(iter_numbers)

        # get the grid geometry
        XC, YC, AngleCS, AngleSN = read_grid_geometry_from_nc(config_dir,L1_model_name)

        Nr = 50

        seconds_per_iter = 300
        unique_year_months = []
        unique_years = []
        year_months = []
        years = []
        for iter_number in iter_numbers:
            year_month, year, month = iter_number_to_year_month(iter_number, seconds_per_iter)
            year_months.append(year_month)
            years.append(str(year))
            if year_month not in unique_year_months:
                unique_year_months.append(year_month)
            if str(year) not in unique_years:
                unique_years.append(str(year))

        if subset in ['EtaN_day_snap','SI_daily_snap','TS_surf_daily_snap','TS_AW_daily_snap']:
            print('        - Found data for the following year-months: '+str(unique_year_months))
        elif subset in ['EtaN_mon_mean','state_3D_mon_mean','state_3D_mon_snap','vel_3D_mon_mean','vel_3D_mon_snap']:
            print('        - Found data for the following years: ' + str(unique_years))
        else:
            raise ValueError('Identify whether this data should be stored monthly or yearly')

        if 'results' not in os.listdir(os.path.join(config_dir,'L1_grid',L1_model_name)):
            os.mkdir(os.path.join(config_dir,'L1_grid',L1_model_name,'results'))
        if subset not in os.listdir(os.path.join(config_dir, 'L1_grid',L1_model_name,'results')):
            os.mkdir(os.path.join(config_dir, 'L1_grid', L1_model_name, 'results', subset))
        output_dir = os.path.join(config_dir, 'L1_grid', L1_model_name, 'results', subset)

        if subset in ['EtaN_day_snap', 'SI_daily_snap', 'TS_surf_daily_snap', 'TS_AW_daily_snap']:
            for year_month in unique_year_months:
                # get the iter bounds
                start_index = year_months.index(year_month)
                end_index = len(year_months) - 1 - year_months[::-1].index(year_month)
                iteration_subset = iter_numbers[start_index:end_index + 1]
                min_iter = np.min(np.array(iteration_subset))
                max_iter = np.max(np.array(iteration_subset))
                output_file = 'L1_' + subset + '.' + year_month + '.' + str(int(min_iter)) + '_' + str(int(max_iter)) + '.nc'
                if output_file not in os.listdir(output_dir):
                    print('            - Creating the file for '+str(year_month))
                    stack_files_to_monthly_nc(output_dir,output_file,subset_folder,subset,
                                      iteration_subset,seconds_per_iter, XC, YC)
                else:
                    print('            - '+output_file+' already created')

        if subset in ['EtaN_mon_mean','state_3D_mon_mean','state_3D_mon_snap','vel_3D_mon_mean','vel_3D_mon_snap']:
            for year in unique_years:
                # get the iter bounds
                start_index = years.index(year)
                end_index = len(years) - 1 - years[::-1].index(year)
                iteration_subset = iter_numbers[start_index:end_index + 1]
                min_iter = np.min(np.array(iteration_subset))
                max_iter = np.max(np.array(iteration_subset))
                output_file = 'L1_' + subset + '.' + year + '.' + str(int(min_iter)) + '_' + str(int(max_iter)) + '.nc'
                if output_file not in os.listdir(output_dir):
                    print('            - Creating the file for '+str(year))
                    stack_files_to_annual_nc(output_dir,output_file,subset_folder,subset,
                                      iteration_subset,seconds_per_iter, XC, YC, Nr)
                else:
                    print('            - '+output_file+' already created')



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
