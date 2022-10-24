import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from datetime import datetime, timedelta
import argparse
import ast
import sys

def iter_number_to_year(iter_number,seconds_per_iter):
    total_seconds = iter_number*seconds_per_iter
    date = datetime(1992,1,1) + timedelta(seconds=total_seconds)
    year = date.year
    #month=date.month
    #year_month = str(year)#+'{:02d}'.format(month)
    return(year)

def date_to_iter_number(date):
    seconds_per_iter = 1200
    seconds_per_iter = 300
    total_seconds = (date-datetime(1992,1,1)).total_seconds()
    iter_number = int(total_seconds/seconds_per_iter)
    return(iter_number)

def iter_number_to_date(iter_number):
    seconds_per_iter = 300
    total_seconds = iter_number*seconds_per_iter
    date = datetime(1992,1,1) + timedelta(seconds=total_seconds)
    return(date)


def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    drF = ds.variables['drF'][:]
    RC = np.cumsum(drF)
    ds.close()
    return(XC, YC, RC)

def stack_files_to_nc(output_dir, output_file, dv_folder, year, N, RC, Nr, iter_numbers, seconds_per_iter, iteration_step):

    # iterations = np.array(iteration_subset)
    # time = np.zeros((len(iteration_subset),))

    var_names = ['THETA','SALT','UVEL','VVEL','AREA','HEFF','HSNOW','UICE','VICE']
    # var_names = ['THETA']
    output_grids = []

    print_debug = False

    if year%4==0:
        n_timesteps = 366 * 24
    else:
        n_timesteps = 365 * 24
    if year==1992:
        n_timesteps = 365 * 24 # 1992 doesnt have day 1

    if year==1992:
        min_year_iter = date_to_iter_number(datetime(year, 1, 2, 0, 30))
    else:
        min_year_iter = date_to_iter_number(datetime(year, 1, 1, 0, 30))
    max_year_iter = date_to_iter_number(datetime(year,12,31,23,30))
    iteration_array = np.arange(min_year_iter,max_year_iter+1,iteration_step)
    time_array = year + np.arange(n_timesteps)/n_timesteps

    if print_debug:
        print('        - This file has iterations '+str(min_year_iter)+' through '+str(max_year_iter) +
              ' (size: '+str(np.size(iteration_array))+' == ' + str(n_timesteps)+')')

    for var_name in var_names:
        if print_debug:
            print('        - Adding data for '+var_name)

        if var_name in ['THETA','SALT','UVEL','VVEL']:
            output_grid = np.zeros((Nr, n_timesteps, N))
        else:
            output_grid = np.zeros((1, n_timesteps, N))

        for iter_number in iter_numbers:
            file_path = os.path.join(dv_folder, 'CTD_mask_'+ var_name + '.' + '{:010d}'.format(iter_number) + '.bin')
            iter_grid = np.fromfile(file_path, '>f4')
            if var_name in ['THETA','SALT','UVEL','VVEL']:
                n_file_timesteps = int(np.size(iter_grid)/(N*Nr))
                iter_grid = np.reshape(iter_grid,(n_file_timesteps,Nr,N))
            else:
                n_file_timesteps = int(np.size(iter_grid) / N)
                iter_grid = np.reshape(iter_grid, (n_file_timesteps, 1, N))
            file_iters = np.arange(iter_number-1,(iter_number-1)+n_file_timesteps*iteration_step,iteration_step)+iteration_step/2
            if not np.min(file_iters)>max_year_iter and not np.max(file_iters)<min_year_iter:
                if print_debug:
                    print('            - Adding data from iter number '+str(iter_number))
                    print('                - This file has iters '+str(int(np.min(file_iters)))+' through '+str(int(np.max(file_iters))))

                if np.min(file_iters)<np.min(iteration_array):
                    min_iter_index = 0
                    min_file_index = np.where(file_iters==iteration_array[0])[0][0]
                else:
                    min_iter_index = np.where(iteration_array==int(np.min(file_iters)))[0][0]
                    min_file_index = 0

                if np.max(file_iters)>np.max(iteration_array):
                    max_iter_index = len(iteration_array)
                    max_file_index = np.where(file_iters==iteration_array[-1])[0][0]+1
                else:
                    max_iter_index = np.where(iteration_array==int(np.max(file_iters)))[0][0]+1
                    max_file_index = len(file_iters)

                if print_debug:
                    print('                - Storing indices '+str(min_file_index)+' to '+str(max_file_index)+
                          ' from file at indices ' + str(min_iter_index) + ' to ' + str(max_iter_index) + ' in the output array')

                for n in range(N):
                    if var_name in ['THETA','SALT','UVEL','VVEL']:
                        for k in range(Nr):
                            output_grid[k, min_iter_index:max_iter_index, n] = iter_grid[min_file_index:max_file_index,k,n]
                    else:
                        output_grid[0, min_iter_index:max_iter_index, n] = iter_grid[min_file_index:max_file_index, 0,n]

        # plt.plot(time_array,output_grid[10,:,1])
        # plt.show()

        output_grids.append(output_grid)

    ds = nc4.Dataset(os.path.join(output_dir,output_file),'w')
    # ds.createDimension('iterations',len(iterations))
    # ds.createDimension('rows',np.shape(XC)[0])
    ds.createDimension('depth', Nr)
    ds.createDimension('time', np.shape(output_grids[0])[1])

    dvar = ds.createVariable('depth', 'f4', ('depth',))
    dvar[:] = RC

    tvar = ds.createVariable('time','f4',('time',))
    ivar = ds.createVariable('iterations', 'f4', ('time',))
    tvar[:] = time_array
    ivar[:] = iteration_array

    for p in range(N):
        grp = ds.createGroup('point_'+'{:02d}'.format(p+1))
        for vn in range(len(var_names)):
            if var_names[vn] in ['THETA','SALT','UVEL','VVEL']:
                var = grp.createVariable(var_names[vn],'f4',('depth','time'))
                var[:, :] = output_grids[vn][:,:,p]
            else:
                var = grp.createVariable(var_names[vn], 'f4', ('time',))
                var[:] = output_grids[vn][0, :, p]

    ds.close()




########################################################################################################################

def stack_data_to_nc(config_dir):

    L1_model_name = 'L1_CE_Greenland'
    # sNx = 180
    # sNy = 180
    # faces = [1,3]
    # face_size_dict = {1: (sNy, 3*sNx), 3: (3*sNy, sNx)}
    #
    # sys.path.insert(1, os.path.join(config_dir, 'L1', L1_model_name, 'utils'))
    # import L1_CE_Greenland_functions as Lf

    # create a list of files
    dv_folder = os.path.join(config_dir, 'L1_grid', L1_model_name, 'run', 'dv','CTD')
    iter_numbers = []
    for file_name in os.listdir(dv_folder):
        if file_name[:3] == 'CTD' and 'THETA' in file_name:
            iter_number = int(file_name.split('.')[-2])
            iter_numbers.append(iter_number)

    iter_numbers = sorted(iter_numbers)

    N = 17

    XC, YC, RC = read_grid_geometry_from_nc(config_dir,L1_model_name)
    Nr = len(RC)
    averaging_period = 3600
    seconds_per_iter = 300
    iteration_step = averaging_period/seconds_per_iter

    unique_years = []
    for iter_number in iter_numbers:
        year = iter_number_to_year(iter_number, seconds_per_iter)
        if year not in unique_years:
            unique_years.append(year)

    print('    - Found data for the following years: '+str(unique_years))

    if 'results' not in os.listdir(os.path.join(config_dir,'L1_grid',L1_model_name)):
        os.mkdir(os.path.join(config_dir,'L1_grid',L1_model_name,'results'))
    if 'CTD' not in os.listdir(os.path.join(config_dir, 'L1_grid',L1_model_name,'results')):
        os.mkdir(os.path.join(config_dir, 'L1_grid', L1_model_name, 'results', 'CTD'))
    output_dir = os.path.join(config_dir, 'L1_grid', L1_model_name, 'results', 'CTD')

    for year in unique_years:
        print('    - Creating the file for '+str(year))
        output_file = 'L1_CE_CTD_profiles_' + str(year)+'.nc'
        stack_files_to_nc(output_dir,output_file, dv_folder, year, N, RC, Nr, iter_numbers, seconds_per_iter, iteration_step)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    stack_data_to_nc(config_dir)
