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
    drF = ds.variables['drF'][:]
    RC = np.cumsum(drF)
    ds.close()
    return(XC, YC, RC)

def stack_files_to_nc(output_dir, output_file, dv_folder, N, RC, Nr, iter_numbers, seconds_per_iter, iteration_step):

    # iterations = np.array(iteration_subset)
    # time = np.zeros((len(iteration_subset),))

    var_names = ['THETA','SALT','UVEL','VVEL']
    output_grids = []


    max_timesteps_per_file = 24*30

    for var_name in var_names:
        n_iters_counted = 0
        output_array = np.zeros((N,len(iter_numbers)*max_timesteps_per_file,Nr))
        print(np.shape(output_array))
        if var_name == var_names[0]:
            iteration_array = np.zeros((len(iter_numbers) * max_timesteps_per_file,))

        for iter_number in iter_numbers:
            file_path = os.path.join(dv_folder, 'CTD_mask_'+ var_name + '.' + '{:010d}'.format(iter_number) + '.bin')
            iter_grid = np.fromfile(file_path, '>f4')
            n_timesteps = int(np.size(iter_grid)/(N*Nr))
            iter_grid = np.reshape(iter_grid,(n_timesteps,Nr,N))
            for n in range(N):
                for timestep in range(n_timesteps):
                    output_array[n, n_iters_counted+timestep, :] = iter_grid[timestep,:,n]
            if var_name == var_names[0]:
                first_iter_number = iter_number-1+iteration_step/2
                file_iter_numbers = np.arange(first_iter_number,first_iter_number+n_timesteps*iteration_step,iteration_step)
                iteration_array[n_iters_counted:n_iters_counted+n_timesteps] = file_iter_numbers

            n_iters_counted += n_timesteps

        print(np.shape(output_array))

        if var_name == var_names[0]:
            nonzero_indices = iteration_array != 0
            iteration_array = iteration_array[nonzero_indices]
        output_array = output_array[:, nonzero_indices, :]

        # print(np.shape(output_array))
        #
        # output_array = output_array[:,:9360,:]
        # if var_name == var_names[0]:
        #     iteration_array = iteration_array[:9360]

        print(np.shape(output_array))

        output_grids.append(output_array)

    #     time[counter] = iter_number*seconds_per_iter
    #     counter+=1

    time_array = iteration_array*seconds_per_iter

    plt.plot(time_array)
    plt.show()

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
            var = grp.createVariable(var_names[vn],'f4',('time','depth'))
            var[:, :] = output_grids[vn][p,:,:]

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
    dv_folder = os.path.join(config_dir, 'L1', L1_model_name, 'run_pleiades', 'dv')
    iter_numbers = []
    for file_name in os.listdir(dv_folder):
        if file_name[:3] == 'CTD' and 'THETA' in file_name:
            iter_number = int(file_name.split('.')[-2])
            iter_numbers.append(iter_number)

    iter_numbers = sorted(iter_numbers)

    print(iter_numbers)

    N = 17

    XC, YC, RC = read_grid_geometry_from_nc(config_dir,L1_model_name)
    Nr = len(RC)
    averaging_period = 3600
    seconds_per_iter = 300
    iteration_step = averaging_period/seconds_per_iter


    # unique_year_months = []
    # year_months = []
    # for iter_number in iter_numbers:
    #     year_month = iter_number_to_year_month(iter_number, seconds_per_iter)
    #     year_months.append(year_month)
    #     if year_month not in unique_year_months:
    #         unique_year_months.append(year_month)
    #
    # print('    - Found data for the following year-months: '+str(unique_year_months))

    if 'results' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name)):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'results'))
    if 'CTD' not in os.listdir(os.path.join(config_dir, 'L1',L1_model_name,'results')):
        os.mkdir(os.path.join(config_dir, 'L1', L1_model_name, 'results', 'CTD'))
    output_dir = os.path.join(config_dir, 'L1', L1_model_name, 'results', 'CTD')

    # for year_month in unique_year_months:
    #     print('        - Creating the file for '+str(year_month))
    #     # get the iter bounds
    #     start_index = year_months.index(year_month)
    #     end_index = len(year_months) - 1 - year_months[::-1].index(year_month)
    #     iteration_subset = iter_numbers[start_index:end_index+1]
    #     min_iter = np.min(np.array(iteration_subset))
    #     max_iter = np.max(np.array(iteration_subset))
    #     output_file = 'L1_'+subset+'.' + year_month+'.' +str(int(min_iter))+'_'+str(int(max_iter))+'.nc'

    output_file = 'CTD_profiles_2002.nc'
    stack_files_to_nc(output_dir,output_file, dv_folder, N, RC, Nr, iter_numbers, seconds_per_iter, iteration_step)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    stack_data_to_nc(config_dir)
