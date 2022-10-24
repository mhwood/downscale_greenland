

import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import ast

########################################################################################################################

def read_grid_geometry(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    drF = ds.variables['drF'][:]
    dxC = ds.variables['dxC'][:, :]
    dyC = ds.variables['dyC'][:, :]
    HFacC = ds.variables['HFacC'][:, :, :]
    HFacS = ds.variables['HFacS'][:, :, :]
    HFacW = ds.variables['HFacW'][:, :, :]
    ds.close()
    return(dxC,dyC,drF,HFacC, HFacS,HFacW)

def calculate_year_time(year):

    if year % 4 == 0:
        time_steps = 366 * 24
    else:
        time_steps = 365 * 24

    time = year + np.arange(time_steps)/float(time_steps)

    # for t in time[:100]:
    #     print(t)

    return(time)

def calculate_fluxes_from_annual_file(config_dir, model_level, model_name, boundary_names, year,
                                      dxC, dyC, drF, HFacC, HFacS, HFacW):

    if year%4==0:
        time_steps = 366*24
    else:
        time_steps = 365*24

    heat_capacity = 3974 # J/kg/C
    density = 1024

    volume_timeseries = []
    heat_timeseries = []
    salt_timeseries = []

    for boundary in boundary_names:

        print('        - Adding data on the '+boundary+' boundary')

        if boundary=='north':
            width = dxC[-1,:-1]
            HFac = HFacS[:,-2,:]
        if boundary=='south':
            width = dxC[0,:-1]
            HFac = HFacS[:, 1, :]
        if boundary=='east':
            width = dyC[:-1,-1]
            HFac = HFacW[:, :, -2]
        if boundary=='west':
            width = dyC[:-1,0]
            HFac = HFacW[:, :, 0]

        if boundary in ['north','south']:
            vel_file_path = os.path.join(config_dir,model_level,model_name,'input','obcs','L1_BC_'+boundary+'_VVEL_'+str(year))
        if boundary in ['west','east']:
            vel_file_path = os.path.join(config_dir,model_level,model_name,'input','obcs','L1_BC_'+boundary+'_UVEL_'+str(year))
        theta_file_path = os.path.join(config_dir,model_level,model_name,'input','obcs','L1_BC_'+boundary+'_THETA_'+str(year))
        salt_file_path = os.path.join(config_dir, model_level, model_name, 'input', 'obcs','L1_BC_' + boundary + '_SALT_' + str(year))

        vel_grid = np.fromfile(vel_file_path,'>f4')
        theta_grid = np.fromfile(theta_file_path, '>f4')
        salt_grid = np.fromfile(salt_file_path, '>f4')
        if boundary in ['north','south']:
            vel_grid = np.reshape(vel_grid,(time_steps,len(drF),np.shape(HFacC)[2]))
            theta_grid = np.reshape(theta_grid, (time_steps, len(drF), np.shape(HFacC)[2]))
            salt_grid = np.reshape(salt_grid, (time_steps, len(drF), np.shape(HFacC)[2]))
        else:
            vel_grid = np.reshape(vel_grid, (time_steps, len(drF), np.shape(HFacC)[1]))
            theta_grid = np.reshape(theta_grid, (time_steps, len(drF), np.shape(HFacC)[1]))
            salt_grid = np.reshape(salt_grid, (time_steps, len(drF), np.shape(HFacC)[1]))

        volume = np.zeros((time_steps,))
        heat = np.zeros((time_steps,))
        salt = np.zeros((time_steps,))

        for time_step in range(time_steps):
            if int(time_step%1000)==0:
                print('          - Calculating timesteps '+str(time_step)+' to '+str(np.min([time_step+1000,time_steps])))
            timestep_flux = 0
            timestep_heatflux = 0
            timestep_saltflux = 0
            for k in range(len(drF)):
                timestep_flux += np.sum(width * vel_grid[time_step,k,:] * drF[k] * HFac[k,:])
                timestep_heatflux += np.sum(width * vel_grid[time_step, k, :] * drF[k] * HFac[k, :] * heat_capacity * density * theta_grid[time_step,k,:])
                timestep_saltflux += np.sum(width * vel_grid[time_step, k, :] * drF[k] * HFac[k, :]  * density * salt_grid[time_step, k,:])
            volume[time_step] = timestep_flux
            heat[time_step] = timestep_heatflux
            salt[time_step] = timestep_saltflux

        if boundary in ['north','east']:
            volume *= -1
            heat *= -1
            salt *= -1

        volume *= 1e-6
        heat *= 1e-12
        salt *= 1e-9

        # plt.subplot(1,3,1)
        # plt.plot(volume)
        # plt.subplot(1, 3, 2)
        # plt.plot(heat)
        # plt.subplot(1, 3, 3)
        # plt.plot(salt)
        # plt.show()

        volume_timeseries.append(volume)
        heat_timeseries.append(heat)
        salt_timeseries.append(salt)

    return(volume_timeseries, heat_timeseries, salt_timeseries)

def store_flux_timeseries_as_nc(config_dir, model_level, model_name, total_time,
                                boundary_names, total_volume_timeseries, total_heat_timeseries, total_salt_timeseries):

    output_file = os.path.join(config_dir,model_level,model_name,'results','analysis',model_name+'_boundary_fluxes.nc')

    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('time',len(total_time))

    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = total_time

    for bn in range(len(boundary_names)):
        boundary_name = boundary_names[bn]
        grp = ds.createGroup(boundary_name)
        fvar = grp.createVariable('volume_flux','f4',('time',))
        fvar[:] = total_volume_timeseries[bn]
        fvar.units = 'Sv'
        hvar = grp.createVariable('heat_flux', 'f4', ('time',))
        hvar[:] = total_heat_timeseries[bn]
        hvar.units = 'TW'
        svar = grp.createVariable('salt_flux', 'f4', ('time',))
        svar[:] = total_salt_timeseries[bn]
        svar.units = 'Mt_per_second'

    ds.note = 'north and east fluxes are multiplied by -1 to represent flux into the domain'

    ds.close()



########################################################################################################################

def calculate_boundary_fluxes(config_dir):

    model_level = 'L1_grid'
    model_name = 'L1_CE_Greenland'
    boundary_names = ['north','west','east','south']

    print('Calculating the boundary fluxes')

    # get the grid geometry
    dxC, dyC, drF, HFacC, HFacS, HFacW = read_grid_geometry(config_dir, model_name)

    # count the number of timesteps which will be in the output timeseries
    years = np.arange(1992,1994).tolist()
    total_timesteps = 0
    for year in years:
        if year%4==0:
            total_timesteps += 366*24
        else:
            total_timesteps += 365*24

    # make bins where all of the annual timeseries will go
    total_time = np.zeros((total_timesteps,))
    total_volume_timeseries = []
    total_heat_timeseries = []
    total_salt_timeseries = []
    for boundary in boundary_names:
        total_volume_timeseries.append(np.zeros((total_timesteps,)))
        total_heat_timeseries.append(np.zeros((total_timesteps,)))
        total_salt_timeseries.append(np.zeros((total_timesteps,)))

    # read in the fluxes for each year
    points_counted = 0
    for year in years:
        print('    - Working on year ' + str(year))

        year_time = calculate_year_time(year)
        total_time[points_counted:points_counted+len(year_time)] = year_time

        volume_timeseries, heat_timeseries, salt_timeseries = \
            calculate_fluxes_from_annual_file(config_dir, model_level, model_name, boundary_names, year,
                                              dxC, dyC, drF, HFacC, HFacS, HFacW)

        for bn in range(len(boundary_names)):
            total_volume_timeseries[bn][points_counted:points_counted + len(year_time)] = volume_timeseries[bn]
            total_heat_timeseries[bn][points_counted:points_counted + len(year_time)] = heat_timeseries[bn]
            total_salt_timeseries[bn][points_counted:points_counted + len(year_time)] = salt_timeseries[bn]

        points_counted += len(year_time)

    # plt.plot(total_time,total_volume_timeseries[0])
    # plt.show()

    store_flux_timeseries_as_nc(config_dir, model_level, model_name, total_time,
                                boundary_names, total_volume_timeseries, total_heat_timeseries, total_salt_timeseries)












if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    calculate_boundary_fluxes(config_dir)