

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
    Depth = ds.variables['Depth'][:, :]
    HFacC = ds.variables['HFacC'][:, :, :]
    HFacS = ds.variables['HFacS'][:, :, :]
    HFacW = ds.variables['HFacW'][:, :, :]
    ds.close()
    return(dxC,dyC,drF,HFacC, HFacS,HFacW, Depth)

def calculate_fluxes_from_annual_file(config_dir, model_level, model_name, boundary_names, year,
                                      dxC, dyC, drF, HFacC, HFacS, HFacW, Depth):

    if year%4==0:
        time_steps = 366*24
    else:
        time_steps = 365*24

    heat_capacity = 3974 # J/kg/C
    density = 1024

    area_grids = []
    dist_vectors = []
    depth_vectors = []

    volume_mean_grids = []
    heat_mean_grids = []
    salt_mean_grids = []

    volume_sum_grids = []
    heat_sum_grids = []
    salt_sum_grids = []

    for boundary in boundary_names:

        print('        - Adding data on the '+boundary+' boundary')

        if boundary=='north':
            width = dxC[-1,:-1]
            HFac = HFacS[:,-2,:]
            depth = Depth[-1,:]
        if boundary=='south':
            width = dxC[1,:-1]
            HFac = HFacS[:, 1, :]
            depth = Depth[1,:]
        if boundary=='east':
            width = dyC[:-1,-1]
            HFac = HFacW[:, :, -2]
            depth = Depth[:,-1]
        if boundary=='west':
            width = dyC[:-1,0]
            HFac = HFacW[:, :, 0]
            depth = Depth[:,0]
        volume_sum_grid = np.zeros_like(HFac)
        heat_sum_grid = np.zeros_like(HFac)
        salt_sum_grid = np.zeros_like(HFac)
        area_grid = np.zeros_like(HFac)

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

        points_counted = 0
        for time_step in range(time_steps):
            if int(time_step%1000)==0:
                print('          - Calculating timesteps '+str(time_step)+' to '+str(np.min([time_step+1000,time_steps])))
            for k in range(len(drF)):
                volume_sum_grid[k, :] += width * vel_grid[time_step,k,:] * drF[k] * HFac[k,:]
                heat_sum_grid[k, :] += width * vel_grid[time_step, k, :] * drF[k] * HFac[k, :] * heat_capacity * density * theta_grid[time_step,k,:] * 1e-6
                salt_sum_grid[k, :] += width * vel_grid[time_step, k, :] * drF[k] * HFac[k, :]  * density * salt_grid[time_step, k,:] * 1e-3
                if time_step==0:
                    area_grid[k,:] = width * drF[k] * HFac[k,:]
            points_counted += 1

        # plt.imshow(area_grid)
        # plt.title(boundary)
        # plt.show()

        if boundary in ['north','east']:
            volume_sum_grid *= -1
            heat_sum_grid *= -1
            salt_sum_grid *= -1

        volume_mean_grid = volume_sum_grid / points_counted
        heat_mean_grid = heat_sum_grid / points_counted
        salt_mean_grid = salt_sum_grid / points_counted

        # plt.subplot(1,3,1)
        # plt.plot(volume)
        # plt.subplot(1, 3, 2)
        # plt.plot(heat)
        # plt.subplot(1, 3, 3)
        # plt.plot(salt)
        # plt.show()
        volume_mean_grids.append(volume_mean_grid)
        volume_sum_grids.append(volume_sum_grid)

        heat_mean_grids.append(heat_mean_grid)
        heat_sum_grids.append(heat_sum_grid)

        salt_mean_grids.append(salt_mean_grid)
        salt_sum_grids.append(salt_sum_grid)

        area_grids.append(area_grid)
        dist_vectors.append(width)
        depth_vectors.append(depth)

    return(area_grids, dist_vectors, depth_vectors, points_counted,
           volume_mean_grids, heat_mean_grids, salt_mean_grids,
           volume_sum_grids, heat_sum_grids, salt_sum_grids)

def store_flux_grids_as_nc(config_dir, model_level, model_name, year, boundary_names,
                           drF, area_grids, dist_vectors, depth_vectors, n_points,
                           volume_mean_grids, heat_mean_grids, salt_mean_grids,
                           volume_sum_grids, heat_sum_grids, salt_sum_grids):

    output_file = os.path.join(config_dir,model_level,model_name,'results','analysis',model_name+'_mean_boundary_flux_'+str(year)+'.nc')

    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('depth', np.shape(volume_mean_grids[0])[0])

    ds.n_points = n_points

    dvar = ds.createVariable('drF','f4',('depth',))
    dvar[:] = drF

    for bn in range(len(boundary_names)):
        boundary_name = boundary_names[bn]
        grp = ds.createGroup(boundary_name)

        grp.createDimension('boundary_points',np.shape(volume_mean_grids[bn])[1])

        fvar = grp.createVariable('volume_flux_sum','f4',('depth','boundary_points'))
        fvar[:,:] = volume_sum_grids[bn]
        fvar = grp.createVariable('volume_flux_mean', 'f4', ('depth', 'boundary_points'))
        fvar[:, :] = volume_mean_grids[bn]

        hvar = grp.createVariable('heat_flux_sum', 'f4', ('depth', 'boundary_points'))
        hvar[:, :] = heat_sum_grids[bn]
        hvar = grp.createVariable('heat_flux_mean', 'f4', ('depth', 'boundary_points'))
        hvar[:, :] = heat_mean_grids[bn]

        svar = grp.createVariable('salt_flux_sum', 'f4', ('depth', 'boundary_points'))
        svar[:, :] = salt_sum_grids[bn]
        svar = grp.createVariable('salt_flux_mean', 'f4', ('depth', 'boundary_points'))
        svar[:, :] = salt_mean_grids[bn]

        avar = grp.createVariable('area', 'f4', ('depth', 'boundary_points'))
        avar[:, :] = area_grids[bn]

        if boundary_name in ['north','south']:
            wvar = grp.createVariable('dxC', 'f4', ('boundary_points',))
        if boundary_name in ['west','east']:
            wvar = grp.createVariable('dyC', 'f4', ('boundary_points',))
        wvar[:] = dist_vectors[bn]

        depvar = grp.createVariable('depth','f4',('boundary_points',))
        depvar[:] = depth_vectors[bn]

    ds.note = 'north and east fluxes are multiplied by -1 to represent flux into the domain'

    ds.close()



########################################################################################################################

def calculate_annual_mean_fluxes(config_dir):

    model_level = 'L1_grid'
    model_name = 'L1_CE_Greenland'
    boundary_names = ['north','west','east','south']

    print('Calculating the mean annual boundary fluxes')

    # get the grid geometry
    dxC, dyC, drF, HFacC, HFacS, HFacW, Depth = read_grid_geometry(config_dir, model_name)

    years = np.arange(1993,2021).tolist()


    # read in the fluxes for each year
    for year in years:
        print('    - Working on year ' + str(year))

        area_grids, dist_vectors, depth_vectors, n_points,\
        volume_mean_grids, heat_mean_grids, salt_mean_grids,\
        volume_sum_grids, heat_sum_grids, salt_sum_grids = \
            calculate_fluxes_from_annual_file(config_dir, model_level, model_name, boundary_names, year,
                                              dxC, dyC, drF, HFacC, HFacS, HFacW, Depth)

        store_flux_grids_as_nc(config_dir, model_level, model_name, year, boundary_names,
                               drF, area_grids, dist_vectors, depth_vectors, n_points,
                               volume_mean_grids, heat_mean_grids, salt_mean_grids,
                               volume_sum_grids, heat_sum_grids, salt_sum_grids)












if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    calculate_annual_mean_fluxes(config_dir)