
import os
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import cmocean.cm as cm


def read_L1_CE_Greenland_boundary_condition(config_dir, L1_model_name, boundary, field_name, sNx, sNy, Nr):

    years = np.arange(1992,2022).tolist()

    N = 0
    for year in years:
        if year%4==0:
            N+=366*24
        else:
            N+=365*24

    if boundary=='west' or boundary=='east':
        hovmoller_grid = np.zeros((2*sNy,N))
    if boundary=='north' or boundary=='south':
        hovmoller_grid = np.zeros((N,3*sNx))
    dec_yrs = np.zeros((N,))

    points_counted = 0
    for year in years:
        print('    - Reading data in year '+str(year))
        if year%4==0:
            n_timesteps = 366 * 24
        else:
            n_timesteps = 365 * 24
        dec_yrs[points_counted:points_counted+n_timesteps] = year + np.arange(n_timesteps)/n_timesteps

        if boundary=='south':
            boundary_file = os.path.join(config_dir,'L1',L1_model_name,'input','obcs',
                                         'L1_BC_'+boundary+'_'+field_name.upper()+'.bin')
            boundary_grid = np.fromfile(boundary_file,'>f4')
            boundary_grid = boundary_grid.reshape((n_timesteps, Nr, sNx*4))
            boundary_grid = boundary_grid[:,:,:sNx*3]

        if boundary=='west':
            boundary_file = os.path.join(config_dir,'L1',L1_model_name,'input','obcs',
                                         'L1_BC_'+boundary+'_'+field_name.upper()+ '_'+str(year))
            boundary_grid = np.fromfile(boundary_file,'>f4')
            boundary_grid = boundary_grid.reshape((n_timesteps, Nr, sNx*4))
            boundary_grid = boundary_grid[:,:,:sNx*2]

        if boundary == 'east':
            boundary_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                         'L1_BC_east_' + field_name.upper() + '_'+str(year))
            boundary_grid = np.fromfile(boundary_file, '>f4')
            boundary_grid = boundary_grid.reshape((n_timesteps, Nr, 4*sNx))
            boundary_grid = boundary_grid[:, :, :sNx]

            boundary_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                         'L1_BC_south_' + field_name.upper() +  '_'+str(year))
            boundary_grid_2 = np.fromfile(boundary_file, '>f4')
            boundary_grid_2 = boundary_grid_2.reshape((n_timesteps, Nr, 4 * sNx))
            boundary_grid_2 = boundary_grid_2[:,:,sNx*3:]

            boundary_grid = np.concatenate([boundary_grid,boundary_grid_2],axis = 2)

        if boundary == 'north':
            boundary_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                         'L1_BC_east_' + field_name.upper() + '_'+str(year))
            boundary_grid = np.fromfile(boundary_file, '>f4')
            boundary_grid = boundary_grid.reshape((n_timesteps, Nr, 4*sNx))
            boundary_grid = boundary_grid[:,:,sNx:]

            hovmoller_grid[points_counted:points_counted+n_timesteps,:] = boundary_grid[:,0,:]

        points_counted+=n_timesteps

    if boundary=='north':
        hovmoller_grid = np.fliplr(hovmoller_grid)

    if boundary=='north' or boundary=='south':
        hovmoller_grid = hovmoller_grid[::24,:]
    dec_yrs = dec_yrs[::24]

    return(dec_yrs,hovmoller_grid)


def plot_vertical_hovmoller_series(config_dir, L1_model_name, boundary, var_name, dec_yrs, hovmoller_grid):

    if var_name == 'AREA':
        vmin = -0.05
        vmax = 1.05
        cmap = cm.ice

    if var_name == 'HEFF':
        vmin = -0.05
        vmax = 3
        cmap = cm.ice

    if var_name == 'HSNOW':
        vmin = -0.05
        vmax = 2
        cmap = cm.ice

    if var_name == 'THETA':
        cmap = cm.thermal
        vmin = min_val
        vmax = max_val

    if var_name == 'SALT':
        cmap = cm.haline
        vmin = min_val
        vmax = max_val

    if var_name == 'UVEL' or var_name=='VVEL':
        cmap = cm.balance
        vmin = -0.15
        vmax = 0.15

    if var_name == 'UICE' or var_name=='VICE':
        cmap = cm.balance
        vmin = -0.5
        vmax = 0.5

    points_along_boundary = np.arange(540)

    N = 365*10+2

    fig = plt.figure(figsize=(12,8))
    plt.style.use('dark_background')

    plt.subplot(1,3,1)
    C = plt.pcolormesh(points_along_boundary,dec_yrs[:N],hovmoller_grid[:N,:],vmin=vmin,vmax=vmax,cmap=cmap, shading='nearest')
    # plt.colorbar(C)
    plt.title('1992-2002')
    plt.yticks(np.arange(1992,2002))
    plt.grid(linestyle='--',alpha=0.4,linewidth=1)

    plt.subplot(1, 3, 2)
    C = plt.pcolormesh(points_along_boundary, dec_yrs[N:2*N], hovmoller_grid[N:2*N,:], vmin=vmin, vmax=vmax, cmap=cmap, shading='nearest')
    # plt.colorbar(C)
    plt.xlabel('Points along boundary')
    plt.title('2002-2012')
    plt.yticks(np.arange(2002, 2012))
    plt.grid(linestyle='--', alpha=0.4, linewidth=1)

    plt.subplot(1, 3, 3)
    C = plt.pcolormesh(points_along_boundary, dec_yrs[2*N:], hovmoller_grid[2*N:,:], vmin=vmin, vmax=vmax, cmap=cmap, shading='nearest')
    # plt.colorbar(C)
    plt.title('2012-2022')
    plt.yticks(np.arange(2012,2022))
    plt.grid(linestyle='--', alpha=0.4, linewidth=1)

    output_file = os.path.join(config_dir, 'L1', L1_model_name, 'plots', 'init_files',
                               L1_model_name + '_BC_' + var_name + '_'+boundary + '_timeseries.png')
    plt.savefig(output_file)
    plt.close(fig)

def plot_L1_CE_Greenland_BC_timeseries(config_dir, L1_model_name, sNx, sNy):

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils', 'plot_creation','init_files'))

    boundary = 'north'

    var_names = ['HEFF']#, 'SALT', 'UVEL', 'VVEL','AREA','HEFF','HSNOW','UICE','VICE']

    for var_name in var_names:
        if var_name in ['THETA','SALT','UVEL','VVEL']:
            Nr = 50
        else:
            Nr = 1

        dec_yrs, hovmoller_grid = read_L1_CE_Greenland_boundary_condition(config_dir, L1_model_name, boundary, var_name,
                                                                 sNx, sNy, Nr)

        if boundary in ['north','south']:
            plot_vertical_hovmoller_series(config_dir, L1_model_name, boundary, var_name, dec_yrs, hovmoller_grid)


        # if var_name in ['THETA', 'SALT', 'UVEL', 'VVEL']:
        #     create_3D_BC_plot(config_dir, L1_model_name, var_name, timestep,
        #                    north_boundary, south_boundary, east_boundary, west_boundary)
        # else:
        #     create_2D_BC_plot(config_dir, L1_model_name, var_name, timestep,
        #                       north_boundary, south_boundary, east_boundary, west_boundary)




    



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_L1_CE_Greenland_init_fields(config_dir)
   

