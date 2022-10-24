
import os
import argparse
import sys
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import cmocean.cm as cm


def read_L1_CE_Greenland_boundary_condition(config_dir, L1_model_name, boundary, field_name,  n_rows_L1, n_cols_L1, Nr, depth_level):

    years = np.arange(1992,2022).tolist()

    N = 0
    for year in years:
        if year%4==0:
            N+=366*24
        else:
            N+=365*24

    dec_yrs = np.zeros((N,))
    points_counted = 0
    for year in years:
        if year%4==0:
            n_timesteps = 366 * 24
        else:
            n_timesteps = 365 * 24
        dec_yrs[points_counted:points_counted+n_timesteps] = year + np.arange(n_timesteps)/n_timesteps
        points_counted += n_timesteps

    # plt.plot(dec_yrs)
    # plt.show()

    if boundary=='west' or boundary=='east':
        hovmoller_grid = np.zeros((n_rows_L1,N))
    if boundary=='north' or boundary=='south':
        hovmoller_grid = np.zeros((N,n_cols_L1))

    # years = np.arange(1992, 1994).tolist()

    points_counted = 0
    for year in years:
        print('    - Reading data in year '+str(year))
        if year%4==0:
            n_timesteps = 366 * 24
        else:
            n_timesteps = 365 * 24

        if boundary=='south' or boundary=='north':
            boundary_file = os.path.join(config_dir,'L1_grid',L1_model_name,'input','obcs',
                                         'L1_BC_'+boundary+'_'+field_name.upper()+'_'+str(year))
            boundary_grid = np.fromfile(boundary_file,'>f4')
            boundary_grid = boundary_grid.reshape((n_timesteps, Nr, n_cols_L1))

            hovmoller_grid[points_counted:points_counted+n_timesteps,:] = boundary_grid[:,depth_level,:]

        points_counted+=n_timesteps

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
        vmin = -2
        vmax = 5

    if var_name == 'SALT':
        cmap = cm.haline
        vmin = min_val
        vmax = max_val

    if var_name == 'UVEL' or var_name=='VVEL':
        cmap = cm.balance
        vmin = -0.05
        vmax = 0.05

    if var_name == 'UICE' or var_name=='VICE':
        cmap = cm.balance
        vmin = -0.5
        vmax = 0.5

    points_along_boundary = np.arange(540)

    hovmoller_grid = np.ma.masked_where(hovmoller_grid==0,hovmoller_grid)

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

    output_file = os.path.join(config_dir, 'L1_grid', L1_model_name, 'plots', 'init_files',
                               L1_model_name + '_BC_' + var_name + '_'+boundary + '_timeseries.png')
    plt.savefig(output_file)
    plt.close(fig)

def plot_L1_CE_Greenland_BC_timeseries(config_dir, L1_model_name, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'L1_grid', 'utils', 'plot_creation','init_files'))

    boundary = 'north'
    depth_level = 19

    grid_file = os.path.join(config_dir, 'nc_grids', L1_model_name + '_grid.nc')
    ds = nc4.Dataset(grid_file)
    XC = ds.variables['XC'][:, :]
    drF = ds.variables['drF'][:]
    n_rows_L1 = np.shape(XC)[0]
    n_cols_L1 = np.shape(XC)[1]
    ds.close()

    var_names = ['THETA','VVEL']#, 'SALT', 'UVEL', 'VVEL','AREA','HEFF','HSNOW','UICE','VICE']

    for var_name in var_names:
        if var_name in ['THETA','SALT','UVEL','VVEL']:
            Nr = 50
        else:
            Nr = 1
            depth_level = 0

        print('Working on plot for '+var_name+' on the '+boundary+' boundary')

        dec_yrs, hovmoller_grid = read_L1_CE_Greenland_boundary_condition(config_dir, L1_model_name, boundary, var_name,
                                                                          n_rows_L1, n_cols_L1, Nr, depth_level)

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
   

