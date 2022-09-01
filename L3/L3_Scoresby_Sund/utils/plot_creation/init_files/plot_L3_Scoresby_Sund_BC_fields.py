
import os
import argparse
import sys
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import cmocean.cm as cm


def read_L3_Scoresby_Sund_boundary_condition(config_dir, L3_model_name, boundary, field_name,
                                            Nr, n_rows_L3, n_cols_L3,  n_timesteps):

    if boundary=='south' or boundary=='north':
        boundary_file = os.path.join(config_dir,'L3',L3_model_name,'input','obcs',
                                     'L3_BC_'+boundary+'_'+field_name.upper()+'.bin')
        boundary_grid = np.fromfile(boundary_file,'>f4')
        boundary_grid = boundary_grid.reshape((n_timesteps, Nr, n_cols_L3))

    if boundary == 'east' or boundary == 'west':
        boundary_file = os.path.join(config_dir, 'L3', L3_model_name, 'input', 'obcs',
                                     'L3_BC_' + boundary + '_' + field_name.upper() + '.bin')
        boundary_grid = np.fromfile(boundary_file, '>f4')
        boundary_grid = boundary_grid.reshape((n_timesteps, Nr, n_rows_L3))


    return(boundary_grid)

def create_2D_BC_plot(config_dir, L3_model_name, var_name, timestep,
                   north_boundary, south_boundary, east_boundary):

    # min_val = np.min([np.min(north_boundary[north_boundary!=0]),
    #                   np.min(south_boundary[south_boundary!=0]),
    #                   np.min(east_boundary[east_boundary!=0]),
    #                   np.min(west_boundary[west_boundary!=0])])
    #
    # max_val = np.max([np.max(north_boundary[north_boundary != 0]),
    #                   np.max(south_boundary[south_boundary != 0]),
    #                   np.max(east_boundary[east_boundary != 0]),
    #                   np.max(west_boundary[west_boundary != 0])])

    print('Plotting '+var_name)

    vmin = -0.5
    vmax = 0.5

    if var_name=='AREA':
        vmin=-0.05
        vmax=1.05
    if var_name=='HEFF':
        vmin = -0.05
        vmax = 3
    if var_name=='HSNOW':
        vmin = -0.05
        vmax = 2

    fig = plt.figure(figsize=(8, 10))
    plt.style.use('dark_background')

    plt.subplot(3, 1, 1)
    C = plt.plot(north_boundary.ravel())
    plt.grid(linestyle='--',alpha=0.5)
    plt.ylabel('North')
    plt.title(var_name+' Boundary Conditions at timestep = '+str(timestep))
    plt.gca().set_ylim([vmin,vmax])

    plt.subplot(3, 1, 2)
    C = plt.plot(south_boundary.ravel())
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('South')
    plt.gca().set_ylim([vmin, vmax])

    plt.subplot(3, 1, 3)
    C = plt.plot(east_boundary.ravel())
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('East')
    plt.xlabel('Points Along Boundary')
    plt.gca().set_ylim([vmin, vmax])

    output_file = os.path.join(config_dir,'L3',L3_model_name,'plots','init_files',L3_model_name+'_BC_'+var_name+'_'+str(timestep)+'.png')
    plt.savefig(output_file)
    plt.close(fig)


def create_3D_BC_plot(config_dir, L3_model_name, var_name, timestep,
                   north_boundary, south_boundary, east_boundary):

    min_val = np.min([np.min(north_boundary[north_boundary!=0]),
                      np.min(south_boundary[south_boundary!=0]),
                      np.min(east_boundary[east_boundary!=0])])

    max_val = np.max([np.max(north_boundary[north_boundary != 0]),
                      np.max(south_boundary[south_boundary != 0]),
                      np.max(east_boundary[east_boundary != 0])])

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

    fig = plt.figure(figsize=(8, 10))
    plt.style.use('dark_background')

    if var_name in ['THETA','SALT','UVEL','VVEL']:
        plt.subplot(3, 1, 1)
        x = np.arange(np.shape(north_boundary)[1])
        d = np.arange(np.shape(north_boundary)[0])
        X, D = np.meshgrid(x,d)
        C = plt.pcolormesh(X, D, north_boundary, cmap=cmap, vmin=vmin, vmax=vmax, shading='nearest')
        plt.colorbar(C)
        plt.gca().invert_yaxis()
        plt.ylabel('North')
        plt.title(var_name+' Boundary Conditions at timestep = '+str(timestep))

        plt.subplot(3, 1, 2)
        x = np.arange(np.shape(south_boundary)[1])
        d = np.arange(np.shape(south_boundary)[0])
        X, D = np.meshgrid(x, d)
        C = plt.pcolormesh(X, D, south_boundary, cmap=cmap, vmin=vmin, vmax=vmax, shading='nearest')
        plt.colorbar(C)
        plt.gca().invert_yaxis()
        plt.ylabel('South')

        plt.subplot(3, 1, 3)
        x = np.arange(np.shape(east_boundary)[1])
        d = np.arange(np.shape(east_boundary)[0])
        X, D = np.meshgrid(x, d)
        C = plt.pcolormesh(X, D, east_boundary, cmap=cmap, vmin=vmin, vmax=vmax, shading='nearest')
        plt.colorbar(C)
        plt.gca().invert_yaxis()
        plt.ylabel('East')
        plt.xlabel('Points Along Boundary')

    output_file = os.path.join(config_dir,'L3',L3_model_name,'plots','init_files',L3_model_name+'_BC_'+var_name+'_'+str(timestep)+'.png')
    plt.savefig(output_file)
    plt.close(fig)

def plot_L3_Scoresby_Sund_BCs(config_dir, L3_model_name, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'L3', 'utils', 'plot_creation','init_files'))

    n_timesteps = 2*24 + 2
    timestep = 31

    grid_file = os.path.join(config_dir, 'nc_grids', L3_model_name + '_grid.nc')
    ds = nc4.Dataset(grid_file)
    XC = ds.variables['XC'][:, :]
    drF = ds.variables['drF'][:]
    n_rows_L3 = np.shape(XC)[0]
    n_cols_L3 = np.shape(XC)[1]
    ds.close()

    var_names = ['THETA', 'SALT', 'UVEL', 'VVEL','AREA','HEFF','HSNOW','UICE','VICE']

    for var_name in var_names:
        if var_name in ['THETA','SALT','UVEL','VVEL']:
            Nr = len(drF)
        else:
            Nr = 1

        north_boundary = read_L3_Scoresby_Sund_boundary_condition(config_dir, L3_model_name, 'north', var_name,
                                                                 Nr, n_rows_L3, n_cols_L3, n_timesteps)
        north_boundary = north_boundary[timestep, :, :]

        south_boundary = read_L3_Scoresby_Sund_boundary_condition(config_dir, L3_model_name, 'south', var_name,
                                                                 Nr, n_rows_L3, n_cols_L3, n_timesteps)
        south_boundary = south_boundary[timestep, :, :]

        east_boundary = read_L3_Scoresby_Sund_boundary_condition(config_dir, L3_model_name, 'east', var_name,
                                                                 Nr, n_rows_L3, n_cols_L3, n_timesteps)
        east_boundary = east_boundary[timestep, :, :]

        if var_name in ['THETA', 'SALT', 'UVEL', 'VVEL']:
            create_3D_BC_plot(config_dir, L3_model_name, var_name, timestep,
                           north_boundary, south_boundary, east_boundary)
        else:
            create_2D_BC_plot(config_dir, L3_model_name, var_name, timestep,
                              north_boundary, south_boundary, east_boundary)




    



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L3, L3, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_L3_Scoresby_Sund_init_fields(config_dir)
   

