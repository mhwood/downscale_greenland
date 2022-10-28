
import os
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import cmocean.cm as cm


def read_L1_W_Greenland_boundary_condition(config_dir, L1_model_name, boundary, field_name, sNx, sNy, Nr, n_timesteps):

    if boundary=='south':

        boundary_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                     'L1_BC_east_' + field_name.upper() + '.bin')
        boundary_grid = np.fromfile(boundary_file, '>f4')
        boundary_grid = boundary_grid.reshape((n_timesteps, Nr, 4 * sNx))
        boundary_grid = boundary_grid[:, :, sNx:]

    if boundary=='west':
        # boundary_file = os.path.join(config_dir,'L1',L1_model_name,'input','obcs',
        #                              'L1_BC_north_'+field_name.upper()+'.bin')
        # boundary_grid = np.fromfile(boundary_file,'>f4')
        # boundary_grid = boundary_grid.reshape((n_timesteps, Nr, sNx*6))
        # boundary_grid = boundary_grid[:,:,sNx*3:]

        boundary_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                     'L1_BC_east_' + field_name.upper() + '.bin')
        boundary_grid_1 = np.fromfile(boundary_file, '>f4')
        boundary_grid_1 = boundary_grid_1.reshape((n_timesteps, Nr, 4 * sNx))
        boundary_grid_1 = boundary_grid_1[:, :, :sNx]

        boundary_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                     'L1_BC_south_' + field_name.upper() + '.bin')
        boundary_grid_2 = np.fromfile(boundary_file, '>f4')
        boundary_grid_2 = boundary_grid_2.reshape((n_timesteps, Nr, sNx * 6))
        boundary_grid_2 = boundary_grid_2[:, :, sNx * 3:]

        boundary_grid = np.concatenate([boundary_grid_1, boundary_grid_2], axis=2)

    if boundary == 'east':
        boundary_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                     'L1_BC_north_' + field_name.upper() + '.bin')
        boundary_grid = np.fromfile(boundary_file, '>f4')
        boundary_grid = boundary_grid.reshape((n_timesteps, Nr, 6*sNx))
        boundary_grid = boundary_grid[:, :, 3*sNx:]

    if boundary == 'north':

        boundary_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'obcs',
                                     'L1_BC_south_' + field_name.upper() + '.bin')
        boundary_grid = np.fromfile(boundary_file, '>f4')
        boundary_grid = boundary_grid.reshape((n_timesteps, Nr, 6 * sNx))
        boundary_grid = boundary_grid[:, :, :3*sNx]

    return(boundary_grid)

def create_2D_BC_plot(config_dir, L1_model_name, var_name, timestep,
                   north_boundary, south_boundary, east_boundary, west_boundary):

    # min_val = np.min([np.min(north_boundary[north_boundary!=0]),
    #                   np.min(south_boundary[south_boundary!=0]),
    #                   np.min(east_boundary[east_boundary!=0]),
    #                   np.min(west_boundary[west_boundary!=0])])
    #
    # max_val = np.max([np.max(north_boundary[north_boundary != 0]),
    #                   np.max(south_boundary[south_boundary != 0]),
    #                   np.max(east_boundary[east_boundary != 0]),
    #                   np.max(west_boundary[west_boundary != 0])])

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

    plt.subplot(4, 1, 1)
    C = plt.plot(north_boundary.ravel())
    plt.grid(linestyle='--',alpha=0.5)
    plt.ylabel('North')
    plt.title(var_name+' Boundary Conditions at timestep = '+str(timestep))
    plt.gca().set_ylim([vmin,vmax])

    plt.subplot(4, 1, 2)
    C = plt.plot(south_boundary.ravel())
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('South')
    plt.gca().set_ylim([vmin, vmax])

    plt.subplot(4, 1, 3)
    C = plt.plot(west_boundary.ravel())
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('West')
    plt.gca().set_ylim([vmin, vmax])

    plt.subplot(4, 1, 4)
    C = plt.plot(east_boundary.ravel())
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('East')
    plt.xlabel('Points Along Boundary')
    plt.gca().set_ylim([vmin, vmax])

    output_file = os.path.join(config_dir,'L1',L1_model_name,'plots','init_files',L1_model_name+'_BC_'+var_name+'_'+str(timestep)+'.png')
    plt.savefig(output_file)
    plt.close(fig)


def create_3D_BC_plot(config_dir, L1_model_name, var_name, timestep,
                   north_boundary, south_boundary, east_boundary, west_boundary):

    min_val = np.min([np.min(north_boundary[north_boundary!=0]),
                      np.min(south_boundary[south_boundary!=0]),
                      np.min(east_boundary[east_boundary!=0]),
                      np.min(west_boundary[west_boundary!=0])])

    max_val = np.max([np.max(north_boundary[north_boundary != 0]),
                      np.max(south_boundary[south_boundary != 0]),
                      np.max(east_boundary[east_boundary != 0]),
                      np.max(west_boundary[west_boundary != 0])])

    # min_val = -3
    # max_val = 8

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
        plt.subplot(4, 1, 1)
        x = np.arange(np.shape(north_boundary)[1])
        d = np.arange(np.shape(north_boundary)[0])
        X, D = np.meshgrid(x,d)
        C = plt.pcolormesh(X, D, north_boundary, cmap=cmap, vmin=vmin, vmax=vmax, shading='nearest')
        plt.colorbar(C)
        plt.gca().invert_yaxis()
        plt.ylabel('North')
        plt.title(var_name+' Boundary Conditions at timestep = '+str(timestep))

        plt.subplot(4, 1, 2)
        x = np.arange(np.shape(south_boundary)[1])
        d = np.arange(np.shape(south_boundary)[0])
        X, D = np.meshgrid(x, d)
        C = plt.pcolormesh(X, D, south_boundary, cmap=cmap, vmin=vmin, vmax=vmax, shading='nearest')
        plt.colorbar(C)
        plt.gca().invert_yaxis()
        plt.ylabel('South')

        plt.subplot(4, 1, 3)
        x = np.arange(np.shape(west_boundary)[1])
        d = np.arange(np.shape(west_boundary)[0])
        X, D = np.meshgrid(x, d)
        C = plt.pcolormesh(X, D, west_boundary, cmap=cmap, vmin=vmin, vmax=vmax, shading='nearest')
        plt.colorbar(C)
        plt.gca().invert_yaxis()
        plt.ylabel('West')

        plt.subplot(4, 1, 4)
        x = np.arange(np.shape(east_boundary)[1])
        d = np.arange(np.shape(east_boundary)[0])
        X, D = np.meshgrid(x, d)
        C = plt.pcolormesh(X, D, east_boundary, cmap=cmap, vmin=vmin, vmax=vmax, shading='nearest')
        plt.colorbar(C)
        plt.gca().invert_yaxis()
        plt.ylabel('East')
        plt.xlabel('Points Along Boundary')

    output_file = os.path.join(config_dir,'L1',L1_model_name,'plots','init_files',L1_model_name+'_BC_'+var_name+'_'+str(timestep)+'.png')
    plt.savefig(output_file)
    plt.close(fig)

def plot_L1_W_Greenland_BCs(config_dir, L1_model_name, sNx, sNy):

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils', 'plot_creation','init_files'))

    n_timesteps = 2*24 + 2
    timestep = 4

    var_names = ['THETA', 'SALT', 'UVEL', 'VVEL','AREA','HEFF','HSNOW','UICE','VICE']
    var_names = ['THETA']

    for var_name in var_names:
        if var_name in ['THETA','SALT','UVEL','VVEL']:
            Nr = 50
        else:
            Nr = 1

        print('    - Reading the north boundary')
        north_boundary = read_L1_W_Greenland_boundary_condition(config_dir, L1_model_name, 'north', var_name,
                                                                 sNx, sNy, Nr, n_timesteps)
        north_boundary = north_boundary[timestep, :, :]

        print('    - Reading the south boundary')
        south_boundary = read_L1_W_Greenland_boundary_condition(config_dir, L1_model_name, 'south', var_name,
                                                                 sNx, sNy, Nr, n_timesteps)
        south_boundary = south_boundary[timestep, :, :]

        print('    - Reading the west boundary')
        west_boundary = read_L1_W_Greenland_boundary_condition(config_dir, L1_model_name, 'west', var_name,
                                                                 sNx, sNy, Nr, n_timesteps)
        west_boundary = west_boundary[timestep, :, :]

        print('    - Reading the east boundary')
        east_boundary = read_L1_W_Greenland_boundary_condition(config_dir, L1_model_name, 'east', var_name,
                                                                 sNx, sNy, Nr, n_timesteps)
        east_boundary = east_boundary[timestep, :, :]

        if var_name in ['THETA', 'SALT', 'UVEL', 'VVEL']:
            create_3D_BC_plot(config_dir, L1_model_name, var_name, timestep,
                           north_boundary, south_boundary, east_boundary, west_boundary)
        else:
            create_2D_BC_plot(config_dir, L1_model_name, var_name, timestep,
                              north_boundary, south_boundary, east_boundary, west_boundary)




    



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_L1_W_Greenland_init_fields(config_dir)
   

