
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from MITgcmutils import mds
from random import randint
import cmocean.cm as cm
import sys
import argparse


def create_mean_exf_timeseries(Lf, config_dir, L1_model_name, sNx, sNy, faces, face_size_dict):

    years = np.arange(1992, 2022)
    var_names = ['ATEMP','AQH','LWDOWN','SWDOWN','PRECIP','RUNOFF','UWIND','VWIND']
    timeseries_set = []

    N = 0
    for year in years:
        if year % 4 == 0:
            N += 366 * 4
        else:
            N += 365 * 4

    input_dir = os.path.join(config_dir,'L1',L1_model_name,'input','exf')

    for var_name in var_names:

        exf_timeseries = np.zeros((N, 2))
        points_counted = 0

        for year in years:
            file_name = 'L1_exf_' + var_name + '_' + str(year)

            if year % 4 == 0:
                n_timesteps = 366 * 4
            else:
                n_timesteps = 365 * 4

            exf_timeseries[points_counted:points_counted + n_timesteps, 0] = year + np.arange(n_timesteps) / n_timesteps

            if file_name in os.listdir(input_dir):
                print('        - Reading ' + file_name)
                file_path = os.path.join(input_dir, file_name)
                grid = np.fromfile(file_path, '>f4')
                grid = np.reshape(grid, (n_timesteps, 6 * sNx, sNx))

                grid = Lf.read_compact_grid_to_stitched_grid(grid, sNx, sNy, faces, face_size_dict, dim=3)

                for g in range(np.shape(grid)[0]):
                    subset = grid[g,:,:]

                    if g==0 and var_name=='ATEMP':
                        mask = subset!=0

                    exf_timeseries[points_counted+g,1] = np.mean(subset[mask])

                del grid

            points_counted += n_timesteps

        timeseries_set.append(exf_timeseries)

    ds = nc4.Dataset(os.path.join(input_dir,L1_model_name+'_mean_exf_timeseries.nc'), 'w')
    ds.createDimension('time',N)
    tvar = ds.createVariable('time','f4',('time',))
    tvar[:] = timeseries_set[0][:,0]
    for vn in range(len(var_names)):
        var = ds.createVariable(var_names[vn],'f4',('time',))
        var[:] = timeseries_set[vn][:,1]
    ds.close()

    return(var_names, timeseries_set)

def read_mean_exf_timeseries(config_dir, L1_model_name):

    var_names = ['ATEMP', 'AQH', 'LWDOWN', 'SWDOWN', 'PRECIP', 'RUNOFF', 'UWIND', 'VWIND']
    timeseries_set = []

    ds = nc4.Dataset(os.path.join(config_dir,'L1',L1_model_name,'input','exf',L1_model_name+'_mean_exf_timeseries.nc'))
    time = ds.variables['time'][:]
    for var_name in var_names:
        var_timeseries = ds.variables[var_name][:]
        timeseries_set.append(np.column_stack([time,var_timeseries]))
    ds.close()

    return (var_names, timeseries_set)


def read_L1_geometry_from_grid(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:, :]
    ds.close()
    rows = np.shape(XC)[0]
    cols = np.shape(XC)[1]
    return(rows, cols)


def read_exf_fields(Lf, exf_dir, var_name, exf_shape, year, timestep,
                    sNx, sNy, faces, face_size_dict, unadjusted):

    if unadjusted:
        if var_name!='RUNOFF':
            file_path = os.path.join(exf_dir,'L1_exf_'+var_name+'_'+str(year))
            print('    Reading '+file_path)
            exf_field = np.fromfile(file_path,'>f4')
            total_timesteps = int(np.size(exf_field)/(exf_shape[0]*exf_shape[1]))
            exf_field = np.reshape(exf_field,(total_timesteps,exf_shape[0],exf_shape[1]))
            exf_field = exf_field[timestep,:,:]
        else:
            file_path = os.path.join(exf_dir, 'L1_exf_' + var_name + '.bin')
            runoff_compact = np.fromfile(file_path,'>f4')
            n_rows = int(np.size(runoff_compact)/(12*sNx))
            runoff_compact = np.reshape(runoff_compact,(12,n_rows,sNx))
            exf_field = Lf.read_compact_grid_to_stitched_grid(runoff_compact, sNx, sNy, faces, face_size_dict, dim=3)
            exf_field = exf_field[timestep,:,:]
    else:
        N = 0
        for face in faces:
            N += face_size_dict[face][0]*face_size_dict[face][1]

        file_path = os.path.join(exf_dir, 'L1_exf_' + var_name + '_'+str(year))
        compact = np.fromfile(file_path, '>f4')
        n_rows = int(N/sNx)
        n_timesteps = int(np.size(compact) / (n_rows * sNx))
        compact = np.reshape(compact, (n_timesteps, n_rows, sNx))
        exf_field = Lf.read_compact_grid_to_stitched_grid(compact, sNx, sNy, faces, face_size_dict, dim=3)
        exf_field = exf_field[timestep, :, :]

    return(exf_field)


def create_exf_plot(config_dir, L1_model_name, exf_shape,
                    sNx, sNy, faces, face_size_dict, unadjusted=False):

    sys.path.insert(1, os.path.join(config_dir, 'L1', L1_model_name, 'utils'))
    import L1_CE_Greenland_functions as Lf

    input_dir = os.path.join(config_dir,'L1', L1_model_name, 'input', 'exf')
    output_dir = os.path.join(config_dir,'L1', L1_model_name, 'plots', 'init_files')

    if L1_model_name+'_mean_exf_timeseries.nc' not in os.listdir(input_dir):
        print('    - Generating the '+L1_model_name+'_mean_exf_timeseries.nc file')
        var_names, timeseries_set = create_mean_exf_timeseries(Lf, config_dir, L1_model_name, sNx, sNy, faces, face_size_dict)
    else:
        print('    - Reading timeseries from the ' + L1_model_name + '_mean_exf_timeseries.nc file')
        var_names, timeseries_set = read_mean_exf_timeseries(config_dir, L1_model_name)


    fig = plt.figure(figsize=(20, 8))
    plt.style.use("dark_background")

    for ff in range(len(var_names)):

        print('    - Plotting '+var_names[ff])

        plt.subplot(4, 2, ff+1)
        cmap = 'viridis'

        timeseries = timeseries_set[ff]

        if var_names[ff] in ['UWIND','VWIND']:
            cmap = cm.balance
            val = np.max(np.abs(timeseries[timeseries[:,1] != 0, 1]))
            vmin = -val
            vmax = val
            units = 'm/s'
        else:
            if np.any(timeseries[:,1]>0):
                vmin = np.min(timeseries[timeseries[:,1] != 0, 1])
                vmax = np.max(timeseries[timeseries[:,1] != 0, 1])
            else:
                vmin=-0.1
                vmax=0.1

        if var_names[ff] in ['AQH']:
            cmap = cm.haline
            units = 'kg/kg'
        if var_names[ff] in ['RUNOFF']:
            cmap = cm.rain
            units = 'm/s'
        if var_names[ff] in ['PRECIP']:
            cmap = cm.rain
            units = 'm/s'
        if var_names[ff] in ['ATEMP']:
            cmap = cm.thermal
            units = 'K'
        if var_names[ff] in ['SWDOWN','LWDOWN']:
            cmap = cm.thermal
            units = 'W/m$^2$'

        C = plt.scatter(timeseries[:,0],timeseries[:,1],s=1,c=timeseries[:,1],cmap=cmap,vmin=vmin,vmax=vmax)

        if ff not in [6,7]:
            plt.gca().set_xticklabels([])

        cbar = plt.colorbar(C)
        cbar.set_label(units)
        plt.ylabel(var_names[ff])

        y_range = vmax-vmin
        plt.gca().set_ylim([vmin-0.05*y_range,vmax+0.05*y_range])

        plt.grid(linestyle='--',linewidth=0.5,alpha=0.3)

    plt.savefig(os.path.join(output_dir, L1_model_name+'_exf_mean_timeseries.png'), bbox_inches='tight')
    plt.close(fig)


########################################################################################################################
# User inputs


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_plot(config_dir)

