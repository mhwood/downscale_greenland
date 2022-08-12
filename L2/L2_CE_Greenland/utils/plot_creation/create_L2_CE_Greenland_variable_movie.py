
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
import cmocean.cm as cm
import shutil
import argparse
import sys


def read_mnc_field(run_dir, field_name, sNx, sNy,
                   ordered_nonblank_tiles, ordered_nonblack_rotations):

    if field_name == 'ETAN':
        file_prefix = 'surfDiag'
    elif field_name in ['SIarea','SIheff','SIhsnow','SIuice','SIvice']:
        file_prefix = 'seaiceDiag'
    elif field_name in ['THETA_AW']:
        file_prefix = 'awDiag'
    else:
        file_prefix = 'dynDiag'

    grid_started = False

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            file_name = file_prefix+'.'+ '{:010d}'.format(2104704) +'.t' + '{:03d}'.format(tile_number) + '.nc'
            for n in range(N):
                if file_name in os.listdir(os.path.join(run_dir, 'mnc_' + '{:04d}'.format(n + 1))):
                    print(file_name)

                    if field_name in ['SPEED']:
                        ds = nc4.Dataset(os.path.join(run_dir, 'mnc_' + '{:04d}'.format(n + 1), file_name))
                        uvel = ds.variables['UVEL'][:, :, :, :]
                        vvel = ds.variables['VVEL'][:, :, :, :]
                        ds.close()
                        uvel = uvel[:, 0, :, :-1]
                        vvel = vvel[:, 0, :-1, :]
                        field = (uvel ** 2 + vvel ** 2) ** 0.5
                    elif field_name in ['VORTICITY']:
                        ds = nc4.Dataset(os.path.join(run_dir, 'mnc_' + '{:04d}'.format(n + 1), file_name))
                        uvel = ds.variables['UVEL'][:, :, :, :]
                        vvel = ds.variables['VVEL'][:, :, :, :]
                        ds.close()
                        uvel = uvel[:, 0, :, :-1]
                        vvel = vvel[:, 0, :-1, :]

                        ds = nc4.Dataset(os.path.join(run_dir, 'mnc_' + '{:04d}'.format(n + 1), 'grid.t' + '{:03d}'.format(tile_number) + '.nc'))
                        YC = ds.variables['YC'][:, :]
                        DXC = ds.variables['dxC'][:, :]
                        DYC = ds.variables['dyC'][:, :]
                        RAZ = ds.variables['rAz'][:, :]
                        DXC = DXC[:, :-1]
                        DYC = DYC[:-1, :]
                        RAZ = RAZ[:-1, :-1]
                        ds.close()

                        Omega = 7.2921e-5
                        f_grid = 2 * Omega * np.sin(np.deg2rad(YC))
                        field = np.zeros((np.shape(uvel)[0],np.shape(uvel)[1],np.shape(uvel)[2]))
                        for timestep in range(np.shape(uvel)[0]):
                            numerator = np.diff(uvel[timestep,:,:] * DXC, axis=0)[:, :-1] + np.diff(vvel[timestep,:,:] * DYC, axis=1)[:-1, :]
                            denominator = RAZ[:-1, :-1]
                            zeta = np.zeros_like(numerator)
                            zeta[denominator != 0] = numerator[denominator != 0] / denominator[denominator != 0]
                            field[timestep, :-1 ,:-1] = zeta / f_grid[:-1, :-1]
                    else:
                        ds = nc4.Dataset(os.path.join(run_dir, 'mnc_' + '{:04d}'.format(n + 1), file_name))
                        if field_name == 'THETA_AW':
                            field = ds.variables['THETA'][:, :, : , :]
                        else:
                            field = ds.variables[field_name][:, :, :, :]
                        ds.close()
                        field = field[:, 0, : ,:]

                    if field_name in ['UVEL','SIuice']:
                        field = field[:,:,:-1]

                    if field_name in ['VVEL','SIvice']:
                        field = field[:,:-1,:]

                    for i in range(ordered_nonblack_rotations[r][c]):
                        field = np.rot90(field,axes=(1,2))

                    if not grid_started:
                        grid_started = True
                        grid = np.zeros((np.shape(field)[0],
                                         sNy * len(ordered_nonblank_tiles),
                                         sNx * len(ordered_nonblank_tiles[0])))

                    grid[:, r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = field

    return (grid)

def compile_panels_to_movie(config_dir,config_name,field_name):
    pwd = os.getcwd()

    panels_dir = os.path.join(config_dir,'L2',config_name,'plots','output',field_name)

    file_name = 'L2_'+field_name+'.mp4'

    os.chdir(panels_dir)
    os.system("ffmpeg -r 5 -i L2_CE_Greenland_"+field_name+"_%04d.png -vcodec mpeg4 -b 1M -y " + file_name)
    os.rename(file_name, os.path.join('..', file_name))

    os.chdir(pwd)


def plot_mnc_fields(config_dir,field_name,remove_old,skip):

    config_name = 'L2_CE_Greenland'

    ordered_nonblank_tiles = [[1,2,3], [4,5,6]]
    ordered_nonblack_rotations = [[0, 0, 0], [0, 0, 0]]

    sNx = 80
    sNy = 120

    run_dir = os.path.join(config_dir, 'L2', config_name, 'run')
    field_grid = read_mnc_field(run_dir, field_name, sNx, sNy,
                                ordered_nonblank_tiles, ordered_nonblack_rotations)

    output_dir = os.path.join(config_dir,'L2',config_name,'plots','output')
    if field_name not in os.listdir(output_dir):
        os.mkdir(os.path.join(output_dir,field_name))
    output_dir = os.path.join(output_dir,field_name)

    if remove_old:
        os.system('rm -rf '+output_dir+'/*')

    if field_name == 'THETA' or field_name == 'THETA_AW':
        vmin = -1.9
        vmax = 3.5
        cmap = cm.thermal
    if field_name == 'THETA_AW':
        vmin = -0
        vmax = 2
        cmap = cm.thermal
    if field_name == 'SALT':
        vmin = 33
        vmax = 35.2
        cmap = cm.haline
    if field_name == 'UVEL' or field_name == 'VVEL':
        vmin = -1
        vmax = 1
        cmap = cm.balance
    if field_name == 'ETAN':
        vmin = -0.4
        vmax = 2
        cmap = 'viridis'
    if field_name == 'SPEED':
        vmin = 0
        vmax = 1
        cmap = cm.tempo_r
    if field_name == 'VORTICITY':
        vmin = -0.25
        vmax = 0.25
        cmap = cm.curl
    if field_name == 'SIarea':
        vmin = 0
        vmax = 1
        cmap = cm.ice

    panel_numbers = np.arange(0,np.shape(field_grid)[0],skip)

    counter = 0
    for i in panel_numbers:

        fig = plt.figure(figsize=(8,6))
        plt.style.use('dark_background')

        print(np.min(field_grid[i,:,:]),np.max(field_grid[i,:,:]))

        C = plt.imshow(field_grid[i,:,:],origin='lower',vmin=vmin,vmax=vmax, cmap=cmap)
        plt.colorbar(C, fraction=0.031, pad=0.04)
        plt.title(field_name+', timestep = '+str(i))

        output_path = os.path.join(output_dir,config_name+'_'+field_name+'_'+'{:04d}'.format(counter)+'.png')
        plt.savefig(output_path)
        plt.close(fig)
        counter += 1

    compile_panels_to_movie(config_dir, config_name, field_name)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L2, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-f", "--field_name", action="store",
                        help="The name of the field to plot.", dest="field_name",
                        type=str, required=True)

    parser.add_argument("-r", "--remove_old", action="store",
                        help="Choose whether to remove old files (1 is true, 0 is false).", dest="remove_old",
                        type=int, required=False, default = 0)

    parser.add_argument("-s", "--skip", action="store",
                        help="Choose how many panels to skip at a time.", dest="skip",
                        type=int, required=False, default=1)

    args = parser.parse_args()
    config_dir = args.config_dir
    field_name = args.field_name
    remove_old = args.remove_old
    skip = args.skip

    if remove_old==0:
        remove_old = False
    else:
        remove_old = True

    plot_mnc_fields(config_dir,field_name,remove_old,skip)
   

