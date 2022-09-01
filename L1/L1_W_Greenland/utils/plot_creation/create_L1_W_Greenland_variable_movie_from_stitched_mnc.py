
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
import shutil
import cmocean.cm as cm
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
            file_name = file_prefix+'.'+ '{:010d}'.format(288) +'.t' + '{:03d}'.format(tile_number) + '.nc'
            for n in range(N):
                if file_name in os.listdir(os.path.join(run_dir, 'mnc_' + '{:04d}'.format(n + 1))):

                    if field_name in ['SPEED','VORTICITY']:
                        ds = nc4.Dataset(os.path.join(run_dir, 'mnc_' + '{:04d}'.format(n + 1), file_name))
                        uvel = ds.variables['UVEL'][:, :, :, :]
                        vvel = ds.variables['VVEL'][:, :, :, :]
                        ds.close()
                        uvel = uvel[:, 0, :, :-1]
                        vvel = vvel[:, 0, :-1, :]
                        field = (uvel ** 2 + vvel ** 2) ** 0.5
                    else:
                        ds = nc4.Dataset(os.path.join(run_dir, 'mnc_' + '{:04d}'.format(n + 1), file_name))
                        field = ds.variables[field_name][:, :, : , :]
                        ds.close()
                        field = field[:, 0 ,: ,:]

                    if field_name=='UVEL':
                        field = field[:,:,:-1]

                    if field_name=='VVEL':
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

    panels_dir = os.path.join(config_dir,'L1',config_name,'plots','output',field_name)

    file_name = 'L1_'+field_name+'.mp4'

    os.chdir(panels_dir)
    os.system("ffmpeg -r 5 -i L1_W_Greenland_"+field_name+"_%04d.png -vcodec mpeg4 -b 1M -y " + file_name)
    os.rename(file_name, os.path.join('..', file_name))

    os.chdir(pwd)


def plot_mnc_fields(config_dir,field_name,remove_old,skip):
    config_name = 'L1_W_Greenland'

    ordered_nonblank_tiles = [[6, 9, 12], [5, 8, 11], [4, 7, 10], [3, 2, 1]]
    ordered_nonblank_rotations = [[1, 1, 1], [1, 1, 1], [1, 1, 1], [2, 2, 2]]

    sNx = 180
    sNy = 180

    run_dir = os.path.join(config_dir, 'L1', config_name, 'run')
    field_grid = read_mnc_field(run_dir, field_name, sNx, sNy,
                                ordered_nonblank_tiles, ordered_nonblank_rotations)

    output_dir = os.path.join(config_dir,'L1',config_name,'plots','output')
    if field_name not in os.listdir(output_dir):
        os.mkdir(os.path.join(output_dir,field_name))
    output_dir = os.path.join(output_dir,field_name)

    if remove_old:
        os.system('rm -rf '+output_dir+'/*')

    if field_name == 'THETA':
        vmin = -1.9
        vmax = 7.5
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
        vmin = -2
        vmax = 0.4
        cmap = 'viridis'
    if field_name == 'SPEED':
        vmin = 0
        vmax = 1
        cmap = 'viridis'
    if field_name == 'SIarea':
        vmin = 0
        vmax = 1
        cmap = cm.ice

    panel_numbers = np.arange(0,np.shape(field_grid)[0],skip)

    counter = 0
    for i in panel_numbers:

        field_subset = field_grid[i,:,:]
        print('Timestep '+str(i)+' data range: '+str(np.min(field_subset[field_subset!=0]))+' to '+str(np.max(field_subset[field_subset!=0])))

        fig = plt.figure(figsize=(8,10))
        plt.style.use('dark_background')

        C = plt.imshow(field_grid[i,:,:],origin='lower',vmin=vmin,vmax=vmax, cmap=cmap)
        plt.colorbar(C, fraction=0.031, pad=0.04)
        plt.title(field_name+', timestep = '+str(i))

        output_path = os.path.join(output_dir,config_name+'_'+field_name+'_'+'{:04d}'.format(counter)+'.png')
        plt.savefig(output_path)
        plt.close(fig)
        counter += 1

    # compile_panels_to_movie(config_dir, config_name, field_name)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
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
   

