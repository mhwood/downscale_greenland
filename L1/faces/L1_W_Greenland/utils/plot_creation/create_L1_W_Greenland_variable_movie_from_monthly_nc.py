
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
import shutil
import cmocean.cm as cm
import argparse
from matplotlib.patches import Rectangle
import sys


def read_field_from_monthly_ncs(config_dir, config_name, field_name):

    results_dir = os.path.join(config_dir, 'L1', config_name, 'results')

    if field_name == 'ETAN':
        file_prefix = 'surfDiag'
    elif field_name in ['Area','SIheff','SIhsnow','SIuice','SIvice']:
        file_prefix = 'seaiceDiag'
    elif field_name in ['Theta_AW','Salt_AW']:
        file_prefix = 'awDiag'
    else:
        file_prefix = 'dynDiag'

    grid_started = False

    year_months = ['199201','199202',
                   '199203','199204',
                   '199205','199206',
                   '199207','199208',
                   '199209','199210',
                   '199211','199212']

    for year_month in year_months:
        for file_name in os.listdir(os.path.join(results_dir,file_prefix)):
            if file_name.split('.')[1]==year_month:
                print('    - Reading from '+file_name)

                if field_name in ['SPEED']:
                    ds = nc4.Dataset(os.path.join(results_dir, file_prefix, file_name))
                    uvel = ds.variables['UVEL'][:, :, :]
                    vvel = ds.variables['VVEL'][:, :, :]
                    ds.close()
                    field = (uvel ** 2 + vvel ** 2) ** 0.5
                elif field_name in ['Vorticity']:
                    ds = nc4.Dataset(os.path.join(results_dir, file_prefix, file_name))
                    uvel = ds.variables['Uvel'][:, :, :]
                    vvel = ds.variables['Vvel'][:, :, :]
                    iter_numbers = ds.variables['iterations'][:]
                    ds.close()

                    grid_file = os.path.join(config_dir,'nc_grids',config_name+'_grid.nc')
                    ds = nc4.Dataset(grid_file)
                    L1_YC = ds.variables['YC'][:,:]
                    L1_DXC = ds.variables['dxC'][:,:]
                    L1_DYC = ds.variables['dyC'][:,:]
                    L1_RAZ = ds.variables['rAz'][:,:]
                    ds.close()

                    # L1_DXC = L1_DXC[:, :-1]
                    # L1_DYC = L1_DYC[:-1, :]
                    # L1_RAZ = L1_RAZ[:-1, :-1]

                    Omega = 7.2921e-5
                    L1_f_grid = 2 * Omega * np.sin(np.deg2rad(L1_YC))

                    field = np.zeros((np.shape(uvel)[0],np.shape(uvel)[1]-1, np.shape(uvel)[2]-1))

                    for i in range(np.shape(field)[0]):
                        numerator = np.diff(uvel[i, :, :] * L1_DXC, axis=0)[:, :-1] + np.diff(vvel[i, :, :] * L1_DYC, axis=1)[:-1, :]
                        denominator = L1_RAZ[:-1, :-1]
                        zeta = np.zeros_like(numerator)
                        zeta[denominator != 0] = numerator[denominator != 0] / denominator[denominator != 0]
                        var_field = zeta / L1_f_grid[:-1, :-1]
                        field[i, :, :] = var_field

                else:
                    ds = nc4.Dataset(os.path.join(results_dir, file_prefix, file_name))
                    field = ds.variables[field_name.split('_')[0]][:, :, :]
                    iter_numbers = ds.variables['iterations'][:]
                    ds.close()
                    field = field[:, : ,:]


                if not grid_started:
                    grid_started = True
                    grid = field
                    all_iter_numbers = np.reshape(iter_numbers,(np.size(iter_numbers),1))
                else:
                    grid = np.concatenate([grid,field],axis=0)
                    all_iter_numbers = np.concatenate([all_iter_numbers,np.reshape(iter_numbers,(np.size(iter_numbers),1))],axis=0)

    return (grid, all_iter_numbers)

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

    field_grid, all_iter_numbers = read_field_from_monthly_ncs(config_dir, config_name, field_name)

    output_dir = os.path.join(config_dir,'L1',config_name,'plots','output')
    # if field_name not in os.listdir(output_dir):
    #     os.mkdir(os.path.join(output_dir,field_name))
    output_dir = os.path.join(output_dir,field_name)

    min_index = 0
    max_index = min_index + (366*24*60*60)/300

    if remove_old:
        os.system('rm -rf '+output_dir+'/*')

    if field_name == 'Theta':
        vmin = -1.9
        vmax = 7.5
        cmap = cm.thermal
    if field_name == 'Theta_AW':
        vmin = -1.9
        vmax = 5
        cmap = cm.thermal
    if field_name == 'Salt':
        vmin = 33
        vmax = 35.2
        cmap = cm.haline
    if field_name == 'Uvel' or field_name == 'Vvel':
        vmin = -1
        vmax = 1
        cmap = cm.balance
    if field_name == 'EtaN':
        vmin = -2
        vmax = 0.4
        cmap = 'viridis'
    if field_name == 'Speed':
        vmin = 0
        vmax = 1
        cmap = 'viridis'
    if field_name == 'Vorticity':
        vmin = -0.5
        vmax = 0.5
        cmap = cm.balance
    if field_name == 'Area':
        vmin = 0
        vmax = 1
        cmap = cm.ice

    panel_numbers = np.arange(0,np.shape(field_grid)[0],skip)

    counter = 0
    for i in panel_numbers:

        field_subset = field_grid[i,:,:]
        print('Timestep '+str(i)+' data range: '+str(np.min(field_subset[field_subset!=0]))+' to '+str(np.max(field_subset[field_subset!=0])))

        fig = plt.figure(figsize=(7,11))
        plt.style.use('dark_background')

        gs2 = GridSpec(15, 12, left=0.05, right=0.95, hspace=0.05)

        ax1 = fig.add_subplot(gs2[:-2,:])
        C = ax1.imshow(field_grid[i,:,:],origin='lower',vmin=vmin,vmax=vmax, cmap=cmap)
        plt.colorbar(C, fraction=0.031, pad=0.04)
        plt.title(field_name+', timestep = '+str(i))

        ax2 = fig.add_subplot(gs2[-1, 2:-1])
        width = (all_iter_numbers[i,0]-min_index)/(max_index-min_index)
        rect = Rectangle((1992,0),width,1,fc='silver',ec='white')
        ax2.add_patch(rect)
        ax2.set_xlim([1992,1993])
        ax2.set_ylim([0,1])
        ax2.set_xticks(np.arange(1992,1993,1/12))
        ax2.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
        ax2.set_yticks([])
        plt.xlabel('Year: 1992')

        output_path = os.path.join(output_dir,config_name+'_'+field_name+'_'+'{:04d}'.format(counter)+'.png')
        plt.savefig(output_path,bbox_inches='tight')
        plt.close(fig)
        counter += 1

    compile_panels_to_movie(config_dir, config_name, field_name)





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
   

