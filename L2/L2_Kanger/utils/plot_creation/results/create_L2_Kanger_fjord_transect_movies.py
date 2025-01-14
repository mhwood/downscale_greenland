
import os
import netCDF4 as nc4
import numpy as np
import argparse
import sys
import cmocean.cm as cm
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import datetime as dt
from datetime import datetime, timedelta
from matplotlib.patches import Rectangle

def read_transect_from_nc(config_dir, experiment, model_name, fjord_name, var_name):

    if experiment=='control':
        results_dir = 'results'
    else:
        results_dir = 'results_'+experiment

    output_file = os.path.join(config_dir, 'L2', model_name, results_dir , 'dv',
                               model_name + '_' + fjord_name + '_Trough_' + var_name +'.nc')
    print('Reading from '+output_file)

    ds = nc4.Dataset(output_file)
    distance = ds.variables['distance'][:]
    depth = ds.variables['depth'][:]
    iterations = ds.variables['iterations'][:]
    # thickness = ds.variables['ice_thickness'][:]
    bathymetry = ds.variables['bathymetry'][:]
    var_grid = ds.variables[field_name.upper()][:, :, :]
    ds.close()

    distance = np.flip(distance)
    bathymetry = np.flip(bathymetry)
    var_grid = np.flip(var_grid, axis=-1)

    return(iterations, depth, distance, var_grid, bathymetry)

def YMD_to_DecYr(year,month,day,hour=0,minute=0,second=0):
    date = dt.datetime(year,month,day,hour,minute,second)
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    decimal_fraction = float(date.toordinal() - start) / year_length
    dec_yr = year+decimal_fraction
    return(dec_yr)

def iter_number_to_date(iter_number):
    seconds_per_iter = 60
    total_seconds = iter_number*seconds_per_iter
    date = datetime(1992,1,1) + timedelta(seconds=total_seconds)
    return(date)

def plot_panel(config_dir, plot_dir,output_file_name, L2_model_name, date_str, depth,
               distance_control, var_grid_control, bathymetry_control,
               distance_melange, var_grid_melange, bathymetry_melange):

    fig = plt.figure(figsize=(10,8))
    plt.style.use('dark_background')

    gs2 = GridSpec(27, 23, left=0.1, right=0.95, top=0.95, bottom=0.05, hspace=0.05)

    vmin = -2
    vmax = 3

    dmin=-0.5
    dmax=0.5

    #############################################################################################################

    ax1 = fig.add_subplot(gs2[:7, :-3])
    ax1.pcolormesh(distance_control, depth, var_grid_control, cmap='turbo', vmin=vmin, vmax=vmax)
    # ax1.plot(distance_control, thickness_control, 'k-')
    ax1.plot(distance_control, bathymetry_control, 'k-')

    bathy_outline_zach = np.vstack([np.column_stack([distance_control, bathymetry_control]),
                                    np.flipud(np.column_stack([distance_control, 2000*np.ones_like(bathymetry_control)]))])
    bathy_poly = Polygon(bathy_outline_zach, facecolor='silver')
    ax1.add_patch(bathy_poly)

    # ice_outline_zach = np.vstack([np.column_stack([distance_control, thickness_control]),
    #                                 np.flipud(np.column_stack([distance_control, -10 * np.ones_like(bathymetry_control)]))])
    # ice_poly = Polygon(ice_outline_zach, facecolor='white')
    # ax1.add_patch(ice_poly)
    ax1.text(np.max(distance_control), 20, ' Zachariae', ha='left', va='top', color='k')

    # ax1.invert_yaxis()
    # max_x_index = np.min(np.where(bathymetry_control == 0)[0] - 1)
    # ax1.set_xlim([distance_control[max_x_index], np.min(distance_control)+150])
    ax1.set_ylim([1000, 0])
    ax1.set_ylabel('Depth (m)')
    ax1.set_title('Along-Fjord Transects')

    #############################################################################################################

    ax2 = fig.add_subplot(gs2[8:14, :-3])
    ax2.pcolormesh(distance_melange, depth, var_grid_melange, cmap='turbo', vmin=vmin, vmax=vmax)
    # ax2.plot(distance_melange, thickness_melange, 'k-')
    ax2.plot(distance_melange, bathymetry_melange, 'k-')

    bathy_outline_melange = np.vstack([np.column_stack([distance_melange, bathymetry_melange]),
                                    np.flipud(np.column_stack([distance_melange, 2000 * np.ones_like(bathymetry_melange)]))])
    bathy_poly = Polygon(bathy_outline_melange, facecolor='silver')
    ax2.add_patch(bathy_poly)

    # ice_outline_melange = np.vstack([np.column_stack([distance_melange, thickness_melange]),
    #                               np.flipud(np.column_stack([distance_melange, -10 * np.ones_like(bathymetry_melange)]))])
    # ice_poly = Polygon(ice_outline_melange, facecolor='white')
    # ax2.add_patch(ice_poly)
    ax2.text(np.max(distance_melange),20,' 79N', ha='left', va='top', color='k')

    # max_x_index = np.min(np.where(bathymetry_melange==0)[0]-1)

    # ax2.invert_yaxis()
    # ax2.set_xlim([distance_melange[max_x_index], np.min(distance_melange)])
    ax2.set_ylim([1000, 0])
    ax2.set_ylabel('Depth (m)')
    ax2.set_xlabel('Distance Along Transect (km)')

    #############################################################################################################

    ax2 = fig.add_subplot(gs2[15:22, :-3])
    ax2.pcolormesh(distance_melange, depth, var_grid_melange - var_grid_control, cmap='seismic', vmin=dmin, vmax=dmax)
    # ax2.plot(distance_melange, thickness_melange, 'k-')
    ax2.plot(distance_melange, bathymetry_melange, 'k-')

    bathy_outline_melange = np.vstack([np.column_stack([distance_melange, bathymetry_melange]),
                                       np.flipud(np.column_stack(
                                           [distance_melange, 2000 * np.ones_like(bathymetry_melange)]))])
    bathy_poly = Polygon(bathy_outline_melange, facecolor='silver')
    ax2.add_patch(bathy_poly)

    # ice_outline_melange = np.vstack([np.column_stack([distance_melange, thickness_melange]),
    #                               np.flipud(np.column_stack([distance_melange, -10 * np.ones_like(bathymetry_melange)]))])
    # ice_poly = Polygon(ice_outline_melange, facecolor='white')
    # ax2.add_patch(ice_poly)
    ax2.text(np.max(distance_melange), 20, ' 79N', ha='left', va='top', color='k')

    # max_x_index = np.min(np.where(bathymetry_melange==0)[0]-1)

    # ax2.invert_yaxis()
    # ax2.set_xlim([distance_melange[max_x_index], np.min(distance_melange)])
    ax2.set_ylim([1000, 0])
    ax2.set_ylabel('Depth (m)')
    ax2.set_xlabel('Distance Along Transect (km)')

    #############################################################################################################

    ax3 = fig.add_subplot(gs2[15:22, -2])
    x = np.array([0,1])
    y = np.linspace(dmin, dmax, 100)
    X, Y = np.meshgrid(x,y)
    ax3.pcolormesh(X,Y,Y,cmap='seismic')
    ax3.yaxis.tick_right()
    ax3.set_ylabel('Temperature Difference ($^{\circ}$C)')
    ax3.yaxis.set_label_position("right")
    ax3.set_xticks([])

    #############################################################################################################

    ax3 = fig.add_subplot(gs2[:13, -2])
    x = np.array([0, 1])
    y = np.linspace(vmin, vmax, 100)
    X, Y = np.meshgrid(x, y)
    ax3.pcolormesh(X, Y, Y, cmap='turbo')
    ax3.yaxis.tick_right()
    ax3.set_ylabel('Temperature ($^{\circ}$C)')
    ax3.yaxis.set_label_position("right")
    ax3.set_xticks([])

    #############################################################################################################

    yr = int(date_str[:4])
    mo = int(date_str[4:6])
    dy = int(date_str[6:8])
    dec_yr = YMD_to_DecYr(yr, mo, dy)

    min_year = 2015
    max_year = 2022
    ax4 = fig.add_subplot(gs2[-1, 1:-3])
    width = (dec_yr - min_year)
    rect = Rectangle((min_year, 0), width, 1, fc='silver', ec='white')
    ax4.add_patch(rect)
    ax4.set_xlim([min_year, max_year])
    ax4.set_ylim([0, 1])
    for i in range(min_year, max_year, 1):
        plt.plot([i, i], [0, 1], 'w-', linewidth=0.5)
    ax4.set_xticks(np.arange(min_year, max_year, 1))
    ax4.set_yticks([])

    ax3 = fig.add_subplot(gs2[-3, 1:-3])
    width = (dec_yr - yr)
    rect = Rectangle((yr, 0), width, 1, fc='silver', ec='white')
    ax3.add_patch(rect)
    ax3.set_xlim([yr, yr + 1])
    ax3.set_ylim([0, 1])
    for i in range(1, 12):
        plt.plot([yr + i / 12, yr + i / 12], [0, 1], 'w-', linewidth=0.5)
    ax3.set_xticks(np.arange(yr + 1 / 24, yr + 1, 1 / 12))
    ax3.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
    ax3.set_yticks([])

    #############################################################################################################

    output_file = os.path.join(config_dir, 'L2', L2_model_name, plot_dir, 'transects',
                 'panels', field_name, output_file_name)

    plt.savefig(output_file)
    plt.close(fig)

def compile_panels_to_movie(config_dir,L2_model_name,field_name, plot_dir):
    pwd = os.getcwd()

    panels_dir = os.path.join(config_dir,'L2',L2_model_name,plot_dir,'transects','panels', field_name)

    # get a list of the file names
    all_file_names = []
    for file_name in os.listdir(panels_dir):
        if file_name[0]!='.' and file_name[-4:]=='.png':
            all_file_names.append(file_name)
    all_file_names = sorted(all_file_names)

    # make a temporary dir where we will link all available images and go there
    panels_dir_tmp = os.path.join(config_dir, 'L2', L2_model_name, plot_dir, 'transects', 'panels_tmp')
    os.mkdir(panels_dir_tmp)
    os.chdir(panels_dir_tmp)

    # link all of the images
    for ff in range(len(all_file_names)):
        os.system('ln -s ' + '../panels/'+field_name+'/'+all_file_names[ff]+' panel_'+'{:05d}'.format(ff)+'.png')

    output_name = L2_model_name+'_'+field_name+'.mp4'

    os.system("ffmpeg -r 5 -i panel_%05d.png -vcodec mpeg4 -b 3M -y " + output_name)
    os.rename(output_name, os.path.join('..', output_name))

    os.chdir(pwd)
    os.system('rm -rf ' + panels_dir_tmp)


def create_movies(config_dir, field_name, print_level):

    L2_model_name = 'L2_Kanger'

    metadata_dict = {'Theta': [-1.9, 5, 'turbo', '$^{\circ}$C', 'Potential Temperature (Surface)'],
                     'Salt': [31.5, 35, cm.haline, 'psu', 'Practical Salinity (Surface)'],  #
                     'Uvel': [-1, 1, cm.balance, 'm/s', 'Zonal Velocity'],  #
                     'Vvel': [-1, 1, cm.balance, 'm/s', 'Meridional Velocity']}
    
    plot_dir = 'plots'


    remove_old = False
    skip = False

    if 'transects' not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,plot_dir)):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,plot_dir,'transects'))
    if 'panels' not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,plot_dir,'transects')):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,plot_dir,'transects','panels'))
    # print(os.listdir(os.path.join(config_dir,'L2',L2_model_name,plot_dir,'transects','panels')))
    # print(field_name in os.listdir(os.path.join(config_dir,'L2',L2_model_name,plot_dir,'transects','panels')))
    if field_name not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,plot_dir,'transects','panels')):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,plot_dir,'transects','panels',field_name))

    iterations, depth, distance_control, var_grid_control, bathymetry_control = \
        read_transect_from_nc(config_dir, 'control', L2_model_name, 'Kanger', field_name.upper())
    distance_control *= 1e-3

    # print(var_grid_control[:, 10, -10])

    _, _, distance_melange, var_grid_melange, bathymetry_melange  = \
        read_transect_from_nc(config_dir, 'melange', L2_model_name, 'Kanger', field_name.upper())
    distance_melange *= 1e-3

    # print(var_grid_melange[:, 10, -10])

    for i in range(111):#len(iterations)-2,len(iterations)-1):
        date = iter_number_to_date(iterations[i])
        date_str = str(date.year)+'{:02d}'.format(date.month)+'{:02d}'.format(date.day)

        output_file_name = L2_model_name+'_Fjord_Transect_'+field_name+'_'+date_str+'.png'

        if output_file_name not in os.listdir(os.path.join(config_dir, 'L2', L2_model_name,
                                                           plot_dir, 'transects',
                                                            'panels', field_name)):
            var_grid_control_timestep = var_grid_control[i, :, :]
            var_grid_melange_timestep = var_grid_melange[i, :, :]
            plot_panel(config_dir, plot_dir, output_file_name, L2_model_name, date_str, depth,
                       distance_control,var_grid_control_timestep, bathymetry_control,
                       distance_melange, var_grid_melange_timestep, bathymetry_melange)

    # print('    - Making the movie for year '+str(year))
    compile_panels_to_movie(config_dir, L2_model_name, field_name, plot_dir)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L2, L2, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-f", "--field_name", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="field_name",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    field_name = args.field_name

    create_movies(config_dir, field_name,  print_level=4)
   

