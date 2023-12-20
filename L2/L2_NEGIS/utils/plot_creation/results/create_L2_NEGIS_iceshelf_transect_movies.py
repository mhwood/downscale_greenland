
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

def read_transect_from_nc(config_dir, model_name, fjord_name, var_name):

    output_file = os.path.join(config_dir, 'L2', model_name, 'results_shelfice', 'dv',
                               model_name + '_' + fjord_name + '_' + var_name + '.nc')
    ds = nc4.Dataset(output_file)
    distance = ds.variables['distance'][:]
    depth = ds.variables['depth'][:]
    iterations = ds.variables['iterations'][:]
    thickness = ds.variables['ice_thickness'][:]
    bathymetry = ds.variables['bathymetry'][:]
    var_grid = ds.variables[field_name.upper()][:, :, :]
    ds.close()

    return(iterations, depth, distance, var_grid, thickness, bathymetry)

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

def plot_panel(config_dir, output_file_name, L2_model_name, date_str, depth,
               distance_Zach, var_grid_Zach, thickness_Zach, bathymetry_Zach,
               distance_79N, var_grid_79N, thickness_79N, bathymetry_79N):

    fig = plt.figure(figsize=(10,8))
    plt.style.use('dark_background')

    gs2 = GridSpec(19, 23, left=0.1, right=0.95, top=0.95, bottom=0.05, hspace=0.05)

    vmin = -2
    vmax = 0.5

    #############################################################################################################

    ax1 = fig.add_subplot(gs2[:7, :-3])
    ax1.pcolormesh(distance_Zach, depth, var_grid_Zach, cmap='turbo', vmin=vmin, vmax=vmax)
    ax1.plot(distance_Zach, thickness_Zach, 'k-')
    ax1.plot(distance_Zach, bathymetry_Zach, 'k-')

    bathy_outline_zach = np.vstack([np.column_stack([distance_Zach, bathymetry_Zach]),
                                    np.flipud(np.column_stack([distance_Zach, 2000*np.ones_like(bathymetry_Zach)]))])
    bathy_poly = Polygon(bathy_outline_zach, facecolor='silver')
    ax1.add_patch(bathy_poly)

    ice_outline_zach = np.vstack([np.column_stack([distance_Zach, thickness_Zach]),
                                    np.flipud(np.column_stack([distance_Zach, -10 * np.ones_like(bathymetry_Zach)]))])
    ice_poly = Polygon(ice_outline_zach, facecolor='white')
    ax1.add_patch(ice_poly)
    ax1.text(np.max(distance_Zach), 20, ' Zachariae', ha='left', va='top', color='k')

    # ax1.invert_yaxis()
    ax1.set_xlim([np.max(distance_Zach), np.max(distance_Zach)-100])
    ax1.set_ylim([1000, 0])
    ax1.set_ylabel('Depth (m)')
    ax1.set_title('Along-Fjord Transects')

    #############################################################################################################

    ax2 = fig.add_subplot(gs2[8:14, :-3])
    ax2.pcolormesh(distance_79N, depth, var_grid_79N, cmap='turbo', vmin=vmin, vmax=vmax)
    ax2.plot(distance_79N, thickness_79N, 'k-')
    ax2.plot(distance_79N, bathymetry_79N, 'k-')

    bathy_outline_79N = np.vstack([np.column_stack([distance_79N, bathymetry_79N]),
                                    np.flipud(np.column_stack([distance_79N, 2000 * np.ones_like(bathymetry_79N)]))])
    bathy_poly = Polygon(bathy_outline_79N, facecolor='silver')
    ax2.add_patch(bathy_poly)

    ice_outline_79N = np.vstack([np.column_stack([distance_79N, thickness_79N]),
                                  np.flipud(np.column_stack([distance_79N, -10 * np.ones_like(bathymetry_79N)]))])
    ice_poly = Polygon(ice_outline_79N, facecolor='white')
    ax2.add_patch(ice_poly)
    ax2.text(np.max(distance_79N),20,' 79N', ha='left', va='top', color='k')

    # ax2.invert_yaxis()
    ax2.set_xlim([np.max(distance_79N), np.max(distance_79N) - 100])
    ax2.set_ylim([1000, 0])
    ax2.set_ylabel('Depth (m)')
    ax2.set_xlabel('Distance Along Transect (km)')

    #############################################################################################################

    ax3 = fig.add_subplot(gs2[:13, -2])
    x = np.array([0,1])
    y = np.linspace(vmin, vmax, 100)
    X, Y = np.meshgrid(x,y)
    ax3.pcolormesh(X,Y,Y,cmap='turbo')
    ax3.yaxis.tick_right()
    ax3.set_ylabel('Temperature ($^{\circ}$C)')
    ax3.yaxis.set_label_position("right")
    ax3.set_xticks([])

    #############################################################################################################

    yr = int(date_str[:4])
    mo = int(date_str[4:6])
    dy = int(date_str[6:8])
    dec_yr = YMD_to_DecYr(yr, mo, dy)

    min_year = 1992
    max_year = 2022
    ax4 = fig.add_subplot(gs2[-1, 1:-3])
    width = (dec_yr - min_year)
    rect = Rectangle((min_year, 0), width, 1, fc='silver', ec='white')
    ax4.add_patch(rect)
    ax4.set_xlim([min_year, max_year])
    ax4.set_ylim([0, 1])
    for i in range(min_year, max_year, 5):
        plt.plot([i, i], [0, 1], 'w-', linewidth=0.5)
    ax4.set_xticks(np.arange(min_year, max_year, 5))
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

    output_file = os.path.join(config_dir, 'L2', L2_model_name, 'plots_shelfice', 'transects',
                 'panels', field_name, output_file_name)

    plt.savefig(output_file)
    plt.close(fig)

def compile_panels_to_movie(config_dir,L2_model_name,field_name):
    pwd = os.getcwd()

    panels_dir = os.path.join(config_dir,'L2',L2_model_name,'plots_shelfice','transects','panels', field_name)

    # get a list of the file names
    all_file_names = []
    for file_name in os.listdir(panels_dir):
        if file_name[0]!='.' and file_name[-4:]=='.png':
            all_file_names.append(file_name)
    all_file_names = sorted(all_file_names)

    # make a temporary dir where we will link all available images and go there
    panels_dir_tmp = os.path.join(config_dir, 'L2', L2_model_name, 'plots_shelfice', 'transects', 'panels_tmp')
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

    L2_model_name = 'L2_NEGIS'

    metadata_dict = {'Theta': [-1.9, 5, 'turbo', '$^{\circ}$C', 'Potential Temperature (Surface)'],
                     'Salt': [31.5, 35, cm.haline, 'psu', 'Practical Salinity (Surface)'],  #
                     'Uvel': [-1, 1, cm.balance, 'm/s', 'Zonal Velocity'],  #
                     'Vvel': [-1, 1, cm.balance, 'm/s', 'Meridional Velocity']}

    remove_old = False
    skip = False

    if 'transects' not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,'plots_shelfice',)):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,'plots_shelfice','transects'))
    if 'panels' not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,'plots_shelfice','transects')):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,'plots_shelfice','transects','panels'))
    if field_name not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,'plots_shelfice','transects','panels')):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,'plots_shelfice','transects','panels',field_name))

    iterations, depth, distance_Zach, var_grid_Zach, thickness_Zach, bathymetry_Zach = \
        read_transect_from_nc(config_dir, L2_model_name, 'Zach', field_name.upper())
    distance_Zach *= 1e-3

    # print(var_grid_Zach[:, 10, -10])

    _, _, distance_79N, var_grid_79N, thickness_79N, bathymetry_79N  = \
        read_transect_from_nc(config_dir, L2_model_name, '79N', field_name.upper())
    distance_79N *= 1e-3

    # print(var_grid_79N[:, 10, -10])

    for i in range(len(iterations)):
        date = iter_number_to_date(iterations[i])
        date_str = str(date.year)+'{:02d}'.format(date.month)+'{:02d}'.format(date.day)
        print(date_str)

        output_file_name = L2_model_name+'_Iceshelf_Transect_'+field_name+'_'+date_str+'.png'
        if output_file_name not in []:#os.listdir(os.path.join(config_dir, 'L2', L2_model_name,
                                                           'plots_shelfice', 'transects',
                                                            'panels', field_name)):
            var_grid_Zach_timestep = var_grid_Zach[i, :, :]
            var_grid_79N_timestep = var_grid_79N[i, :, :]
            plot_panel(config_dir, output_file_name, L2_model_name, date_str, depth,
                       distance_Zach,var_grid_Zach_timestep, thickness_Zach, bathymetry_Zach,
                       distance_79N, var_grid_79N_timestep, thickness_79N, bathymetry_79N)

    compile_panels_to_movie(config_dir, L2_model_name, field_name)






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

    create_movies(config_dir, field_name, print_level=4)
   

