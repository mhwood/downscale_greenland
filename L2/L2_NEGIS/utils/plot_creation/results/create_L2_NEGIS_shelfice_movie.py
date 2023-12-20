
import os
import netCDF4 as nc4
import numpy as np
import argparse
import sys
# import cmocean.cm as cm
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import datetime as dt
from datetime import datetime, timedelta
from matplotlib.patches import Rectangle
from osgeo import gdal

def read_background_imagery(file_path):
    ds = gdal.Open(file_path)
    R = np.array(ds.GetRasterBand(1).ReadAsArray())
    G = np.array(ds.GetRasterBand(2).ReadAsArray())
    B = np.array(ds.GetRasterBand(3).ReadAsArray())
    rows = np.shape(R)[0]
    cols = np.shape(R)[1]
    R = R.reshape((rows,cols, 1))
    G = G.reshape((rows, cols, 1))
    B = B.reshape((rows, cols, 1))
    image = np.concatenate([B,G,R],axis=2)
    brightness_factor = 0.1 # 0 - 1
    image = (np.max(image)-np.max(image)*(brightness_factor))/(np.max(image)-np.min(image))*(image-np.min(image))+np.max(image)*(brightness_factor)
    image[image<0]=0
    image[image>1]=1
    transform = ds.GetGeoTransform()
    extents = [transform[0],transform[0]+transform[1]*np.shape(image)[1],transform[3]+ transform[5] * np.shape(image)[0], transform[3]]
    # x_resolution = transform[1]
    # y_resolution = transform[5]
    return(image,extents)

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

def read_field_from_nc(config_dir, model_name, var_name):

    output_file = os.path.join(config_dir, 'L2', model_name, 'results_iceplume', 'dv',
                               model_name + '_shelfice_' + var_name + '.nc')
    ds = nc4.Dataset(output_file)
    iterations = ds.variables['iterations'][:]
    X = ds.variables['X'][:,:]
    Y = ds.variables['Y'][:,:]
    var_grid = ds.variables[field_name.upper()][:, :, :]
    ds.close()

    return(iterations,X,Y,var_grid)



def plot_panel(config_dir, output_file_name, L2_model_name, date_str,
               X, Y, var_grid_timestep):

    fig = plt.figure(figsize=(8,8))
    plt.style.use('dark_background')

    gs2 = GridSpec(19, 23, left=0.1, right=0.95, top=0.95, bottom=0.05, hspace=0.05)

    vmin = 0#np.min(var_grid_timestep)
    vmax = 50#np.max(var_grid_timestep)

    add_background_imagery = True

    if add_background_imagery:
        # raise ValueError("Need to generate the background imagery first before using this option")
        file_path = "/Volumes/jakobshavn/Research/Projects/Ocean_Modeling/Projects/Downscale_Greenland/MITgcm/" \
                    "configurations/downscale_greenland/L2/L2_NEGIS/plots/L2_NEGIS_MODIS_20220720_3413.tif"
        background_image, extents = read_background_imagery(file_path)

    #############################################################################################################

    extents[0] *= 1e-3
    extents[1] *= 1e-3
    extents[2] *= 1e-3
    extents[3] *= 1e-3

    ax1 = fig.add_subplot(gs2[:13, :-3])
    if add_background_imagery:
        ax1.imshow(background_image, extent=extents, alpha=0.7)

    var_grid_plot = np.ma.masked_where(var_grid_timestep==0, var_grid_timestep)
    ax1.pcolormesh(X/1000, Y/1000, var_grid_plot, cmap='seismic', vmin=vmin, vmax=vmax)

    # ax1.invert_yaxis()
    # ax1.set_xlim([np.max(distance_Zach), np.max(distance_Zach)-100])
    ax1.set_ylim([-1153267/1000, -987941/1000])
    ax1.set_xlim([407537/1000, 579269/1000])
    ax1.set_ylabel('North Distance (km)')
    ax1.set_xlabel('East Distance (km)')
    ax1.set_title('Ice Shelf Melt Rate')

    #############################################################################################################

    ax3 = fig.add_subplot(gs2[:13, -2])
    x = np.array([0,1])
    y = np.linspace(vmin, vmax, 100)
    X, Y = np.meshgrid(x,y)
    ax3.pcolormesh(X,Y,Y,cmap='seismic')
    ax3.yaxis.tick_right()
    ax3.set_ylabel('Melt Rate (m/yr)')
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

    output_file = os.path.join(config_dir, 'L2', L2_model_name, 'plots_iceplume', 'shelfice',
                 'panels', field_name, output_file_name)

    plt.savefig(output_file)
    plt.close(fig)

def compile_panels_to_movie(config_dir,L2_model_name,field_name):
    pwd = os.getcwd()

    panels_dir = os.path.join(config_dir,'L2',L2_model_name,'plots_iceplume','shelfice','panels', field_name)

    # get a list of the file names
    all_file_names = []
    for file_name in os.listdir(panels_dir):
        if file_name[0]!='.' and file_name[-4:]=='.png':
            all_file_names.append(file_name)
    all_file_names = sorted(all_file_names)

    # make a temporary dir where we will link all available images and go there
    panels_dir_tmp = os.path.join(config_dir, 'L2', L2_model_name, 'plots_iceplume', 'shelfice', 'panels_tmp')
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

    # metadata_dict = {'Theta': [-1.9, 5, 'turbo', '$^{\circ}$C', 'Potential Temperature (Surface)'],
    #                  'Salt': [31.5, 35, cm.haline, 'psu', 'Practical Salinity (Surface)'],  #
    #                  'Uvel': [-1, 1, cm.balance, 'm/s', 'Zonal Velocity'],  #
    #                  'Vvel': [-1, 1, cm.balance, 'm/s', 'Meridional Velocity']}
    #
    # remove_old = False
    # skip = False

    if 'shelfice' not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,'plots_iceplume')):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,'plots_iceplume','shelfice'))
    if 'panels' not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,'plots_iceplume','shelfice')):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,'plots_iceplume','shelfice','panels'))
    if field_name not in os.listdir(os.path.join(config_dir,'L2',L2_model_name,'plots_iceplume','shelfice','panels')):
        os.mkdir(os.path.join(config_dir,'L2',L2_model_name,'plots_iceplume','shelfice','panels',field_name))

    iterations, X, Y, var_grid = \
        read_field_from_nc(config_dir, L2_model_name, field_name)

    # kg / m2 / s * (1 m3 / 1000 kg) / * (60 * 60 * 24)
    var_grid *= -1*60*60*24*365/1000

    for i in range(len(iterations)):
        date = iter_number_to_date(iterations[i])
        date_str = str(date.year)+'{:02d}'.format(date.month)+'{:02d}'.format(date.day)
        print(date_str)

        output_file_name = L2_model_name+'_Iceshelf_Transect_shelfice_'+date_str+'.png'
        if output_file_name not in os.listdir(os.path.join(config_dir, 'L2', L2_model_name,
                                                           'plots_iceplume', 'shelfice',
                                                            'panels', field_name)):
            var_grid_timestep = var_grid[i, :, :]
            plot_panel(config_dir, output_file_name, L2_model_name, date_str,
                       X, Y, var_grid_timestep)

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
   

