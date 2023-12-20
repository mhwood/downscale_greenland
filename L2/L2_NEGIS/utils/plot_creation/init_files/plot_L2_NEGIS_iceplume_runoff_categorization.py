
import os
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
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
    transform = ds.GetGeoTransform()
    extents = [transform[0],transform[0]+transform[1]*np.shape(image)[1],transform[3]+ transform[5] * np.shape(image)[0], transform[3]]
    # x_resolution = transform[1]
    # y_resolution = transform[5]
    return(image,extents)

def read_land_ice_points(config_dir, model_name):
    land_file = os.path.join(config_dir,'L2',model_name,'input','iceplume','land_point_categorization.csv')
    land_points = np.genfromtxt(land_file, delimiter=',')

    ice_file = os.path.join(config_dir, 'L2', model_name, 'input', 'iceplume', 'ice_point_categorization.csv')
    ice_points = np.genfromtxt(ice_file, delimiter=',')

    return(land_points, ice_points)

def create_land_points_plot(output_path, config_dir,model_name,land_points):
    print('    - generate the MODIS file first')
    file_path = os.path.join(config_dir, 'L2',model_name, 'plots', model_name+'_MODIS_20220720_3413.tif')
    print(file_path)
    background_image, extents = read_background_imagery(file_path)

    fig = plt.figure(figsize=(12, 12))
    plt.style.use('dark_background')

    plt.imshow(background_image, extent=extents, alpha=0.7)

    plt.plot(land_points[:, 0], land_points[:, 1], 'k.',markersize=7)
    plt.plot(land_points[:, 0], land_points[:, 1], 'r.', markersize=5)

    # plt.gca().set_xticks([])
    # plt.gca().set_yticks([])

    # plt.title(model_name,fontsize=20)

    plt.savefig(output_path, bbox_inches='tight')
    plt.close(fig)


def create_ice_points_plot(output_path, config_dir, model_name,ice_points):

    print('    - generate the MODIS file first')
    file_path = os.path.join(config_dir, 'L2', model_name, 'plots', model_name + '_MODIS_20220720_3413.tif')
    print(file_path)
    background_image, extents = read_background_imagery(file_path)

    fig = plt.figure(figsize=(12, 12))
    plt.style.use('dark_background')

    plt.imshow(background_image, extent=extents, alpha=0.7)

    plt.plot(ice_points[:, 0], ice_points[:, 1], 'k.',markersize=7)
    plt.plot(ice_points[:, 0], ice_points[:, 1], 'y.', markersize=5)

    # plt.gca().set_xticks([])
    # plt.gca().set_yticks([])

    # plt.title(model_name,fontsize=20)

    plt.savefig(output_path, bbox_inches='tight')
    plt.close(fig)


def create_ice_fronts_plot(output_path, config_dir, model_name, fronts):
    print('    - generate the MODIS file first')
    file_path = os.path.join(config_dir, 'L2', model_name, 'plots', model_name + '_MODIS_20220720_3413.tif')
    print(file_path)
    background_image, extents = read_background_imagery(file_path)

    fig = plt.figure(figsize=(12, 12))
    plt.style.use('dark_background')

    # plt.imshow(background_image, extent=extents, alpha=0.7)

    for front in fronts:
        plt.plot(front[:, 0], front[:, 1], '-',color='yellow')

    # plt.gca().set_xticks([])
    # plt.gca().set_yticks([])

    # plt.title(model_name,fontsize=20)

    plt.savefig(output_path, bbox_inches='tight')
    plt.close(fig)


def create_land_points_categories_plot(output_path, config_dir, model_name,land_points):
    print('    - generate the MODIS file first')
    file_path = os.path.join(config_dir, 'L2', model_name, 'plots', model_name + '_MODIS_20220720_3413.tif')
    print(file_path)
    background_image, extents = read_background_imagery(file_path)

    fig = plt.figure(figsize=(12, 12))
    plt.style.use('dark_background')

    plt.imshow(background_image, extent=extents, alpha=0.7)

    land_categories = land_points[:,-1]

    no_glacier_points = land_points[land_categories==0,:]
    plt.plot(no_glacier_points[:, 0], no_glacier_points[:, 1], 'k.',markersize=7)
    plt.plot(no_glacier_points[:, 0], no_glacier_points[:, 1], 'r.', markersize=5)

    glacier_points = land_points[land_categories != 0, :]
    plt.plot(glacier_points[:, 0], glacier_points[:, 1], 'k.', markersize=10)
    plt.plot(glacier_points[:, 0], glacier_points[:, 1], '.', markersize=8, color='deepskyblue')

    # plt.gca().set_xticks([])
    # plt.gca().set_yticks([])

    # plt.title(model_name,fontsize=20)

    plt.savefig(output_path, bbox_inches='tight')
    plt.close(fig)


def create_ice_points_categories_plot(output_path, config_dir, model_name,ice_points):
    print('    - generate the MODIS file first')
    file_path = os.path.join(config_dir, 'L2', model_name, 'plots', model_name + '_MODIS_20220720_3413.tif')
    print(file_path)
    background_image, extents = read_background_imagery(file_path)

    fig = plt.figure(figsize=(12, 12))
    plt.style.use('dark_background')

    plt.imshow(background_image, extent=extents, alpha=0.7)

    ice_categories = ice_points[:,-1]

    no_glacier_points = ice_points[ice_categories == 0, :]
    plt.plot(no_glacier_points[:, 0], no_glacier_points[:, 1], 'k.', markersize=7)
    plt.plot(no_glacier_points[:, 0], no_glacier_points[:, 1], 'y.', markersize=5)

    glacier_points = ice_points[ice_categories != 0, :]
    plt.plot(glacier_points[:, 0], glacier_points[:, 1], 'k.', markersize=10)
    plt.plot(glacier_points[:, 0], glacier_points[:, 1], '.', markersize=8, color='deeppink')

    # plt.gca().set_xticks([])
    # plt.gca().set_yticks([])

    # plt.title(model_name,fontsize=20)

    plt.savefig(output_path, bbox_inches='tight')
    plt.close(fig)


def plot_runoff_field(output_path, config_dir, model_name, X, Y, year, year_day):

    if year%4==0:
        n_days = 366
    else:
        n_days = 365

    file_path = os.path.join(config_dir,'L2',model_name,'input','exf','L2_exf_RUNOFF_Mankoff_'+str(year))
    grid = np.fromfile(file_path,'>f4')
    grid = np.reshape(grid,(n_days,np.shape(X)[0], np.shape(Y)[1]))
    grid = grid[year_day, :, :]

    nonzero_indices = np.where(grid!=0)

    nonzero_values = 1

    # plot_grid = np.zeros_like(grid)
    # plot_grid[nonzero_indices] =
    # plot_grid = np.ma.masked_where(plot_grid>-0.01,plot_grid)

    fig = plt.figure(figsize=(10, 5))
    plt.style.use('dark_background')

    C = plt.scatter(X[nonzero_indices], Y[nonzero_indices], s=3, c=np.log10(grid[nonzero_indices]), cmap='turbo')
    cbar = plt.colorbar(C, ticks=[-12,-10,-8,-6,-4,-2])
    cbar.ax.set_yticklabels(['$10^{-12}$', '$10^{-10}$', '$10^{-8}$', '$10^{-6}$', '$10^{-4}$', '$10^{-2}$'])
    cbar.set_label('Runoff Averaged over Model Cell Area (m/s)')

    plt.gca().set_xlim([np.min(X), np.max(X)])
    plt.gca().set_ylim([np.min(Y), np.max(Y)])

    # C = plt.imshow(plot_grid,extent = [np.min(X),np.max(X),np.min(Y),np.max(Y)],cmap='turbo')
    # plt.colorbar(C, fraction=0.031, pad=0.04)

    plt.savefig(output_path, bbox_inches='tight')
    plt.close(fig)




def plot_L2_NEGIS_runoff(config_dir):

    sys.path.insert(1, os.path.join(config_dir,'L2','utils','init_file_creation'))
    import create_L2_iceplume_files as ip

    runoff_dir = '/Users/mhwood/Documents/Research/Data Repository/Greenland/Runoff'
    model_name = 'L2_NEGIS'

    termpicks_file = '/Users/mhwood/Documents/Research/Data Repository/Greenland/Ice Fronts/TermPicks_V1'
    glacier_IDs = [101, 102]

    print('    - Reading model geometry')
    X, Y, cell_area, Depth, model_boundary, model_boundary_3413 = ip.get_model_grid_boundary(config_dir, model_name)

    # print('    - Reading land and ice points')
    # land_column_names, land_points = ip.get_land_outlet_locations(runoff_dir, model_boundary_3413)
    # ice_column_names, ice_points = ip.get_ice_outlet_locations(runoff_dir, model_boundary_3413)

    land_points, ice_points = read_land_ice_points(config_dir, model_name)

    output_path = os.path.join(config_dir,'L2',model_name,'plots','init_files', 'discharge','Mankoff Land Points.png')
    create_land_points_plot(output_path,config_dir,model_name,land_points)

    output_path = os.path.join(config_dir,'L2',model_name,'plots','init_files', 'discharge','Mankoff Ice Points.png')
    create_ice_points_plot(output_path,config_dir,model_name,ice_points)

    # print('    - Getting fronts')
    # fronts, front_IDs = ip.get_ice_fronts_in_domain(termpicks_file, glacier_IDs)
    #
    # output_path = os.path.join(config_dir,'L2',model_name,'plots','init_files', 'discharge','TermPicks Ice Fronts.png')
    # create_ice_fronts_plot(output_path, config_dir, model_name, fronts)
    #
    # print('    - Categorizing the land points')
    # land_categories = ip.categorize_points(land_points, fronts, front_IDs, x_col=3, y_col=4)
    #
    # print('    - Categorizing the ice points')
    # ice_categories = ip.categorize_points(ice_points, fronts, front_IDs, x_col=9, y_col=10)

    output_path = os.path.join(config_dir,'L2',model_name,'plots','init_files', 'discharge','Mankoff Land Point Categories.png')
    create_land_points_categories_plot(output_path, config_dir, model_name, land_points)

    output_path = os.path.join(config_dir,'L2',model_name,'plots','init_files', 'discharge','Mankoff Ice Point Categories.png')
    create_ice_points_categories_plot(output_path, config_dir, model_name, ice_points)

    # # print('    - Plotting an example discharge timeseries')
    # # output_path = os.path.join(project_dir,'Figures','Ocean','Discharge Strategy','Daugaard Jensen Discharge.png')
    # # plot_subglacial_discharge_for_glacier(output_path, config_dir, model_name,glacierID=118,year=2000)
    #
    # print('    - Plotting an example runoff field')
    # output_path = os.path.join(config_dir,'L2',model_name,'plots','init_files', 'discharge','Runoff Field Example.png')
    # plot_runoff_field(output_path, config_dir, model_name, X, Y, year=2000, year_day=183)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_L2_NEGIS_runoff(config_dir)
   

