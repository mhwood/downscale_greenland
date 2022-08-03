import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
import shutil
import argparse
import sys
import cmocean.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def read_aste_grid(config_dir,var_name, ordered_aste_tiles, ordered_aste_tile_rotations):
    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import aste_functions as af

    aste_dir = '/Users/michwood/Documents/Research/Data Repository/Greenland/Ocean Properties/Models/ASTE'
    aste_XC, aste_YC, aste_AngleCS, aste_AngleSN, aste_hfacC, aste_delR = \
        af.read_aste_grid_geometry(aste_dir, ordered_aste_tiles, ordered_aste_tile_rotations)

    if var_name in ['UVEL', 'VVEL', 'SIuice', 'SIvice']:
        aste_grid = af.read_aste_field_to_stiched_grid(aste_dir, var_name, ordered_aste_tiles,
                                                       ordered_aste_tile_rotations,
                                                       aste_AngleCS, aste_AngleSN, rotate_velocity=True)
    else:
        aste_grid = af.read_aste_field_to_stiched_grid(aste_dir, var_name, ordered_aste_tiles,
                                                       ordered_aste_tile_rotations)
    aste_grid = aste_grid[0, 0, :, :]
    return(aste_XC, aste_YC, aste_grid)

def read_grid_tile_geometry(config_dir,model_name,sNx, sNy,
                            ordered_nonblank_tiles,ordered_nonblack_rotations):

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])
    XC_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    YC_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    AngleCS_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    AngleSN_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    grid_dir = os.path.join(config_dir, 'L05', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
                    AngleCS = ds.variables['AngleCS'][:, :]
                    AngleSN = ds.variables['AngleSN'][:, :]
                    ds.close()

                    for i in range(ordered_nonblack_rotations[r][c]):
                        XC = np.rot90(XC)
                        YC = np.rot90(YC)
                        AngleCS = np.rot90(AngleCS)
                        AngleSN = np.rot90(AngleSN)

                    XC_grid[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = XC
                    YC_grid[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = YC
                    AngleCS_grid[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = AngleCS
                    AngleSN_grid[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = AngleSN

    return(XC_grid, YC_grid, AngleCS_grid, AngleSN_grid)

def rotate_velocity_vectors_to_natural(angle_cos, angle_sin, uvel, vvel):
    zonal_velocity = np.zeros_like(uvel)
    meridional_velocity = np.zeros_like(vvel)
    for k in range(np.shape(uvel)[0]):
        zonal_velocity[k,:,:] = angle_cos * uvel[k,:,:] - angle_sin * vvel[k,:,:]
        meridional_velocity[k,:,:] = angle_sin * uvel[k,:,:] + angle_cos * vvel[k,:,:]
    return (zonal_velocity, meridional_velocity)

def read_mnc_field(config_dir,model_name, field_name, sNx, sNy,
                   ordered_nonblank_tiles, ordered_nonblack_rotations):
    if field_name == 'ETAN':
        file_prefix = 'surfDiag'
    elif field_name in ['SIarea','SIheff','SIhsnow','SIuice', 'SIvice']:
        file_prefix = 'seaiceDiag'
    else:
        file_prefix = 'dynDiag'

    grid_started = False

    XC, YC, AngleCS, AngleSN = read_grid_tile_geometry(config_dir, model_name, sNx, sNy,
                            ordered_nonblank_tiles, ordered_nonblack_rotations)

    run_dir = os.path.join(config_dir, 'L05', model_name, 'run')

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            file_name = file_prefix + '.' + '{:010d}'.format(0) + '.t' + '{:03d}'.format(tile_number) + '.nc'
            for n in range(N):
                if file_name in os.listdir(os.path.join(run_dir, 'mnc_' + '{:04d}'.format(n + 1))):

                    if field_name in ['SPEED', 'VORTICITY','UVEL','VVEL']:
                        ds = nc4.Dataset(os.path.join(run_dir, 'mnc_' + '{:04d}'.format(n + 1), file_name))
                        uvel = ds.variables['UVEL'][:, :, :, :]
                        vvel = ds.variables['VVEL'][:, :, :, :]
                        ds.close()
                        uvel = uvel[:, 0, :, :-1]
                        vvel = vvel[:, 0, :-1, :]
                        if field_name=='SPEED':
                            field = (uvel ** 2 + vvel ** 2) ** 0.5
                        if field_name in ['UVEL','VVEL']:
                            zonal_uvel, meridional_vvel = rotate_velocity_vectors_to_natural(AngleCS[r*sNy:(r+1)*sNy,c*sNx:(c+1)*sNx],
                                                                                             AngleSN[r*sNy:(r+1)*sNy,c*sNx:(c+1)*sNx],
                                                                                             uvel, vvel)
                            if field_name=='UVEL':
                                field = zonal_uvel

                            if field_name=='VVEL':
                                field = meridional_vvel
                    elif field_name in ['SIuice','SIvice']:
                        ds = nc4.Dataset(os.path.join(run_dir, 'mnc_' + '{:04d}'.format(n + 1), file_name))
                        uvel = ds.variables['SIuice'][:, :, :, :]
                        vvel = ds.variables['SIvice'][:, :, :, :]
                        ds.close()
                        uvel = uvel[:, 0, :, :-1]
                        vvel = vvel[:, 0, :-1, :]
                        if field_name in ['SIuice','SIvice']:
                            zonal_uvel, meridional_vvel = rotate_velocity_vectors_to_natural(AngleCS[r*sNy:(r+1)*sNy,c*sNx:(c+1)*sNx],
                                                                                             AngleSN[r*sNy:(r+1)*sNy,c*sNx:(c+1)*sNx],
                                                                                             uvel, vvel)
                            if field_name=='SIuice':
                                field = zonal_uvel

                            if field_name=='SIvice':
                                field = meridional_vvel
                    else:
                        ds = nc4.Dataset(os.path.join(run_dir, 'mnc_' + '{:04d}'.format(n + 1), file_name))
                        field = ds.variables[field_name][:, :, :, :]
                        ds.close()
                        field = field[:, 0, :, :]

                    for i in range(ordered_nonblack_rotations[r][c]):
                        field = np.rot90(field, axes=(1, 2))

                    if not grid_started:
                        grid_started = True
                        grid = np.zeros((np.shape(field)[0],
                                         sNy * len(ordered_nonblank_tiles),
                                         sNx * len(ordered_nonblank_tiles[0])))

                    grid[:, r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = field
    grid = grid[0,:,:]
    return (XC, YC, grid)


def create_plot(config_dir, config_name, var_name,
                aste_XC, aste_YC, aste_grid,
                domain_XC, domain_YC, domain_grid):

    output_dir = os.path.join(config_dir, 'L05', config_name, 'plots', 'ICs')

    if var_name == 'THETA':
        vmin = -1.9
        vmax = 8
        cmap = cm.thermal
        units = 'C'
    if var_name == 'SALT':
        vmin = 33
        vmax = 35.2
        cmap = cm.haline
        units = 'psu'
    if var_name == 'UVEL' or var_name == 'VVEL':
        vmin = -0.5
        vmax = 0.5
        cmap = cm.balance
        units = 'm/s'
    if var_name == 'ETAN':
        vmin = -0.5
        vmax = 0.2
        cmap = 'viridis'
        units = 'm'
    if var_name == 'SPEED':
        vmin = 0
        vmax = 1
        cmap = cm.deep_r
        units = 'm/s'
    if var_name == 'SIarea':
        vmin = 0
        vmax = 1
        cmap = cm.ice
        units = '%'
    if var_name == 'SIheff':
        vmin = 0
        vmax = 2
        cmap = cm.ice
        units = 'm'
    if var_name == 'SIhsnow':
        vmin = 0
        vmax = 1
        cmap = cm.ice
        units = 'm'
    if var_name == 'SIuice' or var_name == 'SIvice':
        vmin = -1
        vmax = 1
        cmap = cm.balance
        units = 'm/s'

    projPC = ccrs.PlateCarree()
    lonW = -42#np.min(domain_XC)
    lonE = np.max(domain_XC)
    latS = np.min(domain_YC)
    latN = np.max(domain_YC)+1
    print('lon',lonW,lonE)
    print('lat', latS, latN)
    cLat = (latN + latS) / 2
    cLon = (lonW + lonE) / 2
    res = '50m'

    fig = plt.figure(figsize=(8, 6))
    plt.style.use('dark_background')

    projStr = ccrs.Stereographic(central_longitude=cLon, central_latitude=cLat)
    fig = plt.figure(figsize=(11, 8.5))
    ax = plt.subplot(1, 1, 1, projection=projStr)

    ax.pcolormesh(aste_XC, aste_YC, aste_grid, vmin=vmin, vmax=vmax, cmap=cmap,
                             transform=ccrs.PlateCarree())

    dataplot = ax.pcolormesh(domain_XC, domain_YC, domain_grid, vmin=vmin, vmax=vmax, cmap=cmap,
                             transform=ccrs.PlateCarree())
    cbar = plt.colorbar(dataplot, orientation='horizontal',fraction=0.054, pad=0.04)
    cbar.set_label(var_name+' ('+units+')')

    if var_name in ['UVEL','VVEL','SIuice','SIvice']:
        ax.plot(domain_XC[:, 0], domain_YC[:, 0], 'k--', linewidth=0.5, transform=ccrs.PlateCarree())
        ax.plot(domain_XC[0, :], domain_YC[0, :], 'k--',linewidth=0.5, transform=ccrs.PlateCarree())
        ax.plot(domain_XC[:, -1], domain_YC[:, -1], 'k--',linewidth=0.5, transform=ccrs.PlateCarree())
        ax.plot(domain_XC[-1, :], domain_YC[-1, :], 'k--',linewidth=0.5, transform=ccrs.PlateCarree())
    else:
        ax.plot(domain_XC[:, 0], domain_YC[:, 0], 'w--', linewidth=0.5, transform=ccrs.PlateCarree())
        ax.plot(domain_XC[0, :], domain_YC[0, :], 'w--', linewidth=0.5, transform=ccrs.PlateCarree())
        ax.plot(domain_XC[:, -1], domain_YC[:, -1], 'w--', linewidth=0.5, transform=ccrs.PlateCarree())
        ax.plot(domain_XC[-1, :], domain_YC[-1, :], 'w--', linewidth=0.5, transform=ccrs.PlateCarree())

    ax.set_title('Downscaled Domain: '+config_name+', Parent Model: ASTE')
    gl = ax.gridlines(draw_labels=True, linewidth=1, color='silver', alpha=0.5, linestyle='--')
    ax.set_extent([lonW, lonE, latS, latN], crs=projPC)
    # ax.coastlines(resolution=res, color='white')

    output_path = os.path.join(output_dir, config_name + '_' + var_name + '.png')
    plt.savefig(output_path)
    plt.close(fig)



def plot_solution_on_ASTE(config_dir, var_name):

    config_name = 'L05_CE_Greenland'

    ordered_aste_tiles = [[27, 5, 6], [27, 14, 11]]
    ordered_aste_tile_rotations = [[1, 0, 0], [2, 3, 3]]  # rotations are counter-clockwise
    aste_XC, aste_YC, aste_grid = read_aste_grid(config_dir,var_name, ordered_aste_tiles, ordered_aste_tile_rotations)

    # plt.imshow(aste_grid,origin='lower')
    # plt.show()

    # aste_XC = aste_YC = aste_grid = []

    sNx = sNy = 90
    ordered_nonblank_tiles = [[1, 2, 3], [6, 5, 4]]
    ordered_nonblack_rotations = [[0, 0, 0], [3, 3, 3]]

    domain_XC, domain_YC, domain_grid = read_mnc_field(config_dir,config_name, var_name, sNx, sNy,
                   ordered_nonblank_tiles, ordered_nonblack_rotations)

    create_plot(config_dir, config_name, var_name,
                aste_XC, aste_YC, aste_grid,
                domain_XC, domain_YC, domain_grid)

    # plt.pcolormesh(aste_XC, aste_YC, aste_grid)
    # plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L05, L05, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-f", "--field_name", action="store",
                        help="The name of the field to plot.", dest="field_name",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    field_name = args.field_name

    plot_solution_on_ASTE(config_dir, field_name)


