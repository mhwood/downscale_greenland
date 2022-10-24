
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
import shutil
import cmocean.cm as cm
import argparse
from matplotlib.patches import Rectangle
from datetime import datetime, timedelta
import sys


def read_field_from_monthly_ncs(results_dir, field_name):

    if field_name in ['SIarea','SIheff','SIhsnow','SIuice','SIvice']:
        file_prefix = 'SI_daily_snap'
    elif field_name in ['Theta','Salt']:
        file_prefix = 'TS_surf_daily_snap'
    elif field_name in ['Theta_AW','Salt_AW']:
        file_prefix = 'TS_AW_daily_snap'
    elif field_name in ['EXFaqh','EXFatemp','EXFpreci','EXFroff','EXFlwdn','EXFswdn','EXFuwind','EXFvwind']:
        file_prefix = 'EXF_day_snap'
    else:
        raise ValueError('Variable not yet implemented')

    grid_started = False

    year_months = []
    for file_name in os.listdir(os.path.join(results_dir, file_prefix)):
        if file_name[-3:]=='.nc' and file_name[0]!='.':
            year_month = file_name.split('.')[1]
            if year_month not in year_months:
                year_months.append(year_month)

    # year_months = ['199201']#,'199202','199203','199204','199205','199206',
    #                #'199207','199208','199209','199210','199211','199212']

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

def field_name_to_plot_metadata(field_name):
    metadata_dict = {'ETAN':[-4,1,'viridis','m','Surface Height Anomaly'],
                        'PHIBOT':[-200,600,cm.haline,'m$^2$/s$^2$','Bottom Pressure Potential Anomaly'],
                        'sIceLoad':[0,5000,cm.ice,'kg/m$^2$','Sea Ice Loading'],
                        'SIarea':[0,1,cm.ice,'m$^2$/m$^2$','Sea Ice Area'],
                        'SIheff':[0,3,cm.ice,'m','Effective Sea Ice Thickness'],
                        'SIhsnow':[0,1,cm.ice,'m','Effective Snow Thickness on Sea Ice'],
                        'ADVxHEFF':[-2000,2000,cm.balance,'m$^2\\cdot$(m/s)','Zonal Advective Flux of Effective Sea Ice Thickness'],
                        'ADVxSNOW':[-400,400,cm.balance,'m$^2\\cdot$(m/s)','Zonal Advective Flux of Effective Snow Thickness on Sea Ice'],
                        'ADVyHEFF':[-2000,2000,cm.balance,'m$^2\\cdot$(m/s)','Meridional Advective Flux of Effective Sea Ice Thickness'],
                        'ADVySNOW':[-400,400,cm.balance,'m$^2\\cdot$(m/s)','Meridional Advective Flux of Effective Snow Thickness on Sea Ice'],
                        'EXFaqh':[0,0.01,cm.tempo,'kg/kg','Surface (2m) Specific Humidity'],
                        'EXFatemp': [250, 290, cm.matter_r, 'Kelvin', 'Surface (2m) Air Temperature'],
                        'EXFempmr': [-1e-7,1e-7, cm.curl, 'm/s', 'Net Upward Freshwater Flux'],
                        'EXFevap': [-5e-8,5e-8, cm.curl, 'm/s', 'Evaporation'],
                        'EXFpreci': [0, 5e-7, cm.rain, 'm/s', 'Precipitation'],
                        'EXFroff': [0, 1e-7, cm.rain, 'm/s', 'River Runoff'],
                        'EXFqnet': [-300, 300, cm.curl, 'W/m$^2$', 'Net Upward Heat Flux'],
                        'EXFhl': [-200, 200, cm.curl, 'W/m$^2$', 'Latent Heat Flux Into the Ocean'],
                        'EXFhs': [-150, 150, cm.curl, 'W/m$^2$', 'Sensible Heat Flux Into the Ocean'],
                        'EXFlwdn': [150, 350, cm.solar, 'W/m$^2$', 'Downward Longwave Radiation'],
                        'EXFlwnet': [-300, 300, cm.curl, 'W/m$^2$', 'Net Upward Longwave Radiation'],
                        'EXFswdn': [30, 300, cm.solar, 'W/m$^2$', 'Downward Shotwave Radiation'],
                        'EXFswnet': [-300, 300, cm.curl, 'W/m$^2$', 'Net Upward Shortwave Radiation'],
                        'EXFuwind': [-15, 15, cm.balance, 'm/s', 'Zonal 10m Wind Speed'],
                        'EXFtaux': [-0.5, 0.5, cm.balance, 'N/m$^2$', 'Zonal Surface Wind Stress'],
                        'EXFvwind': [-15, 15, cm.balance, 'm/s', 'Meridional 10m Wind Speed'],
                        'EXFtauy': [-0.5, 0.5, cm.balance, 'N/m$^2$', 'Meridional Surface Wind Stress'],
                        'oceFWflx': [-0.0005, 0.0005, cm.diff, 'kg/m$^2$/s', 'Net Surface Freshwater Flux Into the Ocean'],
                        'oceQnet': [-2000, 2000, cm.curl, 'W/m$^2$', 'Net Surface Heat Flux Into the Ocean'],
                        'oceQsw': [-300, 300, cm.curl, 'W/m$^2$', 'Net Shortwave Radiation Into the Ocean'],
                        'oceTAUX': [-0.5, 0.5, cm.balance, 'N/m$^2$', 'Zonal Surface Wind Stress'],
                        'oceTAUY': [-0.5, 0.5, cm.balance, 'N/m$^2$', 'Meridional Surface Wind Stress'],
                        'SFLUX': [-0.02, 0.02, cm.tarn, 'g/m^2/s', 'Total Salt Flux'],
                        'TFLUX': [-500,500, cm.curl, 'W/m$^2$', 'Total Heat Flux'],
                        'SIatmFW': [-0.0001, 0.0001, cm.diff_r, 'kg/m^2/s', 'Net Freshwater Flux from Atmosphere and Land'],
                        'SIatmQnt': [-300, 300, cm.curl, 'W/m$^2$', 'Net Atmospheric Heat Flux'],
                        'SIuice': [-1, 1, cm.balance, 'm/s', 'Zonal Sea Ice Velocity'],
                        'SIvice': [-1, 1, cm.balance, 'm/s', 'Meridional Sea Ice Velocity'],
                        'KPPhbl': [0,500, cm.deep, 'm', 'KPP Boundary Layer Depth'],
                        'Theta': [-2, 10, cm.thermal, '$^{\circ}$C', 'Potential Temperature'], #
                        'Theta_AW': [-1, 5, cm.thermal, '$^{\circ}$C', 'Subsurface Potential Temperature'],  #
                        'Salt': [25, 35, cm.haline, 'psu', 'Practical Salinity'], #
                        'ADVr_SLT': [-1e5, 1e5, cm.delta, '(g/kg)$\\cdot$(m$^3$/s)', 'Vertical Advective Flux of Salinity'], #
                        'ADVr_TH': [-25000,25000, cm.curl, '$^{\circ}$C$\\cdot$(m$^3$/s)', 'Vertical Advective Flux of Potential Temperature'], #
                        'ADVx_SLT': [-1e6, 1e6, cm.delta, '(g/kg)$\\cdot$(m$^3$/s)', 'Zonal Advective Flux of Salinity'],
                        'ADVx_TH': [-50000, 50000, cm.curl, '$^{\circ}$C$\\cdot$(m$^3$/s)', 'Zonal Advective Flux of Potential Temperature'],
                        'ADVy_SLT': [-1e6, 1e6, cm.delta, '(g/kg)$\\cdot$(m$^3$/s)', 'Meridional Advective Flux of Salinity'],
                        'ADVy_TH': [-50000, 50000, cm.curl, '$^{\circ}$C$\\cdot$(m$^3$/s)', 'Meridional Advective Flux of Potential Temperature'],
                        'DFrE_SLT': [-1, 1, cm.delta, '(g/kg)$\\cdot$(m$^3$/s)', 'Vertical Diffusive Flux of Salinity (Explicit Part)'],
                        'DFrE_TH': [-1, 1, cm.curl, '$^{\circ}$C$\\cdot$(m$^3$/s)', 'Vertical Diffusive Flux of Potential Temperature (Explicit Part)'],
                        'DFrI_SLT': [-1000,1000, cm.delta, '(g/kg)$\\cdot$(m$^3$/s)', 'Vertical Diffusive Flux of Salinity (Implicit Part)'],
                        'DFrI_TH': [-2000,2000, cm.curl, '$^{\circ}$C$\\cdot$(m$^3$/s)', 'Vertical Diffusive Flux of Potential Temperature (Implicit Part)'],
                        'DFxE_SLT': [-1, 1, cm.delta, '(g/kg)$\\cdot$(m$^3$/s)', 'Zonal Diffusive Flux of Salinity'],
                        'DFxE_TH': [-1, 1, cm.curl, '$^{\circ}$C$\\cdot$(m$^3$/s)', 'Zonal Diffusive Flux of Potential Temperature'],
                        'DFyE_SLT': [-1, 1, cm.delta, '(g/kg)$\\cdot$(m$^3$/s)', 'Meridional Diffusive Flux of Salinity'],
                        'DFyE_TH': [-1, 1, cm.curl, '$^{\circ}$C$\\cdot$(m$^3$/s)', 'Meridional Diffusive Flux of Potential Temperature'],
                        'UVEL': [-1, 1, cm.balance, 'm/s', 'Zonal Velocity'], #
                        'VVEL': [-1, 1, cm.balance, 'm/s', 'Meridional Velocity'], #
                        'WVEL': [-5e-4, 5e-4, cm.balance, 'm/s', 'Vertical Velocity'], #
                        'UVELMASS': [-1, 1, cm.balance, 'm/s', 'Zonal Mass-Weighted Component of Velocity'], #
                        'VVELMASS': [-1, 1, cm.balance, 'm/s', 'Meridional Mass-Weighted Component of Velocity'], #
                        'WVELMASS': [-5e-4, 5e-4, cm.balance, 'm/s', 'Vertical Mass-Weighted Component of Velocity'], #
                        'PHIHYD': [-10, 10, cm.curl, 'm$^2$/s$^2$', 'Hydrostatic Pressure Potential Anomaly'],
                        'PHIHYDcR': [-10, 10, cm.curl, 'm$^2$/s$^2$', 'Hydrostatic Pressure Potential Anomaly at Constant r'],
                        'RHOAnoma': [-10, 10, cm.curl, 'kg/m$^3$', 'Density Anomaly'],
                        'KPPdiffS': [0, 0.5, cm.dense, 'm$^2$/s', 'KPP Vertical Salt Diffusion Coefficient'],
                        'KPPdiffT': [0, 0.5, cm.dense, 'm$^2$/s', 'KPP Vertical Heat Diffusion Coefficient'],
                        'KPPviscA': [0, 0.25, cm.dense, 'm$^2$/s', 'KPP Vertical Eddy Viscosity Coefficient']}

    return(metadata_dict[field_name])

def iter_number_to_date(iter_number):
    seconds_per_iter = 300
    total_seconds = iter_number*seconds_per_iter
    date = datetime(1992,1,1) + timedelta(seconds=total_seconds)
    return(date)

def date_to_iter_number(date):
    seconds_per_iter = 300
    total_seconds = (date-datetime(1992,1,1)).total_seconds()
    iter_number = total_seconds/seconds_per_iter
    return(iter_number)

def compile_panels_to_movie(config_dir,config_name,field_name):
    pwd = os.getcwd()

    panels_dir = os.path.join(config_dir,'L1_grid',config_name,'plots','output',field_name)

    file_name = 'L1_'+field_name+'.mp4'

    os.chdir(panels_dir)
    os.system("ffmpeg -r 5 -i L1_CE_Greenland_"+field_name+"_%04d.png -vcodec mpeg4 -b 1M -y " + file_name)
    os.rename(file_name, os.path.join('..', file_name))

    os.chdir(pwd)

def plot_mnc_fields(config_dir,field_name,remove_old,skip):
    config_name = 'L1_CE_Greenland'

    results_dir = os.path.join(config_dir, 'L1_grid', config_name, 'results')

    field_grid, all_iter_numbers = read_field_from_monthly_ncs(results_dir, field_name)

    output_dir = os.path.join(config_dir,'L1_grid',config_name,'plots','output')
    if field_name not in os.listdir(output_dir):
        os.mkdir(os.path.join(output_dir,field_name))
    output_dir = os.path.join(output_dir,field_name)

    min_index = 432
    max_index = min_index + (366*24*60*60)/300

    if remove_old:
        os.system('rm -rf '+output_dir+'/*')

    metadata = field_name_to_plot_metadata(field_name)
    vmin = metadata[0]
    vmax = metadata[1]
    cmap = metadata[2]

    panel_numbers = np.arange(0,np.shape(field_grid)[0],skip)

    counter = 0
    for i in panel_numbers:

        date = iter_number_to_date(all_iter_numbers[i][0])
        year = date.year
        min_iter = date_to_iter_number(datetime(date.year,1,1))
        max_iter = date_to_iter_number(datetime(date.year+1, 1, 1))

        file_name = config_name+'_'+field_name+'_'+'{:04d}'.format(counter)+'.png'

        if file_name not in os.listdir(output_dir):

            field_subset = field_grid[i,:,:]
            print('Timestep '+str(i)+' data range: '+str(np.min(field_subset[field_subset!=0]))+' to '+str(np.max(field_subset[field_subset!=0])))
            print('    - Iter number: '+str(all_iter_numbers[i][0])+', date: '+str(date))

            fig = plt.figure(figsize=(8,6))
            plt.style.use('dark_background')

            gs2 = GridSpec(15, 12, left=0.05, right=0.95, hspace=0.05)

            ax1 = fig.add_subplot(gs2[:-2,:])
            C = ax1.imshow(field_grid[i,:,:],origin='lower',vmin=vmin,vmax=vmax, cmap=cmap)
            plt.colorbar(C, fraction=0.031, pad=0.04)
            plt.title(field_name+', timestep = '+str(i))
            ax1.text(10, 350, 'Timestep: ' + str(int(all_iter_numbers[i][0])), ha='left', va='top', color='w')

            ax2 = fig.add_subplot(gs2[-1, 2:-1])
            width = (all_iter_numbers[i,0]-min_iter)/(max_iter-min_iter)
            rect = Rectangle((date.year,0),width,1,fc='silver',ec='white')
            ax2.add_patch(rect)
            ax2.set_xlim([date.year,date.year+1])
            ax2.set_ylim([0,1])
            ax2.set_xticks(np.arange(date.year,date.year+1,1/12))
            ax2.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
            ax2.set_yticks([])
            ax2.set_xlabel(str(year))

            output_path = os.path.join(output_dir,file_name)
            plt.savefig(output_path)
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
   

