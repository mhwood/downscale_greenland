
import os
import netCDF4 as nc4
import numpy as np
import argparse
import sys
import cmocean.cm as cm

def create_movies(config_dir, field_name, print_level):

    L1_model_name = 'L1_N_Greenland'

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils','plot_creation','results'))

    TS_fields = ['Theta','Salt']

    TS_AW_fields = ['Theta_AW', 'Salt_AW']

    SI_fields = ['SIarea','SIheff','SIhsnow','SIuice','SIvice']

    etan_fields = ['EtaN']

    if field_name=='All':
        field_names = SI_fields + etan_fields + TS_fields + TS_AW_fields
    else:
        field_names=[field_name]

    metadata_dict = {'EtaN': [-4, 1, 'viridis', 'm', 'Surface Height Anomaly'],
                     'Theta': [-1.9, 5, 'turbo', '$^{\circ}$C', 'Potential Temperature (Surface)'],
                     'Theta_AW': [-1, 5, 'turbo', '$^{\circ}$C', 'Potential Temperature (257 m)'],#
                     'Salt': [31.5, 35, cm.haline, 'psu', 'Practical Salinity (Surface)'],  #
                     'Salt_AW': [32.5, 35, cm.haline, 'psu', 'Practical Salinity (257 m)'],  #
                     'Uvel': [-1, 1, cm.balance, 'm/s', 'Zonal Velocity'],  #
                     'Vvel': [-1, 1, cm.balance, 'm/s', 'Meridional Velocity'],  #
                     'Speed': [0, 0.25, cm.tempo_r, 'm/s', 'Speed (Surface)'],
                     'Speed_AW': [0, 0.25, cm.tempo_r, 'm/s', 'Speed (257 m)'],
                     'Vorticity': [-0.15, 0.15, cm.balance, '$\zeta$/f', 'Vorticity (Surface)'],
                     'Vorticity_AW': [-0.4, 0.4, cm.balance, '$\zeta$/f', 'Subsurface Vorticity (257 m)'],
                     'SIarea': [0, 1, cm.ice, 'm$^2$/m$^2$', 'Sea Ice Concentration'],
                     'SIheff': [0, 4, cm.ice, 'm', 'Sea Ice Thickness'],
                     'SIhsnow': [0, 1, cm.ice, 'm', 'Snow Thickness on Sea Ice'],
                     'SIuice': [-0.5, 0.5, cm.balance, 'm/s', 'Zonal Sea Ice Velocity'],
                     'SIvice': [-0.5, 0.5, cm.balance, 'm/s', 'Meridional Sea Ice Velocity']}

    remove_old = False
    skip = False

    # # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    import create_L1_variable_movies as pvm
    pvm.create_variable_movies(config_dir, L1_model_name, field_names, metadata_dict, remove_old, skip, print_level)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-f", "--field_name", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="field_name",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    field_name = args.field_name

    create_movies(config_dir, field_name, print_level=4)
   

