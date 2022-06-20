
import os
import numpy as np
import matplotlib.pyplot as plt
from MITgcmutils import mds
import argparse

def read_seaice_pickup_file(input_dir):

    seaice_temp_Nr = 7
    global_data, _, global_metadata = mds.rdmds(os.path.join(input_dir, 'pickup_seaice.0000000001'), returnmeta=True)
    print(' Reading from ' + os.path.join(input_dir, 'pickup_seaice.0000000001'))

    has_Nr = {'sitices': True, 'siarea': False, 'siheff': False,
              'sihsnow': False, 'siuice': False, 'sivice': False}

    var_names = []
    row_bounds = []
    var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        if has_Nr[var_name.strip().lower()]:
            end_row = start_row + seaice_temp_Nr
        else:
            end_row = start_row + 1
        var_grid = global_data[start_row:end_row, :, :]
        var_grids.append(var_grid)
        row_bounds.append([start_row, end_row])
        start_row = end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,var_grids,global_metadata)


def create_plot(config_dir):
    model_name = 'L3_Scoresby_Sund'

    if 'plots' not in os.listdir(os.path.join(config_dir,'L3',model_name)):
        os.mkdir(os.path.join(config_dir,'L3',model_name, 'plots'))
    if 'init_files' not in os.listdir(os.path.join(config_dir,'L3',model_name, 'plots')):
        os.mkdir(os.path.join(config_dir,'L3',model_name, 'plots', 'init_files'))

    input_dir = os.path.join(config_dir,'L3',model_name, 'input')
    output_dir = os.path.join(config_dir,'L3',model_name, 'plots', 'init_files')

    var_names, row_bounds, var_grids, global_metadata = read_seaice_pickup_file(input_dir)

    fig = plt.figure(figsize=(12, 8))

    plt.style.use("dark_background")

    counter = 1

    print(len(var_grids), len(var_names))

    for ff in range(len(var_names)):

        plt.subplot(2,3, counter)

        print(' Plotting ' + var_names[ff])  # +' (global min: '+str(np.max(var_grids[ff][var_grids[ff]!=0]))+
        # ', min: '+str(np.min(var_grids[ff][var_grids[ff]!=0]))+')')

        var_grid_subset = var_grids[ff][0, :, :]

        cmap = 'viridis'

        if var_names[ff] in ['siUICE','siVICE']:
            cmap = 'RdBu'
            vmin = -0.15
            vmax = 0.15
        else:
            vmin = np.min(var_grid_subset[var_grids[-3][0, :, :] != 0])
            vmax = np.max(var_grid_subset[var_grids[-3][0, :, :] != 0])

        C = plt.imshow(var_grid_subset, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)

        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

        plt.colorbar(C)
        plt.title(var_names[ff])
        counter += 1

    plt.savefig(os.path.join(output_dir, 'seaice_pickup_file_components.png'), bbox_inches='tight')
    plt.close(fig)


########################################################################################################################
# User inputs


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_plot(config_dir)

