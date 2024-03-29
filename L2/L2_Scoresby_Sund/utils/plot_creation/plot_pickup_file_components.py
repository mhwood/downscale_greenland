
import os
import numpy as np
import matplotlib.pyplot as plt
from MITgcmutils import mds
import argparse

def read_pickup_file(input_dir):

    Nr = 90
    global_data, _, global_metadata = mds.rdmds(os.path.join(input_dir, 'pickup.0000000001'), returnmeta=True)
    print(' Reading from '+os.path.join(input_dir, 'pickup.0000000001'))

    has_Nr = {'uvel': True, 'vvel': True, 'theta': True,
              'salt': True, 'gunm1': True, 'gvnm1': True,
              'gunm2': True, 'gvnm2': True, 'etan': False,
              'detahdt': False, 'etah': False}

    var_names = []
    row_bounds = []
    var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        if has_Nr[var_name.strip().lower()]:
            end_row = start_row + Nr
        else:
            end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]
        var_grids.append(var_grid)
        row_bounds.append([start_row,end_row])
        start_row=end_row
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

    var_names, row_bounds, var_grids, global_metadata = read_pickup_file(input_dir)
    print(var_names)

    fig = plt.figure(figsize=(12, 8))

    plt.style.use("dark_background")

    counter = 1

    for ff in range(len(var_names)):

        plt.subplot(3, 4, counter)

        print(' Plotting ' + var_names[ff])  # +' (global min: '+str(np.max(var_grids[ff][var_grids[ff]!=0]))+
        # ', min: '+str(np.min(var_grids[ff][var_grids[ff]!=0]))+')')

        if np.sum(np.isnan(  var_grids[ff]  ))>0:
            depths,rows,cols = np.where(np.isnan(var_grids[ff]))
            for i in range(20):
                print(depths[i],rows[i],cols[i])

        var_grid_subset = var_grids[ff][0, :, :]

        cmap = 'viridis'

        if var_names[ff] in ['Uvel', 'Vvel', 'GuNm1', 'GuNm2', 'GvNm1', 'GvNm2']:
            cmap = 'RdBu'
            if var_names[ff] in ['Uvel', 'Vvel']:
                vmin = -0.15
                vmax = 0.15
            else:
                vmin = -1.5e-5
                vmax = 1.5e-5
        else:
            if np.any(var_grids[ff][0, :, :])>0:
                vmin = np.min(var_grid_subset[var_grids[ff][0, :, :] != 0])
                vmax = np.max(var_grid_subset[var_grids[ff][0, :, :] != 0])
            else:
                vmin=-0.1
                vmax=0.1

        C = plt.imshow(var_grid_subset, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)

        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

        plt.colorbar(C)
        plt.title(var_names[ff])
        counter += 1

    plt.savefig(os.path.join(output_dir, 'pickup_file_components.png'), bbox_inches='tight')
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

