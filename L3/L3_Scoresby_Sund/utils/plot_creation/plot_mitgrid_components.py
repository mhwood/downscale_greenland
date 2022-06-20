
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import ast

def get_mitgrid_components(file_path,n_rows,n_cols):
    grid_chunk = np.fromfile(file_path,dtype='>f8')
    grid_chunk = np.reshape(grid_chunk,(16,n_rows+1,n_cols+1))
    return(grid_chunk)

def create_plot(config_dir):

    model_name = 'L3_Scoresby_Sund'

    f = open(os.path.join(config_dir, 'domain_sizes.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)
    L_size = size_dict[model_name]
    n_rows = L_size[0]
    n_cols = L_size[1]

    grid_keys = ['XC', 'YC', 'DXF', 'DYF', 'RAC', 'XG', 'YG', 'DXV', 'DYU',
                 'RAZ', 'DXC', 'DYC', 'RAW', 'RAS', 'DXG', 'DYG']

    file_path = os.path.join(config_dir, 'mitgrids', model_name+'.mitgrid')
    entire_grid = np.fromfile(file_path, dtype='>f8')
    entire_grid = np.reshape(entire_grid, (16, n_rows + 1, n_cols + 1))

    fig = plt.figure(figsize=(12, 12))

    plt.style.use('dark_background')

    for i in range(16):
        subset_grid = entire_grid[i,:,:]
        field_name = grid_keys[i]

        if field_name in ['XC', 'YC', 'DXF', 'DYF', 'RAC']:
            subset_grid = subset_grid[:-1,:-1]
        if field_name in ['DXV']:
            subset_grid = subset_grid[:,1:-1]
        if field_name in ['DYU']:
            subset_grid = subset_grid[1:-1,:]
        if field_name in ['RAZ']:
            subset_grid = subset_grid[1:-1,1:-1]
        if field_name in ['DXC','RAW']:
            subset_grid = subset_grid[:-1,1:-1]
        if field_name in ['DYC','RAS']:
            subset_grid = subset_grid[1:-1,:-1]
        if field_name in ['DXG']:
            subset_grid = subset_grid[:,:-1]
        if field_name in ['DYG']:
            subset_grid = subset_grid[:-1,:]

        # vmin = np.min(entire_grid)
        # vmax = np.max(entire_grid)

        plt.subplot(4, 4, i+1)
        C = plt.imshow(subset_grid,origin='lower')#,vmin=vmin,vmax=vmax)
        plt.colorbar(C)
        plt.title(grid_keys[i])

        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

    plt.suptitle('mitgrid file components: '+field_name+' ('+str(n_rows)+' rows by '+str(n_cols)+' columns)')

    output_file = os.path.join(config_dir,'L3',model_name, 'plots', 'mitgrid_components.png')
    plt.savefig(output_file,bbox_inches = 'tight')
    plt.close(fig)




########################################################################################################################
# User inputs


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_plot(config_dir)


