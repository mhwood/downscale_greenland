
import os
import numpy as np
import netCDF4 as nc4
import argparse
import matplotlib.pyplot as plt
import cmocean.cm as cm


def read_geometry_from_grid_nc(config_dir,level_name,model_name):

    grid_path = os.path.join(config_dir,level_name,model_name,'input',model_name+'_grid.nc')

    ds = nc4.Dataset(grid_path)
    Depth = ds.variables['Depth'][:, :]
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    ds.close()

    return(XC, YC, Depth)


def create_nested_plot(config_dir, L1_model_name, L3_model_name):

    output_path = os.path.join(config_dir,'plots','domains', L3_model_name+'_domain.png')

    L1_XC, L1_YC, L1_Depth = read_geometry_from_grid_nc(config_dir, 'L1', L1_model_name)
    L3_XC, L3_YC, L3_Depth = read_geometry_from_grid_nc(config_dir, 'L3', L3_model_name)

    fig = plt.figure(figsize=(8,6))
    plt.style.use('dark_background')

    plt.pcolormesh(L1_XC,L1_YC,L1_Depth,cmap=cm.deep, vmin=0, vmax=3000)

    plt.pcolormesh(L3_XC, L3_YC, L3_Depth, cmap=cm.deep, vmin=0, vmax=3000)

    plt.plot(L3_XC[:, 0], L3_YC[:, 0], 'r-', linewidth=2, label='L3')
    plt.plot(L3_XC[:, -1], L3_YC[:, -1], 'r-', linewidth=2)
    plt.plot(L3_XC[0, :], L3_YC[0, :], 'r-', linewidth=2)
    plt.plot(L3_XC[-1, :], L3_YC[-1, :], 'r-', linewidth=2)

    plt.savefig(output_path)
    plt.close(fig)
    a=1
    a=1






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-L1", "--L1_model_name", action="store",
                        help="The name of the L1 model.", dest="L1_model_name",
                        type=str, required=True)

    parser.add_argument("-L3", "--L3_model_name", action="store",
                        help="The name of the L3 model.", dest="L3_model_name",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    L1_model_name = args.L1_model_name
    L3_model_name = args.L3_model_name

    create_nested_plot(config_dir, L1_model_name, L3_model_name)