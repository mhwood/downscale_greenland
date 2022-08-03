
import os
import numpy as np
import netCDF4 as nc4
import argparse
import matplotlib.pyplot as plt
import cmocean.cm as cm
from PIL import Image, ImageDraw


def read_geometry_from_grid_tiles_nc(config_dir,model_name):

    if model_name == 'L1_CE_Greenland':
        sNx = 180
        sNy = 180
        ordered_tiles = [[1,2,3],[6,5,4]]
        ordered_rotations = [[0,0,0],[3,3,3]]

    XC_grid = np.zeros((sNx*len(ordered_tiles),
                              sNx*len(ordered_tiles[1])))
    YC_grid = np.zeros((sNx * len(ordered_tiles),
                           sNx * len(ordered_tiles[1])))
    Depth_grid = np.zeros((sNx * len(ordered_tiles),
                           sNx * len(ordered_tiles[1])))

    for r in range(len(ordered_tiles)):
        for c in range(len(ordered_tiles[0])):
            tile_number = ordered_tiles[r][c]
            grid_path = os.path.join(config_dir, 'nc_grids', model_name+'_grid_tiles',
                                     model_name + '_grid.t'+'{:03d}'.format(tile_number)+'.nc')

            ds = nc4.Dataset(grid_path)
            Depth = ds.variables['Depth'][:, :]
            XC = ds.variables['XC'][:, :]
            YC = ds.variables['YC'][:, :]
            ds.close()

            for i in range(ordered_rotations[r][c]):
                XC = np.rot90(XC)
                YC = np.rot90(YC)
                Depth = np.rot90(Depth)

            XC_grid[r * sNx:(r + 1) * sNx, c * sNy:(c + 1) * sNy] = XC
            YC_grid[r * sNx:(r + 1) * sNx, c * sNy:(c + 1) * sNy] = YC
            Depth_grid[r * sNx:(r + 1) * sNx, c * sNy:(c + 1) * sNy] = Depth

    return(XC_grid, YC_grid, Depth_grid)

def read_geometry_from_grid_nc(config_dir,model_name):

    grid_path = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')

    ds = nc4.Dataset(grid_path)
    Depth = ds.variables['Depth'][:, :]
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    ds.close()

    return(XC, YC, Depth)

def generate_subdomain_plot(config_dir,model_name,
                            parent_XC,parent_YC,parent_Depth,
                            child_XC, child_YC):

    output_path = os.path.join(config_dir,'plots','bathymetry',model_name+'_bathymetry.png')

    fig = plt.figure(figsize=(12, 12))
    # plt.style.use('dark_background')

    plt.imshow(parent_Depth, cmap=cm.deep, vmin=0, vmax=3000, origin='lower')#,
                   # extent = [np.min(parent_XC), np.max(parent_XC), np.min(parent_YC), np.max(parent_YC)])

    # if len(child_XC)>0:
    #     plt.plot(child_XC[:, 0], child_YC[:, 0], 'r-', linewidth=2)
    #     plt.plot(child_XC[:, -1], child_YC[:, -1], 'r-', linewidth=2)
    #     plt.plot(child_XC[0, :], child_YC[0, :], 'r-', linewidth=2)
    #     plt.plot(child_XC[-1, :], child_YC[-1, :], 'r-', linewidth=2)

    plt.title(model_name,fontsize=20)

    plt.savefig(output_path)
    plt.close(fig)

def combine_panels_to_figure(config_dir, L0_model_name, L1_model_name, L2_model_name, L3_model_name):
    page_width = 1200 + 1200
    page_height = 1200 + 1200

    page = Image.new('RGB', (page_width, page_height), 'black')

    ###########################################
    # Calculate the dimensions

    L0_ulx = 0
    L0_uly = 0
    L1_ulx = 1200
    L1_uly = 0
    L2_ulx = 1200
    L2_uly = 1200
    L3_ulx = 0
    L3_uly = 1200

    # ###########################################
    # # Put all the plots onto one page
    #
    # colorbar_path = os.path.join(config_dir, 'plots', 'multi-level plots', 'panels',
    #                              var_name + '_Colorbar.png')
    # im = Image.open(colorbar_path)
    # page.paste(im, (50, 0))

    L0_file_path = os.path.join(config_dir,'plots','bathymetry',L0_model_name+'_bathymetry.png')
    im = Image.open(L0_file_path)
    page.paste(im, (L0_ulx, L0_uly))

    L1_file_path = os.path.join(config_dir,'plots','bathymetry',L1_model_name+'_bathymetry.png')
    im = Image.open(L1_file_path)
    page.paste(im, (L1_ulx, L1_uly))

    L2_file_path = os.path.join(config_dir,'plots','bathymetry',L2_model_name+'_bathymetry.png')
    im = Image.open(L2_file_path)
    page.paste(im, (L2_ulx, L2_uly))

    L3_file_path = os.path.join(config_dir,'plots','bathymetry',L3_model_name+'_bathymetry.png')
    im = Image.open(L3_file_path)
    page.paste(im, (L3_ulx, L3_uly))

    ###########################################
    # Plot lines to connect all the domains

    draw = ImageDraw.Draw(page)

    L1_ulx_in_L0 = 735
    L1_uly_in_L0 = 385
    L2_ulx_in_L1 = 508
    L2_uly_in_L1 = 657
    L3_urx_in_L2 = 538
    L3_ury_in_L2 = 278

    L1_ulx_in_L1 = 200
    L1_uly_in_L1 = 143
    L2_ulx_in_L2 = 200
    L2_uly_in_L2 = 143
    L3_urx_in_L3 = 1440
    L3_ury_in_L3 = 143

    draw.line([(L0_ulx + L1_ulx_in_L0, L0_uly + L1_uly_in_L0), (L1_ulx + L1_ulx_in_L1, L1_uly + L1_uly_in_L1)],
              fill='white', width=5)
    draw.line([(L1_ulx + L2_ulx_in_L1, L1_uly + L2_uly_in_L1), (L2_ulx + L2_ulx_in_L2, L2_uly + L2_uly_in_L2)],
              fill='white', width=5)
    draw.line([(L2_ulx + L3_urx_in_L2, L2_uly + L3_ury_in_L2), (L3_ulx + L3_urx_in_L3, L3_uly + L3_ury_in_L3)],
              fill='white', width=5)

    draw.line([(L0_ulx + L1_ulx_in_L0, L0_uly + L1_uly_in_L0), (L1_ulx + L1_ulx_in_L1, L1_uly + L1_uly_in_L1)],
              fill='black', width=3)
    draw.line([(L1_ulx + L2_ulx_in_L1, L1_uly + L2_uly_in_L1), (L2_ulx + L2_ulx_in_L2, L2_uly + L2_uly_in_L2)],
              fill='black', width=3)
    draw.line([(L2_ulx + L3_urx_in_L2, L2_uly + L3_ury_in_L2), (L3_ulx + L3_urx_in_L3, L3_uly + L3_ury_in_L3)],
              fill='black', width=3)

    L1_llx_in_L0 = 767
    L1_lly_in_L0 = 562
    L2_urx_in_L1 = 1128
    L2_ury_in_L1 = 661
    L3_lrx_in_L2 = 538
    L3_lry_in_L2 = 895

    L1_llx_in_L1 = 200
    L1_lly_in_L1 = 1068
    L2_urx_in_L2 = 1440
    L2_ury_in_L2 = 144
    L3_lrx_in_L3 = 1440
    L3_lry_in_L3 = 1068

    draw.line([(L0_ulx + L1_llx_in_L0, L0_uly + L1_lly_in_L0), (L1_ulx + L1_llx_in_L1, L1_uly + L1_lly_in_L1)],
              fill='white', width=5)
    draw.line([(L1_ulx + L2_urx_in_L1, L1_uly + L2_ury_in_L1), (L2_ulx + L2_urx_in_L2, L2_uly + L2_ury_in_L2)],
              fill='white', width=5)
    draw.line([(L2_ulx + L3_lrx_in_L2, L2_uly + L3_lry_in_L2), (L3_ulx + L3_lrx_in_L3, L3_uly + L3_lry_in_L3)],
              fill='white', width=5)

    draw.line([(L0_ulx + L1_llx_in_L0, L0_uly + L1_lly_in_L0), (L1_ulx + L1_llx_in_L1, L1_uly + L1_lly_in_L1)],
              fill='black', width=3)
    draw.line([(L1_ulx + L2_urx_in_L1, L1_uly + L2_ury_in_L1), (L2_ulx + L2_urx_in_L2, L2_uly + L2_ury_in_L2)],
              fill='black', width=3)
    draw.line([(L2_ulx + L3_lrx_in_L2, L2_uly + L3_lry_in_L2), (L3_ulx + L3_lrx_in_L3, L3_uly + L3_lry_in_L3)],
              fill='black', width=3)

    ###########################################
    # Output the figure

    output_path = os.path.join(config_dir,'plots','all_domain_bathymetry.png')
    page.save(output_path)

def create_nested_plot(config_dir, L1_model_name, L2_model_name, L3_model_name):

    L1_XC, L1_YC, L1_Depth = read_geometry_from_grid_tiles_nc(config_dir, L1_model_name)
    L2_XC, L2_YC, L2_Depth = read_geometry_from_grid_nc(config_dir, L2_model_name)
    L3_XC, L3_YC, L3_Depth = read_geometry_from_grid_nc(config_dir, L3_model_name)

    generate_subdomain_plot(config_dir, L1_model_name,
                            L1_XC, L1_YC, L1_Depth,
                            L2_XC, L2_YC)

    generate_subdomain_plot(config_dir, L2_model_name,
                            L2_XC, L2_YC, L2_Depth,
                            L3_XC, L3_YC)

    generate_subdomain_plot(config_dir, L3_model_name,
                            L3_XC, L3_YC, L3_Depth,
                            child_XC = [], child_YC = [])

    L0_model_name = 'L0_Global'
    combine_panels_to_figure(config_dir, L0_model_name, L1_model_name, L2_model_name, L3_model_name)







if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-L1", "--L1_model_name", action="store",
                        help="The name of the L1 model.", dest="L1_model_name",
                        type=str, required=True)

    parser.add_argument("-L2", "--L2_model_name", action="store",
                        help="The name of the L2 model.", dest="L2_model_name",
                        type=str, required=True)

    parser.add_argument("-L3", "--L3_model_name", action="store",
                        help="The name of the L3 model.", dest="L3_model_name",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    L1_model_name = args.L1_model_name
    L2_model_name = args.L2_model_name
    L3_model_name = args.L3_model_name

    create_nested_plot(config_dir, L1_model_name, L2_model_name, L3_model_name)