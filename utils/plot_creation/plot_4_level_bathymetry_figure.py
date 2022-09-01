
import os
import numpy as np
import netCDF4 as nc4
import argparse
import matplotlib.pyplot as plt
import cmocean.cm as cm
from PIL import Image, ImageDraw
from osgeo import gdal
from pyproj import Transformer

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

def reproject_points(points,inputCRS,outputCRS,x_column=0,y_column=1):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(points[:, x_column], points[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(points[:, y_column], points[:, x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
        run_test = False
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon = np.copy(points)
    output_polygon[:, x_column] = x2
    output_polygon[:, y_column] = y2
    return output_polygon

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
    print(np.min(image),np.max(image))
    transform = ds.GetGeoTransform()
    extents = [transform[0],transform[0]+transform[1]*np.shape(image)[1],transform[3]+ transform[5] * np.shape(image)[0], transform[3]]
    # x_resolution = transform[1]
    # y_resolution = transform[5]
    return(image,extents)

def generate_L1_subdomain_plot(config_dir, model_name,
                               parent_XC,parent_YC,parent_Depth,
                               child_XC, child_YC):

    add_background_imagery = True
    if add_background_imagery:
        file_path = os.path.join(config_dir,'plots','basemap',model_name+'_MODIS_20220720_row_col.tif')
        print(file_path)
        background_image, extents = read_background_imagery(file_path)

    output_path = os.path.join(config_dir,'plots','bathymetry',model_name+'_bathymetry.png')

    fig = plt.figure(figsize=(12, 12))
    # plt.style.use('dark_background')

    plt.imshow(background_image, extent = extents,alpha=0.7)

    # points = np.column_stack([parent_XC.ravel(), parent_YC.ravel()])
    # points = reproject_points(points, inputCRS=4326, outputCRS=3413)
    # X = np.reshape(points[:, 0], np.shape(parent_XC))
    # Y = np.reshape(points[:, 1], np.shape(parent_YC))
    # masked_depth = np.ma.masked_where(parent_Depth<=0,parent_Depth)
    # plt.pcolormesh(X, Y, masked_depth, shading='nearest',cmap=cm.deep)

    masked_depth = np.ma.masked_where(parent_Depth <= 0, parent_Depth)
    plt.imshow(masked_depth, cmap=cm.deep, vmin=0, vmax=1500, origin='lower')#,
               # extent = [np.min(parent_XC), np.max(parent_XC), np.min(parent_YC), np.max(parent_YC)])

    plt.contour(masked_depth, levels=[500, 1000],colors = 'k',linewidths=0.4)
    # plt.contour(masked_depth, levels=[1500], colors='w', linewidths=0.4)
    plt.contour(parent_Depth, levels=[0], colors='k', linewidths=1)

    interior_shift = 0.5
    outline_color='r'
    # bottom
    plt.plot([extents[0], extents[1]], [extents[2] + interior_shift, extents[2] + interior_shift], 'k-', linewidth=4)
    plt.plot([extents[0], extents[1]], [extents[3] - interior_shift, extents[3] - interior_shift], 'k-', linewidth=4)
    plt.plot([extents[0] + interior_shift, extents[0] + interior_shift], [extents[2], extents[3]], 'k-', linewidth=4)
    plt.plot([extents[1] - interior_shift, extents[1] - interior_shift], [extents[2], extents[3]], 'k-', linewidth=4)

    plt.plot([extents[0], extents[1]], [extents[2] + interior_shift, extents[2] + interior_shift], '-', linewidth=2, color=outline_color)
    plt.plot([extents[0], extents[1]], [extents[3] - interior_shift, extents[3] - interior_shift], '-', linewidth=2, color=outline_color)
    plt.plot([extents[0] + interior_shift, extents[0] + interior_shift], [extents[2], extents[3]], '-', linewidth=2, color=outline_color)
    plt.plot([extents[1] - interior_shift, extents[1] - interior_shift], [extents[2], extents[3]], '-', linewidth=2, color=outline_color)


    left_line = np.zeros((np.shape(child_YC)[0],2))
    for i in range(np.shape(child_YC)[0]):
        dist = ((parent_XC-child_XC[i,0])**2 + (parent_YC-child_YC[i,0])**2)**0.5
        ll_row, ll_col = np.where(dist==np.min(dist))
        left_line[i,0] = ll_col[0]
        left_line[i,1] = ll_row[0]

    right_line = np.zeros((np.shape(child_YC)[0], 2))
    for i in range(np.shape(child_YC)[0]):
        dist = ((parent_XC - child_XC[i, -1]) ** 2 + (parent_YC - child_YC[i, -1]) ** 2) ** 0.5
        ll_row, ll_col = np.where(dist == np.min(dist))
        right_line[i, 0] = ll_col[0]
        right_line[i, 1] = ll_row[0]

    top_line = np.zeros((np.shape(child_YC)[1], 2))
    for i in range(np.shape(child_YC)[1]):
        dist = ((parent_XC - child_XC[-1, i]) ** 2 + (parent_YC - child_YC[-1, i]) ** 2) ** 0.5
        ll_row, ll_col = np.where(dist == np.min(dist))
        top_line[i, 0] = ll_col[0]
        top_line[i, 1] = ll_row[0]

    bottom_line = np.zeros((np.shape(child_YC)[1], 2))
    for i in range(np.shape(child_YC)[1]):
        dist = ((parent_XC - child_XC[0, i]) ** 2 + (parent_YC - child_YC[0, i]) ** 2) ** 0.5
        ll_row, ll_col = np.where(dist == np.min(dist))
        bottom_line[i, 0] = ll_col[0]
        bottom_line[i, 1] = ll_row[0]

    plt.plot(left_line[:, 0], left_line[:, 1], 'k-', linewidth=4)
    plt.plot(right_line[:, 0], right_line[:, 1], 'k-', linewidth=4)
    plt.plot(top_line[:, 0], top_line[:, 1], 'k-', linewidth=4)
    plt.plot(bottom_line[:, 0], bottom_line[:, 1], 'k-', linewidth=4)

    plt.plot(left_line[:, 0], left_line[:, 1], '-', linewidth=2, color='orange')
    plt.plot(right_line[:, 0], right_line[:, 1], '-', linewidth=2, color='orange')
    plt.plot(top_line[:, 0], top_line[:, 1], '-', linewidth=2, color='orange')
    plt.plot(bottom_line[:, 0], bottom_line[:, 1], '-', linewidth=2, color='orange')

    plt.gca().set_xticks([])
    plt.gca().set_yticks([])

    # plt.title(model_name,fontsize=20)

    plt.savefig(output_path, bbox_inches='tight')
    plt.close(fig)


def generate_subdomain_plot(config_dir, model_name,
                            parent_XC, parent_YC, parent_Depth,
                            child_XC, child_YC, interior_shift = 1000, outline_color = 'y'):

    add_background_imagery = True
    if add_background_imagery:
        file_path = os.path.join(config_dir,'plots','basemap',model_name+'_MODIS_20220720_3413.tif')
        background_image, extents = read_background_imagery(file_path)
        print(extents)

    output_path = os.path.join(config_dir,'plots','bathymetry',model_name+'_bathymetry.png')

    fig = plt.figure(figsize=(12, 12))
    # plt.style.use('dark_background')

    plt.imshow(background_image, extent = extents, alpha=0.7)

    points = np.column_stack([parent_XC.ravel(), parent_YC.ravel()])
    points = reproject_points(points, inputCRS=4326, outputCRS=3413)
    X = np.reshape(points[:, 0], np.shape(parent_XC))
    Y = np.reshape(points[:, 1], np.shape(parent_YC))
    masked_depth = np.ma.masked_where(parent_Depth<=0,parent_Depth)
    plt.pcolormesh(X, Y, masked_depth, shading='nearest',cmap=cm.deep, vmin=0, vmax = 1500)

    plt.contour(X, Y, parent_Depth, levels=[500,1000], colors='k', linewidths=0.4)
    plt.contour(X, Y, parent_Depth, levels=[0], colors='k', linewidths=1)

    # bottom
    plt.plot([extents[0],extents[1]],[extents[2]+interior_shift,extents[2]+interior_shift], 'k-', linewidth=4)
    plt.plot([extents[0], extents[1]], [extents[3] - interior_shift, extents[3] - interior_shift], 'k-', linewidth=4)
    plt.plot([extents[0]+ interior_shift, extents[0]+ interior_shift], [extents[2], extents[3]], 'k-', linewidth=4)
    plt.plot([extents[1] - interior_shift, extents[1] - interior_shift], [extents[2], extents[3]], 'k-', linewidth=4)

    plt.plot([extents[0], extents[1]], [extents[2] + interior_shift, extents[2] + interior_shift], '-', linewidth=2,color=outline_color)
    plt.plot([extents[0], extents[1]], [extents[3] - interior_shift, extents[3] - interior_shift], '-', linewidth=2,color=outline_color)
    plt.plot([extents[0] + interior_shift, extents[0] + interior_shift], [extents[2], extents[3]], '-', linewidth=2,color=outline_color)
    plt.plot([extents[1] - interior_shift, extents[1] - interior_shift], [extents[2], extents[3]], '-', linewidth=2,color= outline_color)


    # plt.imshow(parent_Depth, cmap=cm.deep, vmin=0, vmax=3000, origin='lower')#,
                   # extent = [np.min(parent_XC), np.max(parent_XC), np.min(parent_YC), np.max(parent_YC)])

    if len(child_XC)>0:
        points = np.column_stack([child_XC.ravel(), child_YC.ravel()])
        points = reproject_points(points, inputCRS=4326, outputCRS=3413)
        child_X = np.reshape(points[:, 0], np.shape(child_XC))
        child_Y = np.reshape(points[:, 1], np.shape(child_YC))

        # plt.plot(child_X[:, 0], child_Y[:, 0], 'k-', linewidth=4)
        plt.plot(child_X[:, -1], child_Y[:, -1], 'k-', linewidth=4)
        plt.plot(child_X[0, :], child_Y[0, :], 'k-', linewidth=4)
        plt.plot(child_X[-1, :], child_Y[-1, :], 'k-', linewidth=4)

        # plt.plot(child_X[:, 0], child_Y[:, 0], 'y-', linewidth=2)
        plt.plot(child_X[:, -1], child_Y[:, -1], 'y-', linewidth=2)
        plt.plot(child_X[0, :], child_Y[0, :], 'y-', linewidth=2)
        plt.plot(child_X[-1, :], child_Y[-1, :], 'y-', linewidth=2)

    plt.gca().set_xticks([])
    plt.gca().set_yticks([])

    plt.savefig(output_path,bbox_inches='tight')
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

def create_colorbar(output_file):

    vmin = 0
    vmax = 1500
    cmap = cm.deep
    units = 'm'
    aspect = 80

    fig = plt.figure(figsize=(12,4))

    step = (vmax-vmin)/100
    x = np.arange(vmin,vmax+step,step)
    y = np.arange(0,1.1,1)
    X, Y = np.meshgrid(x,y)

    # plt.contourf(x,y,Y,100,cmap=cmap)
    plt.imshow(X,origin='lower',cmap=cmap,extent = [vmin,vmax,0,1],aspect=aspect)
    plt.gca().set_yticklabels([])

    plt.plot([500,500],[0,1],'k-',linewidth=0.6)
    plt.plot([1000, 1000], [0, 1], 'k-', linewidth=0.6)

    # plt.gca().yaxis.tick_right()
    plt.xticks([0,500,1000,1500],fontsize=20)
    plt.gca().set_xticklabels(['0','500','1000','$\geq$1500'])
    plt.gca().set_yticks([])

    # if plot_anomaly:
    #     plt.title(smb_model+' SMB Anomaly',fontsize=20)
    # else:
    #     plt.title(smb_model + ' SMB', fontsize=20)

    plt.title('Depth (' + units + ')', fontsize=20)
    plt.savefig(output_file,bbox_inches='tight')
    plt.close(fig)

def create_nested_plot(config_dir, L1_model_name, L2_model_name, L3_model_name):

    # L1_XC, L1_YC, L1_Depth = read_geometry_from_grid_tiles_nc(config_dir, L1_model_name)
    L1_XC, L1_YC, L1_Depth = read_geometry_from_grid_nc(config_dir, L1_model_name)
    L2_XC, L2_YC, L2_Depth = read_geometry_from_grid_nc(config_dir, L2_model_name)
    L3_XC, L3_YC, L3_Depth = read_geometry_from_grid_nc(config_dir, L3_model_name)

    # generate_L1_subdomain_plot(config_dir, L1_model_name,
    #                         L1_XC, L1_YC, L1_Depth,
    #                         L2_XC, L2_YC)
    #
    # generate_subdomain_plot(config_dir, L2_model_name,
    #                         L2_XC, L2_YC, L2_Depth,
    #                         L3_XC, L3_YC, interior_shift = 1500,
    #                         outline_color = 'orange')
    #
    # generate_subdomain_plot(config_dir,L3_model_name,
    #                         L3_XC, L3_YC, L3_Depth,
    #                         child_XC = [], child_YC = [])


    colorbar_file = os.path.join(config_dir, 'plots', 'bathymetry', 'bathymetry_colorbar.png')
    create_colorbar(colorbar_file)

    # L0_model_name = 'L0_Global'
    # combine_panels_to_figure(config_dir, L0_model_name, L1_model_name, L2_model_name, L3_model_name)







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