
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys

def read_bathy_to_faces(input_dir):

    bathy_file = os.path.join(input_dir,'bathymetry.bin')
    bathy_grid = np.fromfile(bathy_file, '>f4')

    llc = 1080

    grid_faces = {}

    # face 1
    face_1_grid = bathy_grid[:int(5*llc * llc/3)]
    face_1_grid = np.reshape(face_1_grid, (int(5 * llc/3), llc))
    grid_faces[1] = face_1_grid
    points_accumulated = int(5*llc * llc/3)

    # face 3
    face_3_grid = bathy_grid[points_accumulated:points_accumulated+int(llc * llc)]
    face_3_grid = np.reshape(face_3_grid, (llc, llc))
    grid_faces[3] = face_3_grid
    points_accumulated += int(llc * llc)

    # face 4
    face_4_grid = bathy_grid[points_accumulated:points_accumulated+int(2*llc * llc/3)]
    face_4_grid = np.reshape(face_4_grid, (llc, int(2 * llc/3)))
    grid_faces[4] = face_4_grid
    points_accumulated += int(2 * llc * llc / 3)

    # face 3
    face_5_grid = bathy_grid[points_accumulated:points_accumulated+int(5*llc * llc/3)]
    face_5_grid = np.reshape(face_5_grid, (llc, int(5 * llc/3)))
    grid_faces[5] = face_5_grid
    return(grid_faces)

def calculate_grid_corners(bathy_faces,sNx):

    counter = 0
    face_corner_dict = {}

    for face in [1,3,4,5]:
        print('Working on face '+str(face))
        print('    - bathy shape: '+str(np.shape(bathy_faces[face])))
        corner_list = []
        bathy_grid = bathy_faces[face]
        if np.shape(bathy_grid)[0] % sNx !=0:
            raise ValueError('The horizontal axis of face '+str(face)+ ' is not divisible by '+str(sNx))
        elif np.shape(bathy_grid)[1] % sNx != 0:
            raise ValueError('The vertical axis of face ' + str(face) + ' is not divisible by ' + str(sNx))
        else:
            row = 0
            for i in range((np.shape(bathy_grid)[0]//sNx)*(np.shape(bathy_grid)[1] // sNx)):
                counter+=1
                if i!=0 and i*sNx%np.shape(bathy_grid)[1]==0:
                    row+=1
                col = i-row*(np.shape(bathy_grid)[1] // sNx)
                corners = np.array([[col*sNx, row*sNx],
                                    [(col+1)*sNx, (row+1)*sNx]])
                bathy_subset = bathy_grid[row*sNx:(row+1)*sNx,
                                          col*sNx:(col+1)*sNx]
                if np.any(bathy_subset<0):
                    contains_ocean = True
                else:
                    contains_ocean = False
                corner_list.append([corners,counter,contains_ocean])

        face_corner_dict[face] = corner_list

    return(face_corner_dict)

def generate_blank_list(bathy_faces,sNx,ordered_nonblank_tiles):

    blank_list = []
    counter = 0

    for face in [1, 3, 4, 5]:
        bathy_grid = bathy_faces[face]
        if np.shape(bathy_grid)[0] % sNx != 0:
            raise ValueError('The horizontal axis of face ' + str(face) + ' is not divisible by ' + str(sNx))
        elif np.shape(bathy_grid)[1] % sNx != 0:
            raise ValueError('The vertical axis of face ' + str(face) + ' is not divisible by ' + str(sNx))
        else:
            row = 0
            for i in range((np.shape(bathy_grid)[0] // sNx) * (np.shape(bathy_grid)[1] // sNx)):
                counter += 1
                is_blank = True
                for ordered_row in ordered_nonblank_tiles:
                    if counter in ordered_row:
                        is_blank = False
                if is_blank:
                    blank_list.append(counter)
    return(blank_list)

def write_exch2_files(config_dir, config_name, bathy_faces, blank_list):
    exch2_output = ' &W2_EXCH2_PARM01\n'
    exch2_output += '  W2_mapIO   = 1,\n'
    exch2_output += '  preDefTopol = 0,\n'
    # exch2_output += '  dimsFacets = ' + str(np.shape(bathy_faces[1])[1]) + ', '
    # exch2_output += str(np.shape(bathy_faces[1])[0]) + ', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n'
    exch2_output += '  dimsFacets = '+str(np.shape(bathy_faces[1])[1])+', '
    exch2_output += str(np.shape(bathy_faces[1])[0]) + ', 0, 0, '
    exch2_output += str(np.shape(bathy_faces[3])[1]) + ', '
    exch2_output += str(np.shape(bathy_faces[3])[0]) + ', '
    exch2_output += str(np.shape(bathy_faces[4])[1]) + ', '
    exch2_output += str(np.shape(bathy_faces[4])[0]) + ', '
    exch2_output += str(np.shape(bathy_faces[5])[1]) + ', '
    exch2_output += str(np.shape(bathy_faces[5])[0]) + ', 0, 0,\n'
    exch2_output += '  facetEdgeLink(1,1)= 3.4, 0. , 0. , 5.1,\n'
    exch2_output += '  facetEdgeLink(1,2)= 0. , 0. , 0. , 0. ,\n'
    exch2_output += '  facetEdgeLink(1,3)= 5.4, 0. , 4.4, 1.1,\n'
    exch2_output += '  facetEdgeLink(1,4)= 0. , 0. , 0. , 3.3,\n'
    exch2_output += '  facetEdgeLink(1,5)= 1.4, 0. , 0. , 3.1,\n'
    # exch2_output += '  facetEdgeLink(1,1)= 0. , 0. , 0. , 0. ,\n'
    # exch2_output += '  facetEdgeLink(1,2)= 0. , 0. , 0. , 0. ,\n'
    # exch2_output += '  facetEdgeLink(1,3)= 0. , 0. , 0. , 0. ,\n'
    # exch2_output += '  facetEdgeLink(1,4)= 0. , 0. , 0. , 0. ,\n'
    # exch2_output += '  facetEdgeLink(1,5)= 0. , 0. , 0. , 0. ,\n'
    exch2_output += '#\n'
    exch2_output += '  blankList ='
    for counter in blank_list:
        exch2_output += '  '+str(counter)+',\n'
    exch2_output += ' &'

    output_file = os.path.join(config_dir,'L1',config_name,'namelist','data.exch2')
    f = open(output_file,'w')
    f.write(exch2_output)
    f.close()

def create_L1_CE_Greenland_geometry_files(config_dir, config_name):

    llc = 1080
    sNx = 180
    # ordered_nonblank_tiles = [[7,8,9],[10,11,12],[13,14,15]]
    ordered_nonblank_tiles = [[55,56,57],[91,85,79]]

    input_dir = os.path.join(config_dir, 'L1', 'input')

    bathy_faces = read_bathy_to_faces(input_dir)

    face_corner_dict = calculate_grid_corners(bathy_faces, sNx)

    blank_list = generate_blank_list(bathy_faces,sNx,ordered_nonblank_tiles)
    print(blank_list)

    write_exch2_files(config_dir, config_name, bathy_faces, blank_list)

    # # make the point grid in polar coordinates
    # min_x = 543736
    # max_x = 897457
    # new_max_x = max_x+(max_x-min_x)
    # max_y = -1873700
    # min_y = -2092319
    # new_max_y = max_y+2*(max_y-min_y)
    # new_min_y = min_y-2*(max_y-min_y)
    #
    # # define the corners
    # corners = np.array([[min_x, max_y],
    #                     [min_x, min_y],
    #                     [max_x, max_y],
    #                     [max_x, min_y]])
    #
    # # define corner names
    # corner_names = ['UL','LL','UR','LR']
    #
    # # reproject the corners to lat/lon
    # transformer = Transformer.from_crs('EPSG:' + str(3413), 'EPSG:' + str(4326))
    # corner_lats, corner_lons = transformer.transform(corners[:,0].ravel(),corners[:,1].ravel())
    # corners = np.column_stack([corner_lons,corner_lats])
    #
    # # find the locations of the corners in the face tiles
    # input_dir = os.path.join(config_dir,'L1','input')
    # corner_row_col_faces = np.zeros((4,5))
    # for c in range(np.shape(corners)[0]):
    #     min_dist = 1e22
    #     for f in [1,3,4,5]:
    #         XC, YC = read_XC_YC_from_mitgrid(input_dir,f)
    #         dist = ((XC-corners[c,0])**2 + (YC-corners[c,1])**2)**0.5
    #         # print(np.min(dist))
    #         if np.min(dist)<min_dist:
    #             min_dist = np.min(dist)
    #             row, col = np.where(dist==np.min(dist))
    #             row = row[0]
    #             col = col[0]
    #             face = f
    #             corner_row_col_faces[c,:] = [row,col,face,XC[row,col],YC[row,col]]
    #
    # # make arrays of rows, cols, lats, and lons
    # # left side first
    # print(corner_row_col_faces.astype(int))





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-c", "--config_name", action="store",
                        help="The name of the configuration.", dest="config_name",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    config_name = args.config_name

    create_L1_CE_Greenland_geometry_files(config_dir,config_name)
   

