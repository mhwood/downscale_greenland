
import os
import numpy as np
import argparse


def create_mitgrids(config_dir, config_name, ecco_dir,
                    sNx,sNy,ordered_nonblank_tiles,tile_face_index_dict):

    llc = 1080

    face_dimensions = {}
    face_index_extents = {}

    for face in range(1,6):
        min_col = 1e22
        max_col = 0
        min_row = 1e22
        max_row = 0
        face_found = False

        for tile_row in ordered_nonblank_tiles:
            for tile in tile_row:
                f = tile_face_index_dict[tile][0]
                if face==f:
                    face_found = True
                    col = tile_face_index_dict[tile][2]
                    row = tile_face_index_dict[tile][1]
                    if col<min_col:
                        min_col = col
                    if col + sNx > max_col:
                        max_col = col + sNx
                    if row<min_row:
                        min_row = row
                    if row + sNy > max_row:
                        max_row = row + sNy

        if face_found:
            face_dimensions[face] = (max_row-min_row,max_col-min_col)
            face_index_extents[face] = [min_col,max_col,min_row,max_row]

    for face in range(1,6):
        if face in face_dimensions.keys():

            grid_file = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','mitgrid_tiles', 'tile' + '{:03d}'.format(face) + '.mitgrid')
            grid = np.fromfile(grid_file, '>f8')

            output_file = os.path.join(config_dir,'L1',config_name,'input','tile' + '{:03d}'.format(face) + '.mitgrid')

            extents = face_index_extents[face]

            if face==1:
                grid = np.reshape(grid, (16, 3*llc + 1, llc + 1))

            if face==2:
                grid = np.reshape(grid, (16, 3*llc + 1, llc + 1))

            if face==3:
                grid = np.reshape(grid, (16, llc + 1, llc + 1))

            if face==4:
                grid = np.reshape(grid, (16, llc + 1, 3*llc + 1))

            if face==5:
                grid = np.reshape(grid, (16, llc + 1, 3*llc + 1))

            print('  - Creating face '+str(face))
            print('      - Reading in rows '+str(extents[2])+' to '+str(extents[3]+1)+' ('+str(extents[3]-extents[2]+1)+' total)')
            print('      - Reading in cols ' + str(extents[0]) + ' to ' + str(extents[1]+1) + ' (' + str(extents[1] - extents[0]+1) + ' total)')
            grid = grid[:, extents[2]:extents[3]+1, extents[0]:extents[1]+1]

            grid.ravel(order='C').astype('>f8').tofile(output_file)

    dict_file = os.path.join(config_dir, 'L1', config_name, 'namelist', 'face_dimensions.txt')
    f = open(dict_file,'w')
    f.write(str(face_dimensions))
    f.close()





