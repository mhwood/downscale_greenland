
import os
import numpy as np
from scipy.interpolate import griddata
import argparse
import ast
import matplotlib.pyplot as plt

def read_domain_face_geometry_from_mitgrids(config_dir, model_name, face_size_dict):

    XC_faces = {}
    YC_faces = {}

    faces = list(face_size_dict.keys())

    for face in faces:
        grid_file = os.path.join(config_dir,'L1',model_name,'input','tile' + '{:03d}'.format(face) + '.mitgrid')
        grid = np.fromfile(grid_file, '>f8')
        grid = grid.reshape(16,face_size_dict[face][0]+1,face_size_dict[face][1]+1)
        XC_faces[face] = grid[0,:-1,:-1]
        YC_faces[face] = grid[1, :-1, :-1]

    return(faces, XC_faces, YC_faces)

def read_llc_face_geometry(ecco_dir, face, llc=1080):

    grid_file = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'mitgrid_tiles', 'tile' + '{:03d}'.format(face) + '.mitgrid')
    grid = np.fromfile(grid_file, '>f8')

    if face==1 or face==2:
        grid = grid.reshape(16,3*llc+1,llc+1)
    if face==3:
        grid = grid.reshape(16,llc+1,llc+1)
    if face==4 or face==5:
        grid = grid.reshape(16,llc+1,3*llc+1)

    XC = grid[0,:-1,:-1]
    YC = grid[1, :-1, :-1]

    return(XC, YC)

def read_llc_bathymetry_file_to_faces(grid_file,llc=1080):
    grid = np.fromfile(grid_file, '>f4')

    grid_faces = {}

    # face 1
    face_1_grid = grid[:3 * llc * llc]
    face_1_grid = np.reshape(face_1_grid, (3 * llc, llc))
    grid_faces[1] = face_1_grid

    # face 3
    face_3_grid = grid[(6) * llc * llc:(6+1) * llc * llc]
    face_3_grid = np.reshape(face_3_grid, (llc, llc))
    grid_faces[3] = face_3_grid

    # face 4
    face_4_grid = grid[(6+1) * llc * llc:(6+1+3) * llc * llc]
    face_4_grid = np.reshape(face_4_grid, (llc, 3 * llc))
    grid_faces[4] = face_4_grid

    # face 3
    face_5_grid = grid[(6+1+3) * llc * llc:(6+1+6) * llc * llc]
    face_5_grid = np.reshape(face_5_grid, (llc, 3 * llc))
    grid_faces[5] = face_5_grid

    return(grid_faces)

def generate_connected_mask(start_row, start_col, wet_grid):

    if wet_grid[start_row,start_col]==0:
        raise ValueError(' The start row/col location is  dry')

    rows = np.arange(np.shape(wet_grid)[0])
    cols = np.arange(np.shape(wet_grid)[1])
    Cols,Rows = np.meshgrid(cols,rows)

    mask_grid = 1-np.copy(wet_grid)
    mask_grid[start_row,start_col] = 2
    # in the mask, 0 means unverified
    # 1 is verified dry
    # 2 is verified wet

    # plt.imshow(mask_grid)
    # plt.show()

    is_remaining = np.logical_and(mask_grid==0,wet_grid==1)
    n_remaining = np.sum(is_remaining)
    # print(n_remaining)
    continue_iter = True
    for i in range(n_remaining):
        if continue_iter:
            # get the wet rows, cols, and their current mask values
            Wet_Rows = Rows[wet_grid == 1]
            Wet_Cols = Cols[wet_grid == 1]
            Mask_Vals = mask_grid[wet_grid == 1]

            # reduce these to the ones that havent been verified yet
            Wet_Rows = Wet_Rows[Mask_Vals == 0]
            Wet_Cols = Wet_Cols[Mask_Vals == 0]
            Mask_Vals = Mask_Vals[Mask_Vals == 0]

            if len(Mask_Vals)>0:

                # for each row/col, see if its connected to one we've verified is connected
                rows_remaining,cols_remaining = np.where(is_remaining)
                for ri in range(n_remaining):
                    row = rows_remaining[ri]
                    col = cols_remaining[ri]

                    # # this bit allows for diagonal spreading
                    # row_col_dist = ((Wet_Rows.astype(float)-row)**2 + (Wet_Cols.astype(float)-col)**2)**0.5
                    # closest_index = np.argmin(row_col_dist)
                    # if row_col_dist[closest_index]<np.sqrt(2):
                    #     var_grid[row,col] = Wet_Vals[closest_index]

                    # this bit allows for only up/dow/left/right spreading
                    if row<np.shape(wet_grid)[0]-1:
                        if mask_grid[row+1,col] == 2:
                            mask_grid[row,col] = 2
                    if row > 0:
                        if mask_grid[row - 1, col] == 2:
                            mask_grid[row,col] = 2
                    if col<np.shape(wet_grid)[1]-1:
                        if mask_grid[row,col+1] == 2:
                            mask_grid[row,col] = 2
                    if col > 0:
                        if mask_grid[row, col-1] == 2:
                            mask_grid[row,col] = 2


                is_remaining = np.logical_and(mask_grid == 0, wet_grid == 1)
                n_remaining_now = np.sum(is_remaining)

                # plt.subplot(1,2,1)
                # plt.imshow(wet_grid,cmap='Greys_r')
                # plt.subplot(1, 2, 2)
                # plt.imshow(mask_grid)
                # plt.show()

                if n_remaining_now<n_remaining:
                    n_remaining = n_remaining_now
                else:
                    n_remaining = n_remaining_now
                    continue_iter=False
            else:
                continue_iter = False

    return(mask_grid)

def fill_unconnected_areas(interpolated_bathy_faces, sNx, sNy):

    stitched_bathy = np.concatenate([interpolated_bathy_faces[1],np.rot90(interpolated_bathy_faces[3],k=3)],axis=0)

    # C = plt.imshow(stitched_bathy,origin='lower')
    # plt.colorbar(C)
    # plt.show()

    start_row = 150
    start_col = 300

    wet_grid = (stitched_bathy<0).astype(int)
    mask_grid = generate_connected_mask(start_row, start_col, wet_grid)

    # plt.imshow(mask_grid, origin='lower')
    # plt.show()

    stitched_bathy[mask_grid == 0] = 0

    output_faces = {}
    output_faces[1] = stitched_bathy[:sNy,:]
    output_faces[3] = np.rot90(stitched_bathy[sNy:,:])

    return(output_faces)


def create_L1_CE_Greenland_bathymetry(config_dir, model_name, ecco_dir, sNx, sNy, face_size_dict, llc, print_status):

    faces, XC_faces, YC_faces = read_domain_face_geometry_from_mitgrids(config_dir, model_name, face_size_dict)

    if print_status>=1:
        print('    - Reading bathymetry from the ECCO LLC'+str(llc)+' model')
    bathy_file = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','input_init', 'bathy_llc'+str(llc))
    bathy_faces = read_llc_bathymetry_file_to_faces(bathy_file,llc=llc)

    interpolated_bathy_faces = {}

    for face in faces:

        llc_face_XC, llc_face_YC = read_llc_face_geometry(ecco_dir, face)

        XC = XC_faces[face]
        YC = YC_faces[face]

        points = np.column_stack([llc_face_XC.ravel(),llc_face_YC.ravel()])
        values = bathy_faces[face].ravel()

        bathy_grid = griddata(points,values,(XC, YC), method='nearest')
        interpolated_bathy_faces[face] = bathy_grid

        # if not full_grid_started:
        #     full_grid = bathy_grid.reshape((np.size(bathy_grid), 1))
        # else:
        #     full_grid = np.vstack([full_grid, bathy_grid.reshape((np.size(bathy_grid), 1))])
        #
        # plt.subplot(2,2,1)
        # C = plt.imshow(XC,origin='lower')
        # plt.colorbar(C)
        # plt.title('XC')
        #
        # plt.subplot(2, 2, 2)
        # C = plt.imshow(YC, origin='lower')
        # plt.colorbar(C)
        # plt.title('YC')
        #
        # plt.subplot(2, 2, 3)
        # C = plt.imshow(bathy_grid, origin='lower')
        # plt.colorbar(C)
        # plt.title('Bathy (Face '+str(face)+')')
        # plt.show()

    interpolated_bathy_faces = fill_unconnected_areas(interpolated_bathy_faces, sNx, sNy)

    for f in range(len(faces)):
        face = faces[f]
        grid = interpolated_bathy_faces[face]
        n_rows = int((np.shape(grid)[0] * np.shape(grid)[1]) / sNx)
        grid = np.reshape(grid, (n_rows, sNx))
        if face == 1:
            compact_stack = grid
        else:
            compact_stack = np.concatenate([compact_stack, grid], axis=0)

    if print_status>=1:
        print('    - Outputting bathymetry as compact')
    output_file = os.path.join(config_dir, 'L1', model_name, 'input', model_name+'_bathymetry.bin')
    compact_stack.ravel(order='C').astype('>f4').tofile(output_file)


   

