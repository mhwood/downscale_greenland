
import os
import numpy as np
from scipy.interpolate import griddata
import argparse
import ast
import matplotlib.pyplot as plt

def read_domain_face_geometry_from_mitgrids(config_dir, model_name):

    faces_dim_file = os.path.join(config_dir,'L1',model_name,'namelist','face_dimensions.txt')
    f = open(faces_dim_file)
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)

    XC_faces = {}
    YC_faces = {}

    faces = list(size_dict.keys())

    for face in faces:
        grid_file = os.path.join(config_dir,'L1',model_name,'input','tile' + '{:03d}'.format(face) + '.mitgrid')
        grid = np.fromfile(grid_file, '>f8')
        grid = grid.reshape(16,size_dict[face][0]+1,size_dict[face][1]+1)
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

def create_bathymetry(config_dir, model_name, ecco_dir, sNx, sNy, llc, print_status):

    faces, XC_faces, YC_faces = read_domain_face_geometry_from_mitgrids(config_dir, model_name)

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
    output_file = os.path.join(config_dir, 'L1', model_name, 'input', 'bathymetry.bin')
    compact_stack.ravel(order='C').astype('>f4').tofile(output_file)


   

