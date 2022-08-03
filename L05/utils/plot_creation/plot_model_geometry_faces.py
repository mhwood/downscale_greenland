
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from matplotlib.gridspec import GridSpec
import argparse
import sys

def read_domain_face_geometry_from_mitgrids(config_dir, model_name, faces, size_dict):

    XC_faces = {}
    YC_faces = {}

    for face in faces:
        grid_file = os.path.join(config_dir,'L1',model_name,'input','tile' + '{:03d}'.format(face) + '.mitgrid')
        grid = np.fromfile(grid_file, '>f8')
        grid = grid.reshape(16,size_dict[face][0]+1,size_dict[face][1]+1)
        XC_faces[face] = grid[0,:-1,:-1]
        YC_faces[face] = grid[1, :-1, :-1]

    return(XC_faces, YC_faces)

def read_bathy_from_binary(config_dir,model_name,faces,size_dict):

    bathy_file = os.path.join(config_dir,'L1',model_name,'input','bathymetry.bin')
    bathy_grid = np.fromfile(bathy_file,'>f4')

    points_counted = 0
    bathy_faces = {}
    for face in faces:
        n_points = size_dict[face][0]*size_dict[face][1]
        bathy_face = bathy_grid[points_counted:points_counted+n_points]
        bathy_face = bathy_face.reshape((size_dict[face][0], size_dict[face][1]))
        bathy_faces[face] = bathy_face

    return(bathy_faces)

def plot_geometry_faces(config_dir, model_name, faces, size_dict):

    llc = 1080
    sNx = 180

    XC_faces, YC_faces = read_domain_face_geometry_from_mitgrids(config_dir, model_name, faces, size_dict)

    bathy_faces = read_bathy_from_binary(config_dir,model_name,faces,size_dict)

    face_sets = [XC_faces, YC_faces, bathy_faces]
    face_names = ['XC','YC','Bathymetry']

    for s in range(len(face_sets)):

        field_faces = face_sets[s]
        field_name = face_names[s]

        fig = plt.figure(figsize=(8,8))
        plt.style.use('dark_background')

        gs = GridSpec(4,4,left = 0.05, right=0.95)

        for face in range(1,7):
            if face in faces:

                if face==1:
                    ax1 = fig.add_subplot(gs[2:,0])
                    C = ax1.imshow(field_faces[1][:,:],origin='lower')#,
                    plt.colorbar(C)

                if face==3:
                    ax3 = fig.add_subplot(gs[1, 1])
                    ax3.imshow(field_faces[3][:,:],origin='lower')#,
                    plt.colorbar(C)

                if face == 4:
                    ax4 = fig.add_subplot(gs[1, 2])
                    ax4.imshow(field_faces[4][:,:],origin='lower')#,
                    plt.colorbar(C)

                if face == 5:
                    ax5 = fig.add_subplot(gs[0, 2:])
                    ax5.imshow(field_faces[5][:,:],origin='lower')#,
                    plt.colorbar(C)


        output_path = os.path.join(config_dir,'L1',model_name,'plots',model_name+'_pickup_'+field_name+'.png')
        plt.savefig(output_path)
        plt.close(fig)




