
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
import argparse
from pyproj import Transformer
import sys


def read_XC_YC_from_mitgrid(input_dir,face):

    llc = 1080

    if face==1:
        shp = (16,int(5*llc/3)+1,llc+1)
    if face==3:
        shp = (16,llc+1,llc+1)
    if face==4:
        shp = (16,llc+1,int(2*llc/3)+1)
    if face==5:
        shp = (16,llc+1,int(5*llc/3)+1)

    file_path = os.path.join(input_dir,'tile'+'{:03d}'.format(face)+'.mitgrid')
    grid = np.fromfile(file_path,'>f8')
    grid = np.reshape(grid,shp)
    XC = grid[0, :-1, :-1]
    YC = grid[1, :-1, :-1]
    return(XC, YC)

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

def generate_geometry_plot(config_dir,config_name,XC_faces,YC_faces,bathy_faces,face_corner_dict):

    fig = plt.figure(figsize=(8,8))
    plt.style.use('dark_background')

    gs = GridSpec(4,4,left = 0.05, right=0.95)

    ax1 = fig.add_subplot(gs[2:,0])
    ax1.imshow(bathy_faces[1]==0,cmap='Blues_r',origin='lower')#,
               # extent = [np.min(XC_faces[1]), np.max(XC_faces[1]), np.min(YC_faces[1]), np.max(YC_faces[1])])
    # ax1.set_xticks([])
    # ax1.set_yticks([])
    face_corners = face_corner_dict[1]
    for i in range(len(face_corners)):
        rect = Rectangle((face_corners[i][0][0, 0], face_corners[i][0][0, 1]),
                         face_corners[i][0][1, 0] - face_corners[i][0][0, 0],
                         face_corners[i][0][1, 1] - face_corners[i][0][0, 1],
                         fc='none', ec='k')
        ax1.add_patch(rect)
        if face_corners[i][2]:
            plt.text(np.mean([face_corners[i][0][1,0],face_corners[i][0][0,0]]),
                     np.mean([face_corners[i][0][1,1],face_corners[i][0][0,1]]),
                     str(face_corners[i][1]),ha='center',va='center',color='k')

    ax3 = fig.add_subplot(gs[1, 1])
    ax3.imshow(bathy_faces[3]==0,cmap='Blues_r',origin='lower')#,
               # extent = [np.min(XC_faces[3]), np.max(XC_faces[3]), np.min(YC_faces[3]), np.max(YC_faces[3])])
    # ax3.set_xticks([])
    # ax3.set_yticks([])
    face_corners = face_corner_dict[3]
    for i in range(len(face_corners)):
        rect = Rectangle((face_corners[i][0][0, 0], face_corners[i][0][0, 1]),
                         face_corners[i][0][1, 0] - face_corners[i][0][0, 0],
                         face_corners[i][0][1, 1] - face_corners[i][0][0, 1],
                         fc='none', ec='k')
        ax3.add_patch(rect)
        if face_corners[i][2]:
            plt.text(np.mean([face_corners[i][0][1, 0], face_corners[i][0][0, 0]]),
                     np.mean([face_corners[i][0][1, 1], face_corners[i][0][0, 1]]),
                     str(face_corners[i][1]), ha='center', va='center', color='k')

    ax4 = fig.add_subplot(gs[1, 2])
    ax4.imshow(bathy_faces[4]==0,cmap='Blues_r',origin='lower')#,
               # extent = [np.min(XC_faces[4]), np.max(XC_faces[4]), np.min(YC_faces[4]), np.max(YC_faces[4])])
    # ax4.set_xticks([])
    # ax4.set_yticks([])
    face_corners = face_corner_dict[4]
    for i in range(len(face_corners)):
        rect = Rectangle((face_corners[i][0][0, 0], face_corners[i][0][0, 1]),
                         face_corners[i][0][1, 0] - face_corners[i][0][0, 0],
                         face_corners[i][0][1, 1] - face_corners[i][0][0, 1],
                         fc='none', ec='k')
        ax4.add_patch(rect)
        if face_corners[i][2]:
            plt.text(np.mean([face_corners[i][0][1, 0], face_corners[i][0][0, 0]]),
                     np.mean([face_corners[i][0][1, 1], face_corners[i][0][0, 1]]),
                     str(face_corners[i][1]), ha='center', va='center', color='k')

    ax5 = fig.add_subplot(gs[0, 2:])
    ax5.imshow(bathy_faces[5]==0,cmap='Blues_r',origin='lower')#,
               # extent = [np.min(XC_faces[5]), np.max(XC_faces[5]), np.min(YC_faces[5]), np.max(YC_faces[5])])
    # ax5.set_xticks([])
    # ax5.set_yticks([])
    face_corners = face_corner_dict[5]
    for i in range(len(face_corners)):
        rect = Rectangle((face_corners[i][0][0, 0], face_corners[i][0][0, 1]),
                         face_corners[i][0][1, 0] - face_corners[i][0][0, 0],
                         face_corners[i][0][1, 1] - face_corners[i][0][0, 1],
                         fc='none', ec='k')
        ax5.add_patch(rect)
        if face_corners[i][2]:
            plt.text(np.mean([face_corners[i][0][1, 0], face_corners[i][0][0, 0]]),
                     np.mean([face_corners[i][0][1, 1], face_corners[i][0][0, 1]]),
                     str(face_corners[i][1]), ha='center', va='center', color='k')

    output_path = os.path.join(config_dir,'L1',config_name,'plots',config_name+'_domain.png')
    plt.savefig(output_path)
    plt.close(fig)

def create_L1_geometry_plot(config_dir, config_name):

    input_dir = os.path.join(config_dir, 'L1', 'input')

    # make the point grid in polar coordinates
    XC_faces = {}
    YC_faces = {}
    for f in [1,3,4,5]:
        XC, YC = read_XC_YC_from_mitgrid(input_dir,f)
        XC_faces[f] = XC
        YC_faces[f] = YC

    # bathy faces
    bathy_faces = read_bathy_to_faces(input_dir)

    # generate list of proc corners for each face
    sNx = 180
    face_corner_dict = calculate_grid_corners(bathy_faces,sNx)

    generate_geometry_plot(config_dir,config_name,XC_faces,YC_faces,bathy_faces,face_corner_dict)






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

    create_L1_geometry_plot(config_dir, config_name)
   

