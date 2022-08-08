
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import argparse
import cmocean.cm as cm
import ast
import sys


def read_L1_grid_geometry_to_faces(Lf, config_dir, model_name, sNx, sNy, faces, face_size_dict):

    L1_XC_Faces = {}
    L1_YC_Faces = {}
    for face in faces:
        tile_path = os.path.join(config_dir,'L1', model_name, 'input', 'tile'+'{:03d}'.format(face)+ '.mitgrid')
        tile_grid = np.fromfile(tile_path,'>f8')
        tile_grid = tile_grid.reshape((16,face_size_dict[face][0]+1,face_size_dict[face][1]+1))
        L1_XC_Faces[face] = tile_grid[0, :, :]
        L1_YC_Faces[face] = tile_grid[1, :, :]

    bathy_file = os.path.join(config_dir, 'L1', model_name, 'input', 'bathymetry.bin')
    bathy_compact = np.fromfile(bathy_file, '>f4')
    n_rows = int(np.size(bathy_compact) / sNx)
    bathy_compact = np.reshape(bathy_compact, (n_rows, sNx))

    bathy_grid = Lf.read_compact_grid_to_stitched_grid(bathy_compact, sNx, sNy, faces, face_size_dict, dim=2)
    depth = -1 * bathy_grid
    return(L1_XC_Faces, L1_YC_Faces, depth)

def read_mask_to_stitched_grid(Lf, config_dir, model_name, mask_name, sNx, sNy, faces, face_size_dict):
    if mask_name!='surface':
        mask_file = os.path.join(config_dir, 'L1', model_name, 'input', 'dv', 'L2_'+mask_name+'_BC_mask.bin')
    else:
        mask_file = os.path.join(config_dir, 'L1', model_name, 'input', 'dv', 'L2_' + mask_name + '_mask.bin')
    mask_compact = np.fromfile(mask_file, '>f4')
    n_rows = int(np.size(mask_compact) / sNx)
    mask_compact = np.reshape(mask_compact, (n_rows, sNx))

    mask_grid = Lf.read_compact_grid_to_stitched_grid(mask_compact, sNx, sNy, faces, face_size_dict, dim=2)
    mask_grid = np.ma.masked_where(mask_grid==0,mask_grid)

    return(mask_grid)

def read_L2_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    Depth = ds.variables['Depth'][:,:]
    ds.close()
    return(XC, YC, Depth)

def plot_dv_masks(config_dir, L1_model_name, L2_model_name, sNx, sNy, faces, face_size_dict, print_level):
    sys.path.insert(1, os.path.join(config_dir, 'L1', L1_model_name, 'utils'))
    import L1_CE_Greenland_functions as Lf

    L1_XC_Faces, L1_YC_Faces, L1_Depth = read_L1_grid_geometry_to_faces(Lf, config_dir, L1_model_name, sNx, sNy, faces, face_size_dict)
    L2_XC, L2_YC, L2_Depth = read_L2_grid_geometry_from_nc(config_dir, L2_model_name)


    for boundary in ['south','east','north','surface']:
        mask_grid = read_mask_to_stitched_grid(Lf, config_dir, L1_model_name, boundary, sNx, sNy, faces, face_size_dict)
        mask_faces = Lf.read_stitched_grid_to_faces(mask_grid, sNx, sNy, dim=2)

        fig = plt.figure(figsize=(15, 7))
        plt.style.use('dark_background')

        plt.subplot(1, 2, 1)
        plt.contour(L1_Depth,levels=[1],colors='k',linewidths=0.25)
        plt.imshow(L1_Depth,origin='lower',cmap = cm.deep, alpha=0.5)#,vmin=vmin,vmax=vmax)
        C = plt.imshow(mask_grid,origin='lower',cmap='jet')
        plt.colorbar(C)
        plt.title(L1_model_name+' Domain\n'+boundary+' mask')
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

        plt.subplot(1, 2, 2)
        C = plt.pcolormesh(L2_XC, L2_YC, L2_Depth, cmap=cm.deep, shading='nearest', alpha=0.7)  # ,vmin=vmin,vmax=vmax)
        plt.contour(L2_XC, L2_YC, L2_Depth, levels=[1], colors='k', linewidths=0.25)
        plt.colorbar(C)
        for face in faces:
            mask_face = mask_faces[face]
            rows,cols = np.where(mask_face>0)
            x = L1_XC_Faces[face][rows, cols]
            y = L1_YC_Faces[face][rows, cols]
            c = mask_face[rows,cols]
            plt.scatter(x,y,s=3,c=c,vmin=1,vmax=np.max(mask_grid),cmap='jet')

        plt.title(L2_model_name + ' Domain')
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])

        output_file = os.path.join(config_dir, 'L1', L1_model_name, 'plots', 'init_files', L1_model_name+'_dv_mask_'+boundary+'.png')
        plt.savefig(output_file,bbox_inches = 'tight')
        plt.close(fig)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_plot(config_dir)


