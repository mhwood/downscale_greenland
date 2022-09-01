import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import simplegrid as sg
from scipy.interpolate import griddata
from ecco_v4_py.llc_array_conversion import llc_faces_to_compact, llc_compact_to_faces
import argparse
import ast



def read_parent_XC_YC(config_dir,model_name,faces,size_dict):

    XC_faces = {}
    YC_faces = {}
    Bathy_faces = {}

    for face in faces:
        grid_file = os.path.join(config_dir, 'L1', model_name, 'input', 'tile' + '{:03d}'.format(face) + '.mitgrid')
        grid = np.fromfile(grid_file, '>f8')
        grid = grid.reshape(16, size_dict[face][0] + 1, size_dict[face][1] + 1)
        XC_faces[face] = grid[0, :-1, :-1]
        YC_faces[face] = grid[1, :-1, :-1]


    bathy_file = os.path.join(config_dir, 'L1', model_name, 'input', model_name+'_bathymetry.bin')
    bathy_grid = np.fromfile(bathy_file,'>f4')

    points_counted = 0
    for face in faces:
        n_points = size_dict[face][0] * size_dict[face][1]
        bathy_face = bathy_grid[points_counted:points_counted + n_points]
        if face==3:
            bathy_face = bathy_face.reshape((size_dict[face][0], size_dict[face][1]))
            Bathy_faces[face] = bathy_face
        else:
            bathy_face = bathy_face.reshape((size_dict[face][0], size_dict[face][1]))
            Bathy_faces[face] = bathy_face
        points_counted += n_points

    # plt.subplot(2,3,1)
    # C = plt.imshow(XC_faces[1],origin='lower')
    # plt.colorbar(C)
    # plt.title('XC')
    #
    # plt.subplot(2, 3, 2)
    # C = plt.imshow(YC_faces[1], origin='lower')
    # plt.colorbar(C)
    # plt.title('YC')
    #
    # plt.subplot(2, 3, 3)
    # C = plt.imshow(Bathy_faces[1], origin='lower')
    # plt.colorbar(C)
    # plt.title('Bathy')
    #
    # plt.subplot(2, 3, 4)
    # C = plt.imshow(XC_faces[3], origin='lower')
    # plt.colorbar(C)
    #
    # plt.subplot(2, 3, 5)
    # C = plt.imshow(YC_faces[3], origin='lower')
    # plt.colorbar(C)
    #
    # plt.subplot(2, 3, 6)
    # C = plt.imshow(Bathy_faces[3], origin='lower')
    # plt.colorbar(C)
    # plt.title('Bathy')
    #
    # plt.show()

    return(XC_faces, YC_faces, Bathy_faces)

def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    ds.close()
    return(XC, YC)

def create_mask_faces(faces, resolution, XC_faces, YC_faces, Bathy_faces, XC_boundary, YC_boundary):
    mask_faces = {}
    mask_indices = []
    counter = 1
    for face in faces:
        mask_face = np.zeros_like(XC_faces[face])
        bathy_face = Bathy_faces[face]
        for i in range(len(XC_boundary)):
            dist = ((XC_faces[face]-XC_boundary[i])**2 + (YC_faces[face]-YC_boundary[i])**2)**0.5
            rows, cols = np.where(dist<resolution*np.sqrt(2))
            for i in range(len(rows)):
                if mask_face[rows[i],cols[i]]==0 and bathy_face[rows[i],cols[i]]<0:
                    mask_face[rows[i],cols[i]] = counter
                    mask_indices.append([face,rows[i],cols[i]])
                    counter +=1
        mask_faces[face] = mask_face
    mask_indices = np.array(mask_indices)
    return(mask_faces, mask_indices)

def great_circle_distance(lon_ref, lat_ref, Lon, Lat):
    earth_radius = 6371000
    lon_ref_radians = np.radians(lon_ref)
    lat_ref_radians = np.radians(lat_ref)
    lons_radians = np.radians(Lon)
    lats_radians = np.radians(Lat)
    lat_diff = lats_radians - lat_ref_radians
    lon_diff = lons_radians - lon_ref_radians
    d = np.sin(lat_diff * 0.5) ** 2 + np.cos(lat_ref_radians) * np.cos(lats_radians) * np.sin(lon_diff * 0.5) ** 2
    h = 2 * earth_radius * np.arcsin(np.sqrt(d))
    return(h)

def create_CTD_mask_faces(faces, XC_faces, YC_faces, Bathy_faces, ctd_points):
    mask_faces = {}
    mask_indices = []
    counter = 1

    for face in faces:
        mask_face = np.zeros_like(XC_faces[face])
        mask_faces[face] = mask_face

    for p in range(np.shape(ctd_points)[0]):
        dist_err = 1e22
        for face in faces:
            XC_face = XC_faces[face]
            YC_face = YC_faces[face]
            bathy_face = Bathy_faces[face]

            nonzero_XC = XC_face[bathy_face < 0]
            nonzero_YC = YC_face[bathy_face < 0]

            dist = great_circle_distance(ctd_points[p,0],ctd_points[p,1],nonzero_XC,nonzero_YC)
            index = np.argmin(dist)

            x = nonzero_XC[index]
            y = nonzero_YC[index]

            ctd_dist = great_circle_distance(ctd_points[p,0], ctd_points[p,1], x, y)

            if ctd_dist<dist_err:

                row,col = np.where(np.logical_and(XC_face==x,YC_face==y))
                ctd_row = row[0]
                ctd_col = col[0]
                ctd_face = face
                dist_err = ctd_dist

        print('    - CTD point '+str(p+1)+' ('+str(ctd_points[p,0])+', '+str(ctd_points[p,1])+') will be located in face '+str(ctd_face)+' at ('+str(ctd_row)+','+
              str(ctd_col)+') with a depth of '+str(Bathy_faces[ctd_face][ctd_row,ctd_col])+' (dist err = '+str(dist_err)+')')

        mask_faces[ctd_face][ctd_row,ctd_col] = counter
        mask_indices.append([face, ctd_row, ctd_col])
        counter +=1
        #
        # mask_faces[face] = mask_face

    mask_indices = np.array(mask_indices)
    return(mask_faces, mask_indices)

def output_mask_dictionary_to_nc(output_dir,output_file_name,all_mask_dicts,mask_names_list):
    if output_file_name in os.listdir(output_dir):
        os.remove(os.path.join(output_dir,output_file_name))

    ds = nc4.Dataset(os.path.join(output_dir,output_file_name),'w')

    for m in range(len(mask_names_list)):
        grp = ds.createGroup(mask_names_list[m])
        grp.createDimension('n_points', np.shape(all_mask_dicts[m])[0])
        var = grp.createVariable('source_faces', 'i4', ('n_points',))
        var[:] = all_mask_dicts[m][:, 0].astype(int)
        var = grp.createVariable('source_rows', 'i4', ('n_points',))
        var[:] = all_mask_dicts[m][:, 1].astype(int)
        var = grp.createVariable('source_cols', 'i4', ('n_points',))
        var[:] = all_mask_dicts[m][:, 2].astype(int)

    ds.close()


########################################################################################################################

def create_dv_masks(config_dir, L1_model_name, L2_model_name, L3_model_name,
                    sNx, sNy, face_size_dict, print_level):

    # if 'input' not in os.listdir('..'):
    #     os.mkdir(os.path.join('..','input'))
    if 'dv' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'input','dv'))

    ###############################################################################################
    # Read in the grids

    L1_size_dict = face_size_dict
    L1_faces = list(L1_size_dict.keys())

    if print_level>=1:
        print('    - Reading in the L1 face geometry')
    XC_faces, YC_faces, Bathy_faces = read_parent_XC_YC(config_dir,L1_model_name,L1_faces,L1_size_dict)

    # read the extended subset of the XC and YC grids
    L2_XC, L2_YC = read_grid_geometry_from_nc(config_dir, L2_model_name)
    if print_level>=1:
        print('    - The L2 subdomain has shape '+str(np.shape(L2_XC)))

    # read the extended subset of the XC and YC grids
    L3_XC, L3_YC = read_grid_geometry_from_nc(config_dir, L3_model_name)
    if print_level >= 1:
        print('    - The L3 subdomain has shape ' + str(np.shape(L3_XC)))

    ###############################################################################################
    # Create the masks

    resolution = 1/12

    mask_names_list = []
    all_mask_dicts = []

    ctd_points = np.array([[-21.070,69.369], # outer shelf
                           [-21.379,70.091], # shelf, fjord entrance
                           [-22.267,70.220], # points going into fjord
                           [-21.339,70.092],
                           [-23.873,70.246],
                           [-24.889,70.480],
                           [-24.888,71.094], # point before entering narrow part
                           [-25.8428,71.5095], # below here goes to DJG
                           [-27.0630,71.5476],
                           [-28.35017,71.93849], #DJG
                           [-27.1102,70.9003], #next 5 loop counter clockwise around loop
                           [-28.2028,70.5877], #Rolige
                           [-29.1010,70.3950], #Vestfjord
                           [-28.0938,70.4035],
                           [-26.9400,70.4717],
                           [-25.9707,70.3937], # last two are in the southern-most part of the fjord
                           [-27.3367,70.1116]])

    for boundary in ['CTD','L2_north','L2_south','L2_east','L2_surface',
                     'L3_north','L3_south','L3_east','L3_surface']:
        if print_level >= 1:
            print('    - Creating the ' + boundary +' mask')

        if boundary == 'L2_south':
            XC_boundary = L2_XC[0,:]
            YC_boundary = L2_YC[0,:]

        if boundary == 'L2_north':
            XC_boundary = L2_XC[-1,:]
            YC_boundary = L2_YC[-1,:]

        if boundary == 'L2_east':
            XC_boundary = L2_XC[:,-1]
            YC_boundary = L2_YC[:,-1]

        if boundary == 'L2_west':
            XC_boundary = L2_XC[:,0]
            YC_boundary = L2_YC[:,0]

        if boundary == 'L2_surface':
            XC_boundary = L2_XC.ravel()
            YC_boundary = L2_YC.ravel()

        if boundary == 'L3_south':
            XC_boundary = L3_XC[0,:]
            YC_boundary = L3_YC[0,:]

        if boundary == 'L3_north':
            XC_boundary = L3_XC[-1,:]
            YC_boundary = L3_YC[-1,:]

        if boundary == 'L3_east':
            XC_boundary = L3_XC[:,-1]
            YC_boundary = L3_YC[:,-1]

        if boundary == 'L3_west':
            XC_boundary = L3_XC[:,0]
            YC_boundary = L3_YC[:,0]

        if boundary == 'L3_surface':
            XC_boundary = L3_XC.ravel()
            YC_boundary = L3_YC.ravel()


        if boundary!='CTD':
            mask_faces, mask_indices = create_mask_faces(L1_faces, resolution, XC_faces, YC_faces, Bathy_faces, XC_boundary, YC_boundary)
        else:
            mask_faces, mask_indices = create_CTD_mask_faces(L1_faces, XC_faces, YC_faces, Bathy_faces,
                                                            ctd_points)
        all_mask_dicts.append(mask_indices)
        mask_names_list.append(boundary)

        if print_level>=2:
            print('        - The '+boundary+' mask has '+str(np.shape(mask_indices)[0])+' points')

        # plt.subplot(1, 2, 1)
        # C = plt.imshow(mask_faces[1],origin='lower')
        # plt.colorbar(C)
        # plt.subplot(1, 2, 2)
        # C = plt.imshow(mask_faces[3], origin='lower')
        # plt.colorbar(C)
        # plt.show()

        for f in range(len(L1_faces)):
            face = L1_faces[f]
            grid = mask_faces[face]
            n_rows = int((np.shape(grid)[0] * np.shape(grid)[1]) / sNx)
            grid = np.reshape(grid, (n_rows, sNx))
            if face == 1:
                compact_stack = grid
            else:
                compact_stack = np.concatenate([compact_stack, grid], axis=0)

        if boundary in ['L2_west','L2_north','L2_south','L2_east']:
            output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'dv', boundary+'_BC_mask.bin')
        elif boundary=='L2_surface':
            output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'dv', 'L2_surface_mask.bin')
        elif boundary in ['L3_west','L3_north','L3_south','L3_east']:
            output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'dv', boundary+'_BC_mask.bin')
        elif boundary=='L3_surface':
            output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'dv', 'L3_surface_mask.bin')
        elif boundary=='CTD':
            output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'dv', 'CTD_mask.bin')
        else:
            output_file = os.path.join(config_dir, 'L1', L1_model_name, 'input', 'dv', 'L2_'+boundary+'_mask.bin')
        compact_stack.ravel(order='C').astype('>f4').tofile(output_file)

    output_dir = os.path.join(config_dir,'L1',L1_model_name,'input')
    output_file_name = 'L1_dv_mask_reference_dict.nc'
    output_mask_dictionary_to_nc(output_dir,output_file_name,all_mask_dicts,mask_names_list)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_dv_masks(config_dir)
