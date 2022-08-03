
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from MITgcmutils import mds

def get_aste_metadata(metadata_name):
    if metadata_name=='Nr':
        metadata = 50
    elif metadata_name == 'sNx':
        metadata = 90
    elif metadata_name == 'sNy':
        metadata = 90
    elif metadata_name == 'delR':
        metadata = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                              10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                              31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                              93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                              139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                              341.50, 364.50, 387.50, 410.50, 433.50, 456.50])
    else:
        raise ValueError('metadata_name not recognized')
    return(metadata)


#################################################################################################
# grid functions

def read_aste_grid_geometry(aste_dir,ordered_aste_tiles,ordered_aste_tile_rotations):

    aste_Nr = 50
    aste_sNx = 90
    aste_sNy = 90

    aste_XC = np.zeros((aste_sNy*len(ordered_aste_tiles),aste_sNx*len(ordered_aste_tiles[0])))
    aste_YC = np.zeros((aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))
    aste_AngleCS = np.zeros((aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))
    aste_AngleSN = np.zeros((aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))
    aste_hfacC = np.zeros((aste_Nr, aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))

    for r in range(len(ordered_aste_tiles)):
        aste_tile_row = ordered_aste_tiles[r]
        aste_rotation_row = ordered_aste_tile_rotations[r]
        for c in range(len(ordered_aste_tiles[r])):

            # get the hfac grid
            file_name = os.path.join(aste_dir, 'nctiles_grid', 'GRID.' + '{:04d}'.format(aste_tile_row[c])) + '.nc'
            ds = nc4.Dataset(file_name)
            hfac_grid = ds.variables['hFacC'][:, :, :]
            XC = ds.variables['XC'][:, :]
            YC = ds.variables['YC'][:, :]
            AngleCS = ds.variables['AngleCS'][:, :]
            AngleSN = ds.variables['AngleSN'][:, :]
            ds.close()

            # rotate things as necessary
            for n in range(aste_rotation_row[c]):
                hfac_grid = np.rot90(hfac_grid, axes=(1, 2))
                XC = np.rot90(XC)
                YC = np.rot90(YC)
                AngleCS = np.rot90(AngleCS)
                AngleSN = np.rot90(AngleSN)

            # put it into the big grid
            aste_hfacC[:, r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = hfac_grid
            aste_XC[r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = XC
            aste_YC[r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = YC
            aste_AngleCS[r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = AngleCS
            aste_AngleSN[r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = AngleSN

    aste_delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                        10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                        31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                        93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                        139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                        341.50, 364.50, 387.50, 410.50, 433.50, 456.50])

    return(aste_XC,aste_YC,aste_AngleCS,aste_AngleSN,aste_hfacC,aste_delR)

def read_aste_grid_geometry_to_faces(aste_dir):

    aste_Nr = 50
    aste_sNx = 90
    aste_sNy = 90

    aste_XC_faces = {}
    aste_YC_faces = {}
    aste_AngleCS_faces = {}
    aste_AngleSN_faces = {}
    aste_hfacC_faces = {}

    face_dimensions = {1: [5 * aste_sNy, 3 * aste_sNx],
                       3: [3 * aste_sNy, 3 * aste_sNx],
                       4: [3 * aste_sNy, 2 * aste_sNx],
                       5: [3 * aste_sNy, 5 * aste_sNx]}

    for face in [1,3,4,5]:
        dimension = face_dimensions[face]
        aste_XC = np.zeros((dimension[0], dimension[1]))
        aste_YC = np.zeros((dimension[0], dimension[1]))
        aste_AngleCS = np.zeros((dimension[0], dimension[1]))
        aste_AngleSN = np.zeros((dimension[0], dimension[1]))
        aste_hfacC = np.zeros((aste_Nr, dimension[0], dimension[1]))

        aste_XC_faces[face] = aste_XC
        aste_YC_faces[face] = aste_YC
        aste_AngleCS_faces[face] = aste_AngleCS
        aste_AngleSN_faces[face] = aste_AngleSN
        aste_hfacC_faces[face] = aste_hfacC

    aste_tile_to_face_indices = {3: [1,3,0],
                                 5: [1,4,0],
                                 6: [1,4,1],
                                 11: [3,1,0],
                                 12: [3,1,1],
                                 14: [3,2,0],
                                 15: [3,2,1],
                                 24: [5,1,0],
                                 27: [5,2,0]}

    for aste_tile in [3,5,6,11,12,14,15,24,27]:

            # get the grids
            file_name = os.path.join(aste_dir, 'nctiles_grid', 'GRID.' + '{:04d}'.format(aste_tile) + '.nc')
            ds = nc4.Dataset(file_name)
            hfac_grid = ds.variables['hFacC'][:, :, :]
            XC = ds.variables['XC'][:, :]
            YC = ds.variables['YC'][:, :]
            AngleCS = ds.variables['AngleCS'][:, :]
            AngleSN = ds.variables['AngleSN'][:, :]
            ds.close()

            face = aste_tile_to_face_indices[aste_tile][0]
            r = aste_tile_to_face_indices[aste_tile][1]
            c = aste_tile_to_face_indices[aste_tile][2]

            # put it into the big grid
            aste_hfacC_faces[face][:, r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = hfac_grid
            aste_XC_faces[face][r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = XC
            aste_YC_faces[face][r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = YC
            aste_AngleCS_faces[face][r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = AngleCS
            aste_AngleSN_faces[face][r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = AngleSN

    aste_delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                        10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                        31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                        93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                        139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                        341.50, 364.50, 387.50, 410.50, 433.50, 456.50])


    return(aste_XC_faces, aste_YC_faces, aste_AngleCS_faces, aste_AngleSN_faces, aste_hfacC_faces, aste_delR)

def read_aste_grid_geometry_to_subset(aste_dir,ordered_aste_tiles,ordered_aste_tile_rotations, L_XC, L_YC):

    aste_Nr = 50
    aste_sNx = 90
    aste_sNy = 90

    aste_XC = np.zeros((aste_sNy*len(ordered_aste_tiles),aste_sNx*len(ordered_aste_tiles[0])))
    aste_YC = np.zeros((aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))
    aste_AngleCS = np.zeros((aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))
    aste_AngleSN = np.zeros((aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))
    aste_hfacC = np.zeros((aste_Nr, aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))

    for r in range(len(ordered_aste_tiles)):
        aste_tile_row = ordered_aste_tiles[r]
        aste_rotation_row = ordered_aste_tile_rotations[r]
        for c in range(len(ordered_aste_tiles[r])):

            # get the hfac grid
            file_name = os.path.join(aste_dir, 'nctiles_grid', 'GRID.' + '{:04d}'.format(aste_tile_row[c])) + '.nc'
            ds = nc4.Dataset(file_name)
            hfac_grid = ds.variables['hFacC'][:, :, :]
            XC = ds.variables['XC'][:, :]
            YC = ds.variables['YC'][:, :]
            AngleCS = ds.variables['AngleCS'][:, :]
            AngleSN = ds.variables['AngleSN'][:, :]
            ds.close()

            # rotate things as necessary
            for n in range(aste_rotation_row[c]):
                hfac_grid = np.rot90(hfac_grid, axes=(1, 2))
                XC = np.rot90(XC)
                YC = np.rot90(YC)
                AngleCS = np.rot90(AngleCS)
                AngleSN = np.rot90(AngleSN)

            # put it into the big grid
            aste_hfacC[:, r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = hfac_grid
            aste_XC[r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = XC
            aste_YC[r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = YC
            aste_AngleCS[r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = AngleCS
            aste_AngleSN[r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = AngleSN

    # C = plt.imshow(aste_grid[0,:,:],origin='lower')
    # C = plt.imshow(aste_XC, origin='lower')
    # C = plt.imshow(aste_YC, origin='lower')
    # plt.colorbar(C)
    # plt.show()

    ll_dist = ((L_XC[0, 0] - aste_XC) ** 2 + (L_YC[0, 0] - aste_YC) ** 2) ** 0.5
    ll_row, ll_col = np.where(ll_dist == np.min(ll_dist))
    ul_dist = ((L_XC[-1, 0] - aste_XC) ** 2 + (L_YC[-1, 0] - aste_YC) ** 2) ** 0.5
    ul_row, ul_col = np.where(ul_dist == np.min(ul_dist))

    lr_dist = ((L_XC[0, -1] - aste_XC) ** 2 + (L_YC[0, -1] - aste_YC) ** 2) ** 0.5
    lr_row, lr_col = np.where(lr_dist == np.min(lr_dist))
    ur_dist = ((L_XC[-1, -1] - aste_XC) ** 2 + (L_YC[-1, -1] - aste_YC) ** 2) ** 0.5
    ur_row, ur_col = np.where(ur_dist == np.min(ur_dist))

    min_row = np.min([ll_row[0], ul_row[0], lr_row[0], ur_row[0]])
    max_row = np.max([ll_row[0], ul_row[0], lr_row[0], ur_row[0]])
    min_col = np.min([ll_col[0], ul_col[0], lr_col[0], ur_col[0]])
    max_col = np.max([ll_col[0], ul_col[0], lr_col[0], ur_col[0]])

    dist_buffer = 3
    aste_XC = aste_XC[min_row - dist_buffer:max_row + dist_buffer,
              min_col - dist_buffer:max_col + dist_buffer]
    aste_YC = aste_YC[min_row - dist_buffer:max_row + dist_buffer,
              min_col - dist_buffer:max_col + dist_buffer]
    aste_AngleCS = aste_AngleCS[min_row - dist_buffer:max_row + dist_buffer,
              min_col - dist_buffer:max_col + dist_buffer]
    aste_AngleSN = aste_AngleSN[min_row - dist_buffer:max_row + dist_buffer,
              min_col - dist_buffer:max_col + dist_buffer]
    aste_hfacC = aste_hfacC[:, min_row - dist_buffer:max_row + dist_buffer,
                min_col - dist_buffer:max_col + dist_buffer]

    aste_delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                        10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                        31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                        93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                        139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                        341.50, 364.50, 387.50, 410.50, 433.50, 456.50])

    aste_subset_bounds = [min_row,max_row,min_col,max_col]

    return(aste_XC,aste_YC,aste_AngleCS,aste_AngleSN,aste_hfacC,aste_delR,aste_subset_bounds)


##################################################################################################
# pickup function

def read_pickup_file_to_compact(pickup_file_path):

    Nr = 50
    print('      Reading from '+pickup_file_path)
    global_data, _, global_metadata = mds.rdmds(pickup_file_path, returnmeta=True)

    has_Nr = {'uvel': True, 'vvel': True, 'theta': True,
              'salt': True, 'gunm1': True, 'gvnm1': True,
              'gunm2': True, 'gvnm2': True, 'etan': False,
              'detahdt': False, 'etah': False}

    var_names = []
    row_bounds = []
    var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        if has_Nr[var_name.strip().lower()]:
            end_row = start_row + Nr
        else:
            end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]
        var_grids.append(var_grid)
        row_bounds.append([start_row,end_row])
        start_row=end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,var_grids,global_metadata)

def read_tile_from_ASTE_compact(var_grid,tile_number, print_messages=False):

    sNx = 90
    sNy = 90

    if print_messages:
        print('Reading grid for tile number '+str(tile_number))

    # adjust tile number to account for blank cells
    tile_add = 6
    if tile_number>2:
        tile_add+=1
    if tile_number>4:
        tile_add+=1
    if tile_number>23:
        tile_add+=4
    if tile_number>26:
        tile_add+=2
    tile_number += tile_add

    if print_messages:
        print('    - The tile number including blank cells is '+str(tile_number))

    # get the face number
    if tile_number < 16:
        face_number = 1
    if tile_number >= 16 and tile_number < 25:
        face_number = 2
    if tile_number >= 25 and tile_number<31:
        face_number = 3
    if tile_number >= 31:
        face_number = 4

    if print_messages:
        print('    - It is found in face '+str(face_number))

    tiles_per_face = {1:15,2:9,3:6,4:15}
    face_indices_in_compact = {1:[0,15*sNx],
                               2:[15*sNx,(15+9)*sNx],
                               3:[(15+9)*sNx,(15+9+6)*sNx],
                               4:[(15+9+6)*sNx,(15+9+6+15)*sNx]}
    face_dimensions = {1:[5*sNy,3*sNx],
                       2:[3*sNy,3*sNx],
                       3:[3*sNy,2*sNx],
                       4:[3*sNy,5*sNx]}
    face_tile_dimensions = {1: [5, 3],
                           2: [3, 3],
                           3: [3, 2],
                           4: [3, 5]}
    face_first_tile = {1:1,2:16,3:25,4:31}

    indices = face_indices_in_compact[face_number]
    dimensions = face_dimensions[face_number]

    if print_messages:
        print('    - The face will be dimension '+str(dimensions))
        print('    - The face will be read from row indices ' + str(indices)+' from the compact file')

    face_subset = var_grid[:,indices[0]:indices[1],:]
    face_subset = np.reshape(face_subset,(np.shape(face_subset)[0],dimensions[0],dimensions[1]))

    if print_messages:
        print('    - The tile number in the face is '+str((tile_number-face_first_tile[face_number]+1)))
        print('    - The number of tiles before the tile in the face are  '+str(face_tile_dimensions[face_number][1]))

    ll_row = sNy*int(((tile_number-face_first_tile[face_number]) // face_tile_dimensions[face_number][1]))
    ll_col = sNx*((tile_number-face_first_tile[face_number]) % face_tile_dimensions[face_number][1])

    if print_messages:
        print('    - The lower left row in the face is '+str(ll_row))
        print('    - The lower left col in the face is ' + str(ll_col))

    tile_subset = face_subset[:,ll_row:ll_row+sNy,ll_col:ll_col+sNx]

    # plt.imshow(tile_subset[0,:,:],origin='lower')
    # plt.title(str(tile_number-tile_add))
    # plt.show()

    return(tile_subset)

def read_aste_pickup_to_stiched_grid(aste_dir,pickup_file,ordered_aste_tiles,ordered_aste_tile_rotations):

    aste_Nr = 50
    aste_sNx = 90
    aste_sNy = 90

    pickup_file_path = os.path.join(aste_dir,'pickup',pickup_file)
    var_names,row_bounds,compact_var_grids,global_metadata = read_pickup_file_to_compact(pickup_file_path)

    var_grids = []

    for vn in range(len(var_names)):
        compact_var_grid = compact_var_grids[vn]
        if var_names[vn].lower() not in ['etan', 'detahdt', 'etah']:
            aste_grid = np.zeros((aste_Nr, aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))
        else:
            aste_grid = np.zeros((1, aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))

        for r in range(len(ordered_aste_tiles)):
            aste_tile_row = ordered_aste_tiles[r]
            aste_rotation_row = ordered_aste_tile_rotations[r]
            for c in range(len(ordered_aste_tiles[r])):

                # get the variable grid
                var_grid = read_tile_from_ASTE_compact(compact_var_grid, tile_number=aste_tile_row[c])

                # rotate things as necessary
                for n in range(aste_rotation_row[c]):
                    var_grid = np.rot90(var_grid, axes=(1, 2))

                # put it into the big grid
                aste_grid[:, r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = var_grid

        var_grids.append(aste_grid)

    return(var_names, var_grids, global_metadata)

def read_aste_pickup_to_subsetted_stiched_grid(aste_dir,pickup_file,ordered_aste_tiles,ordered_aste_tile_rotations, aste_subset_bounds):

    aste_Nr = 50
    aste_sNx = 90
    aste_sNy = 90

    pickup_file_path = os.path.join(aste_dir,'pickup',pickup_file)
    var_names,row_bounds,compact_var_grids,global_metadata = read_pickup_file_to_compact(pickup_file_path)

    var_grids = []

    dist_buffer = 3
    min_row = aste_subset_bounds[0]
    max_row = aste_subset_bounds[1]
    min_col = aste_subset_bounds[2]
    max_col = aste_subset_bounds[3]

    for vn in range(len(var_names)):
        compact_var_grid = compact_var_grids[vn]
        if var_names[vn].lower() not in ['etan', 'detahdt', 'etah']:
            aste_grid = np.zeros((aste_Nr, aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))
        else:
            aste_grid = np.zeros((1, aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))

        for r in range(len(ordered_aste_tiles)):
            aste_tile_row = ordered_aste_tiles[r]
            aste_rotation_row = ordered_aste_tile_rotations[r]
            for c in range(len(ordered_aste_tiles[r])):

                # get the variable grid
                var_grid = read_tile_from_ASTE_compact(compact_var_grid, tile_number=aste_tile_row[c])

                # rotate things as necessary
                for n in range(aste_rotation_row[c]):
                    var_grid = np.rot90(var_grid, axes=(1, 2))

                # put it into the big grid
                aste_grid[:, r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = var_grid

        # chop it down
        aste_grid = aste_grid[:, min_row - dist_buffer:max_row + dist_buffer,
                    min_col - dist_buffer:max_col + dist_buffer]

        # C = plt.imshow(aste_grid[0,:,:],origin='lower')
        # # C = plt.imshow(aste_XC, origin='lower')
        # # C = plt.imshow(aste_YC, origin='lower')
        # plt.colorbar(C)
        # plt.title(var_names[vn])
        # plt.show()

        var_grids.append(aste_grid)

    return(var_names, var_grids, global_metadata)

def read_aste_pickup_to_faces(aste_dir,pickup_file):

    aste_Nr = 50
    aste_sNx = 90
    aste_sNy = 90

    pickup_file_path = os.path.join(aste_dir,'pickup',pickup_file)
    var_names,row_bounds,compact_var_grids,global_metadata = read_pickup_file_to_compact(pickup_file_path)

    face_indices_in_compact = {1: [0, 15 * aste_sNx],
                               3: [15 * aste_sNx, (15 + 9) * aste_sNx],
                               4: [(15 + 9) * aste_sNx, (15 + 9 + 6) * aste_sNx],
                               5: [(15 + 9 + 6) * aste_sNx, (15 + 9 + 6 + 15) * aste_sNx]}

    face_dimensions = {1: [5 * aste_sNy, 3 * aste_sNx],
                       3: [3 * aste_sNy, 3 * aste_sNx],
                       4: [3 * aste_sNy, 2 * aste_sNx],
                       5: [3 * aste_sNy, 5 * aste_sNx]}

    var_grid_faces = []
    for vn in range(len(compact_var_grids)):
        compact_var_grid = compact_var_grids[vn]
        grid_faces = {}
        for face in [1,3,4,5]:
            indices = face_indices_in_compact[face]
            face_grid = compact_var_grid[:,indices[0]:indices[1],:]
            dimensions = face_dimensions[face]
            face_grid = np.reshape(face_grid,(np.shape(face_grid)[0],dimensions[0],dimensions[1]))
            grid_faces[face] = face_grid

            # if var_names[vn]=='Theta':
            #     plt.imshow(face_grid[0,:,:],origin='lower')
            #     plt.show()

        var_grid_faces.append(grid_faces)

    return(var_names, var_grid_faces, global_metadata)

##################################################################################################
# rotation functions

def rotate_aste_grids_to_natural_grids(var_names, var_grids, aste_AngleCS, aste_AngleSN):

    def rotate_velocity_vectors_to_natural(angle_cos, angle_sin, uvel, vvel):
        zonal_velocity = np.zeros_like(uvel)
        meridional_velocity = np.zeros_like(vvel)
        for k in range(np.shape(uvel)[0]):
            zonal_velocity[k,:,:] = angle_cos * uvel[k,:,:] - angle_sin * vvel[k,:,:]
            meridional_velocity[k,:,:] = angle_sin * uvel[k,:,:] + angle_cos * vvel[k,:,:]
        return (zonal_velocity, meridional_velocity)

    uvel_grid_index = var_names.index('Uvel')
    vvel_grid_index = var_names.index('Vvel')

    zonal_uvel, meridional_vvel = rotate_velocity_vectors_to_natural(aste_AngleCS, aste_AngleSN,
                                                              var_grids[uvel_grid_index], var_grids[vvel_grid_index])

    # plt.subplot(2,2,1)
    # C = plt.imshow(var_grids[uvel_grid_index][0,:,:],origin='lower',cmap='seismic',vmin=-0.4,vmax=0.4)
    # plt.colorbar(C)
    #
    # plt.subplot(2, 2, 2)
    # C = plt.imshow(var_grids[vvel_grid_index][0, :, :], origin='lower', cmap='seismic', vmin=-0.4, vmax=0.4)
    # plt.colorbar(C)
    #
    # plt.subplot(2, 2, 3)
    # C = plt.imshow(zonal_uvel[0, :, :], origin='lower', cmap='seismic', vmin=-0.4, vmax=0.4)
    # plt.colorbar(C)
    #
    # plt.subplot(2, 2, 4)
    # C = plt.imshow(meridional_vvel[0, :, :], origin='lower', cmap='seismic', vmin=-0.4, vmax=0.4)
    # plt.colorbar(C)
    #
    # plt.show()

    gunm1_grid_index = var_names.index('GuNm1')
    gvnm1_grid_index = var_names.index('GvNm1')

    zonal_gunm1, meridional_gvnm1 = rotate_velocity_vectors_to_natural(aste_AngleCS, aste_AngleSN,
                                                                var_grids[gunm1_grid_index], var_grids[gvnm1_grid_index])

    var_grids[uvel_grid_index] = zonal_uvel
    var_grids[vvel_grid_index] = meridional_vvel

    var_grids[gunm1_grid_index] = zonal_gunm1
    var_grids[gvnm1_grid_index] = meridional_gvnm1
    return(var_grids)


