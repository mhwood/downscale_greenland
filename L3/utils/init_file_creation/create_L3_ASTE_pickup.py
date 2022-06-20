
import os
import numpy as np
import netCDF4 as nc4
import ast
import matplotlib.pyplot as plt
from MITgcmutils import mds
import sys

def read_grid_geometry(config_dir,model_name, n_rows, n_cols):
    file_path = os.path.join(config_dir, 'mitgrids', model_name + '.mitgrid')
    entire_grid = np.fromfile(file_path, dtype='>f8')
    entire_grid = np.reshape(entire_grid, (16, n_rows + 1, n_cols + 1))
    XC = entire_grid[0, :, :]
    YC = entire_grid[1, :, :]
    XC = XC[:-1, :-1]
    YC = YC[:-1, :-1]

    delR = np.array([1.00, 1.14, 1.30, 1.49, 1.70,
                         1.93, 2.20, 2.50, 2.84, 3.21,
                         3.63, 4.10, 4.61, 5.18, 5.79,
                         6.47, 7.20, 7.98, 8.83, 9.73,
                         10.69, 11.70, 12.76, 13.87, 15.03,
                         16.22, 17.45, 18.70, 19.97, 21.27,
                         22.56, 23.87, 25.17, 26.46, 27.74,
                         29.00, 30.24, 31.45, 32.65, 33.82,
                         34.97, 36.09, 37.20, 38.29, 39.37,
                         40.45, 41.53, 42.62, 43.73, 44.87,
                         46.05, 47.28, 48.56, 49.93, 51.38,
                         52.93, 54.61, 56.42, 58.38, 60.53,
                         62.87, 65.43, 68.24, 71.33, 74.73,
                         78.47, 82.61, 87.17, 92.21, 97.79,
                         103.96, 110.79, 118.35, 126.73, 136.01,
                         146.30, 157.71, 170.35, 184.37, 199.89,
                         217.09, 236.13, 257.21, 280.50, 306.24,
                         334.64, 365.93, 400.38, 438.23, 479.74, ])

    return(XC, YC, delR)

def read_aste_grid_geometry(aste_dir,ordered_aste_tiles,ordered_aste_tile_rotations, L_XC, L_YC):

    aste_Nr = 50
    aste_sNx = 90
    aste_sNy = 90

    aste_XC = np.zeros((aste_sNy*len(ordered_aste_tiles),aste_sNx*len(ordered_aste_tiles[0])))
    aste_YC = np.zeros((aste_sNy * len(ordered_aste_tiles), aste_sNx * len(ordered_aste_tiles[0])))
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
            ds.close()

            # rotate things as necessary
            for n in range(aste_rotation_row[c]):
                hfac_grid = np.rot90(hfac_grid, axes=(1, 2))
                XC = np.rot90(XC)
                YC = np.rot90(YC)

            # put it into the big grid
            aste_hfacC[:, r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = hfac_grid
            aste_XC[r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = XC
            aste_YC[r * aste_sNy:(r + 1) * aste_sNy, c * aste_sNx:(c + 1) * aste_sNx] = YC

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
    aste_hfacC = aste_hfacC[:, min_row - dist_buffer:max_row + dist_buffer,
                min_col - dist_buffer:max_col + dist_buffer]

    aste_delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                        10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                        31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                        93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                        139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                        341.50, 364.50, 387.50, 410.50, 433.50, 456.50])

    aste_subset_bounds = [min_row,max_row,min_col,max_col]

    return(aste_XC,aste_YC,aste_hfacC,aste_delR,aste_subset_bounds)

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

def read_aste_pickup_to_stiched_grid(aste_dir,pickup_file,ordered_aste_tiles,ordered_aste_tile_rotations, aste_subset_bounds):

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

def read_mask_from_nc(nc_file, hFac='C'):
    ds = nc4.Dataset(nc_file)
    mask = ds.variables['wet_grid_'+hFac][:,:,:]
    ds.close()
    return(mask)

def read_mask_from_grid_nc(nc_file, hFac='C'):
    ds = nc4.Dataset(nc_file)
    mask = ds.variables['HFac'+hFac][:,:,:]
    mask[mask>0]=1
    mask[mask<=0]=0
    ds.close()
    return(mask)

def stack_grids_to_pickup(interp_grids):
    counter = 0
    for grid in interp_grids:

        if counter == 0:
            pickup_grid = grid
        else:
            pickup_grid = np.concatenate([pickup_grid, grid], axis=0)

        counter += 1
    return(pickup_grid)

def write_pickup_file(output_file,dtype,pickup_grid,subset_metadata):

    # output the data subset
    pickup_grid.ravel(order='C').astype(dtype).tofile(output_file+'.data')

    # output the metadata file
    output = " nDims = [   "+str(subset_metadata['ndims'][0])+" ];\n"
    output += " dimList = [\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[2])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[2])+",\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[1])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[1])+"\n"
    output += " ];\n"
    output += " dataprec = [ '"+subset_metadata['dataprec'][0]+"' ];\n"
    output += " nrecords = [   "+str(subset_metadata['nrecords'][0])+" ];\n"
    output += " timeStepNumber = [ "+"{:10d}".format(subset_metadata['timestepnumber'][0])+" ];\n"
    time_interval_exponent = int(np.log10(subset_metadata['timeinterval'][0][0]))
    time_interval_base = subset_metadata['timeinterval'][0][0] / (10 ** time_interval_exponent)
    output += " timeInterval = [  "+"{:.12f}".format(time_interval_base) + "E+" + "{:02d}".format(time_interval_exponent)+  " ];\n"
    output += " nFlds = [   "+str(subset_metadata['nflds'][0])+" ];\n"
    output += " fldList = {\n "
    for var_name in subset_metadata['fldlist']:
        output += "'"+var_name
        for i in range(8-len(var_name)):
            output+= " "
        output+="' "
    output += "\n };"

    f = open(output_file+'.meta','w')
    f.write(output)
    f.close()


def create_L3_ASTE_pickup_file(config_dir,model_name,ordered_aste_tiles,ordered_aste_tile_rotations):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df

    print('    - Creating the pickup file for the '+model_name+' model from ASTE data')

    f = open(os.path.join(config_dir, 'domain_sizes.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)
    L_size = size_dict[model_name]
    n_rows = L_size[0]
    n_cols = L_size[1]

    # step 0: get the model domain
    XC, YC, delR = read_grid_geometry(config_dir,model_name, n_rows, n_cols)

    # step 1: create a stitched grid around the domain using ASTE data
    aste_dir = '/Users/michwood/Documents/Research/Data Repository/Greenland/Ocean Properties/Models/ASTE'

    aste_XC, aste_YC, aste_hfacC, aste_delR, aste_subset_bounds = read_aste_grid_geometry(aste_dir,ordered_aste_tiles,ordered_aste_tile_rotations, XC, YC)
    pickup_file = 'pickup_it55warm.0000000007'
    var_names, var_grids, pickup_metadata = read_aste_pickup_to_stiched_grid(aste_dir, pickup_file,
                                                            ordered_aste_tiles, ordered_aste_tile_rotations,aste_subset_bounds)

    ASTE_wet_cells = np.copy(aste_hfacC)
    ASTE_wet_cells[ASTE_wet_cells>0]=1

    # # plt.pcolormesh(aste_XC,aste_YC,aste_hfacC[0,:,:])
    # plt.pcolormesh(aste_XC, aste_YC, var_grids[2][0,:,:])
    # plt.plot(XC[:, 0], YC[:, 0], 'k-')
    # plt.plot(XC[:, -1], YC[:, -1], 'k-')
    # plt.plot(XC[0, :], YC[0, :], 'k-')
    # plt.plot(XC[-1, :], YC[-1, :], 'k-')
    # plt.show()

    ASTE_grid_on_domain_file = os.path.join(config_dir, 'L3',model_name, 'input', 'ASTE_wetgrid_on_'+model_name+'.nc')
    domain_grid_file = os.path.join(config_dir, 'L3',model_name, 'input', model_name+'_grid.nc')

    print('    - Downscaling the pickup grids')
    interp_grids = []
    output_var_names = []
    for vn in range(len(var_names)):
        var_name = var_names[vn]

        if var_name not in []:  # used for testing
            aste_grid = var_grids[vn]
            print('      - Downscaling ' + var_name)
            print(np.shape(aste_grid))

            if var_name in ['Vvel', 'GvNm1', 'GvNm2']:
                ASTE_wet_cells_on_domain_3D = read_mask_from_nc(ASTE_grid_on_domain_file, hFac='S')
                domain_wet_cells_3D = read_mask_from_grid_nc(domain_grid_file, hFac='S')
                domain_wet_cells_3D = domain_wet_cells_3D[:, 1:, :]
            elif var_name in ['Uvel', 'GuNm1', 'GuNm2']:
                ASTE_wet_cells_on_domain_3D = read_mask_from_nc(ASTE_grid_on_domain_file, hFac='W')
                domain_wet_cells_3D = read_mask_from_grid_nc(domain_grid_file, hFac='W')
                domain_wet_cells_3D = domain_wet_cells_3D[:, :, 1:]
            else:
                ASTE_wet_cells_on_domain_3D = read_mask_from_nc(ASTE_grid_on_domain_file, hFac='C')
                domain_wet_cells_3D = read_mask_from_grid_nc(domain_grid_file, hFac='C')

            mean_vertical_difference = 0
            subset_copy = np.copy(aste_grid)

            if var_name.lower() not in ['etan', 'detahdt', 'etah']:
                aste_grid, ASTE_wet_cells = df.interpolate_var_grid_faces_to_new_depth_levels(
                    aste_grid, ASTE_wet_cells, aste_delR, delR)

            # plt.subplot(1,2,1)
            # plt.imshow(subset_copy[:,10,:])
            # plt.subplot(1, 2, 2)
            # plt.imshow(aste_grid[:, 10, :])
            # plt.show()

            interp_field = df.downscale_3D_field(aste_XC, aste_YC,
                                                 aste_grid, ASTE_wet_cells,
                                                 ASTE_wet_cells_on_domain_3D,
                                                 XC, YC, domain_wet_cells_3D,
                                                 mean_vertical_difference=0, fill_downward=True, remove_zeros=True,
                                                 printing=True)

            # plt.subplot(2, 3, 1)
            # plt.imshow(ASTE_wet_cells[0, :, :], origin='lower')
            # plt.subplot(2, 3, 2)
            # plt.imshow(ASTE_wet_cells_on_domain_3D[0, :, :], origin='lower')
            # plt.subplot(2, 3, 3)
            # plt.imshow(domain_wet_cells_3D[0, :, :], origin='lower')
            # plt.subplot(2, 2, 3)
            # plt.imshow(aste_grid[0, :, :], origin='lower')
            # plt.subplot(2, 2, 4)
            # plt.imshow(interp_field[0, :, :], origin='lower')
            # plt.show()

            # interp_field[domain_wet_cells_3D[:np.shape(interp_field)[0], :, :] == 0] = 0
        else:
            if var_name.lower() not in ['etan', 'detahdt', 'etah']:
                interp_field = np.zeros((len(delR), np.shape(XC)[0], np.shape(XC)[1]))
            else:
                interp_field = np.zeros((1, np.shape(XC)[0], np.shape(XC)[1]))

        interp_grids.append(interp_field)
        output_var_names.append(var_name)
        if var_name=='GuNm1':
            interp_grids.append(interp_field)
            output_var_names.append('GuNm2')
        if var_name=='GvNm1':
            interp_grids.append(interp_field)
            output_var_names.append('GvNm2')

    pickup_grid = stack_grids_to_pickup(interp_grids)
    print(np.shape(pickup_grid))

    output_dir = os.path.join(config_dir, 'L3', model_name, 'input')
    output_file = os.path.join(output_dir, 'pickup.' + '{:010d}'.format(1))
    dtype = '>f8'
    pickup_metadata['timestepnumber'] = [int(1)]
    pickup_metadata['nrecords'] = [int(len(output_var_names) * len(delR) + 3)]
    pickup_metadata['fldlist'] = output_var_names
    write_pickup_file(output_file, dtype, pickup_grid, pickup_metadata)



    # step 2: run the downscale scripts to do the interpolation



   

