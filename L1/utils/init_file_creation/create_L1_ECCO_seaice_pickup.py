
import os
import numpy as np
import netCDF4 as nc4
import ast
import matplotlib.pyplot as plt
from MITgcmutils import mds
from scipy.interpolate import griddata
import sys


def read_grid_tile_geometry(config_dir,model_name,ordered_nonblank_tiles):
    ordered_XC_tiles = []
    ordered_YC_tiles = []
    ordered_AngleCS_tiles = []
    ordered_AngleSN_tiles = []

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L1', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        row_XCs = []
        row_YCs = []
        row_AngleCSs = []
        row_AngleSNs = []
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
                    AngleCS = ds.variables['AngleCS'][:, :]
                    AngleSN = ds.variables['AngleSN'][:, :]
                    ds.close()
                    row_XCs.append(XC)
                    row_YCs.append(YC)
                    row_AngleCSs.append(AngleCS)
                    row_AngleSNs.append(AngleSN)
        ordered_XC_tiles.append(row_XCs)
        ordered_YC_tiles.append(row_YCs)
        ordered_AngleCS_tiles.append(row_AngleCSs)
        ordered_AngleSN_tiles.append(row_AngleSNs)

    return(ordered_XC_tiles, ordered_YC_tiles, ordered_AngleCS_tiles, ordered_AngleSN_tiles)

def interpolate_ecco_wetgrid_to_domain(XC, YC, ecco_XC, ecco_YC, ecco_wet_cells):

    ecco_wet_cells_on_tile_domain = np.zeros((1,np.shape(XC)[0],np.shape(XC)[1]))

    points = np.hstack([np.reshape(ecco_XC, (np.size(ecco_XC), 1)),
                        np.reshape(ecco_YC, (np.size(ecco_YC), 1))])

    values = np.reshape(ecco_wet_cells[0,:,:], (np.size(ecco_wet_cells[0,:,:]), 1))
    grid = griddata(points, values, (XC, YC), method='nearest')[:, :, 0]
    ecco_wet_cells_on_tile_domain[0,:,:] = grid

    return(ecco_wet_cells_on_tile_domain)

def read_tile_wet_cells_from_run_for_grid_dir(grid_dir,var_name,tile_face,n_tiles,tile_number):

    for n in range(n_tiles):
        if 'grid.t' + '{:03d}'.format(tile_number) + '.nc' in os.listdir(
                os.path.join(grid_dir, 'mnc_' + '{:04d}'.format(n + 1))):
            ds = nc4.Dataset(os.path.join(grid_dir, 'mnc_' + '{:04d}'.format(n + 1),
                                          'grid.t' + '{:03d}'.format(tile_number) + '.nc'))
            # if tile_face<3:
            if var_name in ['siVICE']:
                hFac = 'S'
            elif var_name in ['siUICE']:
                hFac = 'W'
            else:
                hFac = 'C'
            # else:
            #     if var_name in ['Vvel', 'GvNm1', 'GvNm2']:
            #         hFac = 'W'
            #     elif var_name in ['Uvel', 'GuNm1', 'GuNm2']:
            #         hFac = 'S'
            #     else:
            #         hFac = 'C'
            wet_cells = np.array(ds.variables['HFac'+hFac][:, :, :])
            ds.close()
            wet_cells[wet_cells>0] = 1
    if hFac == 'S':
        wet_cells = wet_cells[:, :-1, :]
        if tile_number in [1, 2, 3, 4]:
            wet_cells[:, 0, :] = wet_cells[:, 1, :]
    if hFac == 'W':
        wet_cells = wet_cells[:, :, :-1]
        if tile_number in [1]:
            wet_cells[:, :, 0] = wet_cells[:, :, 1]
    return(wet_cells)

def rotate_interpolated_faces_to_domain(var_names,tile_number, tile_face_index_dict,
                                        sNx, sNy,
                                        interp_grid_faces, AngleCS, AngleSN):

    def rotate_velocity_vectors_to_domain(angle_cos, angle_sin, zonal_vel, meridional_vel):
        uvel = np.zeros_like(zonal_vel)
        vvel = np.zeros_like(meridional_vel)
        for k in range(np.shape(uvel)[0]):
            uvel[k, : ,:] = angle_cos * zonal_vel[k,:,:] + angle_sin * meridional_vel[k,:,:]
            vvel[k, : ,:] = -1 * angle_sin * zonal_vel[k,:,:] + angle_cos * meridional_vel[k,:,:]
        return (uvel, vvel)

    face = tile_face_index_dict[tile_number][0]
    min_row = tile_face_index_dict[tile_number][1]
    min_col = tile_face_index_dict[tile_number][2]

    uvel_grid_index = var_names.index('siUICE')
    vvel_grid_index = var_names.index('siVICE')
    natural_uvel = interp_grid_faces[uvel_grid_index][face][:,min_row:min_row+sNy,min_col:min_col+sNx]
    natural_vvel = interp_grid_faces[vvel_grid_index][face][:, min_row:min_row + sNy, min_col:min_col + sNx]
    uvel, vvel = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                   natural_uvel, natural_vvel)
    interp_grid_faces[uvel_grid_index][face][:, min_row:min_row + sNy, min_col:min_col + sNx] = uvel
    interp_grid_faces[vvel_grid_index][face][:, min_row:min_row + sNy, min_col:min_col + sNx] = vvel

    return(interp_grid_faces)

def stack_grids_to_compact_pickup(faces, interp_grids_faces, var_names, compact_tile_size):

    for f in range(len(interp_grids_faces)):
        print('       - Stacking field for '+var_names[f])
        grid_faces = interp_grids_faces[f]
        for i in range(len(faces)):
            face = faces[i]
            grid = grid_faces[face]
            n_rows = int((np.shape(grid)[1]*np.shape(grid)[2])/compact_tile_size)
            grid = np.reshape(grid,(np.shape(grid)[0],n_rows,compact_tile_size))
            if i==0:
                compact_stack = grid
            else:
                compact_stack = np.concatenate([compact_stack, grid], axis=1)

        if f==0:
            pickup_grid = compact_stack
        else:
            pickup_grid = np.concatenate([pickup_grid, compact_stack], axis=0)

    return(pickup_grid)

def write_seaice_pickup_file(output_file,dtype,pickup_grid,subset_metadata):

    # output the data subset
    pickup_grid.ravel(order='C').astype(dtype).tofile(output_file+'.data')

    # output the metadata file
    output = " nDims = [   "+str(subset_metadata['ndims'][0])+" ];\n"
    output += " dimList = [\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[2])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[2])+",\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[1])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[1])+"\n"
    output += " ];\n"
    output += " dataprec = [ '"+subset_metadata['dataprec'][0]+"' ];\n"
    output += " nrecords = [    "+str(subset_metadata['nrecords'][0])+" ];\n"
    output += " timeStepNumber = [ "+"{:10d}".format(subset_metadata['timestepnumber'][0])+" ];\n"
    time_interval_exponent = int(np.log10(subset_metadata['timeinterval'][0][0]))
    time_interval_base = subset_metadata['timeinterval'][0][0] / (10 ** time_interval_exponent)
    output += " timeInterval = [  "+"{:.12f}".format(time_interval_base) + "E+" + "{:02d}".format(time_interval_exponent)+  " ];\n"
    output += " nFlds = [    "+str(subset_metadata['nflds'][0])+" ];\n"
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

def create_L1_ECCO_seaice_pickup_file(config_dir, model_name,
                                      sNx,sNy,ordered_nonblank_tiles,tile_face_index_dict, face_size_dict,
                                      ecco_dir, llc, ordered_ecco_tiles, ordered_ecco_tile_rotations,
                                      parent_model_pickup_iteration, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import ecco_functions as ef

    if print_level >= 1:
        print('    - Creating the pickup file for the ' + model_name + ' model from ECCO data')

    if print_level >= 1:
        print('    - Reading in the L1 tile geometry')
    # step 0: get the model domain
    ordered_XC_tiles, ordered_YC_tiles, ordered_AngleCS_tiles, ordered_AngleSN_tiles = \
        read_grid_tile_geometry(config_dir,model_name,ordered_nonblank_tiles)

    if print_level>=1:
        print('    - Reading in the ECCO tile geometry')
    # step 1: get the ecco faces geometry
    # ecco_AngleCS, ecco_AngleSN, ecco_hfacC,
    ecco_XC, ecco_YC, ecco_AngleCS, ecco_AngleSN, ecco_hfacC, _ = \
        ef.read_ecco_grid_geometry(ecco_dir, llc, ordered_ecco_tiles, ordered_ecco_tile_rotations)

    # C = plt.imshow(ecco_hfacC[0,:,:], origin='lower')
    # plt.colorbar(C)
    # plt.show()

    if print_level >= 1:
        print('    - Reading in the ECCO pickup file')
    pickup_file = 'pickup_seaice.'+'{:010d}'.format(parent_model_pickup_iteration)
    pickup_file_path = os.path.join(config_dir,'L0','run',pickup_file)
    var_names, var_grids, global_metadata = ef.read_ecco_seaice_pickup_to_stiched_grid(pickup_file_path, ordered_ecco_tiles, ordered_ecco_tile_rotations)

    # C = plt.imshow(var_grids[4][0,:,:], origin='lower',vmin=-0.5,vmax=0.5,cmap='seismic')
    # plt.colorbar(C)
    # plt.show()

    if print_level>=1:
        print('    - Rotating oriented fields to natural coordinates')
    var_grids = ef.rotate_ecco_seaice_grids_to_natural_grids(var_names, var_grids, ecco_AngleCS, ecco_AngleSN)

    C = plt.imshow(var_grids[1][0,:,:], origin='lower',vmin=-0.5,vmax=0.5,cmap='seismic')
    plt.colorbar(C)
    plt.show()

    ecco_wet_cells = np.copy(ecco_hfacC)
    ecco_wet_cells[ecco_wet_cells>0]=1

    faces = list(face_size_dict.keys())

    # make some bins where the tiles will be stored
    interp_grid_faces = []
    for vn in range(len(var_names)):
        interp_grid = {}
        depth_levels = 1
        for face in faces:
            interp_grid[face] = np.zeros((depth_levels,
                                          face_size_dict[face][0], face_size_dict[face][1]))
        interp_grid_faces.append(interp_grid)

    run_for_grid_dir = os.path.join(config_dir,'L1',model_name,'run_for_grid')
    n_tiles = len(ordered_nonblank_tiles)*len(ordered_nonblank_tiles[0])

    ####################################################################################################################

    for tile_row in range(len(ordered_nonblank_tiles)):
        for tile_col in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[tile_row][tile_col]
            tile_face = tile_face_index_dict[tile_number][0]
            if print_level >= 2:
                print('        - Downscaling the pickup grids on tile '+str(tile_number))

            XC = ordered_XC_tiles[tile_row][tile_col]
            YC = ordered_YC_tiles[tile_row][tile_col]
            AngleCS = ordered_AngleCS_tiles[tile_row][tile_col]
            AngleSN = ordered_AngleSN_tiles[tile_row][tile_col]

            ecco_wet_cells_on_tile_domain = interpolate_ecco_wetgrid_to_domain(XC, YC, ecco_XC, ecco_YC, ecco_wet_cells)

            for vn in range(len(var_names)):
                var_name = var_names[vn]

                if var_name not in []:  # used for testing
                    ecco_grid = var_grids[var_names.index(var_name)]
                    if print_level >= 3:
                        print('            - Downscaling ' + var_name)

                    domain_wet_cells_3D = read_tile_wet_cells_from_run_for_grid_dir(run_for_grid_dir,var_name,tile_face,n_tiles,tile_number)
                    domain_wet_cells_3D = domain_wet_cells_3D[:1, :, :]

                    # plt.subplot(1, 3, 1)
                    # plt.imshow(ecco_wet_cells[0, :, :], origin='lower')
                    # plt.subplot(1, 3, 2)
                    # plt.imshow(ecco_wet_cells_on_tile_domain[0, :, :], origin='lower')
                    # plt.subplot(1, 3, 3)
                    # plt.imshow(domain_wet_cells_3D[0, :, :], origin='lower')
                    # plt.show()

                    mean_vertical_difference = 0
                    subset_copy = np.copy(ecco_grid)

                    # plt.subplot(1,2,1)
                    # plt.imshow(subset_copy[:,10,:])
                    # plt.subplot(1, 2, 2)
                    # plt.imshow(ecco_grid[:, 10, :])
                    # plt.show()

                    if var_name in ['siAREA','siHEFF','siHSNOW']:
                        remove_zeros = False
                    else:
                        remove_zeros = True

                    interp_field = df.downscale_3D_field_with_zeros(ecco_XC, ecco_YC,
                                                         ecco_grid, ecco_wet_cells,
                                                         ecco_wet_cells_on_tile_domain,
                                                         XC, YC, domain_wet_cells_3D,
                                                         mean_vertical_difference=0, fill_downward=True,
                                                         printing=False)

                    if np.sum(np.isnan(interp_field))>0:
                        print('Setting '+str(np.sum(np.isnan(interp_field)))+' values to 0 in this grid')
                        interp_field[np.isnan(interp_field)] = 0

                    # if var_name=='Theta':
                    #     plt.subplot(2, 3, 1)
                    #     plt.imshow(ecco_wet_cells[0, :, :], origin='lower')
                    #     plt.subplot(2, 3, 2)
                    #     plt.imshow(ecco_wet_cells_on_tile_domain[0, :, :], origin='lower')
                    #     plt.subplot(2, 3, 3)
                    #     plt.imshow(domain_wet_cells_3D[0, :, :], origin='lower')
                    #     plt.subplot(2, 3, 4)
                    #     plt.imshow(ecco_grid[0, :, :], origin='lower')
                    #     plt.subplot(2, 3, 5)
                    #     plt.imshow(interp_field[0, :, :], origin='lower')
                    #     plt.title('Tile: '+str(tile_number))
                    #     plt.show()

                    ##########################################################
                    # put the interpolated field into its face
                    face = tile_face_index_dict[tile_number][0]
                    min_row = tile_face_index_dict[tile_number][1]
                    min_col = tile_face_index_dict[tile_number][2]
                    interp_grid_faces[vn][face][:,min_row:min_row+sNy,min_col:min_col+sNx] = interp_field

            interp_grid_faces = rotate_interpolated_faces_to_domain(var_names, tile_number, tile_face_index_dict,
                                                                    sNx,sNy,
                                                                    interp_grid_faces, AngleCS, AngleSN)

    # for face in faces:
    #     plt.imshow(interp_grid_faces[0][face][0,:,:],origin='lower')
    #     plt.title(var_names[0]+' - face '+str(face))
    #     plt.show()

    print('    - Stacking the interpolated fields into a compact pickup grid')
    pickup_grid = stack_grids_to_compact_pickup(faces,interp_grid_faces,var_names,sNx)

    print('    - Outputting the compact pickup grid to the input directory')
    pickup_metadata = dict(global_metadata)
    output_dir = os.path.join(config_dir, 'L1', model_name, 'input')
    output_file = os.path.join(output_dir, 'pickup_seaice.' + '{:010d}'.format(4*parent_model_pickup_iteration))
    dtype = '>f8'
    pickup_metadata['timestepnumber'] = [4*int(pickup_metadata['timestepnumber'][0])]
    # pickup_metadata['nrecords'] = [np.shape(pickup_grid)[0]]
    # pickup_metadata['fldlist'] = var_names

    write_seaice_pickup_file(output_file, dtype, pickup_grid, pickup_metadata)
