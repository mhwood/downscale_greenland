
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

    delR = np.array([10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
                        10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
                        31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
                        93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
                        139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
                        341.50, 364.50, 387.50, 410.50, 433.50, 456.50])

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L05', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        row_XCs = []
        row_YCs = []
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
                    ds.close()
                    row_XCs.append(XC)
                    row_YCs.append(YC)
        ordered_XC_tiles.append(row_XCs)
        ordered_YC_tiles.append(row_YCs)

    return(ordered_XC_tiles, ordered_YC_tiles, delR)

def interpolate_aste_wetgrid_to_domain(XC, YC, delR, aste_XC, aste_YC, aste_delR, aste_wet_cells):

    Z_bottom = np.cumsum(delR)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])
    Z = (Z_bottom + Z_top) / 2

    aste_Z_bottom = np.cumsum(aste_delR)
    aste_Z_top = np.concatenate([np.array([0]), aste_Z_bottom[:-1]])
    aste_Z = (aste_Z_bottom + aste_Z_top) / 2

    ##########
    # interpolate vertically first

    aste_wet_cells_interpolated = np.zeros((len(delR),np.shape(aste_wet_cells)[1], np.shape(aste_wet_cells)[2]))
    for i in range(len(delR)):
        if np.any(aste_Z>Z[i]):
            aste_Z_subset = aste_Z[aste_Z>Z[i]]
            index = np.argmin(np.abs(aste_Z-np.min(aste_Z_subset)))
        else:
            index = len(aste_delR)-1
        aste_wet_cells_interpolated[i,:,:] = aste_wet_cells[index,:,:]

    ##########
    # then interpolate onto the tile

    aste_wet_cells_on_tile_domain = np.zeros((len(delR),np.shape(XC)[0],np.shape(XC)[1]))

    points = np.hstack([np.reshape(aste_XC, (np.size(aste_XC), 1)),
                        np.reshape(aste_YC, (np.size(aste_YC), 1))])

    for i in range(len(delR)):
        values = np.reshape(aste_wet_cells_interpolated[i,:,:], (np.size(aste_wet_cells_interpolated[i,:,:]), 1))
        grid = griddata(points, values, (XC, YC), method='nearest')[:, :, 0]
        aste_wet_cells_on_tile_domain[i,:,:] = grid

    return(aste_wet_cells_on_tile_domain)

def read_tile_wet_cells_from_run_for_grid_dir(grid_dir,var_name,tile_face,n_tiles,tile_number):

    for n in range(n_tiles):
        if 'grid.t' + '{:03d}'.format(tile_number) + '.nc' in os.listdir(
                os.path.join(grid_dir, 'mnc_' + '{:04d}'.format(n + 1))):
            ds = nc4.Dataset(os.path.join(grid_dir, 'mnc_' + '{:04d}'.format(n + 1),
                                          'grid.t' + '{:03d}'.format(tile_number) + '.nc'))
            if tile_face<3:
                if var_name in ['Vvel', 'GvNm1', 'GvNm2']:
                    hFac = 'S'
                elif var_name in ['Uvel', 'GuNm1', 'GuNm2']:
                    hFac = 'W'
                else:
                    hFac = 'C'
            else:
                if var_name in ['Vvel', 'GvNm1', 'GvNm2']:
                    hFac = 'W'
                elif var_name in ['Uvel', 'GuNm1', 'GuNm2']:
                    hFac = 'S'
                else:
                    hFac = 'C'
            wet_cells = np.array(ds.variables['HFac'+hFac][:, :, :])
            ds.close()
            wet_cells[wet_cells>0] = 1
    if hFac=='S':
        wet_cells = wet_cells[:,1:,:]
    if hFac=='W':
        wet_cells = wet_cells[:,:,:-1]
    return(wet_cells)


def stack_grids_to_compact_pickup(faces, interp_grids_faces, output_var_names, compact_tile_size):

    for f in range(len(interp_grids_faces)):
        print('       - Stacking field for '+output_var_names[f])
        grid_faces = interp_grids_faces[f]
        for face in faces:
            grid = grid_faces[face]
            n_rows = int((np.shape(grid)[1]*np.shape(grid)[2])/compact_tile_size)
            grid = np.reshape(grid,(np.shape(grid)[0],n_rows,compact_tile_size))
            if face==1:
                compact_stack = grid
            else:
                compact_stack = np.concatenate([compact_stack, grid], axis=1)

        if f==0:
            pickup_grid = compact_stack
        else:
            pickup_grid = np.concatenate([pickup_grid, compact_stack], axis=0)

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

def create_L05_ASTE_diffkr_file(config_dir,model_name,sNx,sNy,ordered_nonblank_tiles,tile_face_index_dict,
                               ordered_aste_tiles, ordered_aste_tile_rotations):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import aste_functions as af

    print('    - Creating the pickup file for the '+model_name+' model from ASTE data')

    compact_tile_size = 90

    # step 0: get the model domain
    ordered_XC_tiles, ordered_YC_tiles, delR = \
        read_grid_tile_geometry(config_dir,model_name,ordered_nonblank_tiles)

    # step 1: get the aste faces geometry
    aste_dir = '/Users/michwood/Documents/Research/Data Repository/Greenland/Ocean Properties/Models/ASTE'
    aste_XC, aste_YC, _, _, aste_hfacC, aste_delR = \
        af.read_aste_grid_geometry(aste_dir, ordered_aste_tiles, ordered_aste_tile_rotations)

    var_path = os.path.join(aste_dir,'input','diffkr_i62.bin')
    dims = 3
    diff_kr_faces = af.read_aste_compact_to_faces(var_path,dims,ordered_aste_tiles, ordered_aste_tile_rotations)

    plt.subplot(1, 2, 1)
    C = plt.imshow(diff_kr_faces[1][0,:,:], origin='lower')
    plt.colorbar(C)
    plt.subplot(1, 2, 2)
    C = plt.imshow(diff_kr_faces[3][0, :, :], origin='lower')
    plt.colorbar(C)
    plt.show()

    # var_path = os.path.join(aste_dir, 'input', 'bathy_fill9iU42Ef_noStLA_v1.bin')
    # dims = 2
    # diff_kr_faces = af.read_aste_compact_to_faces(var_path, dims, ordered_aste_tiles, ordered_aste_tile_rotations)
    # plt.subplot(1, 2, 1)
    # C = plt.imshow(diff_kr_faces[1][:, :], origin='lower')
    # plt.colorbar(C)
    # plt.subplot(1, 2, 2)
    # C = plt.imshow(diff_kr_faces[3][:, :], origin='lower')
    # plt.colorbar(C)
    # plt.show()

    # uvel_grid = var_grids[var_names.index('Uvel')]
    # vvel_grid = var_grids[var_names.index('Vvel')]
    # natural_uvel_grid, natural_vvel_grid = af.rotate_aste_grids_to_natural_grids(uvel_grid,vvel_grid, aste_AngleCS, aste_AngleSN)
    # var_grids[var_names.index('Uvel')] = natural_uvel_grid
    # var_grids[var_names.index('Vvel')] = natural_vvel_grid
    #
    # gunm1_grid = var_grids[var_names.index('GuNm1')]
    # gvnm1_grid = var_grids[var_names.index('GvNm1')]
    # natural_gunm1_grid, natural_gvnm1_grid = af.rotate_aste_grids_to_natural_grids(gunm1_grid, gvnm1_grid, aste_AngleCS,
    #                                                                              aste_AngleSN)
    # var_grids[var_names.index('GuNm1')] = natural_gunm1_grid
    # var_grids[var_names.index('GvNm1')] = natural_gvnm1_grid
    #
    # aste_wet_cells = np.copy(aste_hfacC)
    # aste_wet_cells[aste_wet_cells>0]=1
    #
    # # # # plt.pcolormesh(aste_XC,aste_YC,aste_hfacC[0,:,:])
    # # # plt.pcolormesh(aste_XC, aste_YC, var_grids[2][0,:,:])
    # # # plt.plot(XC[:, 0], YC[:, 0], 'k-')
    # # # plt.plot(XC[:, -1], YC[:, -1], 'k-')
    # # # plt.plot(XC[0, :], YC[0, :], 'k-')
    # # # plt.plot(XC[-1, :], YC[-1, :], 'k-')
    # # # plt.show()
    #
    # output_var_names = []
    # for var_name in var_names:
    #    output_var_names.append(var_name)
    #    if var_name=='GuNm1':
    #        output_var_names.append('GuNm2')
    #    if var_name=='GvNm1':
    #        output_var_names.append('GvNm2')
    #
    # faces_dim_file = os.path.join(config_dir, 'L05', model_name, 'namelist', 'face_dimensions.txt')
    # f = open(faces_dim_file)
    # dict_str = f.read()
    # f.close()
    # size_dict = ast.literal_eval(dict_str)
    # faces = list(size_dict.keys())
    #
    # # make some bins where the tiles will be stored
    # interp_grid_faces = []
    # for vn in range(len(output_var_names)):
    #     interp_grid = {}
    #     if output_var_names[vn].lower() not in ['etan', 'detahdt', 'etah']:
    #         depth_levels = len(delR)
    #     else:
    #         depth_levels = 1
    #     for face in faces:
    #         interp_grid[face] = np.zeros((depth_levels,
    #                                       size_dict[face][0], size_dict[face][1]))
    #     interp_grid_faces.append(interp_grid)
    #
    # run_for_grid_dir = os.path.join(config_dir,'L05',model_name,'run_for_grid')
    # n_tiles = len(ordered_nonblank_tiles)*len(ordered_nonblank_tiles[0])
    #
    # ####################################################################################################################
    #
    #
    # for tile_row in range(len(ordered_nonblank_tiles)):
    #     for tile_col in range(len(ordered_nonblank_tiles[0])):
    #         tile_number = ordered_nonblank_tiles[tile_row][tile_col]
    #         tile_face = tile_face_index_dict[tile_number][0]
    #         print('    - Downscaling the pickup grids on tile '+str(tile_number))
    #
    #         XC = ordered_XC_tiles[tile_row][tile_col]
    #         YC = ordered_YC_tiles[tile_row][tile_col]
    #         AngleCS = ordered_AngleCS_tiles[tile_row][tile_col]
    #         AngleSN = ordered_AngleSN_tiles[tile_row][tile_col]
    #
    #         aste_wet_cells_on_tile_domain = interpolate_aste_wetgrid_to_domain(XC, YC, delR,
    #                                                                            aste_XC, aste_YC, aste_delR, aste_wet_cells)
    #
    #         for vn in range(len(output_var_names)):
    #             var_name = output_var_names[vn]
    #
    #             if var_name not in ['Uvel','Vvel',]:  # used for testing
    #                 if var_name=='GuNm2':
    #                     aste_grid = var_grids[var_names.index('GuNm1')]
    #                 elif var_name=='GvNm2':
    #                     aste_grid = var_grids[var_names.index('GvNm1')]
    #                 else:
    #                     aste_grid = var_grids[var_names.index(var_name)]
    #                 print('        - Downscaling ' + var_name)
    #
    #                 domain_wet_cells_3D = read_tile_wet_cells_from_run_for_grid_dir(run_for_grid_dir,var_name,tile_face,n_tiles,tile_number)
    #
    #                 # plt.subplot(1, 3, 1)
    #                 # plt.imshow(aste_wet_cells[0, :, :], origin='lower')
    #                 # plt.subplot(1, 3, 2)
    #                 # plt.imshow(aste_wet_cells_on_tile_domain[0, :, :], origin='lower')
    #                 # plt.subplot(1, 3, 3)
    #                 # plt.imshow(domain_wet_cells_3D[0, :, :], origin='lower')
    #                 # plt.show()
    #
    #                 mean_vertical_difference = 0
    #                 subset_copy = np.copy(aste_grid)
    #                 #
    #                 # if var_name.lower() not in ['etan', 'detahdt', 'etah']:
    #                 #     aste_grid, aste_wet_cells = df.interpolate_var_grid_faces_to_new_depth_levels(
    #                 #         aste_grid, aste_wet_cells, aste_delR, delR)
    #                 # print('     - Skipping the vertical interpolation')
    #
    #                 # plt.subplot(1,2,1)
    #                 # plt.imshow(subset_copy[:,10,:])
    #                 # plt.subplot(1, 2, 2)
    #                 # plt.imshow(aste_grid[:, 10, :])
    #                 # plt.show()
    #
    #                 interp_field = df.downscale_3D_field(aste_XC, aste_YC,
    #                                                      aste_grid, aste_wet_cells,
    #                                                      aste_wet_cells_on_tile_domain,
    #                                                      XC, YC, domain_wet_cells_3D,
    #                                                      mean_vertical_difference=0, fill_downward=True, remove_zeros=True,
    #                                                      printing=False)
    #
    #                 if np.sum(np.isnan(interp_field))>0:
    #                     print('Setting '+str(np.sum(np.isnan(interp_field)))+' values to 0 in this grid')
    #                     interp_field[np.isnan(interp_field)] = 0
    #
    #
    #                 # if var_name=='Theta':
    #                 #     plt.subplot(2, 3, 1)
    #                 #     plt.imshow(aste_wet_cells[0, :, :], origin='lower')
    #                 #     plt.subplot(2, 3, 2)
    #                 #     plt.imshow(aste_wet_cells_on_tile_domain[0, :, :], origin='lower')
    #                 #     plt.subplot(2, 3, 3)
    #                 #     plt.imshow(domain_wet_cells_3D[0, :, :], origin='lower')
    #                 #     plt.subplot(2, 3, 4)
    #                 #     plt.imshow(aste_grid[0, :, :], origin='lower')
    #                 #     plt.subplot(2, 3, 5)
    #                 #     plt.imshow(interp_field[0, :, :], origin='lower')
    #                 #     plt.title('Tile: '+str(tile_number))
    #                 #     plt.show()
    #
    #                 ##########################################################
    #                 # put the interpolated field into its face
    #                 face = tile_face_index_dict[tile_number][0]
    #                 min_row = tile_face_index_dict[tile_number][1]
    #                 min_col = tile_face_index_dict[tile_number][2]
    #                 interp_grid_faces[vn][face][:,min_row:min_row+sNy,min_col:min_col+sNx] = interp_field
    #
    #
    #         interp_grid_faces = rotate_interpolated_faces_to_domain(output_var_names, tile_number, tile_face_index_dict,
    #                                                                 sNx,sNy,
    #                                                                 interp_grid_faces, AngleCS, AngleSN)
    #
    # print('    - Stacking the interpolated fields into a compact pickup grid')
    # pickup_grid = stack_grids_to_compact_pickup(faces,interp_grid_faces,output_var_names,compact_tile_size)
    #
    # print('    - Outputting the compact pickup grid to the input directory')
    # pickup_metadata = dict(global_metadata)
    # output_dir = os.path.join(config_dir, 'L05', model_name, 'input')
    # output_file = os.path.join(output_dir, 'pickup.' + '{:010d}'.format(1))
    # dtype = '>f8'
    # pickup_metadata['timestepnumber'] = [int(1)]
    # pickup_metadata['nflds'] = [11]
    # pickup_metadata['nrecords'] = [np.shape(pickup_grid)[0]]
    # pickup_metadata['fldlist'] = output_var_names
    #
    # write_pickup_file(output_file, dtype, pickup_grid, pickup_metadata)
