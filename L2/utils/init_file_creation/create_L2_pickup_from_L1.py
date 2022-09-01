
import os
#import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from MITgcmutils import mds
import netCDF4 as nc4
import argparse
import ast
import sys

def read_L1_grid_tile_geometry(config_dir,model_name, sNx, sNy, Nr,
                               ordered_nonblank_tiles, ordered_nonblank_rotations,
                               faces, ordered_tiles_faces_dict):

    stitched_XC_grid = np.zeros((sNy*len(ordered_tiles_faces_dict), sNx*len(ordered_nonblank_tiles[1])))
    stitched_YC_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[1])))
    stitched_AngleCS_grid = np.zeros((sNy * len(ordered_tiles_faces_dict), sNx * len(ordered_nonblank_tiles[1])))
    stitched_AngleSN_grid = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[1])))
    stitched_HFac_grid = np.zeros((Nr, sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[1])))

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L1', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
                    AngleCS = ds.variables['AngleCS'][:, :]
                    AngleSN = ds.variables['AngleSN'][:, :]
                    HFac = ds.variables['HFacC'][:,:,:]
                    ds.close()

                    for i in range(ordered_nonblank_rotations[r][c]):
                        XC = np.rot90(XC)
                        YC = np.rot90(YC)
                        AngleCS = np.rot90(AngleCS)
                        AngleSN = np.rot90(AngleSN)
                        HFac = np.rot90(HFac,axes=(1,2))

                    stitched_XC_grid[r * sNy:(r + 1) * sNy,c * sNx: (c + 1) * sNx] = XC
                    stitched_YC_grid[r * sNy:(r + 1) * sNy,c * sNx: (c + 1) * sNx] = YC
                    stitched_AngleCS_grid[r * sNy:(r + 1) * sNy, c * sNx: (c + 1) * sNx] = AngleCS
                    stitched_AngleSN_grid[r * sNy:(r + 1) * sNy, c * sNx: (c + 1) * sNx] = AngleSN
                    stitched_HFac_grid[:, r * sNy:(r + 1) * sNy, c * sNx: (c + 1) * sNx] = HFac

    return(stitched_XC_grid, stitched_YC_grid, stitched_AngleCS_grid, stitched_AngleSN_grid, stitched_HFac_grid)

def read_grid_geometry_from_nc(config_dir,model_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    AngleCS = ds.variables['AngleCS'][:, :]
    AngleSN = ds.variables['AngleSN'][:, :]
    ds.close()
    return(XC,YC,AngleCS,AngleSN)

def read_pickup_file_to_compact(pickup_file_path, Nr):

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

def read_pickup_to_stiched_grid(pickup_file_path, Nr, sNx, sNy,
                                ordered_nonblank_tiles, ordered_nonblank_rotations,
                                faces, ordered_tiles_faces_dict):

    var_names,row_bounds,compact_var_grids,global_metadata = read_pickup_file_to_compact(pickup_file_path, Nr)

    var_grids = []

    for vn in range(len(var_names)):
        compact_var_grid = compact_var_grids[vn]
        if var_names[vn].lower() not in ['etan', 'detahdt', 'etah']:
            grid = np.zeros((Nr, sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
        else:
            grid = np.zeros((1, sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

        face_1 = compact_var_grid[:,:3*sNx,:]
        face_1 = np.reshape(face_1, (np.shape(face_1)[0], sNx, 3 * sNx))
        grid[:,:sNy,:] = face_1

        face_3 = np.rot90(compact_var_grid[:, 3 * sNx:, :],axes=(1,2),k=3)
        grid[:, sNy:, :] = face_3

        var_grids.append(grid)

    return(var_names, var_grids, global_metadata)

def read_grid_mask_from_nc(config_dir,model_name,var_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    if var_name in ['Uvel', 'GuNm1', 'GuNm2']:
        hfac = ds.variables['HFacW'][:,:,:]
        hfac = hfac[:,:,:-1]
    elif var_name in ['Vvel', 'GvNm1', 'GvNm2']:
        hfac = ds.variables['HFacS'][:,:,:]
        hfac = hfac[:, :-1, :]
    else:
        hfac = ds.variables['HFacC'][:, :, :]
    ds.close()
    mask = np.copy(hfac)
    mask[mask>0]=1
    return(mask)

def rotate_grids_to_natural_grids(uvel, vvel, AngleCS, AngleSN):

    def rotate_velocity_vectors_to_natural(angle_cos, angle_sin, uvel, vvel):
        zonal_velocity = np.zeros_like(uvel)
        meridional_velocity = np.zeros_like(vvel)
        for k in range(np.shape(uvel)[0]):
            zonal_velocity[k,:,:] = angle_cos * uvel[k,:,:] - angle_sin * vvel[k,:,:]
            meridional_velocity[k,:,:] = angle_sin * uvel[k,:,:] + angle_cos * vvel[k,:,:]
        return (zonal_velocity, meridional_velocity)

    zonal_uvel, meridional_vvel = rotate_velocity_vectors_to_natural(AngleCS, AngleSN,
                                                                     uvel, vvel)

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

    return(zonal_uvel, meridional_vvel)

def rotate_directional_fields_to_natural(var_grids,var_names,L1_AngleCS, L1_AngleSN):
    uvel_grid = var_grids[var_names.index('Uvel')]
    vvel_grid = var_grids[var_names.index('Vvel')]
    natural_uvel_grid, natural_vvel_grid = rotate_grids_to_natural_grids(uvel_grid, vvel_grid, L1_AngleCS, L1_AngleSN)
    var_grids[var_names.index('Uvel')] = natural_uvel_grid
    var_grids[var_names.index('Vvel')] = natural_vvel_grid

    gunm1_grid = var_grids[var_names.index('GuNm1')]
    gvnm1_grid = var_grids[var_names.index('GvNm1')]
    natural_gunm1_grid, natural_gvnm1_grid = rotate_grids_to_natural_grids(gunm1_grid, gvnm1_grid, L1_AngleCS,
                                                                           L1_AngleSN)
    var_grids[var_names.index('GuNm1')] = natural_gunm1_grid
    var_grids[var_names.index('GvNm1')] = natural_gvnm1_grid

    gunm2_grid = var_grids[var_names.index('GuNm2')]
    gvnm2_grid = var_grids[var_names.index('GvNm2')]
    natural_gunm2_grid, natural_gvnm2_grid = rotate_grids_to_natural_grids(gunm2_grid, gvnm2_grid, L1_AngleCS,
                                                                           L1_AngleSN)
    var_grids[var_names.index('GuNm2')] = natural_gunm2_grid
    var_grids[var_names.index('GvNm2')] = natural_gvnm2_grid
    return(var_grids,var_names)


def rotate_interpolated_grids_to_domain(var_names, var_grids, AngleCS, AngleSN):

    def rotate_velocity_vectors_to_domain(angle_cos, angle_sin, zonal_vel, meridional_vel):
        uvel = np.zeros_like(zonal_vel)
        vvel = np.zeros_like(meridional_vel)
        for k in range(np.shape(uvel)[0]):
            uvel[k, : ,:] = angle_cos * zonal_vel[k,:,:] + angle_sin * meridional_vel[k,:,:]
            vvel[k, : ,:] = -1 * angle_sin * zonal_vel[k,:,:] + angle_cos * meridional_vel[k,:,:]
        return (uvel, vvel)

    uvel_grid_index = var_names.index('Uvel')
    vvel_grid_index = var_names.index('Vvel')
    uvel, vvel = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                   var_grids[uvel_grid_index], var_grids[vvel_grid_index])

    gunm1_grid_index = var_names.index('GuNm1')
    gvnm1_grid_index = var_names.index('GvNm1')
    gunm1, gvnm1 = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                     var_grids[gunm1_grid_index], var_grids[gvnm1_grid_index])

    gunm2_grid_index = var_names.index('GuNm2')
    gvnm2_grid_index = var_names.index('GvNm2')
    gunm2, gvnm2 = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                     var_grids[gunm2_grid_index], var_grids[gvnm2_grid_index])

    var_grids[uvel_grid_index] = uvel
    var_grids[vvel_grid_index] = vvel

    var_grids[gunm1_grid_index] = gunm1
    var_grids[gvnm1_grid_index] = gvnm1

    var_grids[gunm2_grid_index] = gunm2
    var_grids[gvnm2_grid_index] = gvnm2

    return(var_grids)

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

########################################################################################################################


def create_pickup_from_L1(config_dir, L1_model_name, L1_iteration, L2_model_name,
                           sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                           faces, ordered_tiles_faces_dict, print_level):

    if print_level>=1:
        print('    - Creating the pickup file for the ' + L2_model_name + ' model from the '+L1_model_name+' model')

    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import downscale_functions as df

    # this is the dir where the output will be stored
    output_dir = os.path.join(config_dir,'L2',L2_model_name,'input')

    grid_file = os.path.join(config_dir, 'nc_grids', L1_model_name + '_grid.nc')
    ds = nc4.Dataset(grid_file)
    delR_in = np.array(ds.variables['drF'][:])
    ds.close()
    Nr_in = len(delR_in)

    grid_file = os.path.join(config_dir, 'nc_grids', L2_model_name + '_grid.nc')
    ds = nc4.Dataset(grid_file)
    delR_out = np.array(ds.variables['drF'][:])
    ds.close()
    Nr_out = len(delR_out)

    if print_level>=1:
        print('    - Reading in the geometry of the L1 domain')
    L1_XC, L1_YC, L1_AngleCS, L1_AngleSN, L1_Hfac = read_L1_grid_tile_geometry(config_dir,L1_model_name,sNx, sNy, Nr_in,
                                                                        ordered_nonblank_tiles, ordered_nonblank_rotations,
                                                                        faces, ordered_tiles_faces_dict)


    # plt.subplot(2,3,1)
    # C = plt.imshow(L1_XC, origin='lower')
    # plt.colorbar(C)
    # plt.subplot(2, 3, 2)
    # C = plt.imshow(L1_YC, origin='lower')
    # plt.colorbar(C)
    # plt.subplot(2, 3, 3)
    # C = plt.imshow(L1_Hfac[0,:,:], origin='lower')
    # plt.colorbar(C)
    # plt.subplot(2, 3, 4)
    # C = plt.imshow(L1_AngleCS, origin='lower')
    # plt.colorbar(C)
    # plt.subplot(2, 3, 5)
    # C = plt.imshow(L1_AngleSN, origin='lower')
    # plt.colorbar(C)
    # plt.show()

    L2_XC, L2_YC, L2_AngleCS, L2_AngleSN = read_grid_geometry_from_nc(config_dir,L2_model_name)

    pickup_file = 'pickup.'+'{:010d}'.format(L1_iteration)
    pickup_file_path = os.path.join(config_dir,'L1',L1_model_name,'run',pickup_file)
    var_names, var_grids, pickup_metadata = read_pickup_to_stiched_grid(pickup_file_path, Nr_in, sNx, sNy,
                                                                        ordered_nonblank_tiles, ordered_nonblank_rotations,
                                                                        faces, ordered_tiles_faces_dict)

    var_grids, var_names = rotate_directional_fields_to_natural(var_grids, var_names, L1_AngleCS, L1_AngleSN)

    if print_level >= 1:
        print('    - Downscaling the pickup grids')
    interp_grids = []
    output_var_names = []
    for vn in range(len(var_names)):
        var_name = var_names[vn]

        if var_name not in []:  # used for testing
            var_grid = var_grids[vn]
            if print_level >= 2:
                print('        - Downscaling ' + var_name)

            domain_wet_cells_3D = read_grid_mask_from_nc(config_dir,L2_model_name,var_name)

            # mean_vertical_difference = 0
            # subset_copy = np.copy(var_grid)

            L1_wet_cells = np.copy(L1_Hfac)
            L1_wet_cells[L1_wet_cells > 0] = 1

            if var_name.lower() not in ['etan', 'detahdt', 'etah']:
                if Nr_in!=Nr_out:
                    var_grid, L1_wet_cells = df.interpolate_var_grid_faces_to_new_depth_levels(
                        var_grid, L1_wet_cells, delR_in, delR_out)
            else:
                domain_wet_cells_3D = domain_wet_cells_3D[:1,:,:]

            # plt.subplot(1,2,1)
            # plt.imshow(subset_copy[:,10,:])
            # plt.subplot(1, 2, 2)
            # plt.imshow(aste_grid[:, 10, :])
            # plt.show()

            L1_wet_cells_on_domain_3D = np.copy(domain_wet_cells_3D)

            if print_level >= 3:
                print('            - L1_XC shape: '+str(np.shape(L1_XC)))
                print('            - L1_YC shape: ' + str(np.shape(L1_YC)))
                print('            - var_grid shape: ' + str(np.shape(var_grid)))
                print('            - L1_wet_cells shape: ' + str(np.shape(L1_wet_cells)))
                print('            - L2_XC shape: ' + str(np.shape(L2_XC)))
                print('            - L2_YC shape: ' + str(np.shape(L2_YC)))
                print('            - L1_XC shape: ' + str(np.shape(domain_wet_cells_3D)))

            interp_field = df.downscale_3D_field(L1_XC, L1_YC,
                                                 var_grid, L1_wet_cells,
                                                 L1_wet_cells_on_domain_3D,
                                                 L2_XC, L2_YC, domain_wet_cells_3D,
                                                 mean_vertical_difference=0, fill_downward=True, remove_zeros=True,
                                                 printing=True)
            if print_level >= 3:
                print('            - Field output shape: '+str(np.shape(interp_field)))

            if np.sum(np.isnan(interp_field)) > 0:
                print('Setting ' + str(np.sum(np.isnan(interp_field))) + ' values to 0 in this grid')
                interp_field[np.isnan(interp_field)] = 0

            # plt.subplot(2, 3, 1)
            # plt.imshow(L1_wet_cells[0, :, :], origin='lower')
            # plt.subplot(2, 3, 2)
            # plt.imshow(L1_wet_cells_on_domain_3D[0, :, :], origin='lower')
            # plt.subplot(2, 3, 3)
            # plt.imshow(domain_wet_cells_3D[0, :, :], origin='lower')
            # plt.subplot(2, 2, 3)
            # plt.imshow(var_grid[0, :, :], origin='lower')
            # plt.subplot(2, 2, 4)
            # plt.imshow(interp_field[0, :, :], origin='lower')
            # plt.show()

            # interp_field[domain_wet_cells_3D[:np.shape(interp_field)[0], :, :] == 0] = 0
        else:
            if var_name.lower() not in ['etan', 'detahdt', 'etah']:
                interp_field = np.zeros((len(delR_out), np.shape(L2_XC)[0], np.shape(L2_XC)[1]))
            else:
                interp_field = np.zeros((1, np.shape(L2_XC)[0], np.shape(L2_XC)[1]))

        # if var_name.lower() in ['etan','etah']:
        #     interp_field[interp_field!=0]+=2

        interp_grids.append(interp_field)

    # old_grids = [np.copy(interp_grids)[0],np.copy(interp_grids)[1]]

    interp_grids = rotate_interpolated_grids_to_domain(var_names, interp_grids, L2_AngleCS, L2_AngleSN)

    # plt.subplot(2,2,1)
    # plt.imshow(old_grids[0][0,:,:],origin='lower',vmin=-0.5,vmax=0.5,cmap='seismic')
    # plt.subplot(2, 2, 2)
    # plt.imshow(old_grids[1][0, :, :], origin='lower', vmin=-0.5, vmax=0.5, cmap='seismic')
    # plt.subplot(2, 2, 3)
    # plt.imshow(interp_grids[0][0, :, :], origin='lower', vmin=-0.5, vmax=0.5, cmap='seismic')
    # plt.subplot(2, 2, 4)
    # plt.imshow(interp_grids[1][0, :, :], origin='lower', vmin=-0.5, vmax=0.5, cmap='seismic')
    # plt.show()

    pickup_grid = stack_grids_to_pickup(interp_grids)

    output_dir = os.path.join(config_dir, 'L2', L2_model_name, 'input')
    output_file = os.path.join(output_dir, 'pickup.' + '{:010d}'.format(L1_iteration*2))
    dtype = '>f8'
    pickup_metadata['timestepnumber'] = [int(1)]
    pickup_metadata['nFlds'] = [11]
    pickup_metadata['nrecords'] = [int((len(var_names) - 3) * len(delR_out) + 3)]
    pickup_metadata['fldlist'] = var_names
    write_pickup_file(output_file, dtype, pickup_grid, pickup_metadata)



