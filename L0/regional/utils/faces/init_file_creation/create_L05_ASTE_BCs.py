
import os
import numpy as np
import netCDF4 as nc4
import ast
import matplotlib.pyplot as plt
from MITgcmutils import mds
from scipy.interpolate import griddata
import sys


def read_grid_tile_geometry(config_dir,model_name,var_name,ordered_nonblank_tiles,tile_face_index_dict):
    ordered_XC_tiles = []
    ordered_YC_tiles = []
    ordered_AngleCS_tiles = []
    ordered_AngleSN_tiles = []
    ordered_hfac_tiles = []

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
        row_AngleCSs = []
        row_AngleSNs = []
        row_HFacs = []
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            tile_face = tile_face_index_dict[tile_number][0]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
                    AngleCS = ds.variables['AngleCS'][:, :]
                    AngleSN = ds.variables['AngleSN'][:, :]
                    # if tile_face < 3:
                    if var_name in ['VVEL','SIvice']:
                        hFac = 'S'
                    elif var_name in ['UVEL','SIuice']:
                        hFac = 'W'
                    else:
                        hFac = 'C'
                    # else:
                    #     if var_name == 'VVEL':
                    #         hFac = 'W'
                    #     elif var_name == 'UVEL':
                    #         hFac = 'S'
                    #     else:
                    #         hFac = 'C'
                    HFac = ds.variables['HFac' + hFac][:, :, :]
                    if hFac == 'W':
                        HFac = HFac[:,:,:-1]
                        if tile_number==1:
                            HFac[:,:,0] = HFac[:,:,1]
                    if hFac == 'S':
                        HFac = HFac[:,:-1,:]
                        if tile_number in [1,2,3,4]:
                            HFac[:,0,:] = HFac[:,1,:]
                    ds.close()
                    row_XCs.append(XC)
                    row_YCs.append(YC)
                    row_AngleCSs.append(AngleCS)
                    row_AngleSNs.append(AngleSN)
                    row_HFacs.append(HFac)
        ordered_XC_tiles.append(row_XCs)
        ordered_YC_tiles.append(row_YCs)
        ordered_AngleCS_tiles.append(row_AngleCSs)
        ordered_AngleSN_tiles.append(row_AngleSNs)
        ordered_hfac_tiles.append(row_HFacs)

    return(ordered_XC_tiles, ordered_YC_tiles, ordered_AngleCS_tiles, ordered_AngleSN_tiles, ordered_hfac_tiles, delR)

def subset_tile_geometry_to_boundary(boundary, tile_number,ordered_nonblank_tiles,
                                     ordered_XC_tiles, ordered_YC_tiles,
                                     ordered_AngleCS_tiles, ordered_AngleSN_tiles, ordered_hfac_tiles):

    # get the geometry for this tile
    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            if ordered_nonblank_tiles[r][c] == tile_number:
                XC = ordered_XC_tiles[r][c]
                YC = ordered_YC_tiles[r][c]
                AngleCS = ordered_AngleCS_tiles[r][c]
                AngleSN = ordered_AngleSN_tiles[r][c]
                hFac = ordered_hfac_tiles[r][c]

    # subset to the boundary
    if boundary=='south':
        XC_subset = XC[:1, :]
        YC_subset = YC[:1, :]
        AngleCS_subset = AngleCS[:1, :]
        AngleSN_subset = AngleSN[:1, :]
        hFac_subset = hFac[:,:1, :]

    if boundary=='west':
        XC_subset = XC[:, :1]
        YC_subset = YC[:, :1]
        AngleCS_subset = AngleCS[:, :1]
        AngleSN_subset = AngleSN[:, :1]
        hFac_subset = hFac[:,:, :1]

    if boundary=='north':
        XC_subset = XC[-1:, :]
        YC_subset = YC[-1:, :]
        AngleCS_subset = AngleCS[-1:, :]
        AngleSN_subset = AngleSN[-1:, :]
        hFac_subset = hFac[:,-1:, :]

    if boundary=='east':
        XC_subset = XC[:, -1:]
        YC_subset = YC[:, -1:]
        AngleCS_subset = AngleCS[:, -1:]
        AngleSN_subset = AngleSN[:, -1:]
        hFac_subset = hFac[:,:, -1:]

    return(XC_subset, YC_subset, AngleCS_subset, AngleSN_subset, hFac_subset)

def read_grid_tile_geometry_to_boundaries(ordered_nonblank_tiles, compact_tile_size,
                                          ordered_XC_tiles, ordered_YC_tiles,
                                          ordered_AngleCS_tiles, ordered_AngleSN_tiles, ordered_hfac_tiles,
                                          northern_tiles, southern_tiles, eastern_tiles, western_tiles):

    if len(northern_tiles)>0:
        northern_points = np.zeros((len(northern_tiles)*compact_tile_size,2))
        northern_wetgrid = np.zeros((50, len(northern_tiles) * compact_tile_size, ))
        counter = 0
        for tile_number in northern_tiles:
            for r in range(len(ordered_nonblank_tiles)):
                for c in range(len(ordered_nonblank_tiles[0])):
                    if ordered_nonblank_tiles[r][c]==tile_number:
                        northern_points[counter:counter + compact_tile_size, 0] = ordered_XC_tiles[r][c][-1,:]
                        northern_points[counter:counter + compact_tile_size, 1] = ordered_YC_tiles[r][c][-1, :]
                        northern_wetgrid[:,counter:counter + compact_tile_size] = ordered_hfac_tiles[r][c][:,-1,:]
                        northern_wetgrid[northern_wetgrid>0]=1
                        counter += compact_tile_size

        # plt.plot(northern_points[:,0],northern_points[:,1])
        # plt.show()
    else:
        northern_points = []
        northern_wetgrid = []

    if len(southern_tiles)>0:
        southern_points = np.zeros((len(southern_tiles)*compact_tile_size,2))
        southern_wetgrid = np.zeros((50, len(southern_tiles) * compact_tile_size,))
        counter = 0
        for tile_number in southern_tiles:
            for r in range(len(ordered_nonblank_tiles)):
                for c in range(len(ordered_nonblank_tiles[0])):
                    if ordered_nonblank_tiles[r][c]==tile_number:
                        southern_points[counter:counter + compact_tile_size, 0] = ordered_XC_tiles[r][c][0,:]
                        southern_points[counter:counter + compact_tile_size, 1] = ordered_YC_tiles[r][c][0, :]
                        print(np.shape(ordered_hfac_tiles[r][c]))
                        southern_wetgrid[:,counter:counter + compact_tile_size] = ordered_hfac_tiles[r][c][:, 0, :]
                        southern_wetgrid[southern_wetgrid > 0] = 1
                        counter += compact_tile_size

        # plt.plot(southern_points[:,0],southern_points[:,1])
        # plt.show()
    else:
        southern_points = []
        southern_wetgrid = []

    if len(eastern_tiles)>0:
        eastern_points = np.zeros((len(eastern_tiles)*compact_tile_size,2))
        eastern_wetgrid = np.zeros((50, len(eastern_tiles) * compact_tile_size,))
        counter = 0
        for tile_number in eastern_tiles:
            for r in range(len(ordered_nonblank_tiles)):
                for c in range(len(ordered_nonblank_tiles[0])):
                    if ordered_nonblank_tiles[r][c]==tile_number:
                        eastern_points[counter:counter + compact_tile_size, 0] = ordered_XC_tiles[r][c][:,-1]
                        eastern_points[counter:counter + compact_tile_size, 1] = ordered_YC_tiles[r][c][:,-1]
                        eastern_wetgrid[:,counter:counter + compact_tile_size] = ordered_hfac_tiles[r][c][:, :, -1]
                        eastern_wetgrid[eastern_wetgrid > 0] = 1
                        counter += compact_tile_size

        # plt.plot(eastern_points[:,0],eastern_points[:,1])
        # plt.show()
    else:
        eastern_points = []
        eastern_wetgrid = []

    if len(western_tiles)>0:
        western_points = np.zeros((len(western_tiles)*compact_tile_size,2))
        western_wetgrid = np.zeros((50, len(western_tiles) * compact_tile_size,))
        counter = 0
        for tile_number in western_tiles:
            for r in range(len(ordered_nonblank_tiles)):
                for c in range(len(ordered_nonblank_tiles[0])):
                    if ordered_nonblank_tiles[r][c]==tile_number:
                        western_points[counter:counter + compact_tile_size, 0] = ordered_XC_tiles[r][c][:,0]
                        western_points[counter:counter + compact_tile_size, 1] = ordered_YC_tiles[r][c][:,0]
                        western_wetgrid[:,counter:counter + compact_tile_size] = ordered_hfac_tiles[r][c][:, :, 0]
                        western_wetgrid[western_wetgrid > 0] = 1
                        counter += compact_tile_size

        # plt.plot(western_points[:,0],western_points[:,1])
        # plt.show()
    else:
        western_points = []
        western_wetgrid = []


    return(northern_points, southern_points, eastern_points, western_points,
           northern_wetgrid, southern_wetgrid, eastern_wetgrid, western_wetgrid)

def subset_aste_to_boundary_region(XC_subset, YC_subset,
                                   aste_XC, aste_YC, aste_hfacC, aste_grid):

    ll_dist = ((aste_XC-np.min(XC_subset))**2 + (aste_YC-np.min(YC_subset))**2)**0.5
    ll_index_row, ll_index_col = np.where(ll_dist == np.min(ll_dist))

    ul_dist = ((aste_XC - np.min(XC_subset)) ** 2 + (aste_YC - np.max(YC_subset)) ** 2) ** 0.5
    ul_index_row, ul_index_col = np.where(ul_dist == np.min(ul_dist))

    lr_dist = ((aste_XC - np.max(XC_subset)) ** 2 + (aste_YC - np.min(YC_subset)) ** 2) ** 0.5
    lr_index_row, lr_index_col = np.where(lr_dist == np.min(lr_dist))

    ur_dist = ((aste_XC - np.max(XC_subset)) ** 2 + (aste_YC - np.max(YC_subset)) ** 2) ** 0.5
    ur_index_row, ur_index_col = np.where(ur_dist == np.min(ur_dist))


    min_y_index = np.min([ll_index_row[0],ul_index_row[0],lr_index_row[0],ur_index_row[0]])
    max_y_index = np.max([ll_index_row[0], ul_index_row[0], lr_index_row[0], ur_index_row[0]])
    min_x_index = np.min([ll_index_col[0], ul_index_col[0], lr_index_col[0], ur_index_col[0]])
    max_x_index = np.max([ll_index_col[0], ul_index_col[0], lr_index_col[0], ur_index_col[0]])

    buffer = 1
    if min_x_index - buffer>=0:
        min_x_index -= buffer
    if max_x_index + buffer<= np.shape(aste_XC)[1]-1:
        max_x_index += buffer
    if min_y_index - buffer>=0:
        min_y_index -= buffer
    if max_y_index + buffer<= np.shape(aste_XC)[0]-1:
        max_y_index += buffer

    buffer = 1
    if min_x_index - buffer >= 0:
        min_x_index -= buffer
    if max_x_index + buffer <= np.shape(aste_XC)[1] - 1:
        max_x_index += buffer
    if min_y_index - buffer >= 0:
        min_y_index -= buffer
    if max_y_index + buffer <= np.shape(aste_XC)[0] - 1:
        max_y_index += buffer

    aste_XC_subset = aste_XC[min_y_index:max_y_index,min_x_index:max_x_index]
    aste_YC_subset = aste_YC[min_y_index:max_y_index,min_x_index:max_x_index]
    aste_hfacC_subset = aste_hfacC[:, min_y_index:max_y_index, min_x_index:max_x_index]
    aste_grid_subset = aste_grid[:, :, min_y_index:max_y_index, min_x_index:max_x_index]

    # plt.plot(aste_XC_subset,aste_YC_subset,'k.')
    # plt.plot(XC_subset, YC_subset, 'b.')
    # plt.show()

    return(aste_XC_subset, aste_YC_subset, aste_hfacC_subset, aste_grid_subset)


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

def read_tile_wet_cells_from_run_for_grid_dir(grid_dir,n_tiles,tile_number,hFac = 'C'):

    for n in range(n_tiles):
        if 'grid.t' + '{:03d}'.format(tile_number) + '.nc' in os.listdir(
                os.path.join(grid_dir, 'mnc_' + '{:04d}'.format(n + 1))):
            ds = nc4.Dataset(os.path.join(grid_dir, 'mnc_' + '{:04d}'.format(n + 1),
                                          'grid.t' + '{:03d}'.format(tile_number) + '.nc'))
            wet_cells = np.array(ds.variables['HFac'+hFac][:, :, :])
            ds.close()
            wet_cells[wet_cells>0] = 1
    if hFac=='S':
        wet_cells = wet_cells[:,1:,:]
    if hFac=='W':
        wet_cells = wet_cells[:,:,1:]
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

    uvel_grid_index = var_names.index('Uvel')
    vvel_grid_index = var_names.index('Vvel')
    natural_uvel = interp_grid_faces[uvel_grid_index][face][:,min_row:min_row+sNy,min_col:min_col+sNx]
    natural_vvel = interp_grid_faces[vvel_grid_index][face][:, min_row:min_row + sNy, min_col:min_col + sNx]
    uvel, vvel = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                   natural_uvel, natural_vvel)
    interp_grid_faces[uvel_grid_index][face][:, min_row:min_row + sNy, min_col:min_col + sNx] = uvel
    interp_grid_faces[vvel_grid_index][face][:, min_row:min_row + sNy, min_col:min_col + sNx] = vvel

    gunm1_grid_index = var_names.index('GuNm1')
    gvnm1_grid_index = var_names.index('GvNm1')
    natural_gunm1 = interp_grid_faces[gunm1_grid_index][face][:, min_row:min_row + sNy, min_col:min_col + sNx]
    natural_gvnm1 = interp_grid_faces[gvnm1_grid_index][face][:, min_row:min_row + sNy, min_col:min_col + sNx]
    gunm1, gvnm1 = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                     natural_gunm1,natural_gvnm1)
    interp_grid_faces[gunm1_grid_index][face][:, min_row:min_row + sNy, min_col:min_col + sNx] = gunm1
    interp_grid_faces[gvnm1_grid_index][face][:, min_row:min_row + sNy, min_col:min_col + sNx] = gvnm1

    gunm2_grid_index = var_names.index('GuNm2')
    gvnm2_grid_index = var_names.index('GvNm2')
    natural_gunm2 = interp_grid_faces[gunm2_grid_index][face][:, min_row:min_row + sNy, min_col:min_col + sNx]
    natural_gvnm2 = interp_grid_faces[gvnm2_grid_index][face][:, min_row:min_row + sNy, min_col:min_col + sNx]
    gunm2, gvnm2 = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                     natural_gunm2,natural_gvnm2)
    interp_grid_faces[gunm2_grid_index][face][:, min_row:min_row + sNy, min_col:min_col + sNx] = gunm2
    interp_grid_faces[gvnm2_grid_index][face][:, min_row:min_row + sNy, min_col:min_col + sNx] = gvnm2

    return(interp_grid_faces)

def create_L05_BCs(config_dir,model_name,var_name,sNx,sNy,ordered_nonblank_tiles,tile_face_index_dict,
                  ordered_aste_tiles, ordered_aste_tile_rotations,
                  northern_tiles, southern_tiles, eastern_tiles, western_tiles):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df
    import aste_functions as af

    print('    - Creating the BC files for the '+model_name+' model from ASTE data')

    # step 0: get the model domain
    print('    - Reading in the model geometry')
    ordered_XC_tiles, ordered_YC_tiles, ordered_AngleCS_tiles, ordered_AngleSN_tiles, ordered_hfac_tiles, delR = \
        read_grid_tile_geometry(config_dir,model_name,var_name,ordered_nonblank_tiles,tile_face_index_dict)

    # step 1: get the aste faces geometry
    print('    - Reading in the aste geometry')
    aste_dir = '/Users/michwood/Documents/Research/Data Repository/Greenland/Ocean Properties/Models/ASTE'
    aste_XC, aste_YC, aste_AngleCS, aste_AngleSN, aste_hfacC, aste_delR = \
        af.read_aste_grid_geometry(aste_dir, ordered_aste_tiles, ordered_aste_tile_rotations)

    print('    - Reading in the aste grid for variable ' + str(var_name))
    if var_name in ['UVEL','VVEL','SIuice','SIvice']:
        aste_grid = af.read_aste_field_to_stiched_grid(aste_dir, var_name, ordered_aste_tiles,
                                                       ordered_aste_tile_rotations,
                                                       aste_AngleCS, aste_AngleSN, rotate_velocity = True)
    else:
        aste_grid = af.read_aste_field_to_stiched_grid(aste_dir, var_name, ordered_aste_tiles,
                                                       ordered_aste_tile_rotations)

    # plt.imshow(aste_grid[0,0,:,:],origin='lower',cmap='seismic',vmin=-1,vmax=1)
    # plt.show()

    if np.shape(aste_grid)[1]>1:
        Nr = np.size(delR)
    else:
        Nr = 1

    for boundary in ['east','south','west']:
        if boundary == 'south':
            tile_numbers = southern_tiles
        if boundary == 'west':
            tile_numbers = western_tiles
        if boundary == 'north':
            tile_numbers = northern_tiles
        if boundary == 'east':
            tile_numbers = eastern_tiles

        if len(tile_numbers) > 0:

            if boundary in ['north','south']:
                output_grid = np.zeros((np.shape(aste_grid)[0], Nr, sNx*len(tile_numbers)))
            else:
                output_grid = np.zeros((np.shape(aste_grid)[0], Nr, sNy * len(tile_numbers)))

            for tn in range(len(tile_numbers)):
                tile_number = tile_numbers[tn]
                print('      - Downscaling tile '+str(tile_number))

                # subset the geometry to the boundary
                XC_subset, YC_subset, AngleCS_subset, AngleSN_subset, hFac_subset = \
                    subset_tile_geometry_to_boundary(boundary, tile_number,ordered_nonblank_tiles,
                                                     ordered_XC_tiles, ordered_YC_tiles,
                                                     ordered_AngleCS_tiles, ordered_AngleSN_tiles, ordered_hfac_tiles)

                # subset aste to the region around the boundary to run a quicker script
                aste_XC_subset, aste_YC_subset, aste_hfacC_subset, aste_grid_subset = \
                    subset_aste_to_boundary_region(XC_subset, YC_subset,
                                                   aste_XC, aste_YC, aste_hfacC, aste_grid)

                # use the hfacs to make the mask
                aste_wet_cells_subset = np.copy(aste_hfacC_subset)
                aste_wet_cells_subset[aste_wet_cells_subset>0]=1
                wetgrid = np.copy(hFac_subset)
                wetgrid[wetgrid > 0] = 1
                aste_wet_cells_on_domain = np.copy(wetgrid)

                # # interpolate to the new depths of this model
                # if Nr>1:
                #     aste_grid_subset, aste_wet_cells_subset = df.interpolate_var_grid_faces_to_new_depth_levels(
                #         aste_grid_subset, aste_wet_cells_subset, aste_delR, delR)

                for timestep in range(np.shape(aste_grid_subset)[0]):
                    if timestep%50 == 0:
                        print('        - Downscaling timestep '+str(timestep))
                    interp_field = df.downscale_3D_field(aste_XC_subset, aste_YC_subset,
                                                         aste_grid_subset[timestep,:,:,:], aste_wet_cells_subset,
                                                         aste_wet_cells_on_domain,
                                                         XC_subset, YC_subset, wetgrid,
                                                         mean_vertical_difference=0, fill_downward=True, remove_zeros=True,
                                                         printing=False)

                    # plt.subplot(2, 2, 1)
                    # plt.imshow(aste_grid_subset[timestep,:,:,0])
                    # plt.subplot(2, 2, 2)
                    # plt.imshow(aste_wet_cells_subset[:, :, 0])
                    # plt.subplot(2, 2, 3)
                    # plt.imshow(interp_field)
                    # plt.subplot(2, 2, 4)
                    # plt.imshow(wetgrid[:, :, 0])
                    # plt.show()

                    if boundary in ['north','south']:
                        output_grid[timestep,:,tn*sNx:(tn+1)*sNx] = interp_field[:,0,:]
                    if boundary in ['east','west']:
                        output_grid[timestep,:,tn*sNy:(tn+1)*sNy] = interp_field[:,:,0]

            # plt.subplot(2,1,1)
            # plt.imshow(output_grid[10,:,:])
            # plt.subplot(2, 1, 2)
            # plt.imshow(wetgrid[:,:,0])
            # plt.show()

            # plt.plot(output_grid[0,0,:])
            # plt.show()

            # add two more values for interpolation purposes
            output_grid = np.concatenate([output_grid[:1,:,:],output_grid,output_grid[-1:,:,:]])

            if var_name in ['UVEL','VVEL','SIuice','SIvice']:
                output_file = os.path.join(config_dir,'L05',model_name,'input','obcs','L05_BC_'+boundary+'_'+var_name+'_rotated.bin')
            else:
                output_file = os.path.join(config_dir, 'L05', model_name, 'input', 'obcs',
                                           'L05_BC_' + boundary + '_' + var_name + '.bin')

            # this is required because the model still looks for these zeros
            if boundary=='west':
                extra_zeros = np.zeros((np.shape(output_grid)[0],np.shape(output_grid)[1],np.shape(output_grid)[2]*3))
                output_grid = np.concatenate([output_grid,extra_zeros],axis=2)
            print(np.shape(output_grid))

            output_grid.ravel(order='C').astype('>f4').tofile(output_file)



