
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
from pyproj import Transformer
import matplotlib.path as mplPath
import sys


def read_grid_tile_geometry_faces(config_dir,model_name, sNx, sNy, ordered_nonblank_tiles, tile_face_index_dict, face_size_dict):

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])
    Nr = 50

    XC_faces = {}
    YC_faces = {}
    Depth_faces = {}
    hFac_faces = {}

    for face in list(face_size_dict.keys()):
        XC_faces[face] = np.zeros(face_size_dict[face])
        YC_faces[face] = np.zeros(face_size_dict[face])
        Depth_faces[face] = np.zeros(face_size_dict[face])
        hFac_faces[face] = np.zeros((Nr,face_size_dict[face][0],face_size_dict[face][1]))

    # stitched_XC = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    # stitched_YC = np.zeros((sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))
    # stitched_HFac = np.zeros((Nr, sNy * len(ordered_nonblank_tiles), sNx * len(ordered_nonblank_tiles[0])))

    grid_dir = os.path.join(config_dir, 'L1', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            for n in range(N):
                if 'grid.t'+'{:03d}'.format(tile_number)+'.nc' in os.listdir(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1))):
                    ds = nc4.Dataset(os.path.join(grid_dir,'mnc_'+'{:04d}'.format(n+1),'grid.t'+'{:03d}'.format(tile_number)+'.nc'))
                    XC = ds.variables['XC'][:, :]
                    YC = ds.variables['YC'][:, :]
                    Depth = ds.variables['Depth'][:, :]
                    HFac = ds.variables['HFacC'][:, :, :]
                    RC = ds.variables['RC'][:]
                    ds.close()

                    # for i in range(ordered_nonblank_rotations[r][c]):
                    #     XC = np.rot90(XC)
                    #     YC = np.rot90(YC)
                    #     HFac = np.rot90(HFac,axes=(1,2))

                    dest_face = tile_face_index_dict[tile_number][0]
                    dest_row_start = tile_face_index_dict[tile_number][1]
                    dest_col_start = tile_face_index_dict[tile_number][2]

                    # stitched_XC[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = XC
                    # stitched_YC[r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = YC
                    # stitched_HFac[:, r * sNy:(r + 1) * sNy, c * sNx:(c + 1) * sNx] = HFac

                    XC_faces[dest_face][dest_row_start:dest_row_start + sNy, dest_col_start:dest_col_start + sNx] = XC
                    YC_faces[dest_face][dest_row_start:dest_row_start + sNy, dest_col_start:dest_col_start + sNx] = YC
                    Depth_faces[dest_face][dest_row_start:dest_row_start + sNy, dest_col_start:dest_col_start + sNx] = Depth
                    hFac_faces[dest_face][:, dest_row_start:dest_row_start + sNy, dest_col_start:dest_col_start + sNx] = HFac


    # for face in list(face_size_dict.keys()):
    #
    #     plt.subplot(1,3,1)
    #     C = plt.imshow(XC_faces[face],origin='lower')
    #     plt.colorbar(C)
    #
    #     plt.subplot(1, 3, 2)
    #     C = plt.imshow(YC_faces[face], origin='lower')
    #     plt.colorbar(C)
    #
    #     plt.subplot(1, 3, 3)
    #     C = plt.imshow(hFac_faces[face][0,:,:], origin='lower')
    #     plt.colorbar(C)
    #
    #     plt.show()

    return(XC_faces, YC_faces, Depth_faces, hFac_faces, RC)


def reproject_polygon(polygon_array,inputCRS,outputCRS,x_column=0,y_column=1):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    # There seems to be a serious problem with pyproj
    # The x's and y's are mixed up for these transformations
    #       For 4326->3413, you put in (y,x) and get out (x,y)
    #       Foe 3413->4326, you put in (x,y) and get out (y,x)
    # Safest to run check here to ensure things are outputting as expected with future iterations of pyproj

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif inputCRS == 3413 and outputCRS == 4326:
        y2, x2 = transformer.transform(polygon_array[:, x_column], polygon_array[:, y_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    elif str(inputCRS)[:3] == '326' and outputCRS == 3413:
        x2, y2 = transformer.transform(polygon_array[:,y_column], polygon_array[:,x_column])
        x2 = np.array(x2)
        y2 = np.array(y2)
    else:
        raise ValueError('Reprojection with this epsg is not safe - no test for validity has been implemented')

    output_polygon=np.copy(polygon_array)
    output_polygon[:,x_column] = x2
    output_polygon[:,y_column] = y2
    return output_polygon


def read_in_domain_boundary(config_dir, model_name):

    grid_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(grid_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:, :]
    ds.close()

    bottom = np.column_stack([XC[0, :], YC[0, :]])
    right = np.column_stack([XC[:, -1], YC[:, -1]])
    top = np.column_stack([XC[-1, :], YC[-1, :]])
    left = np.column_stack([XC[:,0], YC[:,0]])

    boundary = np.vstack([bottom,right,np.flipud(top),np.flipud(left)])

    boundary = reproject_polygon(boundary, 4326, 3413)

    # plt.plot(boundary[:,0],boundary[:,1])
    # plt.show()

    return(XC, YC, Depth, boundary)


def find_Mankoff_outlets_in_domain(mankoff_dir, domain_boundary):

    file_name = os.path.join(mankoff_dir,'outlets','freshwater','ice','outlets.csv')
    # f = open(file_name)
    # lines = f.read()
    # f.close()
    # lines = lines.split('\n')
    # lines.pop(0)
    outlet_locations = np.genfromtxt(file_name,skip_header=1,delimiter=',')

    points = outlet_locations[:,-2:]
    p = mplPath.Path(domain_boundary)
    inside = p.contains_points(points)

    domain_outlet_locations = outlet_locations[inside,:]
    domain_outlet_locations = domain_outlet_locations[domain_outlet_locations[:,5]<=-10,:]

    # plt.plot(domain_boundary[:, 0], domain_boundary[:, 1])
    # plt.plot(domain_outlet_locations[:,-2],domain_outlet_locations[:,-1],'k.')
    # plt.show()

    outlet_points_xyz = np.column_stack([domain_outlet_locations[:,-4:-2],domain_outlet_locations[:,5]])
    outlet_points_coast_id = domain_outlet_locations[:,0]

    return(outlet_points_xyz, outlet_points_coast_id)


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


def map_Mankoff_outlets_to_model_domain(outlet_points_xyz, XC_faces, YC_faces, Depth_faces, hFac_faces, RC):

    model_fijkxyz = np.zeros((np.shape(outlet_points_xyz)[0],7))

    for i in range(np.shape(outlet_points_xyz)[0]):

        point_err=1e10

        for face in list(XC_faces.keys()):
            domain_X = XC_faces[face]
            domain_Y = YC_faces[face]
            domain_depth = Depth_faces[face]

            xyz = np.column_stack([domain_X.ravel(), domain_Y.ravel(), -1 * domain_depth.ravel()])
            subset_xyz = xyz[xyz[:, 2] <= 0.5 * outlet_points_xyz[i, 2], :]

            distance = great_circle_distance(lon_ref = outlet_points_xyz[i,0], lat_ref = outlet_points_xyz[i,1],
                                             Lon = subset_xyz[:,0], Lat = subset_xyz[:,1])
            index = np.argmin(distance)

            if distance[index]<point_err:
                point_err = distance[index]
                # print('    - Point '+str(i))
                # print('          - Location : '+str(outlet_points_xyz[i,0])+', '+str(outlet_points_xyz[i,1])+' at depth '+str(outlet_points_xyz[i,2])+' m')
                # print('          - Closest point thats deep enough: '+str(subset_xyz[index,0])+', '+str(subset_xyz[index,1])+', at depth '+str(subset_xyz[index,2]))
                # print('                at a distance of '+str(distance[index])+' m')

                row, col = np.where(np.logical_and(domain_X==subset_xyz[index,0],domain_Y==subset_xyz[index,1]))

                k = np.argmin(np.abs(RC-outlet_points_xyz[i,2]))
                #  this ensures that the discharge will be put in a cell that's atleast 75% wet
                for n in range(np.shape(hFac_faces[face])[0]):
                    if hFac_faces[face][k,row[0],col[0]]<0.75 and k!=0:
                        k-=1

                model_fijkxyz[i, 0] = face
                model_fijkxyz[i, 1] = row[0]
                model_fijkxyz[i, 2] = col[0]
                model_fijkxyz[i, 3] = k

                model_fijkxyz[i, 4] = subset_xyz[index,0]
                model_fijkxyz[i, 5] = subset_xyz[index,1]
                model_fijkxyz[i, 6] = RC[k]#subset_xyz[index,2]

    # # plot some lines to see how far away the source and closest model cell are
    # for i in range(np.shape(model_fijkxyz)[0]):
    #     print(model_fijkxyz[i,:])
    #     plt.plot([outlet_points_xyz[i,0],model_ijxyzd[i,2]],[outlet_points_xyz[i,1],model_ijxyzd[i,3]],'k-')
    # plt.show()

    # plt.pcolormesh(domain_X,domain_Y,domain_depth,alpha=0.5)
    # plt.plot(model_ijxyzd[:,2],model_ijxyzd[:,3],'k.')
    # plt.plot(outlet_points_xyz[:,0],outlet_points_xyz[:,1],'g.')
    # plt.show()

    return(model_fijkxyz)


def create_sgd_mask(Lf, config_dir, model_name, model_fijkxyz, hFac_faces, sNx, sNy):

    mask_faces = {}
    for face in list(hFac_faces.keys()):
        mask_faces[face] = np.zeros_like(hFac_faces[face])

    points_with_multiple_sources = 0
    counter = 1
    for i in range(np.shape(model_fijkxyz)[0]):
        face = model_fijkxyz[i,0].astype(int)
        row = model_fijkxyz[i,1].astype(int)
        col = model_fijkxyz[i,2].astype(int)
        k = model_fijkxyz[i,3].astype(int)
        if mask_faces[face][k, row, col]==0:
            mask_faces[face][k, row, col] = counter
            counter+=1
        else:
            points_with_multiple_sources += 1

    print('        - '+str(points_with_multiple_sources)+' will receive fluxes from multiple sources ('+str(np.shape(model_fijkxyz)[0])+' total sources)')

    mask_compact = Lf.read_faces_to_compact(mask_faces, sNx, sNy)

    if 'sgd' not in os.listdir(os.path.join(config_dir, 'L1', model_name,'input')):
        os.mkdir(os.path.join(config_dir, 'L1', model_name,'input','sgd'))

    output_file = os.path.join(config_dir, 'L1', model_name, 'input', 'sgd', 'L1_sgd_mask.bin')
    mask_compact.ravel(order='C').astype('>f4').tofile(output_file)



def create_sgd_timeseries_file(config_dir, model_name,mankoff_dir, model_fijkxyz, outlet_points_coast_id, years):

    # find how many days will be in the output file
    n_days = 0
    for year in years:
        if year%4==0:
            n_days += 366
        else:
            n_days += 365

    # find how many output points are in the output file
    # this looks for duplicates and adds them to an existing line if its found
    counter = 0
    dest_sources = -1*np.ones(np.shape(model_fijkxyz)[0]).astype(int)
    for i in range(len(outlet_points_coast_id)):
        if i==0:
            dest_sources[i]=counter
            counter+=1
        else:
            duplicate_found = False
            for j in range(i): # search in the lines before this one
                if not duplicate_found:
                    if model_fijkxyz[i, 0] == model_fijkxyz[j, 0]:
                        if model_fijkxyz[i, 1] == model_fijkxyz[j, 1]:
                            if model_fijkxyz[i, 2] == model_fijkxyz[j, 2]:
                                if model_fijkxyz[i, 3] == model_fijkxyz[j, 3]:
                                    duplicate_found = True
                                    dest_sources[i]=dest_sources[j]
            if not duplicate_found:
                dest_sources[i] = counter
                counter+=1

    if np.any(dest_sources<0):
        raise ValueError('Some of the outlet points werent assigned sources properly')

    print('        - Found '+str(counter)+' unique points')

    output_grid = np.zeros((n_days, counter))

    days_counted = 0

    for year in years:
        print('        - Adding data from year '+str(year))

        if year%4==0:
            year_days = 366
        else:
            year_days = 365

        discharge_file = os.path.join(mankoff_dir,'runoff','ice','runoff','MAR_'+str(year)+'.nc')
        ds = nc4.Dataset(discharge_file)
        coast_id = ds.variables['station'][:]
        runoff = ds.variables['runoff'][:,:]
        ds.close()

        for i in range(len(outlet_points_coast_id)):
            id = int(outlet_points_coast_id[i])
            id_index = np.where(coast_id==id)[0]
            dst_index = dest_sources[i]
            output_grid[days_counted:days_counted+year_days,dst_index] += runoff[id_index,:].ravel()

        days_counted += year_days

    output_grid[np.isnan(output_grid)] = 0

    # C = plt.pcolormesh(output_grid)
    # plt.colorbar(C)
    # plt.show()

    output_file = os.path.join(config_dir, 'L1', model_name, 'input', 'sgd', 'L1_sgd_flux_volume.bin')
    output_grid.ravel(order='C').astype('>f4').tofile(output_file)



def create_L1_sgd_files(Lf, config_dir, model_name, mankoff_dir,
                        years, sNx, sNy, ordered_nonblank_tiles, tile_face_index_dict, face_size_dict, print_level):

    if print_level >= 1:
        print('    - Creating the subglacial discharge files for the '+model_name+' model from Mankoff data')

    # step 0: get the model domain
    if print_level >= 1:
        print('    - Reading in the model geometry as faces')
    XC_faces, YC_faces, Depth_faces, hFac_faces, RC = read_grid_tile_geometry_faces(config_dir,model_name,sNx,sNy,ordered_nonblank_tiles,
                                                                   tile_face_index_dict, face_size_dict)

    if print_level >= 1:
        print('    - Reading in the model boundary')
    domain_X, domain_Y, domain_depth, domain_boundary = read_in_domain_boundary(config_dir, model_name)

    # step 1: find all of the Mankoff outlets in the domain
    if print_level >= 1:
        print('    - Identifying subglacial discharge outlets within the domain from Mankoff references')
    outlet_points_xyz, outlet_points_coast_id = find_Mankoff_outlets_in_domain(mankoff_dir, domain_boundary)

    # step 2: create a 3D mask for the domain
    # going to try simple approach first:
    # just put the discharge into the closest mostly wet cell - might need to spread the love more if its unstable
    if print_level >= 1:
        print('    - Mapping outlet locations to model grid cells')
    model_fijkxyz = map_Mankoff_outlets_to_model_domain(outlet_points_xyz, XC_faces, YC_faces, Depth_faces, hFac_faces, RC)

    # step 3: create sgd mask
    if print_level >= 1:
        print('    - Creating a mask for use in the sgd package')
    create_sgd_mask(Lf, config_dir, model_name, model_fijkxyz, hFac_faces, sNx, sNy)

    # step 4: create the sgd file
    if print_level >= 1:
        print('    - Creating a file with a timeseries of fluxes for use in the sgd package')
    create_sgd_timeseries_file(config_dir, model_name, mankoff_dir, model_fijkxyz, outlet_points_coast_id, years)




