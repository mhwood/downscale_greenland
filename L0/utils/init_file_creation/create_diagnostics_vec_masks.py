import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import simplegrid as sg
from scipy.interpolate import griddata
from ecco_v4_py.llc_array_conversion import llc_faces_to_compact, llc_compact_to_faces
import argparse
import ast


def read_global_XC_YC(ecco_dir,llc):

    grid_file_dir = os.path.join(ecco_dir,'mitgrid_tiles')
    # grid_file_dir = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files','mitgrid_tiles')
    XC_faces = {}
    YC_faces = {}
    for i in range(1,6):
        if i<3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc, 3*llc)
        if i==3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc, llc)
        if i>3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), 3*llc, llc)
        XC_face = grid_dict['XC'].T
        YC_face = grid_dict['YC'].T
        XC_faces[i] = XC_face
        YC_faces[i] = YC_face

    input_init_dir = os.path.join(ecco_dir, 'input_init')
    bathy_compact = np.fromfile(os.path.join(input_init_dir, 'bathy_llc' + str(llc)), '>f4')
    bathy_compact = np.reshape(bathy_compact, (13 * llc, llc))
    bathy_faces = llc_compact_to_faces(bathy_compact, less_output=True)

    return(XC_faces,YC_faces, bathy_faces)

def read_grid_tile_geometry(config_dir,model_name,ordered_nonblank_tiles,tile_face_index_dict):
    ordered_XC_tiles = []
    ordered_YC_tiles = []

    N = len(ordered_nonblank_tiles) * len(ordered_nonblank_tiles[0])

    grid_dir = os.path.join(config_dir, 'L05', model_name, 'run_for_grid')

    for r in range(len(ordered_nonblank_tiles)):
        row_XCs = []
        row_YCs = []
        for c in range(len(ordered_nonblank_tiles[0])):
            tile_number = ordered_nonblank_tiles[r][c]
            tile_face = tile_face_index_dict[tile_number][0]
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

    return(ordered_XC_tiles, ordered_YC_tiles)

def read_grid_tile_geometry_to_boundaries(ordered_nonblank_tiles, sNx, sNy,
                                          ordered_XC_tiles, ordered_YC_tiles,
                                          northern_tiles, southern_tiles, eastern_tiles, western_tiles):

    if len(northern_tiles)>0:
        northern_points = np.zeros((len(northern_tiles)*sNx,2))
        counter = 0
        for tile_number in northern_tiles:
            for r in range(len(ordered_nonblank_tiles)):
                for c in range(len(ordered_nonblank_tiles[0])):
                    if ordered_nonblank_tiles[r][c]==tile_number:
                        northern_points[counter:counter + sNx, 0] = ordered_XC_tiles[r][c][-1,:]
                        northern_points[counter:counter + sNx, 1] = ordered_YC_tiles[r][c][-1, :]
                        counter += sNx
        # plt.plot(northern_points[:,0],northern_points[:,1])
        # plt.show()
    else:
        northern_points = []

    if len(southern_tiles)>0:
        southern_points = np.zeros((len(southern_tiles)*sNx,2))
        counter = 0
        for tile_number in southern_tiles:
            for r in range(len(ordered_nonblank_tiles)):
                for c in range(len(ordered_nonblank_tiles[0])):
                    if ordered_nonblank_tiles[r][c]==tile_number:
                        southern_points[counter:counter + sNx, 0] = ordered_XC_tiles[r][c][0,:]
                        southern_points[counter:counter + sNx, 1] = ordered_YC_tiles[r][c][0, :]
                        counter += sNx
        # plt.plot(southern_points[:,0],southern_points[:,1])
        # plt.show()
    else:
        southern_points = []

    if len(eastern_tiles)>0:
        eastern_points = np.zeros((len(eastern_tiles)*sNy,2))
        counter = 0
        for tile_number in eastern_tiles:
            for r in range(len(ordered_nonblank_tiles)):
                for c in range(len(ordered_nonblank_tiles[0])):
                    if ordered_nonblank_tiles[r][c]==tile_number:
                        eastern_points[counter:counter + sNy, 0] = ordered_XC_tiles[r][c][:,-1]
                        eastern_points[counter:counter + sNy, 1] = ordered_YC_tiles[r][c][:,-1]
                        counter += sNy
        # plt.plot(eastern_points[:,0],eastern_points[:,1])
        # plt.show()
    else:
        eastern_points = []

    if len(western_tiles)>0:
        western_points = np.zeros((len(western_tiles)*sNy,2))
        counter = 0
        for tile_number in western_tiles:
            for r in range(len(ordered_nonblank_tiles)):
                for c in range(len(ordered_nonblank_tiles[0])):
                    if ordered_nonblank_tiles[r][c]==tile_number:
                        western_points[counter:counter + sNy, 0] = ordered_XC_tiles[r][c][:,0]
                        western_points[counter:counter + sNy, 1] = ordered_YC_tiles[r][c][:,0]
                        counter += sNy
        # plt.plot(western_points[:,0],western_points[:,1])
        # plt.show()
    else:
        western_points = []

    return(northern_points, southern_points, eastern_points, western_points)


def subset_tile_geometry_to_boundary(boundary, tile_number,ordered_nonblank_tiles,
                                     ordered_XC_tiles, ordered_YC_tiles):

    # get the geometry for this tile
    for r in range(len(ordered_nonblank_tiles)):
        for c in range(len(ordered_nonblank_tiles[0])):
            if ordered_nonblank_tiles[r][c] == tile_number:
                XC = ordered_XC_tiles[r][c]
                YC = ordered_YC_tiles[r][c]

    # subset to the boundary
    if boundary=='south':
        XC_subset = XC[:1, :]
        YC_subset = YC[:1, :]

    if boundary=='west':
        XC_subset = XC[:, :1]
        YC_subset = YC[:, :1]

    if boundary=='north':
        XC_subset = XC[-1:, :]
        YC_subset = YC[-1:, :]

    if boundary=='east':
        XC_subset = XC[:, -1:]
        YC_subset = YC[:, -1:]

    return(XC_subset, YC_subset)

def create_boundary_masks(llc, resolution, XC_faces, YC_faces, bathy_faces,
                         northern_points, southern_points, eastern_points, western_points):

    all_masks = []
    mask_dicts = []

    for points in [northern_points,southern_points,eastern_points,western_points]:



        if len(points)>0:
            mask_faces = {}
            mask_faces[1] = np.zeros((3 * llc, llc))
            mask_faces[2] = np.zeros((3 * llc, llc))
            mask_faces[3] = np.zeros((llc, llc))
            mask_faces[4] = np.zeros((llc,3 * llc))
            mask_faces[5] = np.zeros((llc,3 * llc))

            mask_dict = np.zeros((llc * llc * 5, 3))

            counter = 1
            for face in [1,2,3,4,5]:
                XC = XC_faces[face]
                YC = YC_faces[face]
                bathy = bathy_faces[face]
                for n in range(np.shape(points)[0]):
                    x = points[n,0]
                    y = points[n,1]
                    dist = ((XC-x)**2 + (YC-y)**2)**0.5
                    rows, cols = np.where(dist <= resolution)
                    if len(rows)>0:
                        for i in range(len(rows)):
                            row = rows[i]
                            col = cols[i]
                            if bathy[row,col]<0:
                                if mask_faces[face][row,col]==0:
                                    mask_faces[face][row,col] = counter
                                    counter += 1
                                    mask_dict[counter-1,0] = face
                                    mask_dict[counter-1,1] = row
                                    mask_dict[counter-1,2] = col
            # print(counter-1)
            mask_dict = mask_dict[:np.sum(mask_dict[:, 0] != 0), :]

            all_masks.append(mask_faces)
            mask_dicts.append(mask_dict)
        else:
            all_masks.append({})
            mask_dicts.append([])

    return(all_masks, mask_dicts)

def output_mask_dictionary_to_nc(output_dir,output_file_name,all_mask_dicts,mask_names):
    if output_file_name in os.listdir(output_dir):
        os.remove(os.path.join(output_dir,output_file_name))

    ds = nc4.Dataset(os.path.join(output_dir,output_file_name),'w')

    for m in range(len(mask_names)):
        print(mask_names[m],np.shape(all_mask_dicts[m]))
        if len(all_mask_dicts[m])>0:
            grp = ds.createGroup(mask_names[m])
            grp.createDimension('n_points', np.shape(all_mask_dicts[m])[0])
            var = grp.createVariable('source_faces', 'i4', ('n_points',))
            var[:] = all_mask_dicts[m][:, 0].astype(int)
            var = grp.createVariable('source_rows', 'i4', ('n_points',))
            var[:] = all_mask_dicts[m][:, 1].astype(int)
            var = grp.createVariable('source_cols', 'i4', ('n_points',))
            var[:] = all_mask_dicts[m][:, 2].astype(int)

    ds.close()



########################################################################################################################

def create_dv_masks(config_dir, ecco_path,print_status):

    if print_status:
        print('Creating the diagnostics_vec masks to use in the L0 domain')

    llc = 270
    resolution = 1/3

    level = 'L05'
    model_name = 'L05_CE_Greenland'
    sNx = 90
    sNy = 90
    tile_face_index_dict = {1: [1, 0, 0],
                            2: [1, 0, sNx],
                            3: [1, 0, 2 * sNx],
                            4: [3, 0, 0],
                            5: [3, sNy, 0],
                            6: [3, 2 * sNy, 0]}
    ordered_nonblank_tiles = [[1, 2, 3], [6, 5, 4]]
    northern_tiles = []
    southern_tiles = [1, 2, 3, 4]
    eastern_tiles = [3, 4, 5, 6]
    western_tiles = [1]

    ###############################################################################################
    # Read in the grids

    if print_status:
        print('    Reading in the L0 domain files')

    # read the mitgrids to faces
    ecco_dir = os.path.join(ecco_path,'LLC'+str(llc)+'_Files')
    XC_faces, YC_faces, bathy_faces = read_global_XC_YC(ecco_dir, llc)

    # step 0: get the model domain
    print('    - Reading in the model geometry')
    ordered_XC_tiles, ordered_YC_tiles = read_grid_tile_geometry(config_dir, model_name, ordered_nonblank_tiles, tile_face_index_dict)

    ###############################################################################################
    # Find the boundary points

    northern_points, southern_points, eastern_points, western_points = \
        read_grid_tile_geometry_to_boundaries(ordered_nonblank_tiles, sNx, sNy,
                                              ordered_XC_tiles, ordered_YC_tiles,
                                              northern_tiles, southern_tiles, eastern_tiles, western_tiles)

    ###############################################################################################
    # Mask masks of the boundary points in the faces

    mask_names = ['north','south','east','west']
    all_masks, mask_dicts = create_boundary_masks(llc, resolution, XC_faces, YC_faces, bathy_faces,
                                   northern_points, southern_points, eastern_points, western_points)

    ###############################################################################################
    # Convert the masks to compact and output as a binary file

    for m in range(len(mask_names)):
        mask_faces = all_masks[m]
        if 1 in mask_faces.keys():
            print('   - Outputting the '+mask_names[m]+' mask to binary')
            compact_mask = llc_faces_to_compact(mask_faces)
            output_file = os.path.join(config_dir, 'L0', 'input', mask_names[m]+'_boundary_mask.bin')
            compact_mask.ravel('C').astype('>f4').tofile(output_file)

    output_dir = os.path.join(config_dir, 'L0', 'input')
    output_file_name = 'L0_dv_mask_reference_dict.nc'
    output_mask_dictionary_to_nc(output_dir,output_file_name,mask_dicts,mask_names)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The name of the configuration.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=True)

    parser.add_argument("-p", "--print_status", action="store",
                        help="Print status of routine (1 for True, 0 for False).", dest="print_status",
                        type=int, required=False, default=1)

    args = parser.parse_args()
    ecco_path = args.ecco_path
    config_dir = args.config_dir
    print_status = args.print_status

    if print_status>0:
        print_status=True
    else:
        print_status=False

    create_dv_masks(config_dir,ecco_path,print_status)
