
import os
import argparse
import sys
import numpy as np
import ast

def create_L1_SE_Greenland_files(config_dir):

    L1_model_name = 'L1_SE_Greenland'
    L2_model_name = 'L2_Sermilik'

    ##########################################################################################################
    # L1 model metadata

    parent_model = 'L0'
    parent_model_pickup_iteration = 72
    parent_llc = 270
    llc = 1080

    f = open(os.path.join(config_dir, 'L1', L1_model_name, 'namelist', L1_model_name + '_geometry.dict'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)

    sNx = 180
    sNy = 180

    face_size_dict = {1: (2*sNy, 2*sNx), 5: (1*sNy, 2*sNx)}

    global_tile_face_index_dict = {
                                    # top left corner of face 3
                                    103: [1, 17 * sNy, 0],
                                    104: [1, 17 * sNy, sNx],
                                    97:  [1, 16 * sNy, 0],
                                    98:  [1, 16 * sNy, sNx],

                                    # top left corner of face 5
                                    451: [5, 5 * sNy, 0],
                                    452: [5, 5 * sNy, sNx]}

    global_ordered_nonblank_tiles = [[452, 97, 98], [451,103,104]]

    tile_face_index_dict = {1: [1, 0, 0],
                            2: [1, 0, 1*sNx],
                            3: [1, sNy, 0],
                            4: [1, sNy, sNx],
                            5: [5, 0, 0],
                            6: [5, 0, sNx]}

    ordered_tiles = [[6, 1, 2], [5, 3, 4]]
    ordered_tile_rotations = [[1, 0, 0], [1, 0, 0]]

    # ordered_ecco_tiles = [[109, 25, 26], [109, 61, 58]]
    ordered_ecco_tiles = [[110, 22, 23], [109, 25, 26], [61, 61, 58]]
    ordered_ecco_tile_rotations = [[1, 0, 0], [1, 0, 0], [2, 3, 3]]  # rotations are counter-clockwise
    ecco_dir = '/Users/michwood/Documents/Research/Projects/Ocean_Modeling/ECCO'


    # ordered_aste_tiles = [[27, 5, 6], [27, 14, 11]]
    # ordered_aste_tile_rotations = [[1, 0, 0], [2, 3, 3]]  # rotations are counter-clockwise

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils', 'init_file_creation'))
    sys.path.insert(1, os.path.join(config_dir, 'L1', L1_model_name, 'utils', 'init_file_creation'))

    ##########################################################################################################

    print_level = 4

    steps = [6]

    # step 1: make the grid
    if 1 in steps:
        print('Step 1: Creating the mitgrid file for the ' + L1_model_name + ' model')
        import create_L1_SE_Greenland_mitgrids as cm
        cm.create_L1_SE_Greenland_mitgrid_files(config_dir, L1_model_name, ecco_dir,
                                                sNx,sNy,global_ordered_nonblank_tiles,global_tile_face_index_dict,
                                                print_level)

    # step 2: make the bathymetry
    if 2 in steps:
        print('Step 2: Creating the bathymetry file for the ' + L1_model_name + ' model')
        import create_L1_SE_Greenland_bathymetry as cb
        cb.create_L1_SE_Greenland_bathymetry(config_dir, L1_model_name, ecco_dir, sNx, sNy, face_size_dict, llc, print_level)

    #########################################################################################################
    # After the bathymetry is made, the *_for_grid model should be run to make the grid
    # This bit checks that the grid file has been made
    # if np.any(np.array(steps)>2):
    #     grid_path = os.path.join(config_dir, 'nc_grids', L1_model_name+'_grid.nc')
    #     if not os.path.exists(grid_path):
    #         sys.exit("Need to make the grid netcdf file for this model for reference\n"
    #                  "   - Run the *for_grid model\n"
    #                  "   - Run the stitch_L1_SE_Greenland_nc_grid_files_for_ref.py script")

    # note: there is no diff_kr in ECCOv5

    # step 3: make the initial conditions
    if 3 in steps:
        print('Step 3: Creating the pickup (initial conditions) file for the ' + L1_model_name + ' model')
        import create_L1_SE_Greenland_pickup as cp
        cp.create_pickup_file(config_dir, L1_model_name,
                              sNx, sNy, ordered_tiles, tile_face_index_dict, face_size_dict,
                              ecco_dir, parent_llc, ordered_ecco_tiles, ordered_ecco_tile_rotations,
                              parent_model_pickup_iteration, print_level)

    # step 4: make the seaice initial conditions
    if 4 in steps:
        print('Step 4: Creating the seaice pickup (initial conditions) file for the ' + L1_model_name + ' model')
        import create_L1_SE_Greenland_seaice_pickup as cp
        cp.create_seaice_pickup_file(config_dir, L1_model_name,
                                     sNx,sNy,ordered_tiles,tile_face_index_dict, face_size_dict,
                                     ecco_dir, parent_llc, ordered_ecco_tiles, ordered_ecco_tile_rotations,
                                     parent_model_pickup_iteration, print_level)

    # step 5: make the external forcing conditions
    if 5 in steps:
        print('Step 5: Creating the external forcing conditions for the ' + L1_model_name + ' model')
        import create_L1_SE_Greenland_exf as ce
        ce.create_exf_files(config_dir, L1_model_name, ecco_dir, parent_llc, ordered_ecco_tiles, ordered_ecco_tile_rotations,
                            sNx, sNy, ordered_tiles, ordered_tile_rotations, print_level)

    # step 6: make the boundary conditions
    if 6 in steps:
        print('Step 6: Creating the boundary conditions for the ' + L1_model_name + ' model')
        import create_L1_SE_Greenland_BCs as cbc
        cbc.create_BCs(config_dir, L1_model_name,
                       sNx, sNy, ordered_tiles, tile_face_index_dict,
                       ecco_dir, parent_llc, print_level)

    # step 7: make the dv masks for constructing the L2 model
    if 7 in steps:
        print('Step 7: Creating the diagnostics_vec masks for the ' + L1_model_name + ' model')
        import create_L1_SE_Greenland_dv_masks as cdv
        cdv.create_dv_masks(config_dir, L1_model_name, L2_model_name, L3_model_name, sNx, sNy, face_size_dict, print_level)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_L1_SE_Greenland_files(config_dir)
   

