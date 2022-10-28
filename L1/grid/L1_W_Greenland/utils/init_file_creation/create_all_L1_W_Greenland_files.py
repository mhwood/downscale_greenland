
import os
import argparse
import sys
import numpy as np

def create_L1_W_Greenland_files(config_dir):

    L1_model_name = 'L1_W_Greenland'
    L3_model_name = 'L3_Disko_Bay'

    ##########################################################################################################
    # L1 model metadata

    parent_model = 'L0'
    parent_model_pickup_iteration = 73
    parent_llc = 270
    llc = 1080

    sNx = 180
    sNy = 180

    n_rows = 720
    n_cols = 540

    global_tile_face_index_dict = {
        # top row in face 3
        247: [3, 5 * sNy, 0],
        248: [3, 5 * sNy, sNx],
        249: [3, 5 * sNy, 2 * sNx],
        # first 3 cells in row 4 in face 5
        415: [5, 3 * sNy, 0],
        416: [5, 3 * sNy, sNx],
        417: [5, 3 * sNy, 2 * sNx],
        # first 3 cells in row 5 in face 5
        433: [5, 4 * sNy, 0],
        434: [5, 4 * sNy, sNx],
        435: [5, 4 * sNy, 2 * sNx],
        # first 3 cells in row 6 in face 5
        451: [5, 5 * sNy, 0],
        452: [5, 5 * sNy, sNx],
        453: [5, 5 * sNy, 2 * sNx]}

    global_ordered_nonblank_tiles = [[417, 435, 453], [416, 434, 452], [415, 433, 451], [247, 248, 249]]

    ordered_ecco_tiles = [[101, 110], [100, 109], [62, 61]]
    ordered_ecco_tile_rotations = [[1, 1], [1, 1], [2, 2]]  # rotations are counter-clockwise
    ecco_dir = '/Users/michwood/Documents/Research/Projects/Ocean_Modeling/ECCO'

    mankoff_dir = '/Users/michwood/Documents/Research/Data Repository/Greenland/Runoff/Mankoff_liquid'

    sys.path.insert(1, os.path.join(config_dir, 'L1_grid', 'utils', 'init_file_creation'))
    sys.path.insert(1, os.path.join(config_dir, 'L1_grid', L1_model_name, 'utils', 'init_file_creation'))

    ##########################################################################################################

    print_level = 4

    steps = [9]

    # step 1: make the grid
    if 1 in steps:
        print('Step 1: Creating the mitgrid file for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_mitgrid as cm
        cm.create_L1_W_Greenland_mitgrid_file(config_dir, L1_model_name, ecco_dir,
                                               sNx, sNy, global_ordered_nonblank_tiles, global_tile_face_index_dict,
                                               print_level)

    # step 2: make the bathymetry
    if 2 in steps:
        print('Step 2: Creating the bathymetry file for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_bathymetry as cb
        cb.create_L1_W_Greenland_bathymetry(config_dir, L1_model_name, ecco_dir, n_rows, n_cols, llc, print_level)

    #########################################################################################################
    # After the bathymetry is made, the *_for_grid model should be run to make the grid
    # This bit checks that the grid file has been made
    # if np.any(np.array(steps)>2):
    #     grid_path = os.path.join(config_dir, 'nc_grids', L1_model_name+'_grid.nc')
    #     if not os.path.exists(grid_path):
    #         sys.exit("Need to make the grid netcdf file for this model for reference\n"
    #                  "   - Run the *for_grid model\n"
    #                  "   - Run the stitch_L1_W_Greenland_nc_grid_files_for_ref.py script")

    # note: there is no diff_kr in ECCOv5

    # step 3: make the initial conditions
    if 3 in steps:
        print('Step 3: Creating the pickup (initial conditions) file for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_pickup as cp
        cp.create_pickup_file(config_dir, L1_model_name,
                              ecco_dir, parent_llc, ordered_ecco_tiles, ordered_ecco_tile_rotations,
                              parent_model_pickup_iteration, print_level)

    # step 4: make the seaice initial conditions
    if 4 in steps:
        print('Step 4: Creating the seaice pickup (initial conditions) file for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_seaice_pickup as cp
        cp.create_seaice_pickup_file(config_dir, L1_model_name,
                                     ecco_dir, parent_llc, ordered_ecco_tiles, ordered_ecco_tile_rotations,
                                     parent_model_pickup_iteration, print_level)

    # step 5: make the seaice initial conditions
    if 5 in steps:
        print('Step 5: Creating the ptracer pickup (initial conditions) file for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_ptracer_pickup as cp
        cp.create_ptracer_pickup_file(config_dir, L1_model_name,
                                     ecco_dir, parent_llc, ordered_ecco_tiles, ordered_ecco_tile_rotations,
                                     parent_model_pickup_iteration, print_level)

    # step 6: make the darwmin initial conditions
    if 6 in steps:
        print('Step 6: Creating the darwin pickup (initial conditions) file for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_darwin_pickup as cp
        cp.create_darwin_pickup_file(config_dir, L1_model_name,
                                      ecco_dir, parent_llc, ordered_ecco_tiles, ordered_ecco_tile_rotations,
                                      parent_model_pickup_iteration, print_level)

    # step 7: make the external forcing conditions
    if 7 in steps:
        print('Step 7: Creating the external forcing conditions for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_exf as ce
        ce.create_exf_files(config_dir, L1_model_name, ecco_dir, mankoff_dir, parent_llc, print_level)

    # step 8: make the boundary conditions
    if 8 in steps:
        print('Step 8: Creating the boundary conditions for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_BCs as cbc
        cbc.create_BCs(config_dir, L1_model_name, ecco_dir, parent_llc, print_level)

    # step 9: make the dv masks for constructing the L2 model
    if 9 in steps:
        print('Step 9: Creating the diagnostics_vec masks for the ' + L1_model_name + ' model')
        import create_L1_W_Greenland_dv_masks as cdv
        cdv.create_dv_masks(config_dir, L1_model_name, L3_model_name, print_level)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_L1_W_Greenland_files(config_dir)
   

