
import os
import argparse
import sys
import numpy as np

def create_L2_CE_Greenland_files(config_dir):

    ####################################################################################################################
    # model metadata

    L2_model_name = 'L2_CE_Greenland'

    parent_model = 'L1_CE_Greenland'
    parent_model_pickup_iteration = 1077984

    L3_model_name = 'L3_Scoresby_Sund'

    central_wet_row = 120
    central_wet_col = 239
    hFacMinDr = 1
    hFacMin = 0.1
    delR = np.array([1.00, 1.14])

    sNx = 180
    sNy = 180
    ordered_nonblank_tiles = [[1, 2, 3], [6, 5, 4]]
    ordered_nonblank_rotations = [[0, 0, 0], [3, 3, 3]]

    faces = [1, 3]
    ordered_tiles_faces_dict = {1: [[1, 2, 3]],
                                3: [[4], [5], [6]]}

    ####################################################################################################################

    sys.path.insert(1, os.path.join(config_dir, 'L2', 'utils', 'init_file_creation'))
    sys.path.insert(1, os.path.join(config_dir, 'L2', L2_model_name, 'utils', 'init_file_creation'))

    print_level = 3

    steps = [5]

    # step 1: make the grid
    if 1 in steps:
        print('Step 1: Creating the mitgrid file for the ' + L2_model_name + ' model')
        import create_L2_CE_Greenland_mitgrid as cm
        cm.create_mitgrid(config_dir, print_level)

    # step 2: make the bathymetry
    if 2 in steps:
        print('Step 2: Creating the bathymetry file for the ' + L2_model_name + ' model')
        import create_L2_CE_Greenland_bathymetry as cb
        cb.create_bathy_file(config_dir, L2_model_name, central_wet_row, central_wet_col,
                             hFacMinDr, hFacMin, delR, print_level)

    #########################################################################################################
    # After the bathymetry is made, the *_for_grid model should be run to make the grid
    # This bit checks that the grid file has been made
    if np.any(np.array(steps)>2):
        grid_path = os.path.join(config_dir, 'nc_grids', L2_model_name+'_grid.nc')
        if not os.path.exists(grid_path):
            sys.exit("Need to make the grid netcdf file for this model for reference\n"
                     "   - ln -s the bathymetry and tile into input_for_grid\n"
                     "   - Run the *for_grid model\n"
                     "   - ln -s the nc_grid to config_dir/nc_grids")

    # step 3: make the initial conditions
    if 3 in steps:
        print('Step 3: Creating the pickup (initial conditions) file for the ' + L2_model_name + ' model')
        import create_L2_CE_Greenland_pickup as cp
        cp.create_pickup_file(config_dir, parent_model, parent_model_pickup_iteration, L2_model_name,
                              sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                              faces, ordered_tiles_faces_dict, print_level)

    # step 4: make the seaice initial conditions
    if 4 in steps:
        print('Step 4: Creating the seaice pickup (initial conditions) file for the ' + L2_model_name + ' model')
        import create_L2_CE_Greenland_seaice_pickup as cp
        cp.create_seaice_pickup_file(config_dir, parent_model, parent_model_pickup_iteration, L2_model_name,
                                     sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                                     faces, ordered_tiles_faces_dict, print_level)

    # step 5: make the external forcing conditions
    if 5 in steps:
        print('Step 5: Creating the external forcing conditions for the ' + L2_model_name + ' model')
        import create_L2_CE_Greenland_exf as ce
        ce.create_exf_files(config_dir, L2_model_name, parent_model,
                            sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                            faces, ordered_tiles_faces_dict, print_level)

    # step 6: make the boundary conditions
    if 6 in steps:
        print('Step 6: Creating the boundary conditions for the ' + L2_model_name + ' model')
        import create_L2_CE_Greenland_BCs as cbc
        cbc.create_BCs(config_dir, L2_model_name, parent_model,
                       sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                       faces, ordered_tiles_faces_dict, print_level)

    # step 7: make the diagnostics_vec masks
    if 7 in steps:
        print('Step 7: Creating the diagnostics_vec masks for the ' + L2_model_name + ' model')
        import create_L2_CE_Greenland_dv_masks as cdvm
        cdvm.create_dv_masks(config_dir, L2_model_name, L3_model_name, print_level)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_L2_CE_Greenland_files(config_dir)
   

