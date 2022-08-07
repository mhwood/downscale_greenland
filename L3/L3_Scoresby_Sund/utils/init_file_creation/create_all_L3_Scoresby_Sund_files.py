
import os
import argparse
import sys
import numpy as np

def create_L3_Scoresby_Sund_files(config_dir):

    L3_model_name = 'L3_Scoresby_Sund'

    parent_model = 'L2_CE_Greenland'
    parent_model_pickup_iteration = 8064

    sys.path.insert(1, os.path.join(config_dir, 'L3', 'utils', 'init_file_creation'))
    sys.path.insert(1, os.path.join(config_dir, 'L3', L3_model_name, 'utils', 'init_file_creation'))

    print_level = 3

    steps = [4]

    # step 1: make the grid
    if 1 in steps:
        print('Step 1: Creating the stretched mitgrid file for the ' + L3_model_name + ' model')
        import create_L3_Scoresby_Sund_mitgrid_stretched as cm
        cm.create_mitgrid(config_dir, print_level)

    # step 2: make the bathymetry
    if 2 in steps:
        print('Step 2: Creating the bathymetry file for the ' + L3_model_name + ' model')
        import create_L3_Scoresby_Sund_bathymetry as cb
        cb.create_bathy_file(config_dir, print_level)

    #########################################################################################################
    # After the bathymetry is made, the *_for_grid model should be run to make the grid
    # This bit checks that the grid file has been made
    if np.any(np.array(steps)>2):
        grid_path = os.path.join(config_dir, 'nc_grids', L3_model_name+'_grid.nc')
        if not os.path.exists(grid_path):
            sys.exit("Need to make the grid netcdf file for this model for reference\n"
                     "   - Run the stitch_L3_Scoresby_Sund_nc_grid_files_for_ref.py script with the -f 1 flag (and -r 150 -c 210)\n"
                     "   - Run the *for_grid model\n"
                     "   - Run the stitch_L3_Scoresby_Sund_nc_grid_files_for_ref.py script")

    # step 3: make the diff_kr file

    # step 4: make the initial conditions
    if 4 in steps:
        print('Step 4: Creating the pickup (initial conditions) file for the ' + L3_model_name + ' model')
        import create_L3_Scoresby_Sund_pickup as cp
        cp.create_pickup_file(config_dir, L3_model_name,
                              parent_model, parent_model_pickup_iteration, print_level)

    # step 5: make the seaice initial conditions
    if 5 in steps:
        print('Step 5: Creating the seaice pickup (initial conditions) file for the ' + L3_model_name + ' model')
        import create_L3_Scoresby_Sund_seaice_pickup as cp
        cp.create_seaice_pickup_file(config_dir, L3_model_name,
                                     parent_model, parent_model_pickup_iteration, print_level)

    # step 6: make the external forcing conditions
    if 6 in steps:
        print('Step 6: Creating the external forcing conditions for the ' + L3_model_name + ' model')
        import create_L3_Scoresby_Sund_exf as ce
        ce.create_exf_files(config_dir, L3_model_name, parent_model, print_level)

    # step 7: make the boundary conditions
    if 7 in steps:
        print('Step 7: Creating the boundary conditions for the ' + L3_model_name + ' model')
        import create_L3_Scoresby_Sund_BCs as cbc
        cbc.create_BCs(config_dir, L3_model_name, parent_model, print_level)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_L3_Scoresby_Sund_files(config_dir)
   

