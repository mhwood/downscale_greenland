
import os
import argparse
import sys
import numpy as np

def create_L2_Kanger_files(config_dir):

    L2_model_name = 'L2_Kanger'

    parent_model_level = 'L1'
    parent_model_name = 'L1_CE_Greenland'
    parent_model_pickup_iteration = 2419488

    sys.path.insert(1, os.path.join(config_dir, 'L2', 'utils', 'init_file_creation'))
    sys.path.insert(1, os.path.join(config_dir, 'L2', L2_model_name, 'utils', 'init_file_creation'))

    runoff_dir = '/Users/mhwood/Documents/Research/Data Repository/Greenland/Runoff'
    termpicks_file = '/Users/mhwood/Documents/Research/Data Repository/Greenland/Ice Fronts/TermPicks_V1'
    glacier_IDs = [144, 145, 146, 147, 149, 152, 156, 157, 158, 160, 161, 162, 163, 164, 167, 170, 171, 172]
    glacier_names = []

    print_level = 5

    # steps = [3,4,8,9]
    # steps = [6]
    steps = [5, 7]
    steps = [7]

    # step 1: make the grid
    if 1 in steps:
        print('Step 1: Creating the mitgrid file for the ' + L2_model_name + ' model')
        import create_L2_Kanger_mitgrid as cm
        cm.create_mitgrid(config_dir, print_level)

    # step 2: make the bathymetry
    if 2 in steps:
        print('Step 2: Creating the bathymetry file for the ' + L2_model_name + ' model')
        import create_L2_Kanger_bathymetry as cb
        cb.create_bathy_file(config_dir, print_level)

    #########################################################################################################
    # After the bathymetry is made, the *_for_grid model should be run to make the grid
    # This bit checks that the grid file has been made
    if np.any(np.array(steps)>2):
        grid_path = os.path.join(config_dir, 'nc_grids', L2_model_name+'_grid.nc')
        if not os.path.exists(grid_path):
            sys.exit("Need to make the grid netcdf file for this model for reference\n"
                     "   - Run the stitch_L2_Kanger_nc_grid_files_for_ref.py script with the -f 1 flag (and -r 150 -c 210)\n"
                     "   - Run the *for_grid model\n"
                     "   - Run the stitch_L2_Kanger_nc_grid_files_for_ref.py script")

    # step 3: make the initial conditions
    if 3 in steps:
        print('Step 3: Creating the pickup (initial conditions) file for the ' + L2_model_name + ' model')
        import create_L2_Kanger_pickup as cp
        cp.create_pickup_file(config_dir, L2_model_name,
                              parent_model_level, parent_model_name, parent_model_pickup_iteration, print_level)

    # step 4: make the seaice initial conditions
    if 4 in steps:
        print('Step 4: Creating the seaice pickup (initial conditions) file for the ' + L2_model_name + ' model')
        import create_L2_Kanger_seaice_pickup as cp
        cp.create_seaice_pickup_file(config_dir, L2_model_name,
                                     parent_model_level, parent_model_name, parent_model_pickup_iteration, print_level)

    # step 5: make the external forcing conditions
    if 5 in steps:
        print('Step 5: Creating the external forcing conditions for the ' + L2_model_name + ' model')
        # import create_L2_Kanger_exf as ce
        # ce.create_exf_files(config_dir, L2_model_name, parent_model_level, parent_model_name, print_level)
        import create_L2_iceplume_files as cip
        cip.create_L2_iceplume_files(config_dir, L2_model_name, runoff_dir, termpicks_file, glacier_IDs,
                                     years=np.arange(2014, 2022).tolist(), print_level=print_level, runoff_only=True)

    # step 6: make the boundary conditions
    if 6 in steps:
        print('Step 6: Creating the boundary conditions for the ' + L2_model_name + ' model')
        import create_L2_Kanger_BCs as cbc
        cbc.create_BCs(config_dir, parent_model_name, L2_model_name, print_level)

    # step 7: make the iceplume files
    if 7 in steps:
        print('Step 7: Creating the iceplume files for the ' + L2_model_name + ' model')
        import create_L2_Kanger_iceplume as cip
        cip.create_L2_iceplume_files(config_dir, runoff_dir, termpicks_file, glacier_IDs, print_level)

    # step 8: make the dv masks
    if 8 in steps:
        print('Step 8: Creating the diagnostic_vec masks for the ' + L2_model_name + ' model')
        import create_L2_Kanger_dv_masks as cdv
        cdv.create_dv_masks(config_dir, L2_model_name, print_level)

    # step 9: make the shelfice masks
    if 9 in steps:
        print('Step 9: Creating the shelfice files for the ' + L2_model_name + ' model')
        import create_L2_Kanger_shelfice as cs
        cs.create_shelfice_files(config_dir, L2_model_name, print_level)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L2 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_L2_Kanger_files(config_dir)
   

