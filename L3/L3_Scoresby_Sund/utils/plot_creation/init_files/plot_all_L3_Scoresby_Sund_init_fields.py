
import os
import argparse
import sys
import numpy as np

def plot_L3_Scoresby_Sund_init_fields(config_dir):

    L3_model_name = 'L3_Scoresby_Sund'
    pickup_iteration = 16128

    sys.path.insert(1, os.path.join(config_dir, 'L3', L3_model_name, 'utils', 'plot_creation','init_files'))

    print_level = 3

    steps = [6]

    # step 1: plot the mitgrid components
    if 1 in steps:
        print('Step 1: Plotting the mitgrid components for the ' + L3_model_name + ' model')
        import plot_L3_Scoresby_Sund_mitgrid_components as pmc
        pmc.plot_L3_Scoresby_Sund_mitgrid_components(config_dir)

    # step 2: plot the bathymetry and wetgrid
    if 2 in steps:
        print('Step 2: Plotting the bathymetry for the ' + L3_model_name + ' model')
        import plot_L3_Scoresby_Sund_bathymetry as pb
        pb.plot_L3_Scoresby_Sund_bathymetry(config_dir)

    # # step 3: make the diff_kr file

    # step 4: plot the initial conditions
    if 4 in steps:
        print('Step 4: Plotting the pickup (initial conditions) fields for the ' + L3_model_name + ' model')
        import plot_L3_Scoresby_Sund_pickup_fields as cp
        cp.plot_L3_Scoresby_Sund_pickup(config_dir, pickup_iteration)

    # step 5: plot the seaice initial conditions
    if 5 in steps:
        print('Step 5: Plotting the seaice pickup (initial conditions) fields for the ' + L3_model_name + ' model')
        import plot_L3_Scoresby_Sund_seaice_pickup_fields as csp
        csp.plot_L3_Scoresby_Sund_seaice_pickup(config_dir, pickup_iteration)

    # step 6: plot the external forcing fields at a random time step
    if 6 in steps:
        print('Step 6: Plotting the external forcing conditions for the ' + L3_model_name + ' model at a random timestep')
        import plot_L3_Scoresby_Sund_exf_fields as cef
        cef.plot_L3_Scoresby_Sund_exfs(config_dir)

    # # step 6: make the external forcing conditions
    # if 6 in steps:
    #     print('Step 6: Creating the external forcing conditions for the ' + L3_model_name + ' model')
    #     import create_L3_Scoresby_Sund_exf as ce
    #     ce.create_exf_files(config_dir, L3_model_name, parent_model, print_level)
    #
    # # step 7: make the boundary conditions
    # if 7 in steps:
    #     print('Step 7: Creating the boundary conditions for the ' + L3_model_name + ' model')
    #     import create_L3_Scoresby_Sund_BCs as cbc
    #     cbc.create_BCs(config_dir, L3_model_name, parent_model, print_level)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_L3_Scoresby_Sund_init_fields(config_dir)
   

