
import os
import argparse
import sys


def create_boundary_condition_files(config_dir):

    sys.path.insert(1, os.path.join(config_dir,'L3','utils','init_file_creation'))
    import create_L3_ASTE_tracer_BCs as tbc
    import create_L3_ASTE_vector_BCs as vbc

    model_name = 'L3_Scoresby_Sund'
    ordered_aste_tiles = [[14,5]]
    ordered_aste_tile_rotations = [[2,3]] # rotations are counter-clockwise

    tracer_var_names = ['THETA','SALT']
    boundaries = ['north','south','east']

    for var_name in tracer_var_names:
        for boundary in boundaries:
            tbc.create_L3_ASTE_tracer_boundary_condition(config_dir,model_name,
                                                  var_name,boundary,
                                                  ordered_aste_tiles,ordered_aste_tile_rotations)

    # note: uvel and vvel get done at the same time
    vector_var_name_sets = [['UVEL','VVEL']] # U should be first
    boundaries = ['north','south','east']

    for var_set in vector_var_name_sets:
        for boundary in boundaries:
            vbc.create_L3_ASTE_vector_boundary_condition(config_dir,model_name,
                                                  var_set[0],var_set[1],boundary,
                                                  ordered_aste_tiles,ordered_aste_tile_rotations)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_boundary_condition_files(config_dir)
   

