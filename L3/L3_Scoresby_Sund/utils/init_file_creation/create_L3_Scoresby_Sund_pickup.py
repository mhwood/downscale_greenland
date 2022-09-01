
import os
import argparse
import sys
import ast


def create_pickup_file(config_dir, L3_model_name, parent_model_level, parent_model_name, parent_model_pickup_iteration, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L3','utils','init_file_creation'))

    if parent_model_name == 'ASTE':
        import create_L3_ASTE_pickup as cap

        ordered_aste_tiles = [[14,5]]
        ordered_aste_tile_rotations = [[2,3]] # rotations are counter-clockwise

        cap.create_L3_ASTE_pickup_file(config_dir,L3_model_name,
                                       ordered_aste_tiles,ordered_aste_tile_rotations)

    if parent_model_level == 'L1':
        import create_L3_pickup_from_L1 as c31

        f = open(os.path.join(config_dir,'L1',parent_model_name,'namelist',parent_model_name+'_geometry.dict'))
        dict_str = f.read()
        f.close()
        size_dict = ast.literal_eval(dict_str)
        sNx = size_dict['sNx']
        sNy = size_dict['sNy']
        ordered_nonblank_tiles = size_dict['ordered_nonblank_tiles']
        ordered_nonblank_rotations = size_dict['ordered_nonblank_rotations']
        faces = size_dict['faces']
        ordered_tiles_faces_dict = size_dict['ordered_tiles_faces_dict']

        c31.create_L3_pickup_file(config_dir, parent_model_name, parent_model_pickup_iteration, L3_model_name,
                           sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                           faces, ordered_tiles_faces_dict, print_level)

    if parent_model_level == 'L2':
        import create_L3_pickup_from_L2 as c32
        c32.create_L3_pickup_file(config_dir, L3_model_name,
                                  parent_model, parent_model_pickup_iteration, print_level)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

