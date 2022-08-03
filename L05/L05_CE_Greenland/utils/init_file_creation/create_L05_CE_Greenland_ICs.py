
import os
import argparse
import numpy as np
import sys


def create_pickup_file(config_dir):

    sys.path.insert(1, os.path.join(config_dir,'L05','utils','init_file_creation'))
    import create_L05_ASTE_ICs as caic
    import rotate_ICs_to_domain_orientation as rot

    model_name = 'L05_CE_Greenland'

    sNx = 90
    sNy = 90

    faces = [1,3]
    face_shapes = [(sNy,3*sNx),(3*sNy,sNx)]
    tile_face_index_dict = {1: [1, 0, 0, 0],
                            2: [1, 0, sNx, 0],
                            3: [1, 0, 2 * sNx, 0],
                            4: [3, 0, 0, 1],
                            5: [3, sNy, 0, 1],
                            6: [3, 2 * sNy, 0, 1]}

    ordered_nonblank_tiles = [[1, 2, 3], [6, 5, 4]]
    ordered_nonblank_rotations = [[0, 0, 0], [1,1,1]]

    ordered_aste_tiles = [[27, 5, 6], [27, 14, 11]]
    ordered_aste_tile_rotations = [[1, 0, 0], [2, 3, 3]]  # rotations are counter-clockwise

    var_names = ['UVEL','VVEL','ETAN','SIarea','SIheff','SIhsnow','SIuice','SIvice'] #'THETA','SALT',
    var_names = ['UVEL','VVEL','SIuice','SIvice']

    # for var_name in var_names:
    #     caic.create_L05_ICs(config_dir,model_name,var_name,
    #                        sNx,sNy, faces,face_shapes,
    #                        ordered_nonblank_tiles,tile_face_index_dict,
    #                        ordered_aste_tiles, ordered_aste_tile_rotations)

    rot.rotate_ICs(config_dir,model_name,sNx,sNy,faces, face_shapes,
                   ordered_nonblank_tiles,ordered_nonblank_rotations,
                   tile_face_index_dict)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L05, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

