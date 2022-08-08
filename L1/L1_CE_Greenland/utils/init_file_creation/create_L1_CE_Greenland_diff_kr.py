
import os
import argparse
import sys


def create_pickup_file(config_dir):

    sys.path.insert(1, os.path.join(config_dir,'L1','utils','init_file_creation'))
    import create_L1_ASTE_pickup as cap

    model_name = 'L1_CE_Greenland'
    llc = 1080

    sNx = 180
    sNy = 180

    tile_face_index_dict = {1: [1, 0, 0],
                            2: [1, 0, sNx],
                            3: [1, 0, 2 * sNx],
                            4: [3, 0, 0],
                            5: [3, sNy, 0],
                            6: [3, 2 * sNy, 0]}

    ordered_nonblank_tiles = [[1, 2, 3], [6, 5, 4]]

    ordered_aste_tiles = [[5,6],[14,11]]
    ordered_aste_tile_rotations = [[0,0],[3,3]] # rotations are counter-clockwise

    cap.create_L1_ASTE_pickup_file(config_dir,model_name,
                                   sNx,sNy,ordered_nonblank_tiles,tile_face_index_dict,
                                   ordered_aste_tiles, ordered_aste_tile_rotations)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

