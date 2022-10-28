
import os
import argparse
import sys


def create_pickup_file(config_dir):

    sys.path.insert(1, os.path.join(config_dir,'L05','utils','init_file_creation'))
    import create_L05_ECCO_pickup as cap

    model_name = 'L05_CE_Greenland'

    sNx = 90
    sNy = 90

    tile_face_index_dict = {1: [1, 0, 0],
                            2: [1, 0, sNx],
                            3: [1, 0, 2 * sNx],
                            4: [3, 0, 0],
                            5: [3, sNy, 0],
                            6: [3, 2 * sNy, 0]}

    ordered_nonblank_tiles = [[1, 2, 3], [6, 5, 4]]

    ###################################################################
    # these are the ecco functions

    ordered_ecco_tiles = [[109, 25, 26], [109, 61, 58]]
    ordered_ecco_tile_rotations = [[1, 0, 0], [2, 3, 3]]  # rotations are counter-clockwise

    cap.create_L05_ECCO_pickup_file(config_dir, model_name,
                                    sNx, sNy, ordered_nonblank_tiles, tile_face_index_dict,
                                    ordered_ecco_tiles, ordered_ecco_tile_rotations)

    # ###################################################################
    # # these are the aste functions
    #
    # ordered_aste_tiles = [[27, 5, 6], [27, 14, 11]]
    # ordered_aste_tile_rotations = [[1, 0, 0], [2, 3, 3]]  # rotations are counter-clockwise
    #
    # cap.create_L05_ASTE_pickup_file(config_dir,model_name,
    #                                sNx,sNy,ordered_nonblank_tiles,tile_face_index_dict,
    #                                ordered_aste_tiles, ordered_aste_tile_rotations)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L05, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

