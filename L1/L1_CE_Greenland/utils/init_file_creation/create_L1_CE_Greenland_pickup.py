
import os
import argparse
import sys


def create_pickup_file(config_dir, model_name,
                       sNx,sNy,ordered_nonblank_tiles,tile_face_index_dict, face_size_dict,
                       ecco_dir, llc, ordered_ecco_tiles, ordered_ecco_tile_rotations,
                       parent_model_pickup_iteration, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L1','utils','init_file_creation'))

    source = 'ECCO'
    import create_L1_ECCO_pickup as cep

    ###################################################################
    # these are the ecco functions

    cep.create_L1_ECCO_pickup_file(config_dir, model_name,
                                   sNx, sNy, ordered_nonblank_tiles, tile_face_index_dict, face_size_dict,
                                   ecco_dir, llc, ordered_ecco_tiles, ordered_ecco_tile_rotations,
                                   parent_model_pickup_iteration, print_level)

    # ###################################################################
    # # these are the aste functions
    # import create_L1_ASTE_pickup as cap
    #
    # ordered_aste_tiles = [[27, 5, 6], [27, 14, 11]]
    # ordered_aste_tile_rotations = [[1, 0, 0], [2, 3, 3]]  # rotations are counter-clockwise
    #
    # cap.create_L1_ASTE_pickup_file(config_dir,model_name,
    #                                sNx,sNy,ordered_nonblank_tiles,tile_face_index_dict,
    #                                ordered_aste_tiles, ordered_aste_tile_rotations)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

