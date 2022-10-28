
import os
import argparse
import numpy as np
import sys


def create_pickup_file(config_dir):

    sys.path.insert(1, os.path.join(config_dir,'L05','utils','init_file_creation'))
    source = 'ECCO'

    if source == 'ASTE':
        import create_L05_ASTE_BCs as cabc
        import rotate_BCs_to_domain_orientation as rot
    if source == 'ECCO':
        import create_L05_ECCO_BCs_ref as bcref
        import create_L05_ECCO_BCs_from_ref as cebc
        import combine_and_rotate_L05_daily_BC_files as crd

    model_name = 'L05_CE_Greenland'

    sNx = 90
    sNy = 90
    Nr = 50

    tile_face_index_dict = {1: [1, 0, 0],
                            2: [1, 0, sNx],
                            3: [1, 0, 2 * sNx],
                            4: [3, 0, 0],
                            5: [3, sNy, 0],
                            6: [3, 2 * sNy, 0]}

    ordered_nonblank_tiles = [[1, 2, 3], [6, 5, 4]]

    if source == 'ASTE':
        ordered_aste_tiles = [[5,6],[14,11]]
        ordered_aste_tile_rotations = [[0,0],[3,3]] # rotations are counter-clockwise
    if source == 'ECCO':
        ordered_ecco_tiles = [[109, 25, 26], [109, 61, 58]]
        ordered_ecco_tile_rotations = [[1, 0, 0], [2, 3, 3]]  # rotations are counter-clockwise

    # tile, row, col
    northern_tiles = []
    southern_tiles = [1,2,3,4]
    eastern_tiles = [3,4,5,6]
    western_tiles = [1,0,0,0] # put the 0's to have 0'd BCs where model expects it

    start_year = 2002
    final_year = 2002
    start_month = 1
    final_month = 2
    start_day = 1
    final_day = 28

    var_names = ['UVEL','VVEL','THETA','SALT','AREA','HEFF','HSNOW','UICE', 'VICE']
    # var_names = ['THETA']
    # var_names = ['UICE','VICE']
    # var_names = ['AREA']
    # var_names = []

    if source=='ECCO':
        bcref.create_bc_fields_reference_dict(config_dir)

    ####################################################################################################################
    # This part masks the BCs

    for var_name in var_names:
        if source == 'ASTE':
            ##################################################################################
            # this is for ASTE
            cabc.create_L05_BCs(config_dir,model_name,var_name,
                               sNx,sNy,ordered_nonblank_tiles,tile_face_index_dict,
                               ordered_aste_tiles, ordered_aste_tile_rotations,
                               northern_tiles, southern_tiles, eastern_tiles, western_tiles)
            ##################################################################################
        # this is for ECCO
        if source == 'ECCO':
            cebc.create_L05_BCs(config_dir, model_name, var_name,
                                Nr, sNx, sNy, ordered_nonblank_tiles, tile_face_index_dict,
                                northern_tiles, southern_tiles, eastern_tiles, western_tiles,
                                start_year, final_year, start_month,final_month, start_day, final_day)

    ####################################################################################################################
    # This part combines/rotates things

    if source == 'ASTE':
        rot.rotate_BCs(config_dir,model_name,sNx,sNy,ordered_nonblank_tiles,tile_face_index_dict,
                   northern_tiles, southern_tiles, eastern_tiles, western_tiles)
    if source == 'ECCO':
        for proc_id in range(27):
            crd.combine_L05_daily_BC_files(config_dir, model_name, Nr, sNx, sNy,ordered_nonblank_tiles,
                                           northern_tiles, southern_tiles, eastern_tiles, western_tiles,
                                           proc_id, start_year, final_year, start_month, final_month, start_day, final_day)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L05, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

