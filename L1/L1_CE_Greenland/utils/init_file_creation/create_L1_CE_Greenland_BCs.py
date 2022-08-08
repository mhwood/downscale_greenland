
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys


def create_BCs(config_dir, L1_model_name,
               sNx, sNy, ordered_nonblank_tiles, tile_face_index_dict,
               ecco_dir, llc, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils','init_file_creation'))

    start_year = 2002
    start_month = 1
    start_day = 1

    final_year = 2002
    final_month = 1
    final_day = 3

    # tile, row, col
    northern_tiles = []
    southern_tiles = [1, 2, 3, 4]
    eastern_tiles = [3, 4, 5, 6]
    western_tiles = [1, 0, 0, 0]  # put the 0's to have 0'd BCs where model expects it

    Nr = 50

    ###################################################################################
    # The BC fields are created in 3 steps

    # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    # import create_L1_BC_field_ref as ebcr
    # ebcr.create_L1_BC_ref_file(config_dir, L1_model_name, print_level)

    var_names = ['THETA','SALT','UVEL','VVEL']
    var_names = ['AREA','HEFF','HSNOW','UICE','VICE']

    # step 2: using the reference dict, organize downscaled BC into daily files
    import create_L1_ECCO_BCs_from_ref as cef
    for var_name in var_names:
        cef.create_L1_BCs(config_dir, L1_model_name, var_name,
                          Nr, sNx, sNy, ordered_nonblank_tiles, tile_face_index_dict,
                          ecco_dir, llc,
                          northern_tiles, southern_tiles, eastern_tiles, western_tiles,
                          start_year, final_year, start_month, final_month, start_day, final_day, print_level)

    # step 3: combine all of the BC fields into a single file
    import combine_and_rotate_L1_daily_BC_files as com
    for var_name in var_names:
        com.combine_L1_daily_BC_files(config_dir, L1_model_name, var_name,
                                      Nr, sNx, sNy, ordered_nonblank_tiles,
                                      northern_tiles, southern_tiles, eastern_tiles, western_tiles,
                                      start_year, final_year, start_month, final_month, start_day, final_day,
                                      print_level)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-i", "--proc_id", action="store",
                        help="The id of the process to run.", dest="proc_id",
                        type=int, required=True)

    parser.add_argument("-S", "--start_year", action="store",
                        help="The start year.", dest="start_year",
                        type=int, required=True)

    parser.add_argument("-s", "--start_month", action="store",
                        help="The start month.", dest="start_month",
                        type=int, required=True)

    parser.add_argument("-sd", "--start_day", action="store",
                        help="The start day.", dest="start_day",
                        type=int, required=True)

    parser.add_argument("-F", "--final_year", action="store",
                        help="The final year.", dest="final_year",
                        type=int, required=True)

    parser.add_argument("-f", "--final_month", action="store",
                        help="The final ymonth.", dest="final_month",
                        type=int, required=True)

    parser.add_argument("-fd", "--final_day", action="store",
                        help="The final day.", dest="final_day",
                        type=int, required=True)


    args = parser.parse_args()
    config_dir = args.config_dir
    proc_id = args.proc_id
    start_year = args.start_year
    final_year = args.final_year
    start_month = args.start_month
    final_month = args.final_month
    start_day = args.start_day
    final_day = args.final_day

    create_BCs(config_dir,proc_id,
               start_year, final_year, start_month,
               final_month, start_day, final_day)
   

