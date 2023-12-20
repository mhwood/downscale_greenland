
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys


def create_BCs(config_dir, L1_model_name,
               ecco_dir, llc, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'L1_grid', 'utils','init_file_creation'))

    start_year = 2008
    start_month = 3
    start_day = 10

    final_year = 2008
    final_month = 3
    final_day = 10

    averaging_period = 3600
    seconds_per_iter = 1200
    n_timesteps_per_day = 24

    # this references the old faces boundary geometry files to the
    # "normal" grid geometry
    boundary_read_group_dict = {'south':['south'],
                                'west':['west'],
                                'north':['east'],
                                'east':['east', 'south']}

    ###################################################################################
    # The BC fields are created in 3 steps

    # # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    # import create_L1_BC_field_ref as ebcr
    # ebcr.create_L1_BC_ref_file(config_dir, L1_model_name, averaging_period, seconds_per_iter, print_level)

    var_names = ['THETA','SALT','UVEL','VVEL','AREA','HEFF','HSNOW','UICE','VICE']
    var_names = ['VVEL']

    # step 2: using the reference dict, organize downscaled BC into daily files
    import create_L1_ECCO_BCs_from_ref as cef
    for var_name in var_names:
        cef.create_L1_BCs(config_dir, L1_model_name, var_name,
                          ecco_dir, llc, n_timesteps_per_day, boundary_read_group_dict,
                          start_year, final_year, start_month, final_month, start_day, final_day, print_level)

    # # step 3: combine all of the BC fields into a single file
    # import combine_and_rotate_L1_daily_BC_files_annual as com
    #
    # for var_name in var_names:
    #     com.combine_L1_daily_BC_files_annual(config_dir, L1_model_name, var_name,
    #                                          start_year, final_year, n_timesteps_per_day,
    #                                          print_level)



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
   

