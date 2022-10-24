
import os
import argparse
import numpy as np
import sys

def create_exf_files(config_dir, L1_model_name, ecco_dir, mankoff_dir, llc, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L1_grid','utils','init_file_creation'))

    if 'exf' not in os.listdir(os.path.join(config_dir,'L1_grid',L1_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L1_grid',L1_model_name,'input','exf'))

    start_year = 1992
    start_month = 1
    start_day = 1

    final_year = 1992
    final_month = 12
    final_day = 30

    ###################################################################################################################
    # these are the fields output from diagnostics vec (after ECCO corrections have been applied)

    # # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    # import create_L1_exf_field_ref as cer
    # cer.create_L1_exf_ref_file(config_dir, L1_model_name, print_level, read_darwin = True)

    # step 2: make monthly stacks of the exf files
    var_names = ['ATEMP', 'AQH', 'LWDOWN', 'SWDOWN', 'UWIND', 'VWIND', 'PRECIP','RUNOFF','ATMOSCO2','IRONDUST']

    # import create_L1_ECCO_exfs_from_ref as cexf
    # for var_name in var_names:
    #     cexf.create_L1_exfs(config_dir, L1_model_name, var_name, ecco_dir, llc,
    #                         start_year, final_year, start_month, final_month, start_day, final_day, print_level, read_darwin = True)

    import combine_and_rotate_L1_daily_exf_files_annual as com
    for var_name in var_names:
        com.combine_L1_daily_exf_files_annual(config_dir, L1_model_name, var_name,
                                              start_year, final_year, print_level)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

