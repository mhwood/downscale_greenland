
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys


def create_exf(config_dir, proc_id,
               start_year, final_year, start_month,
               final_month, start_day, final_day):

    sys.path.insert(1, os.path.join(config_dir,'L2', 'utils','init_file_creation'))
    import combine_L2_daily_exf_files as cex

    L2_model_name = 'L2_CE_Greenland'

    n_rows = 240
    n_cols = 240

    print('Combining the exf fields for the '+L2_model_name+' model')

    # pass to general function to generate mitgrid
    cex.combine_L2_daily_exf_files(config_dir, L2_model_name, n_rows, n_cols, proc_id,
                                   start_year, final_year, start_month,
                                   final_month, start_day, final_day)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
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

    create_exf(config_dir,proc_id,
               start_year, final_year, start_month,
               final_month, start_day, final_day)
   

