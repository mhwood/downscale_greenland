
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys


def create_BCs(config_dir, proc_id,
               start_year, final_year, start_month,
               final_month, start_day, final_day):

    sys.path.insert(1, os.path.join(config_dir,'L2', 'utils','init_file_creation'))
    import create_L2_daily_exf_from_ref as cex

    L1_model_name = 'L05_CE_Greenland'
    L2_model_name = 'L2_CE_Greenland'

    sNx = 90
    sNy = 90
    ordered_nonblank_tiles = [[1, 2, 3], [6, 5, 4]]
    ordered_nonblank_rotations = [[0, 0, 0], [3, 3, 3]]

    faces = [1,3]
    ordered_tiles_faces_dict = {1:[[1,2,3]],
                                3:[[4],[5],[6]]}

    print('Creating the BCs for the '+L2_model_name+' model')

    # pass to general function to generate mitgrid
    cex.create_exf_fields_via_interpolation(config_dir, L1_model_name, L2_model_name, proc_id,
                                            sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                                            faces, ordered_tiles_faces_dict,
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

    create_BCs(config_dir,proc_id,
               start_year, final_year, start_month,
               final_month, start_day, final_day)
   

