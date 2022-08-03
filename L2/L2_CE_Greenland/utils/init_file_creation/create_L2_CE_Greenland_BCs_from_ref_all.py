
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys
import subprocess


def create_all_BCs(config_dir,
                   start_year, final_year, start_month,
                   final_month, start_day, final_day):

    path = os.path.join(config_dir,'L2','L2_CE_Greenland','utils','init_file_creation','create_L2_CE_Greenland_BCs_from_ref.py')

    for proc in range(12,27):
        command = 'python3 '+path+' -d '+config_dir+' -i '+str(proc)+' -S '+str(start_year)+' -s '+str(start_month)+' -sd '+str(start_day)+\
                  ' -F '+str(final_year)+' -f '+str(final_month)+' -fd '+str(final_day)
        print(command)
        os.system(command)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

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
    start_year = args.start_year
    final_year = args.final_year
    start_month = args.start_month
    final_month = args.final_month
    start_day = args.start_day
    final_day = args.final_day

    create_all_BCs(config_dir,
                   start_year, final_year, start_month,
                   final_month, start_day, final_day)
   

