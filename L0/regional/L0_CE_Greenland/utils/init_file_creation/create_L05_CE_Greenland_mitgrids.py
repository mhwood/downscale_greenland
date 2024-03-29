
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyproj import Transformer
import sys

def create_L1_CE_Greenland_mitgrid_files(config_dir, ecco_dir):
    sys.path.insert(1, os.path.join(config_dir, 'L05', 'utils', 'init_file_creation'))
    import create_L05_mitgrids as cm

    model_name = 'L05_CE_Greenland'

    sNx = 90
    sNy = 90

    tile_face_index_dict = {103: [1, 17 * sNy, 0],
                            104: [1, 17 * sNy, sNx],
                            105: [1, 17 * sNy, 2 * sNx],
                            235: [3, 3 * sNy, 0],
                            241: [3, 4 * sNy, 0],
                            247: [3, 5 * sNy, 0]}

    ordered_nonblank_tiles = [[103, 104, 105], [247, 241, 235]]

    cm.create_mitgrids(config_dir, model_name, ecco_dir,
                       sNx,sNy,ordered_nonblank_tiles,tile_face_index_dict)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-e", "--ecco_dir", action="store",
                        help="The directory where the ECCO mitgrid files are stored.", dest="ecco_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    ecco_dir = args.ecco_dir

    create_L1_CE_Greenland_mitgrid_files(config_dir,ecco_dir)
   

