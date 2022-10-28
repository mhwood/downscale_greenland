
import os
import argparse
import numpy as np
import sys


def create_exf_files(config_dir):

    sys.path.insert(1, os.path.join(config_dir,'L05','utils','init_file_creation'))
    import create_L05_ECCO_exfs as cexf

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
    ordered_nonblank_rotations = [[0, 0, 0], [3, 3, 3]]

    var_names = ['ATEMP','AQH','LWDOWN','SWDOWN','UWIND','VWIND','PRECIP']#,'RUNOFF']
    file_prefixes = {'ATEMP':'EIG_tmp2m_degC_plus_ECCO_v4r1_ctrl',
                     'AQH':'EIG_spfh2m_plus_ECCO_v4r1_ctrl',
                     'LWDOWN':'EIG_dlw_plus_ECCO_v4r1_ctrl',
                     'SWDOWN':'EIG_dsw_plus_ECCO_v4r1_ctrl',
                     'UWIND':'EIG_u10m',
                     'VWIND':'EIG_v10m',
                     'PRECIP':'EIG_rain_plus_ECCO_v4r1_ctrl',
                     'RUNOFF':'runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin'}

    years = [2001]

    for vn in range(len(var_names)):
        var_name = var_names[vn]
        file_prefix = file_prefixes[var_name]
        for year in years:
            if var_name == 'RUNOFF':
                cexf.create_L05_exfs(config_dir,model_name,var_name,file_prefix,year,
                                     sNx,sNy,ordered_nonblank_tiles,ordered_nonblank_rotations,
                                     is_runoff = True)
            else:
                cexf.create_L05_exfs(config_dir, model_name, var_name, file_prefix, year,
                                     sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L05, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_exf_files(config_dir)
   

