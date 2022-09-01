
import os
import argparse
import sys

def create_exf_files(config_dir, L1_model_name, ecco_dir, llc, ordered_ecco_tiles, ordered_ecco_tile_rotations,
                     sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L1','utils','init_file_creation'))
    import create_L1_ECCO_exfs as cexf

    if 'exf' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'input','exf'))

    sys.path.insert(1, os.path.join(config_dir, 'L1', L1_model_name, 'utils'))
    import L1_W_Greenland_functions as Lf

    years = [1992]

    var_names = ['ATEMP', 'AQH', 'LWDOWN', 'SWDOWN', 'UWIND', 'VWIND', 'PRECIP','RUNOFF']
    file_prefixes = {'ATEMP': 'EIG_tmp2m_degC_plus_ECCO_v4r1_ctrl',
                     'AQH': 'EIG_spfh2m_plus_ECCO_v4r1_ctrl',
                     'LWDOWN': 'EIG_dlw_plus_ECCO_v4r1_ctrl',
                     'SWDOWN': 'EIG_dsw_plus_ECCO_v4r1_ctrl',
                     'UWIND': 'EIG_u10m',
                     'VWIND': 'EIG_v10m',
                     'PRECIP': 'EIG_rain_plus_ECCO_v4r1_ctrl',
                     'RUNOFF': 'runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin'}

    for vn in range(len(var_names)):
        var_name = var_names[vn]
        file_prefix = file_prefixes[var_name]
        if var_name == 'RUNOFF':
            year = -12 # means nothing for runoff, just a filler
            cexf.create_L1_exfs(Lf, config_dir, L1_model_name, var_name, ecco_dir, llc, ordered_ecco_tiles,
                                ordered_ecco_tile_rotations, file_prefix, year,
                                sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations, print_level,
                                is_runoff=True)
        else:
            for year in years:
                cexf.create_L1_exfs(Lf, config_dir, L1_model_name, var_name, ecco_dir, llc, ordered_ecco_tiles, ordered_ecco_tile_rotations, file_prefix, year,
                                     sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations, print_level)







if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

