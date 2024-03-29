
import os
import argparse
import numpy as np
import sys

def create_exf_files(config_dir, L1_model_name, ecco_dir, mankoff_dir, llc, ordered_ecco_tiles, ordered_ecco_tile_rotations,
                     sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations, tile_face_index_dict, face_size_dict, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L1','utils','init_file_creation'))

    if 'exf' not in os.listdir(os.path.join(config_dir,'L1',L1_model_name,'input')):
        os.mkdir(os.path.join(config_dir,'L1',L1_model_name,'input','exf'))

    sys.path.insert(1, os.path.join(config_dir, 'L1', L1_model_name, 'utils'))
    import L1_CE_Greenland_functions as Lf

    start_year = 2017
    start_month = 1
    start_day = 1

    final_year = 2018
    final_month = 12
    final_day = 31

    years = np.arange(start_year,final_year+1).tolist()

    ###################################################################################################################
    # these are the fields output from diagnostics vec (after ECCO corrections have been applied)

    # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    import create_L1_exf_field_ref as cer
    cer.create_L1_exf_ref_file(config_dir, L1_model_name, print_level)

    # step 2: make monthly stacks of the exf files
    var_names = ['ATEMP', 'AQH', 'LWDOWN', 'SWDOWN', 'UWIND', 'VWIND', 'PRECIP','RUNOFF']

    import create_L1_ECCO_exfs_from_ref as cexf
    for var_name in var_names:
        cexf.create_L1_exfs(Lf, config_dir, L1_model_name, var_name, ecco_dir, llc, ordered_ecco_tiles,
                            ordered_ecco_tile_rotations,
                            start_year, final_year, start_month, final_month, start_day, final_day,
                            sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations, tile_face_index_dict, face_size_dict, print_level)

    import combine_and_rotate_L1_daily_exf_files_annual as com
    for var_name in var_names:
        com.combine_L1_daily_exf_files_annual(Lf, config_dir, L1_model_name, var_name,
                                          sNx, sNy, ordered_nonblank_tiles,ordered_nonblank_rotations,
                                          start_year, final_year,
                                          print_level)

    # ####################################################################################################################
    # # these are the files for the iceplume files
    # import create_L1_iceplume_files as cipf
    #
    # cipf.create_L1_iceplume_files(Lf, config_dir, L1_model_name, mankoff_dir,
    #                               years, sNx, sNy, ordered_nonblank_tiles, tile_face_index_dict, face_size_dict, print_level)

    # ####################################################################################################################
    # # these are the files for the subglacial discharge files
    # import create_L1_sgd as csgd
    #
    # csgd.create_L1_sgd_files(Lf, config_dir, L1_model_name, mankoff_dir,
    #                          years, sNx, sNy, ordered_nonblank_tiles, tile_face_index_dict, face_size_dict, print_level)

    # ####################################################################################################################
    # # these are the "normal" exf fields
    # import create_L1_ECCO_exfs_unadjusted as cexf
    #
    # var_names = ['ATEMP', 'AQH', 'LWDOWN', 'SWDOWN', 'UWIND', 'VWIND', 'PRECIP','RUNOFF']
    # file_prefixes = {'ATEMP': 'EIG_tmp2m_degC_plus_ECCO_v4r1_ctrl',
    #                  'AQH': 'EIG_spfh2m_plus_ECCO_v4r1_ctrl',
    #                  'LWDOWN': 'EIG_dlw_plus_ECCO_v4r1_ctrl',
    #                  'SWDOWN': 'EIG_dsw_plus_ECCO_v4r1_ctrl',
    #                  'UWIND': 'EIG_u10m',
    #                  'VWIND': 'EIG_v10m',
    #                  'PRECIP': 'EIG_rain_plus_ECCO_v4r1_ctrl',
    #                  'RUNOFF': 'runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin'}
    #
    # for vn in range(len(var_names)):
    #     var_name = var_names[vn]
    #     file_prefix = file_prefixes[var_name]
    #     if var_name == 'RUNOFF':
    #         year = -12 # means nothing for runoff, just a filler
    #         cexf.create_L1_exfs(Lf, config_dir, L1_model_name, var_name, ecco_dir, llc, ordered_ecco_tiles,
    #                             ordered_ecco_tile_rotations, file_prefix, year,
    #                             sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations, print_level,
    #                             is_runoff=True)
    #     else:
    #         for year in years:
    #             cexf.create_L1_exfs(Lf, config_dir, L1_model_name, var_name, ecco_dir, llc, ordered_ecco_tiles,
    #                                 ordered_ecco_tile_rotations, file_prefix, year,
    #                                 sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations, print_level)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

