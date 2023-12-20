
import os
import argparse
import sys
import ast
import numpy as np

def create_exf_files(config_dir, L3_model_name, parent_model_level, parent_model_name, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L3','utils','init_file_creation'))

    start_year = 2005
    final_year = 2009

    ###################################################################################
    # The exf fields are created in 3 steps

    if parent_model_level == 'L1_grid':

        ############################################################################################
        # load all of the L1 model stuff

        ############################################################################################
        # interpolate from annual L1 input files as a grid
        proc_ids = np.arange(7)

        import create_L3_daily_exf_files_from_L1_input_grid as cef
        for proc_id in proc_ids:
            cef.create_exf_fields_via_interpolation(config_dir, L3_model_name, parent_model_name, proc_id,
                                                    start_year, final_year, print_level)

        # combine all of the exf fields into a single file
        import combine_and_rotate_L3_daily_exf_files as com
        for proc_id in proc_ids:
            com.combine_and_rotate_L3_daily_exf_files(config_dir, L3_model_name, proc_id, start_year, final_year, print_level)


        # ############################################################################################
        # # these codes just interpolate from annual L1 input files as faces
        # proc_ids = np.arange(4,8)
        #
        # f = open(os.path.join(config_dir, 'L1', parent_model_name, 'namelist', parent_model_name + '_geometry.dict'))
        # dict_str = f.read()
        # f.close()
        # size_dict = ast.literal_eval(dict_str)
        # sNx = size_dict['sNx']
        # sNy = size_dict['sNy']
        # ordered_nonblank_tiles = size_dict['ordered_nonblank_tiles']
        # ordered_nonblank_rotations = size_dict['ordered_nonblank_rotations']
        # faces = size_dict['faces']
        # ordered_tiles_faces_dict = size_dict['ordered_tiles_faces_dict']
        #
        # face_size_dict = {}
        # for key in list(ordered_tiles_faces_dict.keys()):
        #     rows = len(ordered_tiles_faces_dict[key])
        #     cols = len(ordered_tiles_faces_dict[key][0])
        #     face_size_dict[key] = (rows * sNy, cols * sNx)
        #
        # sys.path.insert(1, os.path.join(config_dir, 'L1', parent_model_name, 'utils'))
        # import L1_CE_Greenland_functions as L1f
        #
        # import create_L3_daily_exf_files_from_L1_input_faces as cef
        # for proc_id in proc_ids:
        #     cef.create_exf_fields_via_interpolation(L1f, config_dir, L3_model_name, parent_model_name, proc_id,
        #                                             sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
        #                                             faces, face_size_dict, ordered_tiles_faces_dict,
        #                                             start_year, final_year, print_level)


        ############################################################################################
        # these codes use the dv output from L1

        # # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
        # import create_L3_exf_field_ref as efr
        # efr.create_L3_exf_ref_file(config_dir, L3_model_name, parent_model_level, parent_model_name, print_level)
        #
        # # step 2: using the reference dict, organize downscaled exf into daily files
        #
        # proc_ids = np.arange(8)
        #
        # import create_L3_daily_exf_fields_from_L1_ref as cef
        # for proc_id in proc_ids:
        #     cef.create_exf_fields_via_interpolation(config_dir, L3_model_name, parent_model_name, proc_id,
        #                                             sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
        #                                             faces, ordered_tiles_faces_dict,
        #                                             start_year, final_year, start_month,
        #                                             final_month, start_day, final_day, print_level)
        #
        # # step 3: combine all of the exf fields into a single file
        # import combine_and_rotate_L3_daily_exf_files as com
        # for proc_id in proc_ids:
        #     com.combine_and_rotate_L3_daily_exf_files(config_dir, L3_model_name, proc_id,
        #                                    start_year, final_year, start_month, final_month, start_day, final_day, print_level)



    if parent_model_level == 'L2':
        import create_L3_daily_exf_fields_from_ref as cef
        for proc_id in range(1):
            cef.create_exf_fields_via_interpolation(config_dir, L3_model_name, parent_model, proc_id,
                                                    start_year, final_year,
                                                    start_month, final_month,
                                                    start_day, final_day, print_level)

        # step 3: combine all of the exf fields into a single file
        import combine_L3_daily_exf_files as com
        for proc_id in range(1):
            com.combine_L3_daily_exf_files(config_dir, L3_model_name, proc_id,
                                           start_year, final_year, start_month, final_month, start_day, final_day, print_level)







if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

