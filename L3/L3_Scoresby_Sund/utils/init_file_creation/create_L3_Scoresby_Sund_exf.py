
import os
import argparse
import sys
import ast
import numpy as np

def create_exf_files(config_dir, L3_model_name, parent_model_level, parent_model_name, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L3','utils','init_file_creation'))

    start_year = 2002
    start_month = 1
    start_day = 1

    final_year = 2002
    final_month = 2
    final_day = 28

    ###################################################################################
    # The exf fields are created in 3 steps

    # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    import create_L3_exf_field_ref as efr
    efr.create_L3_exf_ref_file(config_dir, L3_model_name, parent_model_level, parent_model_name, print_level)


    # step 2: using the reference dict, organize downscaled exf into daily files

    if parent_model_level == 'L1':

        f = open(os.path.join(config_dir, 'L1', parent_model_name, 'namelist', parent_model_name + '_geometry.dict'))
        dict_str = f.read()
        f.close()
        size_dict = ast.literal_eval(dict_str)
        sNx = size_dict['sNx']
        sNy = size_dict['sNy']
        ordered_nonblank_tiles = size_dict['ordered_nonblank_tiles']
        ordered_nonblank_rotations = size_dict['ordered_nonblank_rotations']
        faces = size_dict['faces']
        ordered_tiles_faces_dict = size_dict['ordered_tiles_faces_dict']

        proc_ids = np.arange(8)

        import create_L3_daily_exf_fields_from_L1_ref as cef
        for proc_id in proc_ids:
            cef.create_exf_fields_via_interpolation(config_dir, L3_model_name, parent_model_name, proc_id,
                                                    sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                                                    faces, ordered_tiles_faces_dict,
                                                    start_year, final_year, start_month,
                                                    final_month, start_day, final_day, print_level)

        # step 3: combine all of the exf fields into a single file
        import combine_and_rotate_L3_daily_exf_files as com
        for proc_id in proc_ids:
            com.combine_and_rotate_L3_daily_exf_files(config_dir, L3_model_name, proc_id,
                                           start_year, final_year, start_month, final_month, start_day, final_day, print_level)



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
   

