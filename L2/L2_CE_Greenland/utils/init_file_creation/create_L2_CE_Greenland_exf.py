
import os
import argparse
import sys

def create_exf_files(config_dir, L2_model_name, parent_model,
                     sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
                     faces, ordered_tiles_faces_dict, print_level):

    sys.path.insert(1, os.path.join(config_dir,'L2', 'utils','init_file_creation'))

    start_year = 2002
    start_month = 1
    start_day = 1

    final_year = 2002
    final_month = 1
    final_day = 3

    ###################################################################################
    # The exf fields are created in 3 steps

    # # step 1: make a reference whereby the diagnostics_vec files are organized in a dictionary
    # import create_L2_exf_field_ref as efr
    # efr.create_L2_exf_ref_file(config_dir, L2_model_name, parent_model, print_level)

    # step 2: using the reference dict, organize downscaled exf into daily files
    # import create_L2_daily_exf_from_ref as cef
    # for proc_id in range(4,6):
    #     cef.create_exf_fields_via_interpolation(config_dir, L2_model_name, parent_model, proc_id,
    #                                             sNx, sNy, ordered_nonblank_tiles, ordered_nonblank_rotations,
    #                                             faces, ordered_tiles_faces_dict,
    #                                             start_year, final_year,
    #                                             start_month, final_month,
    #                                             start_day, final_day, print_level)

    # step 3: combine all of the exf fields into a single file
    import combine_and_rotate_L2_daily_exf_files as com
    for proc_id in range(4,6):
        com.combine_and_rotate_L2_daily_exf_files(config_dir, L2_model_name, proc_id,
                                       start_year, final_year, start_month, final_month, start_day, final_day,
                                       print_level)







if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L1 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_pickup_file(config_dir)
   

