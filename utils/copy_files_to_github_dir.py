

import argparse
import shutil
import os


def copy_files(config_dir, github_dir):

    level_names = ['L1','L2','L3','L05','L0']

    for level_name in level_names:
        if level_name not in os.listdir(github_dir):
            os.mkdir(os.path.join(github_dir,level_name))

    level_name_model_dict = {'L1':['L1_CE_Greenland'],
                             'L2':['L2_CE_Greenland'],
                             'L3':['L3_Scoresby_Sund'],
                             'L05':['L05_CE_Greenland'],
                             'L0':[]}
    subdirs = ['code','code_for_grid','namelist','namelist_for_grid']
    utils_subdirs = ['init_file_creation','plot_creation']

    for level_name in level_names:

        #################################################################
        # copy the model files first
        model_names = level_name_model_dict[level_name]

        for model_name in model_names:
            if model_name not in os.listdir(os.path.join(github_dir,level_name)):
                os.mkdir(os.path.join(github_dir,level_name,model_name))

            #####################################
            # copy the code and namelist files
            for dir_name in subdirs:
                if dir_name in os.listdir(os.path.join(config_dir,level_name,model_name)):
                    if dir_name not in os.listdir(os.path.join(github_dir,level_name,model_name)):
                        os.mkdir(os.path.join(github_dir,level_name,model_name,dir_name))
                    for file_name in os.listdir(os.path.join(config_dir,level_name,model_name,dir_name)):
                        if file_name[0]!='.' and file_name[-3:]!='.nc':
                            shutil.copyfile(os.path.join(config_dir,level_name,model_name,dir_name,file_name),
                                            os.path.join(github_dir,level_name,model_name,dir_name,file_name))

            #####################################
            # copy the utils files
            if 'utils' not in os.listdir(os.path.join(github_dir, level_name, model_name)):
                os.mkdir(os.path.join(github_dir, level_name, model_name, 'utils'))
            for dir_name in utils_subdirs:
                if dir_name in os.listdir(os.path.join(config_dir, level_name, model_name, 'utils')):
                    if dir_name not in os.listdir(os.path.join(github_dir, level_name, model_name,'utils')):
                        os.mkdir(os.path.join(github_dir, level_name, model_name, 'utils', dir_name))
                    for file_name in os.listdir(os.path.join(config_dir, level_name, model_name, 'utils',dir_name)):
                        if file_name[0] != '.' and file_name[-3:]=='.py' and file_name[0] != '_':
                            shutil.copyfile(os.path.join(config_dir, level_name, model_name, 'utils',dir_name, file_name),
                                            os.path.join(github_dir, level_name, model_name, 'utils',dir_name, file_name))
            for file_name in os.listdir(os.path.join(config_dir, level_name, model_name, 'utils')):
                if file_name[-3:]=='.sh' or file_name[-3:]=='.py':
                    shutil.copyfile(os.path.join(config_dir, level_name, model_name, 'utils', file_name),
                                    os.path.join(github_dir, level_name, model_name, 'utils', file_name))

        #####################################
        # copy the general level utils files
        if 'utils' in os.listdir(os.path.join(config_dir, level_name)):
            if 'utils' not in os.listdir(os.path.join(github_dir, level_name)):
                os.mkdir(os.path.join(github_dir, level_name, 'utils'))
            for dir_name in utils_subdirs:
                if dir_name in os.listdir(os.path.join(config_dir, level_name, 'utils')):
                    if dir_name not in os.listdir(os.path.join(github_dir, level_name, 'utils')):
                        os.mkdir(os.path.join(github_dir, level_name, 'utils', dir_name))
                    for file_name in os.listdir(os.path.join(config_dir, level_name, 'utils', dir_name)):
                        if file_name[0] != '.' and file_name[0] != '_':
                            shutil.copyfile(
                                os.path.join(config_dir, level_name, 'utils', dir_name, file_name),
                                os.path.join(github_dir, level_name, 'utils', dir_name, file_name))
            for file_name in os.listdir(os.path.join(config_dir, level_name, 'utils')):
                if file_name[-3:] == '.sh' or file_name[-3:] == '.py':
                    shutil.copyfile(os.path.join(config_dir, level_name, 'utils', file_name),
                                    os.path.join(github_dir, level_name, 'utils', file_name))

    ####################################################################################################################
    # copy the general code and namelist files
    if 'utils' not in os.listdir(github_dir):
        os.mkdir(os.path.join(github_dir, 'utils'))
    for dir_name in utils_subdirs:
        if dir_name in os.listdir(os.path.join(config_dir, 'utils')):
            if dir_name not in os.listdir(os.path.join(github_dir, 'utils')):
                os.mkdir(os.path.join(github_dir, 'utils', dir_name))
            for file_name in os.listdir(os.path.join(config_dir, 'utils', dir_name)):
                if file_name[0] != '.' and file_name[0] != '_':
                    shutil.copyfile(
                        os.path.join(config_dir,  'utils', dir_name, file_name),
                        os.path.join(github_dir,  'utils', dir_name, file_name))
    for file_name in os.listdir(os.path.join(config_dir,  'utils')):
        if file_name[-3:] == '.sh' or file_name[-3:] == '.py' and file_name[0] != '_':
            shutil.copyfile(os.path.join(config_dir, 'utils', file_name),
                            os.path.join(github_dir, 'utils', file_name))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-g", "--github_dir", action="store",
                        help="The directory which will be pushed to Github.", dest="github_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir
    github_dir = args.github_dir

    copy_files(config_dir, github_dir)