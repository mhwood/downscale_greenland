

import argparse
import shutil
import os


def copy_files(config_dir, github_dir):

    if 'L3' not in os.listdir(github_dir):
        os.mkdir(os.path.join(github_dir,'L3'))

    L3_models = ['L3_Scoresby_Sund']
    subdirs = ['code','code_for_grid','namelist_for_grid']
    utils_subdirs = ['init_file_creation','plot_creation']

    for model_name in L3_models:
        if model_name not in os.listdir(os.path.join(github_dir,'L3')):
            os.mkdir(os.path.join(github_dir,'L3',model_name))

        #####################################
        # copy the code and namelist files
        for dir_name in subdirs:
            if dir_name in os.listdir(os.path.join(config_dir,'L3',model_name)):
                if dir_name not in os.listdir(os.path.join(github_dir,'L3',model_name)):
                    os.mkdir(os.path.join(github_dir,'L3',model_name,dir_name))
                for file_name in os.listdir(os.path.join(config_dir,'L3',model_name,dir_name)):
                    if file_name[0]!='.':
                        shutil.copyfile(os.path.join(config_dir,'L3',model_name,dir_name,file_name),
                                        os.path.join(github_dir,'L3',model_name,dir_name,file_name))

        #####################################
        # copy the code and namelist files
        for dir_name in utils_subdirs:
            if dir_name in os.listdir(os.path.join(config_dir, 'L3', model_name, 'utils')):
                if dir_name not in os.listdir(os.path.join(github_dir, 'L3', model_name,'utils')):
                    os.mkdir(os.path.join(github_dir, 'L3', model_name, 'utils', dir_name))
                for file_name in os.listdir(os.path.join(config_dir, 'L3', model_name, 'utils',dir_name)):
                    if file_name[0] != '.':
                        shutil.copyfile(os.path.join(config_dir, 'L3', model_name, 'utils',dir_name, file_name),
                                        os.path.join(github_dir, 'L3', model_name, 'utils',dir_name, file_name))
        for file_name in os.listdir(os.path.join(config_dir, 'L3', model_name, 'utils')):
            if file_name[-3:]=='.sh' or file_name[-3:]=='.py':
                shutil.copyfile(os.path.join(config_dir, 'L3', model_name, 'utils', file_name),
                                os.path.join(github_dir, 'L3', model_name, 'utils', file_name))

    #####################################
    # copy the L3 code and namelist files
    if 'utils' not in os.listdir(os.path.join(github_dir, 'L3')):
        os.mkdir(os.path.join(github_dir, 'L3', 'utils'))
    for dir_name in utils_subdirs:
        if dir_name in os.listdir(os.path.join(config_dir, 'L3', 'utils')):
            if dir_name not in os.listdir(os.path.join(github_dir, 'L3', 'utils')):
                os.mkdir(os.path.join(github_dir, 'L3', 'utils', dir_name))
            for file_name in os.listdir(os.path.join(config_dir, 'L3', 'utils', dir_name)):
                if file_name[0] != '.' and file_name[0] != '_':
                    shutil.copyfile(
                        os.path.join(config_dir, 'L3', 'utils', dir_name, file_name),
                        os.path.join(github_dir, 'L3', 'utils', dir_name, file_name))
    for file_name in os.listdir(os.path.join(config_dir, 'L3', 'utils')):
        if file_name[-3:] == '.sh' or file_name[-3:] == '.py':
            shutil.copyfile(os.path.join(config_dir, 'L3', 'utils', file_name),
                            os.path.join(github_dir, 'L3', 'utils', file_name))

    #####################################
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
        if file_name[-3:] == '.sh' or file_name[-3:] == '.py':
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