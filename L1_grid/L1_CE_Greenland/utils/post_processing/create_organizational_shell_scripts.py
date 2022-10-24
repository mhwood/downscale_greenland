
import os
import numpy as np
import argparse

def create_diag_tar_scripts(config_dir, var_names, prefixes):

    tar_dir = config_dir+'/L1_grid/L1_CE_Greenland/results/tar_diag_scripts'

    for prefix in prefixes:

        output_file = 'tar_L1_diag_files_'+prefix+'00000.sh'
        output = ''

        for var_name in var_names:
            command = 'echo "tarring '+var_name+' files for set '+prefix+'"'
            output += command + '\n'

            command = 'cd ../../run/diags/'+var_name
            # print(command)
            output += command + '\n'

            command = 'tar -czvf ../../../results/diags/'+var_name+'/'+var_name+'.'+prefix+'00000.tar.gz '+var_name+'.'+prefix+'*'
            # print(command)
            output += command + '\n'

            command = 'cd ../../../results/tar_diag_scripts/'
            # print(command)
            output += command + '\n'

        f = open(os.path.join(tar_dir, output_file), 'w')
        f.write(output)
        f.close()

def create_dv_tar_scripts(config_dir, mask_names, prefixes,dv_var_names, dv_surf_names):

    tar_dir = config_dir+'/L1_grid/L1_CE_Greenland/results/tar_dv_scripts'

    for prefix in prefixes:

        output_file = 'tar_L1_dv_files_'+prefix+'00000.sh'
        output = ''

        for mask_name in mask_names:
            command = 'echo "tarring '+mask_name+' files for set '+prefix+'"'
            output += command + '\n'

            command = 'cd ../../run/dv'
            # print(command)
            output += command + '\n'

            if 'surface' in mask_name:
                for var_name in dv_surf_names:
                    command = 'tar -czvf ../../results/dv/'+mask_name+'/'+mask_name+'_'+var_name+'.'+prefix+'00000.tar.gz '+mask_name+'*'+var_name+'.'+prefix+'*'
                    # print(command)
                    output += command + '\n'
            else:
                for var_name in dv_var_names:
                    command = 'tar -czvf ../../results/dv/' + mask_name + '/' + mask_name+'_'+var_name + '.' + prefix + '00000.tar.gz ' + mask_name+'*'+var_name + '.' + prefix + '*'
                    # print(command)
                    output += command + '\n'

            command = 'cd ../../results/tar_dv_scripts/'
            # print(command)
            output += command + '\n'

        f = open(os.path.join(tar_dir, output_file), 'w')
        f.write(output)
        f.close()

def create_lou_mv_diag_scripts(config_dir, var_names, prefixes):

    lou_dir = config_dir+'/L1_grid/L1_CE_Greenland/results/lou_mv_diag_scripts'

    path_name = '/nobackup/mwood7/Greenland/Nested_Models/MITgcm/configurations/downscale_greenland/L1_grid/L1_CE_Greenland/results/diags/'

    for prefix in prefixes:
        output_file = 'mv_diag_files_from_nobackupp_'+prefix+'.sh'
        output = ''
        for var_name in var_names:
            command = 'mv '+path_name+'/'+var_name+'/'+var_name+'.'+prefix+'00000.tar.gz ../diags/'+var_name+'/'
            output+=command+'\n'

        f = open(os.path.join(lou_dir,output_file),'w')
        f.write(output)
        f.close()

def create_lou_mv_dv_scripts(config_dir, mask_names, prefixes, dv_var_names, dv_surf_names):

    lou_dir = config_dir+'/L1_grid/L1_CE_Greenland/results/lou_mv_dv_scripts'

    path_name = '/nobackup/mwood7/Greenland/Nested_Models/MITgcm/configurations/downscale_greenland/L1_grid/L1_CE_Greenland/results/dv/'

    for prefix in prefixes:
        output_file = 'mv_dv_files_from_nobackupp_'+prefix+'.sh'
        output = ''
        for mask_name in mask_names:
            if 'surface' in mask_name:
                for var_name in dv_surf_names:
                    command = 'mv '+path_name+'/'+mask_name+'/'+mask_name+'_'+var_name + '.' + prefix + '00000.tar.gz ../dv/'+mask_name
                    output += command + '\n'
            else:
                for var_name in dv_var_names:
                    command = 'mv ' + path_name + '/' + mask_name + '/' + mask_name+'_'+var_name + '.' + prefix + '00000.tar.gz ../dv/' + mask_name
                    output+=command+'\n'

        f = open(os.path.join(lou_dir,output_file),'w')
        f.write(output)
        f.close()

def create_diag_scp_scripts(config_dir, var_names, prefixes):

    diag_dir = config_dir+'/L1_grid/L1_CE_Greenland/run/diags/scp_scripts'

    path_name = 'mwood7@pfe22.nas.nasa.gov:/home4/mwood7/nobackupp12/mwood7/Greenland/Nested_Models/MITgcm/' \
                'configurations/downscale_greenland/L1_grid/L1_CE_Greenland/results/diags'

    for prefix in prefixes:
        output_file = 'scp_L1_CE_'+prefix+'00000.sh'
        output = ''
        for var_name in var_names:
            command = 'sup scp  '+path_name+'/'+var_name+'/'+var_name+'.'+prefix+'00000.tar.gz ../'
            output+=command+'\n'

        f = open(os.path.join(diag_dir,output_file),'w')
        f.write(output)
        f.close()

def create_diag_untar_scripts(config_dir, var_names, prefixes):

    diag_dir = config_dir + '/L1_grid/L1_CE_Greenland/run/diags/untar_scripts'

    for prefix in prefixes:
        output_file = 'untar_L1_CE_' + prefix + '00000.sh'
        output = ''

        command = 'echo "untarring files for set ' + prefix + '"'
        output += command + '\n'

        for var_name in var_names:
            command = 'echo "untarring ' + var_name + '"'
            output += command + '\n'

            command = 'cd ..'
            output += command + '\n'

            command = 'tar -xf ' + var_name + '.' + prefix + '00000.tar.gz -C '+var_name
            output += command + '\n'

            command = 'rm ' + var_name + '.' + prefix + '00000.tar.gz'
            output += command + '\n'

            command = 'cd untar_scripts'
            output += command + '\n'

        f = open(os.path.join(diag_dir, output_file), 'w')
        f.write(output)
        f.close()

def create_dv_scp_scripts(config_dir, mask_names, prefixes, dv_var_names, dv_surf_names):

    diag_dir = config_dir+'/L1_grid/L1_CE_Greenland/run/dv/scp_dv_scripts'

    path_name = 'mwood7@lfe6.nas.nasa.gov:/u/mwood7/Greenland/L1_grid/L1_CE_Greenland/results/dv'

    for prefix in prefixes:
        output_file = 'scp_L1_CE_'+prefix+'00000.sh'
        output = ''
        for mask_name in mask_names:
            command = 'sup scp  '+path_name+'/'+mask_name+'/'+mask_name+'*.'+prefix+'00000.tar.gz ../'
            output+=command+'\n'

        f = open(os.path.join(diag_dir,output_file),'w')
        f.write(output)
        f.close()

def create_dv_untar_scripts(config_dir, mask_names, prefixes, dv_var_names, dv_surf_names):

    diag_dir = config_dir + '/L1_grid/L1_CE_Greenland/run/dv/untar_dv_scripts'

    for prefix in prefixes:
        output_file = 'untar_L1_CE_' + prefix + '00000.sh'
        output = ''

        command = 'echo "untarring dv files for set ' + prefix + '"'
        output += command + '\n'

        for mask_name in mask_names:

            command = 'echo "untarring ' + mask_name + ' files"'
            output += command + '\n'

            command = 'cd ..'
            output += command + '\n'

            if 'surface' in mask_name:
                for var_name in dv_surf_names:
                    command = 'echo "untarring ' + var_name + '"'
                    output += command + '\n'
                    command = 'tar -xf ' + mask_name + '_' + var_name + '.' + prefix + '00000.tar.gz -C ' + mask_name
                    output += command + '\n'
                    command = 'rm ' + mask_name + '_' + var_name + '.' + prefix + '00000.tar.gz'
                    output += command + '\n'
            else:
                for var_name in dv_var_names:
                    command = 'echo "untarring ' + var_name + '"'
                    output += command + '\n'
                    command = 'tar -xf ' + mask_name + '_' + var_name + '.' + prefix + '00000.tar.gz -C ' + mask_name
                    output += command + '\n'
                    command = 'rm ' + mask_name + '_' + var_name + '.' + prefix + '00000.tar.gz'
                    output += command + '\n'

            command = 'cd untar_dv_scripts'
            output += command + '\n'

        f = open(os.path.join(diag_dir, output_file), 'w')
        f.write(output)
        f.close()


def create_shell_scripts(config_dir):

    var_names = ['EtaN_day_snap', 'EtaN_mon_mean', 'SI_daily_snap', 'TS_AW_daily_snap', 'TS_surf_daily_snap',
                 'state_3D_mon_mean', 'state_3D_mon_snap', 'vel_3D_mon_mean', 'vel_3D_mon_snap']

    dv_var_names = ['AREA','HEFF','HSNOW','SALT','THETA','UICE','UVEL','VICE','VVEL']
    dv_surf_names = ['ATEMP','AQH','LWDOWN','SWDOWN','RUNOFF','PRECIP','UWIND','VWIND']

    mask_names = ['CTD', 'L3_east', 'L3_north', 'L3_south', 'L3_surface']

    prefixes = []
    for index in range(32):
        prefixes.append('{:05d}'.format(index))
    # prefixes = ['00000']

    print('Creating shell scripts to ...')

    print('    ... tar diag files on Pleiades')
    create_diag_tar_scripts(config_dir, var_names, prefixes)

    print('    ... tar dv files on Pleiades')
    create_dv_tar_scripts(config_dir, mask_names, prefixes, dv_var_names, dv_surf_names)

    print('    ... mv diag files from Pleiades to Lou')
    create_lou_mv_diag_scripts(config_dir, var_names, prefixes)

    print('    ... mv dv files from Pleiades to Lou')
    create_lou_mv_dv_scripts(config_dir, mask_names, prefixes, dv_var_names, dv_surf_names)

    print('    ... scp diag files from Pleiades to local drive')
    create_diag_scp_scripts(config_dir, var_names, prefixes)

    print('    ... untar diag files on local drive')
    create_diag_untar_scripts(config_dir, var_names, prefixes)

    print('    ... scp dv files from Lou to local drive')
    create_dv_scp_scripts(config_dir, mask_names, prefixes, dv_var_names, dv_surf_names)

    print('    ... untar diag files on local drive')
    create_dv_untar_scripts(config_dir, mask_names, prefixes, dv_var_names, dv_surf_names)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_shell_scripts(config_dir)