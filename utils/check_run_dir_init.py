

import os
import argparse

def check_mandatory_files(run_dir):

    errors_detected = 0
    error_statements = []

    mandatory_files = ['data','data.pkg','eedata','mitgcmuv']

    for file_name in mandatory_files:
        if file_name in os.listdir(run_dir):
            print('        - Found '+file_name)
        else:
            errors_detected += 1
            error_statements.append('Missing file: '+file_name)
            print(' - Missing file: '+file_name)

    return(errors_detected, error_statements)


def check_package_files(run_dir):

    errors_detected = 0
    error_statements = []

    f = open(os.path.join(run_dir,'data.pkg'))
    lines = f.read()
    f.close()

    lines = lines.split('\n')

    pkgs = []
    for line in lines:
        if '.TRUE.' in line:
            line = line.split()
            if line[0]!='#':
                for part in line:
                    if "use" in part:
                        pkg=part.split("=")[0][3:].lower()
                        pkgs.append(pkg)

    for pkg in pkgs:
        if 'data.'+pkg in os.listdir(run_dir):
            print('        - '+pkg+' is activated and data.'+pkg+' is found')
        else:
            errors_detected += 1
            error_statements.append('        - ' + pkg + ' is activated and but data.' + pkg + ' is not found')
            print('        - ' + pkg + ' is activated but data.' + pkg + ' is NOT FOUND')

    return(errors_detected, error_statements, pkgs)


def check_bathymetry_file(run_dir):

    f = open(os.path.join(run_dir,'data'))
    lines = f.read()
    f.close()
    lines = lines.split('\n')

    errors_detected = 0
    error_statements = []

    for line in lines:
        if 'bathyFile' in line:
            line = line.split()
            if "=" in line[0]:
                bathy_file = line[0].split("=")[1][1:-2]
            else:
                bathy_file = line[2][1:-2]
            print('        - Bathymetry file identified in data: '+bathy_file)
            if bathy_file in os.listdir(run_dir):
                print('        - '+bathy_file+' file found')
            else:
                errors_detected += 1
                error_statements.append('        - '+bathy_file+' file not found')
                print('        - '+bathy_file+' file not found')

    return(errors_detected, error_statements)


def check_obcs_files(run_dir):

    f = open(os.path.join(run_dir, 'data.obcs'))
    lines = f.read()
    f.close()
    lines = lines.split('\n')

    errors_detected = 0
    error_statements = []

    file_names = []

    for line in lines:
        line=line.strip()
        if len(line)>0:
            if line[0]!='#':
                if 'file' in line.lower():
                    line = line.split()
                    if '=' in line[0]:
                        file_name = line[0].split('=')[2][1:-2]
                    else:
                        file_name = line[2][1:-2]
                    file_names.append(file_name)

    for file_name in file_names:
        if '/' in file_name:
            dir_name = file_name.split('/')[0]
            bc_name = file_name.split('/')[1]
            error_statements.append('        - obcs directory detected: '+dir_name)
            if dir_name not in os.listdir(run_dir):
                errors_detected += 1
                print('        - searching for '+bc_name+' in ' + dir_name + ' but directory not found')
            else:
                if file_name in os.listdir(os.path.join(run_dir,dir_name)):
                    print('        - '+file_name+' found')
                else:
                    errors_detected +=1
                    print('        - ' + file_name + ' not found')
        else:
            if file_name in os.listdir(run_dir):
                print('        - ' + file_name + ' found')
            else:
                errors_detected += 1
                print('        - ' + file_name + ' not found')

    return(errors_detected,error_statements)



def check_run_dir_init(run_dir):
    print('Checking the run directory for simple run-time issues')

    print('    - Checking mandatory files are present')
    errors_detected, error_statements = check_mandatory_files(run_dir)

    if errors_detected == 0:
        print('    - Checking package files are present')
        errors_detected, error_statements, pkgs = check_package_files(run_dir)

    if errors_detected == 0:
        print('    - Checking bathymetry file is present')
        errors_detected, error_statements = check_bathymetry_file(run_dir)

    if errors_detected == 0 and 'obcs' in pkgs:
        print('    - Checking obcs files')
        errors_detected, error_statements = check_obcs_files(run_dir)

    print('Number of errors found: '+str(errors_detected))

    # check exf files
    # if add mass, make sure it is compiled, select_fluid >= 1, add_mass is in pickup




# testing
run_dir = '/Users/michwood/Documents/Research/Projects/Ocean_Modeling/Projects/Downscale_LLC540_tides_Pacific/MITgcm/' \
          'configurations/pacific_llc540_tides_forcing/L1_1080/run_2'
check_run_dir_init(run_dir)

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#
#     parser.add_argument("-d", "--run_dir", action="store",
#                         help="The directory where the configuration will be run.", dest="run_dir",
#                         type=str, required=True)
#
#     args = parser.parse_args()
#     run_dir = args.run_dir
#
#     check_run_dir_init(run_dir)