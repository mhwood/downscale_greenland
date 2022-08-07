
import os
from datetime import datetime
import numpy as np
import netCDF4 as nc4
import argparse
import sys


def find_dv_files_to_read(config_dir, model_name, averaging_period, seconds_per_iter,
                          var_name, start_iter, final_iter):

    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import time_functions as tf

    iters_per_output = averaging_period/seconds_per_iter

    nc_dict_file = os.path.join(config_dir, 'L05', model_name, 'namelist', 'L05_dv_mask_reference_dict.nc')
    ds = nc4.Dataset(nc_dict_file)
    grp = ds.groups['surface']
    faces = grp.variables['source_faces'][:]
    N = len(faces)
    ds.close()

    # loop through the dv files and see if they have any of the data we want, keep a list
    file_names = []
    file_iters = []
    iters_per_files = []
    # dv_dir = os.path.join(config_dir, 'L0_540', 'run', 'dv', mask_name + '_i0', var_name)
    dv_dir = os.path.join(config_dir, 'L05', model_name, 'run', 'dv')
    for file_name in os.listdir(dv_dir):
        if var_name in file_name:
            file_iter = int(file_name.split('.')[-2])
            iters_per_file = int(np.size(np.fromfile(os.path.join(dv_dir,file_name),'>f4'))/(N))
            iter_midpoints = tf.dv_file_name_iter_to_iter_midpoints(file_iter, iters_per_output, iters_per_file)
            if not np.min(iter_midpoints)>final_iter+iters_per_output and not np.max(iter_midpoints)<start_iter-iters_per_output:
                file_iters.append(file_iter)
                file_names.append(file_name)
                iters_per_files.append(iters_per_file)

    sorted_indices = sorted(range(len(file_iters)), key=lambda k: file_iters[k])
    file_names_sorted = []
    file_iters_sorted = []
    iters_per_files_sorted = []
    for index in sorted_indices:
        file_names_sorted.append(file_names[index])
        file_iters_sorted.append(file_iters[index])
        iters_per_files_sorted.append(iters_per_files[index])

    # make dicts of the index and corresponding iters which will be read from these files
    iter_midpoint_dict = {}
    for i in range(len(file_iters_sorted)):
        file_iter = file_iters_sorted[i]
        file_name = file_names_sorted[i]
        iters_per_file = iters_per_files_sorted[i]
        iter_midpoints = tf.dv_file_name_iter_to_iter_midpoints(file_iter, iters_per_output, iters_per_file)
        iter_midpoint_dict[file_name] = iter_midpoints
        if file_name == 'L2_surface_mask_ATEMP.0000006481.bin':
            iter_midpoints+=1

    return(file_names, iter_midpoint_dict)

def create_destination_file_list(config_dir, averaging_period, seconds_per_iter,
                                 var_name, start_iter, final_iter, file_names, iter_midpoint_dict):
    sys.path.insert(1, os.path.join(config_dir, 'utils','init_file_creation'))
    import time_functions as tf

    start_seconds = 0
    iters_per_output = averaging_period / seconds_per_iter

    # create a list of monthly bounds (in seconds)
    date_bounds = []
    for year in range(2002, 2003):
        for month in range(1, 13):
            if month in [1, 3, 5, 7, 8, 10, 12]:
                nDays = 31
            elif month in [4, 6, 9, 11]:
                nDays = 30
            else:
                if year % 4 == 0:
                    nDays = 29
                else:
                    nDays = 28
            for day in range(1, nDays):
                date_bounds.append([datetime(year, month, day), datetime(year, month, day + 1)])
            if month != 12:
                date_bounds.append([datetime(year, month, day + 1), datetime(year, month + 1, 1)])
        date_bounds.append([datetime(year, 12, 31), datetime(year + 1, 1, 1)])

    dest_files = []
    dest_file_iter_bounds = {}
    # convert these to iteration numbers
    for date_set in date_bounds:
        date_0 = (date_set[0] - datetime(2002, 1, 1)).total_seconds()
        date_1 = (date_set[1] - datetime(2002, 1, 1)).total_seconds()
        iter_0 = tf.elapsed_seconds_to_iters(start_seconds, seconds_per_iter, date_0)
        iter_1 = tf.elapsed_seconds_to_iters(start_seconds, seconds_per_iter, date_1)
        if not iter_0>final_iter+1 and not iter_1<start_iter-1:
            dest_file = 'L2_exf_'+str(var_name)+'.'+str(date_set[0].year)+'{:02d}'.format(date_set[0].month)+'{:02d}'.format(date_set[0].day)+'.bin'
            dest_files.append(dest_file)
            dest_file_iter_bounds[dest_file] = [iter_0,iter_1]

    source_file_read_dict = {}
    source_file_read_index_sets = {}

    # # make dicts of the index and corresponding iters which will be read from these files
    for dest_file in dest_files:
        print('Creating reference dictionary for '+dest_file)

        source_file_read_dict[dest_file] = []
        source_file_read_index_sets[dest_file] = []

        dest_file_iter_0 = dest_file_iter_bounds[dest_file][0]
        dest_file_iter_1 = dest_file_iter_bounds[dest_file][1]
        if dest_file_iter_0<start_iter:
            dest_file_iter_0=start_iter
        if dest_file_iter_1 > final_iter:
            dest_file_iter_1 = final_iter

        dest_file_iters = np.arange(dest_file_iter_0 + 1,dest_file_iter_1+1,iters_per_output)
        dest_file_iters += iters_per_output / 2
        print('      - Iteration range: '+str(np.min(dest_file_iters))+' to '+str(np.max(dest_file_iters))+' ('+str(np.max(dest_file_iters)-np.min(dest_file_iters))+')')

        source_file_names = []

        # get a list of files to use for each iteration
        for dest_iter in dest_file_iters:
            file_to_use = ''
            for file_name in file_names:
                iter_midpoints = iter_midpoint_dict[file_name]
                if dest_iter >= np.min(iter_midpoints) and dest_iter <= np.max(iter_midpoints):
                    if file_to_use == '':
                        file_to_use = file_name
                    else:
                        if int(file_to_use.split('.')[-2]) > int(file_name.split('.')[-2]):
                            file_to_use = file_name
            if file_to_use!='':
                source_file_names.append(file_to_use)

        if len(source_file_names)>0:

            # use the list of source file names to make a dict that lists which iters to read from it
            unique_list = list(set(source_file_names))
            print('      - Reading from file(s) '+', '.join(unique_list))
            for file_name in unique_list:
                source_file_read_dict[dest_file].append(file_name)
            source_file_read_dict[dest_file] = sorted(source_file_read_dict[dest_file])

            # loop through the source file names to identify which indices will be read from each one
            for file_name in source_file_read_dict[dest_file]:
                min_dest_iter_index = source_file_names.index(file_name)
                max_dest_iter_index = len(source_file_names) - 1 - source_file_names[::-1].index(file_name)
                min_dest_iter = dest_file_iters[min_dest_iter_index]
                max_dest_iter = dest_file_iters[max_dest_iter_index]
                # print('        - The file ' + file_name+ ' will cover iters '+str(min_dest_iter)+' to '+str(max_dest_iter))
                source_file_iters = iter_midpoint_dict[file_name]
                # print('            - This file has iters '+str(np.min(source_file_iters))+' to '+str(np.max(source_file_iters)))
                try:
                    index_0 = list(source_file_iters).index(min_dest_iter)
                    index_1 = list(source_file_iters).index(max_dest_iter)+1
                    if file_name == source_file_read_dict[dest_file][-1]:
                        index_1 -= 1
                    # print('            - This corresponds to indices '+str(index_0)+' ('+str(source_file_iters[index_0]) +') through '+str(index_1-1)+' ('+str(source_file_iters[index_1-1]) +')')
                    source_file_read_index_sets[dest_file].append([index_0, index_1])
                except:
                    a=1
                    # print('            - Didnt find the correct indices in this file')


    return(dest_files, source_file_read_dict, source_file_read_index_sets)

########################################################################################################################


def create_exf_fields_reference_dict(config_dir, L05_model_name, averaging_period, seconds_per_iter):

    var_name = 'ATEMP'
    start_iter = 0
    final_iter = 525600

    iter_count_check = 4#int(averaging_period/seconds_per_iter)

    print('Creating the exf field reference to cover iterations ' + str(start_iter) + ' to ' + str(final_iter))

    # calculate the file name iterations to read from
    #     along with the iterations these files cover
    file_names, iter_midpoint_dict = find_dv_files_to_read(config_dir, L05_model_name, averaging_period, seconds_per_iter,
                                                           var_name, start_iter, final_iter)
    print('  Source file summary:')
    output = '{\n'
    for file_name in file_names:
        print('    - The file ' + file_name + ' has iterations ' + str(np.min(iter_midpoint_dict[file_name])) + \
              ' to ' + str(np.max(iter_midpoint_dict[file_name])))
        output += ' \''+file_name.split('.')[-2]+'\': '+'['+str(np.min(iter_midpoint_dict[file_name]))+', '+str(np.max(iter_midpoint_dict[file_name]))+'],\n'
    output += '}'

    output_file = os.path.join(config_dir,'L05',L05_model_name,'run','dv','daily_exf_source_ref.txt')
    f = open(output_file,'w')
    f.write(output)
    f.close()

    # calculate the destination file iteration bounds
    #     along with the source files data they will contain, and which inidces of those source files they will contain
    dest_files, source_file_read_dict, source_file_read_index_sets = create_destination_file_list(config_dir,
                                                                                                  averaging_period,
                                                                                                  seconds_per_iter,
                                                                                                  var_name,
                                                                                                  start_iter, final_iter,
                                                                                                  file_names,
                                                                                                  iter_midpoint_dict)

    print(dest_files)
    print('  Destination file summary:')
    output = '{\n'
    for file_name in dest_files:
        print('    - The file ' + file_name + ' will be created from the following data:')
        output+=' \''+file_name.split('.')[-2]+'\': ['
        source_files = source_file_read_dict[file_name]
        index_sets = source_file_read_index_sets[file_name]
        iter_count = 0
        add_line = ''
        for s in range(len(source_files)):
            if len(index_sets)>0:
                print('         - From ' + source_files[s] + ', will read indices ' + str(
                    index_sets[s][0]) + ' though ' + str(index_sets[s][1]))
                iter_count += index_sets[s][1] - index_sets[s][0] + 1

                add_line += '[\'' + source_files[s].split('.')[-2] + '\', ' + '[' + str(index_sets[s][0]) + ', ' + str(index_sets[s][1]) + ']], '
        if int(file_name.split('.')[-2][-2:]) in [1,3,5,7,8,10,12]:
            if int(iter_count/iter_count_check)==1:
                output+=add_line[:-2]
            else:
                output+=''
        elif int(file_name.split('.')[-2][-2:]) in [4,6,9,11]:
            if int(iter_count / iter_count_check) == 1:
                output += add_line[:-2]
            else:
                output += ''
        else:
            if int(file_name.split('.')[-2][:4])%4==0:
                if int(iter_count / iter_count_check) == 1:
                    output += add_line[:-2]
                else:
                    output += ''
            else:
                if int(iter_count / iter_count_check) == 1:
                    output += add_line[:-2]
                else:
                    output += ''
        print('         - Total iterations for this file: '+str(iter_count)+' (= '+str(iter_count/iter_count_check)+' days)')
        output+='],\n'
    output += '}'

    output_file = os.path.join(config_dir,'L05',L05_model_name,'run','dv','daily_exf_dest_ref.txt')
    f = open(output_file,'w')
    f.write(output)
    f.close()


