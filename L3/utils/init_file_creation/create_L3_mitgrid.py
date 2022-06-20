
import os
import simplegrid as sg
import numpy as np
import argparse
import ast

def create_L3_mitgrid_file(config_dir,model_name,Lat_C, Lon_C,Lat_G, Lon_G):

    print('    - Using simplegrid to create an mitgrid file for the '+model_name+' model')

    Lon_C = np.flipud(Lon_C)
    Lat_C = np.flipud(Lat_C)
    Lon_G = np.flipud(Lon_G)
    Lat_G = np.flipud(Lat_G)

    mitgrid_matrices = dict()
    mitgrid_matrices['XG'] = Lon_G
    mitgrid_matrices['YG'] = Lat_G
    mitgrid_matrices['XC'] = Lon_C
    mitgrid_matrices['YC'] = Lat_C


    if 'mitgrids' not in os.listdir(config_dir):
        os.mkdir(os.path.join(config_dir,'mitgrids'))

    output_file = os.path.join(config_dir,'mitgrids',model_name+'.mitgrid')

    factor = 1
    mitgrid, n_rows, n_cols = \
        sg.regrid.regrid(mitgrid_matrices=mitgrid_matrices, \
                         lon_subscale=factor, lat_subscale=factor, \
                         lon1=Lon_G[0, 0], lat1=Lat_G[0, 0],
                         lon2=Lon_G[-1, -1], lat2=Lat_G[-1, -1],
                         verbose=False, outfile=output_file)

    # mitgrid, n_rows, n_cols = sg.mkgrid.mkgrid(
    #     lon1=Lon_G[0, 0], lat1=Lat_G[0, 0],
    #     lon2=Lon_G[-1, -1], lat2=Lat_G[-1, -1],
    #     lon_subscale=110, lat_subscale=180)

    sg.gridio.write_mitgridfile(output_file, mitgrid, n_rows, n_cols)


    if 'domain_sizes.txt' in os.listdir(config_dir):
        f = open(os.path.join(config_dir, 'domain_sizes.txt'))
        dict_str = f.read()
        f.close()
        size_dict = ast.literal_eval(dict_str)
    else:
        size_dict = {}

    size_dict[model_name] = (n_cols, n_rows)
    f = open(os.path.join(config_dir, 'domain_sizes.txt'),'w')
    f.write(str(size_dict))
    f.close()

   

