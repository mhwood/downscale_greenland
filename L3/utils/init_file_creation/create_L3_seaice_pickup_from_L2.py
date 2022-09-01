
import os
import numpy as np
import netCDF4 as nc4
import ast
import matplotlib.pyplot as plt
from MITgcmutils import mds
import sys

def read_grid_geometry_from_nc(config_dir, model_name):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = ds.variables['XC'][:,:]
    YC = ds.variables['YC'][:,:]
    delR = ds.variables['drF'][:]
    ds.close()
    return(XC, YC, delR)

def read_wetgrid_from_nc(config_dir, model_name, hFac):
    file_path = os.path.join(config_dir, 'nc_grids', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    hFac_grid = ds.variables['HFac'+hFac][:,:]
    ds.close()

    wet_grid = np.copy(hFac_grid[:1,:,:])
    wet_grid[wet_grid>0] = 1

    if hFac=='S':
        wet_grid = wet_grid[:,:-1,:]
    if hFac=='W':
        wet_grid = wet_grid[:,:,:-1]

    return(wet_grid)

def read_seaice_pickup_file(seaice_pickup_file_path):

    global_data, _, global_metadata = mds.rdmds(seaice_pickup_file_path, returnmeta=True)

    var_names = []
    row_bounds = []
    var_grids = []

    start_row = 0
    for var_name in global_metadata['fldlist']:
        end_row = start_row + 1
        var_grid = global_data[start_row:end_row,:,:]
        var_grids.append(var_grid)
        row_bounds.append([start_row,end_row])
        start_row=end_row
        var_names.append(var_name.strip())

    return(var_names,row_bounds,var_grids,global_metadata)

def stack_grids_to_pickup(interp_grids):
    counter = 0
    for grid in interp_grids:

        if counter == 0:
            pickup_grid = grid
        else:
            pickup_grid = np.concatenate([pickup_grid, grid], axis=0)

        counter += 1
    return(pickup_grid)

def write_seaice_pickup_file(output_file,dtype,pickup_grid,subset_metadata):

    # output the data subset
    pickup_grid.ravel(order='C').astype(dtype).tofile(output_file+'.data')

    # output the metadata file
    output = " nDims = [   "+str(subset_metadata['ndims'][0])+" ];\n"
    output += " dimList = [\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[2])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[2])+",\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[1])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[1])+"\n"
    output += " ];\n"
    output += " dataprec = [ '"+subset_metadata['dataprec'][0]+"' ];\n"
    output += " nrecords = [    "+str(subset_metadata['nrecords'][0])+" ];\n"
    output += " timeStepNumber = [ "+"{:10d}".format(subset_metadata['timestepnumber'][0])+" ];\n"
    time_interval_exponent = int(np.log10(subset_metadata['timeinterval'][0][0]))
    time_interval_base = subset_metadata['timeinterval'][0][0] / (10 ** time_interval_exponent)
    output += " timeInterval = [  "+"{:.12f}".format(time_interval_base) + "E+" + "{:02d}".format(time_interval_exponent)+  " ];\n"
    output += " nFlds = [    "+str(subset_metadata['nflds'][0])+" ];\n"
    output += " fldList = {\n "
    for var_name in subset_metadata['fldlist']:
        output += "'"+var_name
        for i in range(8-len(var_name)):
            output+= " "
        output+="' "
    output += "\n };"

    f = open(output_file+'.meta','w')
    f.write(output)
    f.close()

def create_L3_seaice_pickup_file(config_dir, L3_model_name, L2_model_name, iter_number, print_level):

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df

    if print_level >= 1:
        print('    - Creating the seaice pickup file for the ' + L3_model_name + ' model from L2 output')

    # step 0: get the model domain
    L2_XC, L2_YC, L2_delR = read_grid_geometry_from_nc(config_dir, L2_model_name)
    L3_XC, L3_YC, L3_delR = read_grid_geometry_from_nc(config_dir, L3_model_name)

    if print_level >= 1:
        print('    - Reading in variables from the L2 seaice pickup file')
    seaice_pickup_file_path = os.path.join(config_dir,'L2',L2_model_name,'run','pickup_seaice.'+'{:010d}'.format(iter_number))
    var_names, _, var_grids, pickup_metadata = read_seaice_pickup_file(seaice_pickup_file_path)

    # for i in range(len(pickup_var_names)):
    #     C = plt.imshow(pickup_pacific_grids[i][0, :, :])
    #     plt.title(pickup_var_names[i])
    #     plt.colorbar(C)
    #     plt.show()


    min_row = np.min([np.argmin(np.abs(L2_YC[:, 0] - np.min(L3_YC))),
                      np.argmin(np.abs(L2_YC[:, -1] - np.min(L3_YC)))])
    max_row = np.max([np.argmin(np.abs(L2_YC[:, 0] - np.max(L3_YC))),
                      np.argmin(np.abs(L2_YC[:, -1] - np.max(L3_YC)))])

    min_col = np.min([np.argmin(np.abs(L2_XC[0, :] - np.min(L3_XC))),
                      np.argmin(np.abs(L2_XC[-1, :] - np.min(L3_XC)))])
    max_col = np.max([np.argmin(np.abs(L2_XC[0, :] - np.max(L3_XC))),
                      np.argmin(np.abs(L2_XC[-1, :] - np.max(L3_XC)))])

    buffer = 3
    for i in range(buffer):
        if min_row>0:
            min_row -= 1
    for i in range(buffer):
        if max_row<np.shape(L3_XC)[0]-1:
            max_row += 1
    for i in range(buffer):
        if min_col>0:
            min_col -= 1
    for i in range(buffer):
        if max_col<np.shape(L3_XC)[1]-1:
            max_col += 1

    L2_XC = L2_XC[min_row:max_row, :]
    L2_YC = L2_YC[min_row:max_row, :]
    L2_XC = L2_XC[:, min_col:max_col]
    L2_YC = L2_YC[:, min_col:max_col]

    if print_level >= 1:
        print('    - Downscaling the pickup grids')
    interp_grids = []
    for vn in range(len(var_names)):
        var_name = var_names[vn]

        # if var_name in ['EtaN']:
        L2_var = var_grids[vn]
        if print_level >= 2:
            print('        - Downscaling ' + var_name)

        if var_name in ['siVICE']:
            L2_wet_cells = read_wetgrid_from_nc(config_dir, L2_model_name, hFac='S')
            # L2_wet_cells_on_L3 = read_mask_from_nc(os.path.join('..', 'input', 'L2_wetgrid_on_L3.nc'), hFac='S')
            L3_wet_cells = read_wetgrid_from_nc(config_dir, L3_model_name, hFac='S')
        elif var_name in ['siUICE']:
            L2_wet_cells = read_wetgrid_from_nc(config_dir, L2_model_name, hFac='W')
            # L2_wet_cells_on_L3 = read_mask_from_nc(os.path.join('..', 'input', 'L2_wetgrid_on_L3.nc'), hFac='W')
            L3_wet_cells = read_wetgrid_from_nc(config_dir, L3_model_name, hFac='W')
        else:
            L2_wet_cells = read_wetgrid_from_nc(config_dir, L2_model_name, hFac='C')
            # L2_wet_cells_on_L3 = read_mask_from_nc(os.path.join('..', 'input', 'L2_wetgrid_on_L3.nc'), hFac='C')
            L3_wet_cells = read_wetgrid_from_nc(config_dir, L3_model_name, hFac='C')

        L2_wet_cells = L2_wet_cells[:, min_row:max_row, :]
        L2_wet_cells = L2_wet_cells[:, :, min_col:max_col]

        L2_var = L2_var[:, min_row:max_row, :]
        L2_var = L2_var[:, :, min_col:max_col]

        L2_wet_cells_on_L3 = np.copy(L3_wet_cells)

        # plt.imshow(L2_wet_cells[0,:,:])
        # plt.show()

        if print_level >= 3:
            print('        - Variable shapes:')
            print('            - L2_XC: ' + str(np.shape(L2_XC)))
            print('            - L2_YC: ' + str(np.shape(L2_YC)))
            print('            - L2_var: ' + str(np.shape(L2_var)))
            print('            - L2_wet_cells: ' + str(np.shape(L2_wet_cells)))
            # print('            L2_wet_cells_on_L3: ' + str(np.shape(L2_wet_cells_on_L3)))
            print('            - L3_XC: ' + str(np.shape(L3_XC)))
            print('            - L3_YC: ' + str(np.shape(L3_YC)))
            print('            - L3_wet_cells ' + str(np.shape(L3_wet_cells)))

        if var_name in ['siAREA','siHEFF','siHSNOW']:
            remove_zeros = False
        else:
            remove_zeros = True

        if print_level >= 4:
            printing = True
        else:
            printing = False

        interp_field = df.downscale_3D_field(L2_XC, L2_YC, L2_var,
                                          L2_wet_cells, L2_wet_cells_on_L3,
                                          L3_XC, L3_YC, L3_wet_cells,
                                          printing=printing, remove_zeros = remove_zeros)

        # plt.subplot(1,2,1)
        # plt.imshow(L2_var[0, :,:],origin='lower')
        # plt.subplot(1,2,2)
        # plt.imshow(interp_field[0,:,:],origin='lower')
        # plt.show()

        # else:
        #     if var_name.lower() not in ['etan','detahdt','etah']:
        #         interp_field = np.zeros((Nr,np.shape(L3_XC)[0],np.shape(L3_XC)[1]))
        #     else:
        #         interp_field = np.zeros((1,np.shape(L3_XC)[0],np.shape(L3_XC)[1]))

        interp_grids.append(interp_field)

    pickup_grid = stack_grids_to_pickup(interp_grids)

    output_dir = os.path.join(config_dir,'L3', L3_model_name, 'input')
    output_file = os.path.join(output_dir, 'pickup_seaice.'+'{:010d}'.format(2*iter_number))
    pickup_metadata['timestepnumber'] = [2 * iter_number]
    dtype = '>f8'
    write_seaice_pickup_file(output_file, dtype, pickup_grid, pickup_metadata)




   

