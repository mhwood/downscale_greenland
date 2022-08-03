
import os
import numpy as np
import netCDF4 as nc4
import ast
import matplotlib.pyplot as plt
from MITgcmutils import mds
from pyproj import Transformer
import sys
import argparse

def read_grid_geometry(config_dir,model_name):

    file_path = os.path.join(config_dir, 'L3', model_name, 'input', model_name + '_grid.nc')
    ds = nc4.Dataset(file_path)
    XC = np.array(ds.variables['XC'][:, :])
    YC = np.array(ds.variables['YC'][:, :])
    delR = np.array(ds.variables['drF'][:])
    ds.close()

    return(XC, YC, delR)

def read_first_BC_field(config_dir,model_name,var_name,n_rows,n_cols,Nr):

    north_file = os.path.join(config_dir,'L3',model_name,'input','obcs',model_name+'_ASTE_'+var_name+'_north_BC.bin')
    north_grid = np.fromfile(north_file,'>f4')
    n_timesteps = int(np.size(north_grid)/(n_cols*Nr))
    north_grid = np.reshape(north_grid,(n_timesteps,Nr,n_cols))
    north_grid = north_grid[0,:,:]
    # north_grid[np.isnan(north_grid)] = 0
    # north_grid[north_grid>30] = 30

    # C = plt.imshow(north_grid)
    # plt.colorbar(C)
    # plt.show()

    south_file = os.path.join(config_dir, 'L3', model_name, 'input', 'obcs',
                              model_name + '_ASTE_' + var_name + '_south_BC.bin')
    south_grid = np.fromfile(south_file,'>f4')
    n_timesteps = int(np.size(south_grid) / (n_cols * Nr))
    south_grid = np.reshape(south_grid, (n_timesteps, Nr, n_cols))
    south_grid = south_grid[0, :, :]
    # south_grid[np.isnan(south_grid)] = 0
    # south_grid[south_grid > 30] = 30

    # C = plt.imshow(south_grid)
    # plt.colorbar(C)
    # plt.show()

    east_file = os.path.join(config_dir, 'L3', model_name, 'input', 'obcs',
                              model_name + '_ASTE_' + var_name + '_east_BC.bin')
    east_grid = np.fromfile(east_file,'>f4')
    n_timesteps = int(np.size(east_grid) / (n_rows * Nr))
    east_grid = np.reshape(east_grid, (n_timesteps, Nr, n_rows))
    east_grid = east_grid[0, :, :]
    # east_grid[np.isnan(east_grid)]=0
    # east_grid[east_grid>30] = 30

    # C = plt.imshow(east_grid)
    # plt.colorbar(C)
    # plt.show()

    return(north_grid,south_grid,east_grid)

def read_mask_from_grid_nc(nc_file, hFac='C'):
    ds = nc4.Dataset(nc_file)
    mask = ds.variables['HFac'+hFac][:,:,:]
    mask[mask>0]=1
    mask[mask<=0]=0
    ds.close()
    return(mask)

def calculate_mean_profile(north_profile,south_profile,east_profile):

    mean_profile = np.zeros((np.shape(north_profile)[0],))
    for k in range(np.shape(north_profile)[0]):
        boundary_values = np.vstack([np.reshape(north_profile[k,:],(np.size(north_profile[k,:]),1)),
                                     np.reshape(south_profile[k,:],(np.size(south_profile[k,:]),1)),
                                     np.reshape(east_profile[k,:],(np.size(east_profile[k,:]),1))])
        if np.any(boundary_values!=0) and k<36:
            boundary_values = boundary_values[boundary_values!=0]
            mean_profile[k] = np.mean(boundary_values)
        else:
            mean_profile[k]=mean_profile[k-1]

    plt.plot(mean_profile)
    plt.show()

    return(mean_profile)

def interpolate_weighted_pickup_grid(interp_field,XC_3413, YC_3413,domain_wet_cells_3D,
                                     mean_profile,north_profile,south_profile,east_profile):

    distance_threshold = 50000

    # south_profile[south_profile!=0] = 1
    # north_profile[north_profile != 0] = 1
    # east_profile[east_profile != 0] = 1

    for k in range(np.shape(domain_wet_cells_3D)[0]):
        rows, cols = np.where(domain_wet_cells_3D[k,:,:]!=0)
        if len(rows)>0:
            boundary_values = np.vstack([np.column_stack([XC_3413[-1,:].ravel(),YC_3413[-1,:].ravel(),north_profile[k,:].ravel()]),
                                         np.column_stack([XC_3413[0,:].ravel(),YC_3413[0,:].ravel(),south_profile[k,:].ravel()]),
                                         np.column_stack([XC_3413[:,-1].ravel(),YC_3413[:,-1].ravel(),east_profile[k,:].ravel()])])
            boundary_values = boundary_values[boundary_values[:,2]!=0,:]
            for ri in range(len(rows)):
                x = XC_3413[rows[ri],cols[ri]]
                y = YC_3413[rows[ri],cols[ri]]
                dist = ((boundary_values[:,0]-x)**2 + (boundary_values[:,1]-y)**2)**0.5
                if np.any(dist<distance_threshold):
                    boundary_values_subset = boundary_values[dist<distance_threshold,2]
                    dist_subset = dist[dist<distance_threshold]
                    dist_weights = 1/((1+dist_subset)**2)
                    interpolated_val = np.sum(boundary_values_subset*dist_weights)/np.sum(dist_weights)

                    weight = np.min(dist)/distance_threshold

                    assigned_val = (1-weight)*interpolated_val + weight*mean_profile[k]

                    interp_field[k, rows[ri], cols[ri]] = assigned_val

                else:
                    interp_field[k, rows[ri], cols[ri]] = mean_profile[k]

    # C = plt.imshow(interp_field[10,:,:],origin='lower')
    # plt.colorbar(C)
    # plt.show()
    return(interp_field)



def rotate_interpolated_grids_to_domain(var_names, var_grids, AngleCS, AngleSN):

    def rotate_velocity_vectors_to_domain(angle_cos, angle_sin, zonal_vel, meridional_vel):
        uvel = np.zeros_like(zonal_vel)
        vvel = np.zeros_like(meridional_vel)
        for k in range(np.shape(uvel)[0]):
            uvel[k, : ,:] = angle_cos * zonal_vel[k,:,:] + angle_sin * meridional_vel[k,:,:]
            vvel[k, : ,:] = -1 * angle_sin * zonal_vel[k,:,:] + angle_cos * meridional_vel[k,:,:]
        return (uvel, vvel)

    uvel_grid_index = var_names.index('Uvel')
    vvel_grid_index = var_names.index('Vvel')
    uvel, vvel = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                   var_grids[uvel_grid_index], var_grids[vvel_grid_index])

    gunm1_grid_index = var_names.index('GuNm1')
    gvnm1_grid_index = var_names.index('GvNm1')
    gunm1, gvnm1 = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                     var_grids[gunm1_grid_index], var_grids[gvnm1_grid_index])

    gunm2_grid_index = var_names.index('GuNm2')
    gvnm2_grid_index = var_names.index('GvNm2')
    gunm2, gvnm2 = rotate_velocity_vectors_to_domain(AngleCS, AngleSN,
                                                     var_grids[gunm2_grid_index], var_grids[gvnm2_grid_index])

    var_grids[uvel_grid_index] = uvel
    var_grids[vvel_grid_index] = vvel

    var_grids[gunm1_grid_index] = gunm1
    var_grids[gvnm1_grid_index] = gvnm1

    var_grids[gunm2_grid_index] = gunm2
    var_grids[gvnm2_grid_index] = gvnm2

    return(var_grids)

def stack_grids_to_pickup(interp_grids):
    counter = 0
    for grid in interp_grids:
        print(np.shape(grid))

        if counter == 0:
            pickup_grid = grid
        else:
            pickup_grid = np.concatenate([pickup_grid, grid], axis=0)

        counter += 1
    return(pickup_grid)

def write_pickup_file(output_file,dtype,pickup_grid,subset_metadata):

    # output the data subset
    pickup_grid.ravel(order='C').astype(dtype).tofile(output_file+'.data')

    # output the metadata file
    output = " nDims = [   "+str(subset_metadata['ndims'][0])+" ];\n"
    output += " dimList = [\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[2])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[2])+",\n"
    output += " "+"{:5d}".format(np.shape(pickup_grid)[1])+",    1,"+"{:5d}".format(np.shape(pickup_grid)[1])+"\n"
    output += " ];\n"
    output += " dataprec = [ '"+subset_metadata['dataprec'][0]+"' ];\n"
    output += " nrecords = [   "+str(subset_metadata['nrecords'][0])+" ];\n"
    output += " timeStepNumber = [ "+"{:10d}".format(subset_metadata['timestepnumber'][0])+" ];\n"
    time_interval_exponent = int(np.log10(subset_metadata['timeinterval'][0][0]))
    time_interval_base = subset_metadata['timeinterval'][0][0] / (10 ** time_interval_exponent)
    output += " timeInterval = [  "+"{:.12f}".format(time_interval_base) + "E+" + "{:02d}".format(time_interval_exponent)+  " ];\n"
    output += " nFlds = [   "+str(subset_metadata['nflds'][0])+" ];\n"
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


def create_L3_ASTE_pickup_file(config_dir):

    model_name = 'L3_Scoresby_Sund'

    sys.path.insert(1, os.path.join(config_dir, 'utils', 'init_file_creation'))
    import downscale_functions as df

    print('    - Creating the pickup file for the '+model_name+' model from ASTE data')


    f = open(os.path.join(config_dir, 'domain_sizes.txt'))
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)
    L_size = size_dict[model_name]
    n_rows = L_size[0]
    n_cols = L_size[1]
    Nr = 90

    # step 0: get the model domain
    XC, YC, delR = read_grid_geometry(config_dir,model_name)

    transformer = Transformer.from_crs('EPSG:' + str(4326), 'EPSG:' + str(3413))
    XC_3413, YC_3413 = transformer.transform(YC.ravel(), XC.ravel())
    YC_3413 = np.reshape(YC_3413, np.shape(XC))
    XC_3413 = np.reshape(XC_3413, np.shape(XC))

    var_names = ['Uvel', 'Vvel', 'Theta', 'Salt', 'GuNm1', 'GuNm2', 'GvNm1', 'GvNm2', 'EtaN', 'dEtaHdt', 'EtaH']

    # plt.imshow(XC_3413)
    # # plt.plot(XC[:, 0], YC[:, 0], 'k-')
    # # plt.plot(XC[:, -1], YC[:, -1], 'k-')
    # # plt.plot(XC[0, :], YC[0, :], 'k-')
    # # plt.plot(XC[-1, :], YC[-1, :], 'k-')
    # plt.show()

    domain_grid_file = os.path.join(config_dir, 'L3',model_name, 'input', model_name+'_grid.nc')

    print('    - Downscaling the pickup grids')
    interp_grids = []
    output_var_names = []
    for vn in range(len(var_names)):
        var_name = var_names[vn]

        if var_name.lower() not in ['etan', 'detahdt', 'etah']:
            interp_field = np.zeros((len(delR), np.shape(XC)[0], np.shape(XC)[1]))
        else:
            interp_field = np.zeros((1, np.shape(XC)[0], np.shape(XC)[1]))

        if var_name in ['Theta','Salt']:  # used for testing

            print('      - Creating the ' + var_name+' pickup field')
            if var_name in ['Vvel', 'GvNm1', 'GvNm2']:
                domain_wet_cells_3D = read_mask_from_grid_nc(domain_grid_file, hFac='S')
                domain_wet_cells_3D = domain_wet_cells_3D[:, :-1, :]
            elif var_name in ['Uvel', 'GuNm1', 'GuNm2']:
                domain_wet_cells_3D = read_mask_from_grid_nc(domain_grid_file, hFac='W')
                domain_wet_cells_3D = domain_wet_cells_3D[:, :, :-1]
            else:
                domain_wet_cells_3D = read_mask_from_grid_nc(domain_grid_file, hFac='C')

            # plt.subplot(1,2,1)
            # plt.imshow(subset_copy[:,10,:])
            # plt.subplot(1, 2, 2)
            # plt.imshow(aste_grid[:, 10, :])
            # plt.show()

            north_profile,south_profile,east_profile = read_first_BC_field(config_dir, model_name, var_name, n_rows, n_cols, Nr)


            if var_name in ['Theta','Salt']:
                mean_profile = calculate_mean_profile(north_profile,south_profile,east_profile)
            else:
                mean_profile = np.zeros((Nr, ))

            interp_field = interpolate_weighted_pickup_grid(interp_field,XC_3413, YC_3413, domain_wet_cells_3D,
                                             mean_profile,north_profile,south_profile,east_profile)

        interp_grids.append(interp_field)
        output_var_names.append(var_name)

    pickup_grid = stack_grids_to_pickup(interp_grids)

    output_dir = os.path.join(config_dir, 'L3', model_name, 'input')
    output_file = os.path.join(output_dir, 'pickup.' + '{:010d}'.format(1))
    dtype = '>f8'

    pickup_metadata = {'ndims': [2],
                       'dataprec': ['float64'],
                       'nrecords': [int((len(output_var_names)-3) * len(delR) + 3)],
                       'nflds': [11],
                       'fldlist': output_var_names,
                       'dimlist': [90, 4050],
                       'timestepnumber': [1],
                       'timeinterval': [[2678400.0]]}

    write_pickup_file(output_file, dtype, pickup_grid, pickup_metadata)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_L3_ASTE_pickup_file(config_dir)