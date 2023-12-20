
import os
import numpy as np
import netCDF4 as nc4
from pyproj import Transformer
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import ast

def read_grid_geometry_from_nc(config_dir,model_name):
    nc_file = os.path.join(config_dir,'nc_grids',model_name+'_grid.nc')
    ds = nc4.Dataset(nc_file)
    XC = ds.variables['XC'][:, :]
    YC = ds.variables['YC'][:, :]
    Depth = ds.variables['Depth'][:, :]
    drF = ds.variables['drF'][:]
    HFacC = ds.variables['HFacC'][:, :, :]
    return(XC, YC, Depth, HFacC, drF)

def bulkmodjmd95(s,theta,p):
    """ Compute bulk modulus
    """

    # coefficients in pressure coordinates for
    # 3. secant bulk modulus K of fresh water at p = 0
    eosJMDCKFw = [1.965933e+04,
                  1.444304e+02,
                  - 1.706103e+00,
                  9.648704e-03,
                  - 4.190253e-05,
                  ]
    # 4. secant bulk modulus K of sea water at p = 0
    eosJMDCKSw = [5.284855e+01,
                  - 3.101089e-01,
                  6.283263e-03,
                  - 5.084188e-05,
                  3.886640e-01,
                  9.085835e-03,
                  - 4.619924e-04,
                  ]
    # 5. secant bulk modulus K of sea water at p
    eosJMDCKP = [3.186519e+00,
                 2.212276e-02,
                 - 2.984642e-04,
                 1.956415e-06,
                 6.704388e-03,
                 - 1.847318e-04,
                 2.059331e-07,
                 1.480266e-04,
                 2.102898e-04,
                 - 1.202016e-05,
                 1.394680e-07,
                 - 2.040237e-06,
                 6.128773e-08,
                 6.207323e-10,
                 ]

    # make sure arguments are floating point
    s = np.asfarray(s)
    t = np.asfarray(theta)
    p = np.asfarray(p)

    t2 = t*t
    t3 = t2*t
    t4 = t3*t

#    if np.any(s<0):
#        sys.stderr.write('negative salinity values! setting to nan\n')
#       the sqrt will take care of this
#        if s.ndim > 0:
#            s[s<0] = np.nan
#        else:
#            s = np.nan

    s3o2 = s*np.sqrt(s)

    #p = pressure(i,j,k,bi,bj)*SItoBar
    p2 = p*p
    # secant bulk modulus of fresh water at the surface
    bulkmod = ( eosJMDCKFw[0]
              + eosJMDCKFw[1]*t
              + eosJMDCKFw[2]*t2
              + eosJMDCKFw[3]*t3
              + eosJMDCKFw[4]*t4
              )
    # secant bulk modulus of sea water at the surface
    bulkmod = ( bulkmod
              + s*(      eosJMDCKSw[0]
                       + eosJMDCKSw[1]*t
                       + eosJMDCKSw[2]*t2
                       + eosJMDCKSw[3]*t3
                       )
              + s3o2*(   eosJMDCKSw[4]
                       + eosJMDCKSw[5]*t
                       + eosJMDCKSw[6]*t2
                       )
               )
    # secant bulk modulus of sea water at pressure p
    bulkmod = ( bulkmod
              + p*(   eosJMDCKP[0]
                    + eosJMDCKP[1]*t
                    + eosJMDCKP[2]*t2
                    + eosJMDCKP[3]*t3
                  )
              + p*s*(   eosJMDCKP[4]
                      + eosJMDCKP[5]*t
                      + eosJMDCKP[6]*t2
                    )
              + p*s3o2*eosJMDCKP[7]
              + p2*(   eosJMDCKP[8]
                     + eosJMDCKP[9]*t
                     + eosJMDCKP[10]*t2
                   )
              + p2*s*(  eosJMDCKP[11]
                      + eosJMDCKP[12]*t
                      + eosJMDCKP[13]*t2
                     )
               )

    return bulkmod

def densjmd95(s, theta, p):
    """
    Computes in-situ density of sea water

    Density of Sea Water using Jackett and McDougall 1995 (JAOT 12)
    polynomial (modified UNESCO polynomial).

    Parameters
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    theta : array_like
        potential temperature [degree C (IPTS-68)];
        same shape as s
    p : array_like
        pressure [dbar]; broadcastable to shape of s

    Returns
    -------
    dens : array
        density [kg/m^3]

    Example
    -------
    >>> densjmd95(35.5, 3., 3000.)
    1041.83267

    Notes
    -----
    AUTHOR:  Martin Losch 2002-08-09  (mlosch@mit.edu)

    Jackett and McDougall, 1995, JAOT 12(4), pp. 381-388
    """

    eosJMDCFw = [999.842594,6.793952e-02,
                 -9.095290e-03,
                 1.001685e-04,
                 -1.120083e-06,
                 6.536332e-09]
    # 2. density of sea water at p = 0
    eosJMDCSw = [8.244930e-01,
                 -4.089900e-03,
                 7.643800e-05,
                 -8.246700e-07,
                 5.387500e-09,
                 -5.724660e-03,
                 1.022700e-04,
                 -1.654600e-06,
                 4.831400e-04]

    # make sure arguments are floating point
    s = np.asfarray(s)
    t = np.asfarray(theta)
    p = np.asfarray(p)

    # convert pressure to bar
    p = .1 * p

    t2 = t * t
    t3 = t2 * t
    t4 = t3 * t

    if np.any(s < 0):
        sys.stderr.write('negative salinity values! setting to nan\n')
    #       the sqrt will take care of this
    #        if s.ndim > 0:
    #            s[s<0] = np.nan
    #        else:
    #            s = np.nan

    s3o2 = s * np.sqrt(s)

    # density of freshwater at the surface
    rho = (eosJMDCFw[0]
           + eosJMDCFw[1] * t
           + eosJMDCFw[2] * t2
           + eosJMDCFw[3] * t3
           + eosJMDCFw[4] * t4
           + eosJMDCFw[5] * t4 * t
           )
    # density of sea water at the surface
    rho = (rho
           + s * (
                   eosJMDCSw[0]
                   + eosJMDCSw[1] * t
                   + eosJMDCSw[2] * t2
                   + eosJMDCSw[3] * t3
                   + eosJMDCSw[4] * t4
           )
           + s3o2 * (
                   eosJMDCSw[5]
                   + eosJMDCSw[6] * t
                   + eosJMDCSw[7] * t2
           )
           + eosJMDCSw[8] * s * s
           )

    rho = rho / (1. - p / bulkmodjmd95(s, t, p))

    return rho

def density_equation(k, s, t, p, rhoConst, equation_of_state):
    if equation_of_state=='linear':
        drho = rhoConst * (1 - talpha * (t[k] - tref[k]) + sbeta * (s[k] - sref[k])) - rhoConst
    elif equation_of_state=='jmd95':
        drho = densjmd95(s[k], t[k], p[k]) - rhoConst
    else:
        raise ValueError('EOS not implemented')
    return(drho)

def reproject_array(X_array, Y_array, inputCRS,outputCRS):

    transformer = Transformer.from_crs('EPSG:' + str(inputCRS), 'EPSG:' + str(outputCRS))

    x = np.ravel(X_array)
    y = np.ravel(Y_array)

    if inputCRS == 4326 and outputCRS == 3413:
        x2, y2 = transformer.transform(y, x)
        x2 = np.array(x2)
        y2 = np.array(y2)

    X_out = np.reshape(x2,np.shape(X_array))
    Y_out = np.reshape(y2,np.shape(Y_array))

    # plt.subplot(2, 2, 1)
    # C = plt.imshow(XC, origin='lower')
    # plt.colorbar(C)
    # plt.subplot(2, 2, 2)
    # C = plt.imshow(X, origin='lower')
    # plt.colorbar(C)
    # plt.subplot(2, 2, 3)
    # C = plt.imshow(YC, origin='lower')
    # plt.colorbar(C)
    # plt.subplot(2, 2, 4)
    # C = plt.imshow(Y, origin='lower')
    # plt.colorbar(C)
    # plt.show()

    return(X_out, Y_out)

def create_mask_subset(mask_file, mask_file_subset, extents, years):

    #####################
    # open the main file and subset to extents

    print('    - Reading in the time array')
    ds = nc4.Dataset(mask_file)
    time = ds.variables['time'][:]
    ds.close()

    date_indices = []
    for index in range(len(time)):
        days = time[index]
        full_date = dt.datetime(1900,1,1) + dt.timedelta(days = days)
        if full_date.month==1:
            for yr in years:
                if full_date.year==yr:
                    date_indices.append(index)

    date_indices = np.array(date_indices)
    # years = years[:2]

    for di in range(len(date_indices)):

        print('      - Reading in the big mask arrays for index '+str(date_indices[di]))
        ds = nc4.Dataset(mask_file)
        if di==0:
            x = ds.variables['x'][:]
            y = ds.variables['y'][:]
            time = ds.variables['time'][:]
            rock = ds.variables['rock'][:, :]
        ice = ds.variables['ice'][date_indices[di], :, :]
        ds.close()

        if di == 0:
            time = time[date_indices]

            min_x_index = np.argmin(np.abs(x-extents[0]))
            max_x_index = np.argmin(np.abs(x-extents[2]))
            max_y_index = np.argmin(np.abs(y-extents[1]))
            min_y_index = np.argmin(np.abs(y-extents[3]))

            x_subset = x[min_x_index:max_x_index]
            y_subset = y[min_y_index:max_y_index]
            rock_subset = rock[min_y_index:max_y_index, min_x_index:max_x_index]

        ice_subset = ice[min_y_index:max_y_index, min_x_index:max_x_index]

        if di==0:
            full_data = np.reshape(ice_subset,(1,np.shape(ice_subset)[0], np.shape(ice_subset)[1]))
        else:
            full_data = np.concatenate([full_data, np.reshape(ice_subset,(1,np.shape(ice_subset)[0], np.shape(ice_subset)[1]))], axis=0)

    ##################################
    # write out to a small file

    print('    - Writing out the small file')

    ds = nc4.Dataset(mask_file_subset,'w')

    ds.createDimension('x',len(x_subset))
    ds.createDimension('y',len(y_subset))
    ds.createDimension('time',len(time))

    xvar = ds.createVariable('x','f4',('x',))
    xvar[:] = x_subset

    yvar = ds.createVariable('y','f4',('y',))
    yvar[:] = y_subset

    tvar = ds.createVariable('years','f4',('time',))
    tvar[:] = years

    ivar = ds.createVariable('ice','i2',('time','y','x'))
    ivar[:, : ,:] = full_data

    rvar = ds.createVariable('rock','i2',('y','x'))
    rvar[:,:] = rock_subset

    ds.close()

def read_mask_subset(mask_file, mask_file_subset, extents, years):

    if not os.path.exists(mask_file_subset):
        create_mask_subset(mask_file, mask_file_subset, extents, years)

    ds = nc4.Dataset(mask_file_subset)
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    time = ds.variables['years'][:]
    rock = ds.variables['rock'][:, :]
    ice = ds.variables['ice'][:, :, :]
    ds.close()

    return(x, y, time, rock, ice)

def read_ice_thickness_grids(output_folder, year):

    file_list = ['79n_icedraft.nc',
                 'zachariae_icedraft.nc']

    for file_name in file_list:
        ds = nc4.Dataset(os.path.join(output_folder,'Ice Thickness',file_name))
        time = ds.variables['time'][:]
        x = ds.variables['x'][:]
        y = ds.variables['y'][:]
        thickness = ds.variables['icedraft'][:, :, :]
        ds.close()

        year_index = np.argmin(np.abs(time-year))
        thickness = thickness[year_index, :, :]

        X, Y = np.meshgrid(x,y)
        points = np.column_stack([X.ravel(),Y.ravel()])

        if file_name == file_list[0]:
            all_points = points
            all_thickness = np.reshape(thickness,(np.size(thickness),1))
        else:
            all_points = np.vstack([all_points, points])
            all_thickness = np.vstack([all_thickness, np.reshape(thickness, (np.size(thickness), 1))])

    nonzero_indices = all_thickness.ravel()>0
    all_points = all_points[nonzero_indices,:]
    all_thickness = all_thickness[nonzero_indices]

    return(all_points, all_thickness)

def generate_connected_mask(start_row, start_col, ice_mask):

    if ice_mask[start_row,start_col]==0:
        mask_grid = np.zeros_like(ice_mask)
    else:
        rows = np.arange(np.shape(ice_mask)[0])
        cols = np.arange(np.shape(ice_mask)[1])
        Cols,Rows = np.meshgrid(cols,rows)

        mask_grid = 1-np.copy(ice_mask)
        mask_grid[start_row,start_col] = 2
        # in the mask, 0 means unverified
        # 1 is verified not shelfice
        # 2 is verified shelfice

        # plt.imshow(mask_grid)
        # plt.show()

        is_remaining = np.logical_and(mask_grid==0,ice_mask==1)
        n_remaining = np.sum(is_remaining)
        # print(n_remaining)
        continue_iter = True
        for i in range(n_remaining):
            if continue_iter:
                # get the wet rows, cols, and their current mask values
                Wet_Rows = Rows[ice_mask == 1]
                Wet_Cols = Cols[ice_mask == 1]
                Mask_Vals = mask_grid[ice_mask == 1]

                # reduce these to the ones that havent been verified yet
                Wet_Rows = Wet_Rows[Mask_Vals == 0]
                Wet_Cols = Wet_Cols[Mask_Vals == 0]
                Mask_Vals = Mask_Vals[Mask_Vals == 0]

                if len(Mask_Vals)>0:

                    # for each row/col, see if its connected to one we've verified is connected
                    rows_remaining,cols_remaining = np.where(is_remaining)
                    for ri in range(n_remaining):
                        row = rows_remaining[ri]
                        col = cols_remaining[ri]

                        # this bit allows for only up/dow/left/right spreading
                        if row<np.shape(ice_mask)[0]-1:
                            if mask_grid[row+1,col] == 2:
                                mask_grid[row,col] = 2
                        if row > 0:
                            if mask_grid[row - 1, col] == 2:
                                mask_grid[row,col] = 2
                        if col<np.shape(ice_mask)[1]-1:
                            if mask_grid[row,col+1] == 2:
                                mask_grid[row,col] = 2
                        if col > 0:
                            if mask_grid[row, col-1] == 2:
                                mask_grid[row,col] = 2


                    is_remaining = np.logical_and(mask_grid == 0, ice_mask == 1)
                    n_remaining_now = np.sum(is_remaining)

                    # plt.subplot(1,2,1)
                    # plt.imshow(ice_mask,cmap='Greys_r')
                    # plt.subplot(1, 2, 2)
                    # plt.imshow(mask_grid)
                    # plt.show()

                    if n_remaining_now<n_remaining:
                        n_remaining = n_remaining_now
                    else:
                        n_remaining = n_remaining_now
                        continue_iter=False
                else:
                    continue_iter = False

    return(mask_grid)

def compute_annual_mask(X, Y, ice_mask_on_domain, year):
    # print('        - Spreading mask otuward from glacier points')
    glacier_points = [[489269, -1090645], [458469, -1036779]]
    glacier_names = ['Zachariae', '79N']

    full_ice_mask = np.zeros_like(ice_mask_on_domain)

    for g in range(len(glacier_names)):
        # print('        - Spreading in mask for ' + glacier_names[g])
        point = glacier_points[g]
        dist = ((X - point[0]) ** 2 + (Y - point[1]) ** 2) ** 0.5
        start_row, start_col = np.where(dist == np.min(dist))
        start_row = start_row[0]
        start_col = start_col[0]

        mask = generate_connected_mask(start_row, start_col, ice_mask_on_domain)
        # print(start_row, start_col)

        glacier_mask = np.copy(ice_mask_on_domain)
        glacier_mask[mask == 0] = 0

        full_ice_mask[glacier_mask == 1] = 1

    if year<=2014:
        full_ice_mask[367:372,176:182] = 1

    # C = plt.imshow(full_ice_mask, origin='lower')
    # plt.colorbar(C)
    # plt.show()
    return(full_ice_mask)

def create_icetopo_grid_from_bedmachine():

    if print_level >= 2:
        print('        - Interpolating from the BedMachine dataset')

    ice_topo = np.zeros_like(XC)

    bedMachine_file = '/Users/mhwood/Documents/Research/Data Repository/Greenland/Bathymetry/BedMachineGreenland-v5.nc'

    ds = nc4.Dataset(bedMachine_file)
    x = ds.variables['x'][:]
    y = ds.variables['y'][:]
    mask = ds.variables['mask'][:, :]
    surface = ds.variables['surface'][:, :]
    thickness = ds.variables['thickness'][:, :]
    ds.close()

    min_x_index = np.argmin(np.abs(x - np.min(XC)))
    max_x_index = np.argmin(np.abs(x - np.max(XC)))
    max_y_index = np.argmin(np.abs(y - np.min(YC)))
    min_y_index = np.argmin(np.abs(y - np.max(YC)))

    x = x[min_x_index:max_x_index]
    y = y[min_y_index:max_y_index]
    mask = mask[min_y_index:max_y_index, min_x_index:max_x_index]
    surface = surface[min_y_index:max_y_index, min_x_index:max_x_index]
    thickness = thickness[min_y_index:max_y_index, min_x_index:max_x_index]

    # this is a manual nearest neighbor method
    bedmachine_resolution = 150
    for i in range(np.shape(XC)[0]):
        for j in range(np.shape(YC)[1]):
            x_index = np.argmin(np.abs(x - XC[i, j]))
            y_index = np.argmin(np.abs(y - YC[i, j]))
            if mask[y_index, x_index] == 3 and Depth[i, j] > 0:
                ice_depth = surface[y_index, x_index] - thickness[y_index, x_index]
                if ice_depth < -1 * Depth[i, j]:
                    ice_depth = -1 * Depth[i, j]
                ice_topo[i, j] = ice_depth

def create_surface_load_file(icetopo, hFacC, Nr, dz, tref, sref, rhoConst, equation_of_state = 'jmd95'):

    # define the vertical grid
    zgp1 = np.cumsum(dz)
    zgp1 = np.insert(zgp1, 0, 0)
    zc = 0.5 * (zgp1[:-1] + zgp1[1:])
    zg = zgp1[:-1]

    # define the reference vertical profiles
    gravity = 9.81
    t = tref
    s = sref
    dzm = 0.5 * np.diff(zc)
    dzm = np.insert(dzm, 0, np.abs(zg[0] - zc[0]))
    dzp = 0.5 * np.diff(zc)
    dzp = np.insert(dzp, -1, np.abs(zg[-1] - zc[-1]))
    p = np.abs(zc) * gravity * rhoConst * 1e-4

    # compute the phiHydC and phiHydF
    dp = p
    kp = -1
    phiHydC = np.zeros((Nr,))
    phiHydF = np.zeros((Nr + 1,))
    for iteration in range(50):
        p0 = p
        kp = kp + 1
        for k in range(Nr):
            drho = density_equation(k, s, t, p, rhoConst, equation_of_state)
            phiHydC[k] = phiHydF[k] + dzm[k] * gravity * drho / rhoConst
            phiHydF[k + 1] = phiHydC[k] + dzp[k] * gravity * drho / rhoConst
        dp = p - p0
        RMS = np.sqrt(np.mean(np.sum(dp ** 2)))

        # plt.subplot(1, 2, 1)
        # plt.plot(phiHydC, zg)
        # plt.gca().invert_yaxis()
        # plt.title(str(iteration))
        # plt.subplot(1, 2, 2)
        # plt.plot(phiHydF, zgp1)
        # plt.gca().invert_yaxis()
        # plt.title(RMS)
        # plt.show()

    mask = np.sum(hFacC,axis=0)
    mask[mask>0]=1
    phi0surf = np.zeros((np.shape(hFacC)[1],np.shape(hFacC)[2]))
    for i in range(np.shape(hFacC)[2]):
        for j in range(np.shape(hFacC)[1]):
            k = np.argmin(np.abs(zg)<np.abs(icetopo[j,i]))
            kp1 = np.min([k+1,Nr])
            drloc = 1 - hFacC[k,j,i]
            dphi = phiHydF[kp1] - phiHydF[k]
            phi0surf[j,i] = (phiHydF[k]+drloc*dphi)*rhoConst*mask[j,i]

    return(phi0surf)

def create_icetopo_grid_from_millan(config_dir, model_name, N_Greenland_folder, X, Y, hFacC, Nr, dz, Depth, tRef, sRef,
                                    rhoConst, equation_of_state = 'jmd95'):

    extents = [np.min(X), np.min(Y), np.max(X), np.max(Y)]
    years = np.arange(1992, 2022)

    print('    - Reading in the Greene ice mask')
    mask_file = ''
    mask_file_subset = os.path.join(N_Greenland_folder,'Mask','greenland_ice_masks_1972-2022_v1_NEGIS_subset.nc')
    mask_x, mask_y, mask_years, mask_rock, mask_ice = read_mask_subset(mask_file, mask_file_subset, extents, years)
    mask_X, mask_Y = np.meshgrid(mask_x, mask_y)

    for year in years:

        index = np.argmin(np.abs(mask_years - year))
        ice_mask = mask_ice[index, :, :]

        # interpolate mask to the domain
        print('        - Interpolating mask on domain for time ' + str(int(mask_years[index])))
        ice_mask_on_domain = griddata(np.column_stack([mask_X.ravel(), mask_Y.ravel()]), ice_mask.ravel(), (X, Y),
                                      method='nearest')
        ice_mask_on_domain[Depth <= 0] = 0

        # plt.subplot(1,2,1)
        # C = plt.imshow(Depth,origin='lower')
        # plt.colorbar(C)
        # plt.subplot(1, 2, 2)
        # C = plt.imshow(ice_mask_on_domain, origin='lower')
        # plt.colorbar(C)
        # plt.show()

        # find connected regions on the actual shelves
        print('               - Computing the annual mask')
        annual_mask = compute_annual_mask(X, Y, ice_mask_on_domain, year)

        # connect h
        # C = plt.imshow(annual_mask, origin='lower')
        # plt.colorbar(C)
        # plt.show()

        print('               - Searching for miscellaneous side points')
        min_row = 310
        max_row = 560
        min_col = 40
        max_col = 290
        land = np.array(Depth <= 0)
        ocean = np.array(Depth > 0)

        land_and_ice = np.logical_or(land == 1, annual_mask == 1)
        # C = plt.imshow(land_and_ice, origin='lower')
        # plt.colorbar(C)
        # plt.show()
        for row in range(min_row, max_row):
            for col in range(min_col, max_col):
                if ocean[row, col] == 1:
                    if land_and_ice[row + 1, col] and land_and_ice[row - 1, col] and land_and_ice[row, col + 1] and \
                            land_and_ice[row, col - 1]:
                        if annual_mask[row + 1, col] == 1 or annual_mask[row - 1, col] == 1 or annual_mask[
                            row, col + 1] == 1 or annual_mask[row, col - 1] == 1:
                            annual_mask[row, col] = 1
                            # print('                    - filled row '+str(row)+' and col '+str(col))

        land_and_ice = np.logical_or(land == 1, annual_mask == 1)
        # C = plt.imshow(land_and_ice, origin='lower')
        # plt.colorbar(C)
        # plt.show()
        for row in range(min_row, max_row):
            for col in range(min_col, max_col):
                if ocean[row, col] == 1:
                    if land_and_ice[row + 1, col] and land_and_ice[row - 1, col] and land_and_ice[row, col + 1] and \
                            land_and_ice[row, col - 1]:
                        if annual_mask[row + 1, col] == 1 or annual_mask[row - 1, col] == 1 or annual_mask[
                            row, col + 1] == 1 or annual_mask[row, col - 1] == 1:
                            annual_mask[row, col] = 1
                            # print('                    - filled row '+str(row)+' and col '+str(col))

        # C = plt.imshow(annual_mask, origin='lower')
        # plt.colorbar(C)
        # plt.show()

        # compute thickness on the shelf for a given year
        thickness_points, thickness = read_ice_thickness_grids(N_Greenland_folder, year)
        annual_thickness = griddata(thickness_points, thickness.ravel(), (X, Y), method='nearest')
        annual_thickness[annual_mask == 0] = 0

        # need ice draft as a negative number for shelfice
        annual_thickness *= -1

        # C = plt.imshow(annual_thickness, origin='lower')
        # plt.colorbar(C)
        # plt.show()

        # if print_level >= 2:
        #     print('        - Computing the surface load field')
        phi0surf = create_surface_load_file(annual_thickness, hFacC, Nr, dz, tRef, sRef, rhoConst=rhoConst)

        # plt.subplot(2, 2, 1)
        # C = plt.imshow(Depth, origin='lower')
        # plt.colorbar(C)
        # plt.title('Depth')
        #
        # plt.subplot(2, 2, 2)
        # C = plt.imshow(annual_thickness, origin='lower')
        # plt.colorbar(C)
        # plt.title('Ice Topography')
        #
        # plt.subplot(2, 2, 3)
        # C = plt.imshow(annual_mask, origin='lower')
        # plt.colorbar(C)
        # plt.title('Mask')
        #
        # plt.subplot(2, 2, 4)
        # C = plt.imshow(phi0surf, origin='lower')
        # plt.colorbar(C)
        # plt.title('Phi0Surf')
        #
        # plt.show()

        shelfice_dir = os.path.join(config_dir, 'L2', model_name, 'input', 'shelfice')

        # output the grid to binary
        output_path = os.path.join(shelfice_dir, 'ice_thickness_' + str(year))
        annual_thickness.ravel('C').astype('>f4').tofile(output_path)

        output_file = os.path.join(shelfice_dir, 'ice_phi0surf_' + str(year))
        phi0surf.ravel('C').astype('>f4').tofile(output_file)



def create_shelfice_files(config_dir, model_name, print_level=1):

    level_name = 'L2'
    N_Greenland_folder = '/Users/mhwood/Documents/Research/Projects/North Greenland/Data'

    if print_level>=1:
        print('    - Generating the shelfice files for the ' + model_name + ' model from BedMachine')

    Lon_C, Lat_C, Depth, hFacC, drF = read_grid_geometry_from_nc(config_dir, model_name)
    Nr = len(drF)

    # plt.subplot(1,2,1)
    # C = plt.imshow(Lon_C)
    # plt.colorbar(C)
    # plt.subplot(1, 2, 2)
    # C = plt.imshow(Lat_C)
    # plt.colorbar(C)
    # plt.show()

    # x2, y2 = transformer.transform(polygon_array[:, y_column], polygon_array[:, x_column])
    transformer = Transformer.from_crs('EPSG:' + str(4326), 'EPSG:' + str(3413))
    XC, YC = transformer.transform(Lat_C.ravel(), Lon_C.ravel())

    XC = XC.reshape(np.shape(Lon_C))
    YC = YC.reshape(np.shape(Lat_C))

    ###############################################################

    tRef = np.array([18.89, 18.89, 18.89, 18.89, 18.89, 18.87,
                     18.85, 18.82, 18.80, 18.73, 18.65, 18.57,
                     18.40, 18.22, 18.00, 17.74, 17.44, 17.12,
                     16.76, 16.39, 15.98, 15.55, 15.08, 14.59,
                     14.07, 13.53, 12.99, 12.47, 11.97, 11.49,
                     11.02, 10.57, 10.12, 9.71, 9.27, 8.88,
                     8.46, 8.09, 7.71, 7.37, 7.03, 6.72,
                     6.42, 6.13, 5.86, 5.59, 5.34, 5.09,
                     4.87, 4.65, 4.45, 4.26, 4.08, 3.91,
                     3.75, 3.60, 3.47, 3.33, 3.20, 3.08,
                     2.96,  2.84,  2.73,  2.62,  2.51,  2.42, 2.32,
                     2.23,  2.14,  2.06,  1.98,  1.90,
                     1.81,  1.73,  1.65,  1.57,  1.49,  1.41,
                     1.33,  1.24,  1.15,  1.06,  0.98,  0.94,
                     0.91,  0.92,  0.98,  0.98,  0.98,  0.98])
    tRef = -1.9*np.ones((Nr,))

    sRef = np.array([34.84, 34.84, 34.84, 34.84, 34.84, 34.84,
                     34.85, 34.85, 34.85, 34.86, 34.87, 34.88,
                     34.89, 34.90, 34.92, 34.94, 34.96, 34.98,
                     35.00, 35.02, 35.04, 35.06, 35.07, 35.07,
                     35.07, 35.05, 35.03, 35.01, 34.98, 34.95,
                     34.92, 34.89, 34.85, 34.82, 34.79, 34.76,
                     34.73, 34.71, 34.68, 34.66, 34.64, 34.62,
                     34.61, 34.60, 34.59, 34.59, 34.58, 34.58,
                     34.59, 34.59, 34.60, 34.60, 34.61, 34.62,
                     34.63, 34.64, 34.65, 34.66, 34.67, 34.68,
                     34.69, 34.70, 34.71, 34.71, 34.72, 34.72, 34.73,
                     34.73, 34.74, 34.74, 34.74, 34.74,
                     34.75, 34.74, 34.74, 34.74, 34.74, 34.74,
                     34.74, 34.74, 34.73, 34.73, 34.73, 34.73,
                     34.73, 34.72, 34.72, 34.72, 34.72, 34.72])
    sRef = 34 * np.ones((Nr,))

    rhoConst = 1030

    create_icetopo_grid_from_millan(config_dir, model_name, N_Greenland_folder, XC, YC, hFacC,
                                    Nr, drF, Depth, tRef, sRef, rhoConst)


    # mask = np.zeros_like(ice_topo)
    # counter = 1
    # for j in range(np.shape(hFacC)[1]):
    #     for i in range(np.shape(hFacC)[2]):
    #         if ice_topo[j, i] != 0:
    #             mask[j, i] = counter
    #             counter += 1
    #         if phi0surf[j,i]>0:
    #             print(phi0surf[j,i], Depth[j,i], ice_topo[j,i])
    #
    # # output_file = '../input/shelf_mask.bin'
    # # mask.ravel('C').astype('>f4').tofile(output_file)
