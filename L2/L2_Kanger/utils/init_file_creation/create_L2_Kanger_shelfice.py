
import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
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

def read_ice_elevation_grid(config_dir, XC, YC):

    grid = np.fromfile(os.path.join(config_dir,'L2','L2_Kanger',
                                    'input','shelfice','melange_elevation.bin'), '>f4')
    grid = np.reshape(grid, np.shape(XC))

    return(grid)

def compute_ice_draft(elevation_grid, Depth):

    rho_w = 1024
    rho_ice = 917

    # (h_w + h_a)*rho_ice = h_w*rho_w
    # h_a*rho_ice = h_w*(rho_w-rho_ice)
    # h_a*rho_ice/(rho_w-rho_ice) = h_w

    ice_draft = elevation_grid*rho_ice/(rho_w-rho_ice)

    for i in range(np.shape(Depth)[0]):
        for j in range(np.shape(Depth)[1]):
            if elevation_grid[i,j]>0:
                if Depth[i,j]<ice_draft[i,j]:
                    ice_draft[i,j] = Depth[i,j]-20

    ice_draft *= -1

    return(ice_draft)

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

def plot_output_grids(Depth, hFacC, elevation_grid, ice_draft, surface_load, dv_mask):

    rows, cols = np.where(elevation_grid!=0)

    min_row = np.min(rows)-5
    max_row = np.max(rows)+6

    min_col = np.min(cols)-5
    max_col = np.max(cols)+5

    print(min_row, max_row)
    print(min_col, max_col)

    plt.subplot(2, 4, 1)
    C = plt.pcolormesh(hFacC[0, min_row:max_row, min_col:max_col])
    plt.colorbar(C)
    plt.title('hFacC')

    plt.subplot(2, 4, 2)
    C = plt.pcolormesh(hFacC[0, min_row:max_row, min_col:max_col])
    mask = elevation_grid[min_row:max_row, min_col:max_col] != 0
    mask = np.ma.masked_where(mask==0, mask)
    plt.pcolormesh(mask,cmap='Greys')
    plt.colorbar(C)
    plt.title('hFac with ice extent')

    plt.subplot(2, 4, 3)
    C = plt.pcolormesh(elevation_grid[min_row:max_row, min_col:max_col])
    plt.colorbar(C)
    plt.title('Ice Elevation')

    plt.subplot(2, 4, 4)
    C = plt.pcolormesh(ice_draft[min_row:max_row, min_col:max_col])
    plt.colorbar(C)
    plt.title('Ice Draft')

    plt.subplot(2, 4, 5)
    C = plt.pcolormesh(ice_draft[min_row:max_row, min_col:max_col] + Depth[min_row:max_row, min_col:max_col])
    plt.colorbar(C)
    plt.title('Bathymetry - Ice Draft')

    plt.subplot(2, 4, 6)
    C = plt.pcolormesh(surface_load[min_row:max_row, min_col:max_col])
    plt.colorbar(C)
    plt.title('Surface Load')

    plt.subplot(2, 4, 7)
    C = plt.pcolormesh(dv_mask[min_row:max_row, min_col:max_col])
    plt.colorbar(C)
    plt.title('dv mask')

    plt.show()

def create_shelfice_files(config_dir, model_name, print_level=1):

    if print_level>=1:
        print('    - Generating the shelfice files for the ' + model_name + ' model from BedMachine')

    XC, YC, Depth, hFacC, drF = read_grid_geometry_from_nc(config_dir, model_name)
    Nr = len(drF)

    # plt.subplot(1,2,1)
    # C = plt.imshow(Lon_C)
    # plt.colorbar(C)
    # plt.subplot(1, 2, 2)
    # C = plt.imshow(Lat_C)
    # plt.colorbar(C)
    # plt.show()


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

    elevation_grid = read_ice_elevation_grid(config_dir, XC, YC)
    elevation_grid[elevation_grid<0]=0

    elevation_grid[533+17,146+32] = elevation_grid[533+16,146+32]

    ice_draft = compute_ice_draft(elevation_grid, Depth)

    ice_draft_file = os.path.join(config_dir,'L2','L2_Kanger','input','shelfice','melange_draft.bin')
    ice_draft.ravel('C').astype('>f4').tofile(ice_draft_file)

    surface_load = create_surface_load_file(ice_draft, hFacC, Nr, drF, tRef, sRef, rhoConst)

    ice_load_file = os.path.join(config_dir, 'L2', 'L2_Kanger', 'input', 'shelfice', 'melange_load.bin')
    surface_load.ravel('C').astype('>f4').tofile(ice_load_file)

    dv_mask = np.zeros_like(ice_draft)
    counter = 1
    for i in range(np.shape(dv_mask)[0]):
        for j in range(np.shape(dv_mask)[1]):
            if ice_draft[i,j]!=0:
                dv_mask[i, j] = counter
                counter+=1
    print(counter)

    dv_mask_file = os.path.join(config_dir, 'L2', 'L2_Kanger', 'input', 'shelfice', 'melange_dv_mask.bin')
    dv_mask.ravel('C').astype('>f4').tofile(dv_mask_file)

    plot_output_grids(Depth, hFacC, elevation_grid, ice_draft, surface_load, dv_mask)
