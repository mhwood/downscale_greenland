# Building the L0 model

This configuration is nearly identical to that described in the [LLC270 experiment]([0/tides_exp](https://github.com/MITgcm-contrib/llc_hires/tree/master/llc_270)) - the only major change is that the [diagnostics_vec](https://github.com/mhwood/diagnostics_vec) package is added to output pertinent variables for the boundary of the L1 domain. 

Here, we outline the steps to run this model on the Pleiades computing cluster. If the model is to be run on different machines, small changes might be necessary. Note also that these instructions have a redundancy in that the model is compiled twice - the steps here are written as they were implemented for the successful run of the model. 

## Building the LLC270 ECCOv5 model
To start, follow all of the instructions in the LLC270 [readme](https://github.com/MITgcm-contrib/llc_hires/blob/master/llc_270/readme.txt) up until the job submission (last line with qsub). For convenience, the list of steps from this directory have been copied here (wuth comments removed):
```
git clone https://github.com/MITgcm/MITgcm.git
cd MITgcm
git checkout checkpoint64x
cd ..
svn checkout https://github.com/MITgcm-contrib/llc_hires/trunk/llc_270
# For the following requests you need your Earthdata username and WebDAV password (different from Earthdata password)
# Find it at :https://ecco.jpl.nasa.gov/drive
wget -r -nH -np --user=USERNAME --ask-password https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/input_forcing
wget -r -nH -np --user=USERNAME --ask-password https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/input_ecco
wget -r -nH -np --user=USERNAME --ask-password https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/input_init
wget -r -nH -np --user=USERNAME --ask-password https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/XX
mv drive/files/Version5/Alpha/input_forcing llc_270/
mv drive/files/Version5/Alpha/input_ecco    llc_270/
mv drive/files/Version5/Alpha/input_init    llc_270/
mv drive/files/Version5/Alpha/XX            llc_270/
rm -r drive/

cd MITgcm
mkdir build run
cd build
module purge
module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of ../../llc_270/code_ad/linux_amd64_ifort+mpi_ice_nas \
-mo ../../llc_270/code_ad
make depend
make -j 16
 
cd ../run
mkdir diags tapes
cp ../../llc_270/input_ad/* .
ln -s ../build/mitgcmuv .
ln -s ../../llc_270/input_forcing era_xx
ln -s ../../llc_270/input_ecco/* .
ln -s ../../llc_270/input_init/* .
ln -s ../../llc_270/XX/* .
 ```


## Adding diagnostics_vec to MITgcm
There are 3 main steps to add diagnostics_vec to the model

### Step 1: Add the `diagnostics_vec` package files to MITgcm
The diagnostics_vec package files can be easily added to MITgcm using a Python utilities function in the diagnostics_vec directory:
```
git clone https://github.com/mhwood/diagnostics_vec.git
cd ../../ # back to the head directory
cd diagnostics_vec/utils/
copy_pkg_files_to_MITgcm.py -m ../../MITgcm
cd ../../MITgcm/build
```

### Step 2: Add/Edit complile time files:
There are 7 files which need to be added or edited manually.

1. `PARAMS.h`: This file is inside the `model/inc` directory. Here, we will add two lines for the `diagnostics_vec` package. First define the `useDiagnostics_vec` as a LOGICAL and then put it in the name list. Its easiest to find where the usual `diagnostics` package is defined, and then add lines accordingly. For example:
```
1069      LOGICAL useDiagnostics
1070      LOGICAL useDiagnostics_vec                                  # new line added
1071      LOGICAL useREGRID
```
and
```
1088     &        useDiagnostics, useREGRID, useLayers, useMNC,
1089     &        useDiagnostics_vec,                                  # new line added
1090     &        useRunClock, useEMBED_FILES,
```

2. `packages_boot.F`: This file is inside the `model/src` directory. For this file, we will similarly add the `diagnostics_vec` package in the same location where the `diagnostics` package is include. This occurs in three places:
```
90     &          useDiagnostics,
91     &          useDiagnostics_vec,                  # added line
92     &          useREGRID,
```
and
```
156      useDiagnostics  =.FALSE.
157      useDiagnostics_vec  =.FALSE.                         # added line
158      useREGRID       =.FALSE.
```
and
```
390 #ifdef ALLOW_DIAGNOSTICS
391       CALL PACKAGES_PRINT_MSG( useDiagnostics,'Diagnostics', ' ' )
392 #endif
393 #ifdef ALLOW_DIAGNOSTICS_VEC                                                   # added line
394       CALL PACKAGES_PRINT_MSG( useDiagnostics_vec,                             # added line
395      &                         'Diagnostics_vec', ' ' )                        # added line
396 #endif                                                                         # added line
397 #ifdef ALLOW_REGRID
398       CALL PACKAGES_PRINT_MSG( useREGRID,     'REGRID',      ' ' )
399 #endif
```

3. `packages_init_fixed.F`: This file is inside the `${MOD540}/code` directory. For this file, we will add the `diagnostics_vec` package near the end of the file before the `ALLOW_DIAGNOSTIC` block:
```
660 #endif /* ALLOW_CTRL */
661 
662 #ifdef ALLOW_DIAGNOSTICS_VEC                                                   # added line
663       IF ( useDiagnostics_vec ) THEN                                           # added line
664 # ifdef ALLOW_DEBUG                                                            # added line
665         IF (debugMode)                                                         # added line
666      & CALL DEBUG_CALL('DIAGNOSTICS_VEC_INIT_FIXED',myThid)                    # added line
667 # endif                                                                        # added line
668         CALL DIAGNOSTICS_VEC_INIT_FIXED( myThid )                              # added line
669       ENDIF                                                                    # added line
670 #endif                                                                         # added line
671 
672 #ifdef ALLOW_DIAGNOSTICS
```

4. `packages_readparms.F`: This file is inside the `model/src` directory. For this file, we will add the `diagnostics_vec` package after the block for the `diagnostics` package:
```
363 #endif /* ALLOW_DIAGNOSTICS */
364 
365 #ifdef ALLOW_DIAGNOSTICS_VEC
366       CALL DIAGNOSTICS_VEC_READPARMS( myThid )
367 #endif /* ALLOW_DIAGNOSTICS_VEC */
```

5. `packages_check.F`: This file is inside the `model/src` directory. For this file, we will add the `diagnostics_vec` package after the block for the `diagnostics` package:
```
440 #ifdef ALLOW_DIAGNOSTICS_VEC
441       IF (useDiagnostics_vec) CALL DIAGNOSTICS_VEC_CHECK( myThid )
442 #else
443       IF (useDiagnostics_vec)
444      &   CALL PACKAGES_ERROR_MSG( 'Diagnostics_vec', ' ', myThid )
445 #endif
```

6. `do_the_model_io.F`: This file is inside the `model/src` directory. For this file, we will add the `diagnostics_vec` package before the block for the `diagnostics` package:
```
#ifdef ALLOW_DIAGNOSTICS_VEC
      IF ( useDiagnostics_vec )
     &     CALL DIAGNOSTICS_VEC_OUTPUT( myTime, myIter, myThid )
#endif
```

7. `packages.conf`: This file is inside the `code` directory. For this file, we simply add a line for `diagnostics_vec`.

8. `CPP_OPTIONS.h`: This file is inside the `code` directory. For this file, we will add three lines so that `diagnostics_vec` can access external forcing variables (not available when `ecco` is used for some reason I can't figure out...):
```
18 C-- Forcing code options:
19 
20 #define ALLOW_ATM_TEMP                                 # added line
21 #define ALLOW_DOWNWARD_RADIATION                       # added line
22 #define ALLOW_RUNOFF                                   # added line
```

9. `DIAGNOSTICS_VEC_SIZE.h`: This is a new file, provided in this repository, which we will copy to the `$code` directory.
Add this file from the `code` directory in this repository to the `code` directory of the model.

After the package is added and code modification files are edited, the model can be rebuilt using the same commands in LLC270 [readme](https://github.com/MITgcm-contrib/llc_hires/blob/master/llc_270/readme.txt) which are copied here for convenience:
```
 module purge
 module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
 ../tools/genmake2 -of ../../llc_270/code_ad/linux_amd64_ifort+mpi_ice_nas \
 -mo ../../llc_270/code_ad
 make depend
 make -j 16
 cd ..
 ```
 
 

## Run the model
Now that the model is built, it is nearly ready to run. To run this model, follow the steps provided in the [Step 1.2: Running the L0](https://github.com/mhwood/downscaled_east_pacific/blob/main/instructions/step_12_running_L0.md) instructions.
