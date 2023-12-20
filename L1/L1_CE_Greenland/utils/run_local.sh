cd ../
rm run_local/*
rm -r run_local/dv
rm -r run_local/diags
cd run_local
mkdir diags
mkdir diags/EtaN_day_snap
mkdir diags/EtaN_mon_mean
mkdir diags/SI_daily_snap
mkdir diags/TS_surf_daily_snap
mkdir diags/TS_AW_daily_snap
mkdir diags/state_3D_mon_snap
mkdir diags/vel_3D_mon_snap
mkdir diags/state_3D_mon_mean
mkdir diags/vel_3D_mon_mean
ln -s ../namelist/* .
cp -r ../input/dv .
ln -s ../input/* .
ln -s ../build/mitgcmuv .
mpirun -np 6 ./mitgcmuv
