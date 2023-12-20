cd ../
rm -r run_local/diags
rm run_local/*
cd run_local
mkdir diags
mkdir diags/EtaN_day_snap
mkdir diags/EtaN_mon_mean
mkdir diags/SI_daily_snap
mkdir diags/TS_surf_daily_snap
mkdir diags/vel_surf_daily_snap
mkdir diags/TS_AW_daily_snap
mkdir diags/state_3D_mon_snap
mkdir diags/vel_3D_mon_snap
mkdir diags/state_3D_mon_mean
mkdir diags/vel_3D_mon_mean
ln -s ../namelist/* .
ln -s ../input/* .
ln -s ../build/mitgcmuv .
mpirun -np 9 ./mitgcmuv
