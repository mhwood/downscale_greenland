cd ../
rm -r run_for_grid/mnc_*
rm run_for_grid/*
cd run_for_grid
ln -s ../input_for_grid/* .
ln -s ../namelist_for_grid/* .
ln -s ../build_for_grid/mitgcmuv .
mpirun --oversubscribe -np 16 ./mitgcmuv