cd ../
mkdir run_for_grid
rm -r run_for_grid/mnc_0001
rm run_for_grid/*
cd run_for_grid
ln -s ../input_for_grid/* .
ln -s ../build_for_grid/mitgcmuv .
mpirun -np 1 ./mitgcmuv
