cd ../
mkdir run
rm -r run/mnc_*
rm run/*
cd run
ln -s ../namelist/* .
ln -s ../input/* .
ln -s ../build/mitgcmuv .
mpirun -np 6 ./mitgcmuv
