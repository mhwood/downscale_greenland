cd ../
# mkdir run
rm -rf run/mnc_*
rm -r run/dv
rm run/*
cd run
cp -r ../input/dv .
ln -s ../input/* .
ln -s ../namelist/* .
ln -s ../build/mitgcmuv .
mpirun -np 6 ./mitgcmuv
