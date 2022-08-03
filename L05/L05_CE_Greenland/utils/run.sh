cd ../
#mkdir run
rm run/*
rm -r run/dv
rm -rf run/mnc*
cd run
ln -s ../namelist/* .
cp -r ../input/dv .
ln -s ../input/* .
ln -s ../build/mitgcmuv .
mpirun -np 6 ./mitgcmuv
