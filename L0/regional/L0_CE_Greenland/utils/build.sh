mkdir ../build
cd ../build
../../../../../tools/genmake2 -mods ../code -mpi -optfile ../../../../../tools/build_options/darwin_amd64_gfortran_mw
make depend
make
cd ..