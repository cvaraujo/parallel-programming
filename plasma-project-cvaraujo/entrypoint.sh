#!/bin/sh -l

# Compile PLASMA
export CC=clang
export CXX=clang++

mkdir build
cd build
cmake ..
make -j24

# Execute PLASMA within OmpCluster container
./plasmatest ssysv --dim=256 --nb=64 --tol=100
# TODO: The command is currently testing the execution of the dgeadd
# application. You need to change the name to test the application that was
# attributed to you.
