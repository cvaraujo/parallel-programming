#!/bin/sh -l

run_tests() {
  # Run serial
  build/serial tests/$1.in \
      1>serial.$1.out \
      2>serial.$1.time
  echo "Serial stdout:"
  cat serial.$1.out
  echo "Serial stderr:"
  cat serial.$1.time

  # Run parallel
  build/parallel tests/$1.in \
      1>parallel.$1.out \
      2>parallel.$1.time
  echo "Parallel stdout:"
  cat parallel.$1.out
  echo "Parallel stderr:"
  cat parallel.$1.time

  # Compare results
  diff -u serial.$1.out parallel.$1.out
}

# Compile the lab
export CC=clang
export CXX=clang++

mkdir -p build
cd build
cmake ..
make -j12

# Execute the Lab within OmpCluster container
cd ..
run_tests 1
run_tests 2
