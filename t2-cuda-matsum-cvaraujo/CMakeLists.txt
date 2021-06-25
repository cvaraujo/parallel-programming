cmake_minimum_required(VERSION 3.9)

project(matsum
    LANGUAGES C CXX CUDA)

find_package(CUDA REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(serial   src/matsum-serial.c)
add_executable(parallel src/matsum-parallel.cu)

target_link_libraries(serial
    PUBLIC OpenMP::OpenMP_C)

target_link_libraries(parallel
    PUBLIC OpenMP::OpenMP_C)

target_compile_options(parallel
    PRIVATE "-arch=sm_35")

set_property(
    TARGET   parallel
    PROPERTY CUDA_ARCHITECTURES OFF)
