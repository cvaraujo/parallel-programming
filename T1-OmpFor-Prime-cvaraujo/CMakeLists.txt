cmake_minimum_required(VERSION 3.9)

project(prime)

find_package(OpenMP REQUIRED)

set(CMAKE_C_STANDARD 99)

add_executable(serial   src/prime-serial.c)
add_executable(parallel src/prime-parallel.c)

target_link_libraries(serial
    PUBLIC OpenMP::OpenMP_C)

target_link_libraries(parallel
    PUBLIC OpenMP::OpenMP_C
    m)
