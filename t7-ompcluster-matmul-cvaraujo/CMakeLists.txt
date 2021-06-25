cmake_minimum_required(VERSION 3.9)

project(matmul LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)


find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

option(USE_OMPC "Build for OmpCluster" ON)
# Ensure we are using clang
if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    if(USE_OMPC)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp-targets=x86_64-pc-linux-gnu")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-targets=x86_64-pc-linux-gnu")
    else()
        message(WARNING "OmpCluster offloading manually disable")
    endif()
else()
    message(WARNING "OmpCluster is not compatible with `${CMAKE_C_COMPILER_ID}`!")
endif()


add_executable(serial   src/matmul-serial.c)
add_executable(parallel src/matmul-parallel.c)
