cmake_minimum_required(VERSION 3.9)
project(untitled C)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_C_STANDARD 99)



add_executable(untitled main.c)