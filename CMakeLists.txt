cmake_minimum_required(VERSION 3.6)
project(cs267_hw2)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS -fopenmp)

set(COMMON_SOURCE
        common.cpp
        common.h)

set(SERIAL_SOURCE
        serial.cpp)

add_executable(serial ${COMMON_SOURCE} ${SERIAL_SOURCE})