cmake_minimum_required(VERSION 3.27)
project(simplePGA)

set(CMAKE_C_STANDARD 11)
#set(CMAKE_CXX_STANDARD 17)

set (SOURCES_FILES main.c
        src/genetic_alg_BE.c
        src/Graph.c
        src/communesFunctions.c
        src/genetic_alg_DVTC.c)
add_executable(simplePGA ${SOURCES_FILES})

