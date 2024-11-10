cmake_minimum_required(VERSION 3.0)

project(atomcomp)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_C_FLAGS "-O3")
set(CMAKE_CXX_FLAGS "-O3")

set(SRC_DIR  ${PROJECT_SOURCE_DIR})
set(SRC_FILES   ${SRC_DIR}/main.c
                ${SRC_DIR}/cutils.c
                ${SRC_DIR}/crand.cpp
                ${SRC_DIR}/logcomp.cpp
                ${SRC_DIR}/gencomp.c)

add_executable(atomcomp ${SRC_FILES})
