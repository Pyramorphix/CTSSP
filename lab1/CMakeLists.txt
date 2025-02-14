cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

project(gmsh_laba)

add_subdirectory(stage0)
add_subdirectory(stage1)
add_subdirectory(stage2)

