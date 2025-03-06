cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11")

project(vtk_laba)

add_executable(main src/main.cpp)


add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/usr/local/lib/libgmsh.so"
    INTERFACE_INCLUDE_DIRECTORIES "/usr/local/include/"
)
set(GMSH_LIBRARIES gmsh)


find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

target_link_libraries(main ${VTK_LIBRARIES})
target_link_libraries(main ${GMSH_LIBRARIES})
