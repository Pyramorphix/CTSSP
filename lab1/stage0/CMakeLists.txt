add_executable(cube cube.cpp)
add_executable(circle circle.cpp)
add_executable(cylinder cylinder.cpp)


add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/usr/lib/libgmsh.so"
    INTERFACE_INCLUDE_DIRECTORIES "/usr/include"
)
set(GMSH_LIBRARIES gmsh)


target_link_libraries(cube ${GMSH_LIBRARIES})
target_link_libraries(circle ${GMSH_LIBRARIES})
target_link_libraries(cylinder ${GMSH_LIBRARIES})
