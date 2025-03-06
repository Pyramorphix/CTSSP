add_executable(model model.cpp)


add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/usr/local/lib/libgmsh.so"
    INTERFACE_INCLUDE_DIRECTORIES "/usr/local/include"
)
set(GMSH_LIBRARIES gmsh)


target_link_libraries(model ${GMSH_LIBRARIES})
