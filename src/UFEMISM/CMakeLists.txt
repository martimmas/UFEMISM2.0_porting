# Collect LADDIE/UFEMISM source files
file( GLOB_RECURSE LADDIE_sources_raw  "${PROJECT_SOURCE_DIR}/src/LADDIE/*.f90")
list(REMOVE_ITEM LADDIE_sources_raw "${PROJECT_SOURCE_DIR}/src/LADDIE/main/LADDIE_program.f90")

file( GLOB_RECURSE UFEMISM_sources_raw "${CMAKE_CURRENT_SOURCE_DIR}/*.f90")
list(REMOVE_ITEM UFEMISM_sources_raw "${CMAKE_CURRENT_SOURCE_DIR}/main/UFEMISM_program.f90")

# UFEMISM program
add_executable( UFEMISM_program ${LADDIE_sources_raw} ${UFEMISM_sources_raw} "${CMAKE_CURRENT_SOURCE_DIR}/main/UFEMISM_program.f90")

target_link_libraries( UFEMISM_program PUBLIC UPSY_static_library)
link_to_dependencies(UFEMISM_program)
add_compiler_flags(UFEMISM_program)
