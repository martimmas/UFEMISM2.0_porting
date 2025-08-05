# Collect UPSY/LADDIE/UFEMISM source files
file( GLOB_RECURSE UPSY_sources_raw    "${PROJECT_SOURCE_DIR}/src/UPSY/*.f90")
list(REMOVE_ITEM UPSY_sources_raw "${PROJECT_SOURCE_DIR}/src/UPSY/UPSY_unit_test_program.f90")
list(REMOVE_ITEM UPSY_sources_raw "${PROJECT_SOURCE_DIR}/src/UPSY/UPSY_component_test_program.f90")

file( GLOB_RECURSE LADDIE_sources_raw  "${CMAKE_CURRENT_SOURCE_DIR}/*.f90")
list(REMOVE_ITEM LADDIE_sources_raw "${CMAKE_CURRENT_SOURCE_DIR}/main/LADDIE_program.f90")

file( GLOB_RECURSE UFEMISM_sources_raw "${PROJECT_SOURCE_DIR}/src/UFEMISM/*.f90")
list(REMOVE_ITEM UFEMISM_sources_raw "${PROJECT_SOURCE_DIR}/src/UFEMISM/main/UFEMISM_program.f90")

# LADDIE program
add_executable( LADDIE_program ${UPSY_sources_raw} ${LADDIE_sources_raw} ${UFEMISM_sources_raw} "${CMAKE_CURRENT_SOURCE_DIR}/main/LADDIE_program.f90")

link_to_dependencies(LADDIE_program)
add_compiler_flags(LADDIE_program)
