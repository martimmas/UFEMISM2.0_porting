# Collect UPSY/LADDIE/UFEMISM source files
file( GLOB_RECURSE UPSY_sources_raw    "${PROJECT_SOURCE_DIR}/src/UPSY/*.f90")

# Create static library for UPSY
add_library( UPSY_static_library STATIC ${UPSY_sources_raw})
link_to_dependencies( UPSY_static_library)
add_compiler_flags( UPSY_static_library)

target_include_directories( UPSY_static_library PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Initialise an empty list of executables
set( all_programs "")

# UPSY programs
add_executable( UPSY_unit_test_program      "${CMAKE_CURRENT_SOURCE_DIR}/UPSY_unit_test_program.f90")
add_executable( UPSY_component_test_program "${CMAKE_CURRENT_SOURCE_DIR}/UPSY_component_test_program.f90")
list( APPEND all_programs UPSY_unit_test_program)
list( APPEND all_programs UPSY_component_test_program)

file( GLOB_RECURSE LADDIE_sources_raw  "${PROJECT_SOURCE_DIR}/src/LADDIE/*.f90")
list(REMOVE_ITEM LADDIE_sources_raw "${PROJECT_SOURCE_DIR}/src/LADDIE/main/LADDIE_program.f90")

file( GLOB_RECURSE UFEMISM_sources_raw "${CMAKE_CURRENT_SOURCE_DIR}/*.f90")
list(REMOVE_ITEM UFEMISM_sources_raw "${CMAKE_CURRENT_SOURCE_DIR}/main/UFEMISM_program.f90")

# UFEMISM programs
add_executable( UFEMISM_program ${UPSY_sources_raw} ${LADDIE_sources_raw} ${UFEMISM_sources_raw} "${CMAKE_CURRENT_SOURCE_DIR}/main/UFEMISM_program.f90")

# Link all executables to the UPSY static library and the external dependencies,
# and set the compiler flags
foreach(prog IN LISTS all_programs)
    target_link_libraries(${prog} PUBLIC UPSY_static_library)
    link_to_dependencies(${prog})
    add_compiler_flags(${prog})
endforeach()
