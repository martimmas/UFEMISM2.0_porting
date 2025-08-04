cmake_minimum_required(VERSION 3.10)
project( UPSY_models LANGUAGES Fortran)

set( CMAKE_Fortran_STANDARD 2018)

include( cmake/link_to_dependencies.cmake)
include( cmake/add_compiler_flags.cmake)

add_subdirectory(src/UPSY)
add_subdirectory(src/UFEMISM)
