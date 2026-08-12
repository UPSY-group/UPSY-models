cmake_minimum_required(VERSION 3.10)
project( UPSY_models LANGUAGES Fortran C)

set( CMAKE_Fortran_STANDARD 2018)

option(BUILD_UFEMISM "Build the UFEMISM model" OFF)
option(BUILD_LADDIE "Build the LADDIE model" OFF)

include( cmake/link_to_dependencies.cmake)
include( cmake/add_compiler_flags.cmake)
include( cmake/fortran_comment_sanitizer.cmake)

add_subdirectory(src/UPSY)

if(BUILD_UFEMISM)
    add_subdirectory(src/UFEMISM)
endif(BUILD_UFEMISM)

if(BUILD_LADDIE)
    add_subdirectory(src/LADDIE)
endif(BUILD_LADDIE)
