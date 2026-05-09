#  _________________________________________________________
#  Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
#  Author: Varga, Steven <steven@vargaconsulting.ca>
#  _________________________________________________________

set(name armadillo)

set(libname lib${name})
project(${libname} LANGUAGES CXX)

set(current_version_dir "v12.6.7")

# Save and restore BUILD_SHARED_LIBS to avoid affecting parent project
set(_saved_build_shared_libs ${BUILD_SHARED_LIBS})
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

add_subdirectory(${current_version_dir})

set(BUILD_SHARED_LIBS ${_saved_build_shared_libs} CACHE BOOL "" FORCE)

add_library(${libname} ALIAS armadillo)
