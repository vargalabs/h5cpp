#  _________________________________________________________
#  Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
#  Author: Varga, Steven <steven@vargaconsulting.ca>
#  _________________________________________________________

set(name ublas)

set(libname lib${name})
project(${libname} LANGUAGES CXX)

set(current_version_dir "v1.74.0")
add_library(${libname} INTERFACE)
target_include_directories(${libname} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir})
