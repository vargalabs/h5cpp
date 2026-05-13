#  _________________________________________________________
#  Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
#  Author: Varga, Steven <steven@vargaconsulting.ca>
#  _________________________________________________________

set(name xtensor)

set(libname lib${name})
project(${libname} LANGUAGES CXX)

set(current_version_dir "0.25.0")
add_library(${libname} INTERFACE)
add_library(xtensor::xtensor ALIAS ${libname})
target_include_directories(${libname} INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir}
)
target_link_libraries(${libname} INTERFACE libxtl)
