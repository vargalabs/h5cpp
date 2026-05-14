#  _________________________________________________________
#  Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
#  Author: Varga, Steven <steven@vargaconsulting.ca>
#  _________________________________________________________

set(name xtensor-blas)

set(libname lib${name})
project(${libname} LANGUAGES CXX)

set(current_version_dir "0.23.0")
add_library(${libname} INTERFACE)
add_library(xtensor-blas::xtensor-blas ALIAS ${libname})
target_include_directories(${libname} INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir}
)
find_package(BLAS)
find_file(XFLENS_CXXBLAS_CXX
    NAMES "xflens/cxxblas/cxxblas.cxx"
    HINTS "${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir}"
    NO_DEFAULT_PATH
)
if(BLAS_FOUND AND XFLENS_CXXBLAS_CXX)
    set(XFLENS_FOUND TRUE CACHE INTERNAL "xflens bundled headers present")
    target_link_libraries(${libname} INTERFACE libxtensor ${BLAS_LIBRARIES})
elseif(BLAS_FOUND)
    message(STATUS "xtensor-blas: xflens headers not found — BLAS/LAPACK operations unavailable")
    set(XFLENS_FOUND FALSE CACHE INTERNAL "xflens bundled headers not found")
    target_link_libraries(${libname} INTERFACE libxtensor)
else()
    set(XFLENS_FOUND FALSE CACHE INTERNAL "xflens bundled headers not found")
    target_link_libraries(${libname} INTERFACE libxtensor)
endif()
