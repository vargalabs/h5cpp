#  _________________________________________________________
#  Copyright (c) 2018-2020 Steven Varga, Toronto,ON Canada
#  Author: Varga, Steven <steven@vargaconsulting.ca>
#  _________________________________________________________

set(name armadillo)

set(libname lib${name})
project(${libname} LANGUAGES CXX)

set(current_version_dir "v15.2.6")

# Armadillo is used in header-only mode (ARMA_DONT_USE_WRAPPER).
# The wrapper library (libarmadillo) is not needed for h5cpp's use case,
# which only requires basic linear algebra object introspection
# (memptr(), n_elem, n_cols, n_rows, constructors).
# Consumers must link BLAS/LAPACK directly when needed.
add_library(${libname} INTERFACE)
add_library(armadillo::armadillo ALIAS ${libname})
target_include_directories(${libname} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/${current_version_dir}/include)
target_compile_definitions(${libname} INTERFACE ARMA_DONT_USE_WRAPPER)
