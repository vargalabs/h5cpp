#
# \file CMakeLists.txt
# \brief main cmake configuration file
# \author Bogdan Cristea
#
# -------------------------------------------------------------------------
#
# Copyright (C) 1995-2012  (see AUTHORS file for a list of contributors)
#
# This file is part of IT++ - a C++ library of mathematical, signal
# processing, speech processing, and communications classes and functions.
#
# IT++ is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along
# with IT++.  If not, see <http://www.gnu.org/licenses/>.
#
# -------------------------------------------------------------------------

cmake_minimum_required ( VERSION 2.8.6 )

set ( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/ ${CMAKE_MODULE_PATH} )

project(ITPP)

#check for BLAS library
find_package ( BLAS )
if (BLAS_FOUND)
  if (BLA_VENDOR MATCHES "^ACML")
    set (HAVE_ACML 1)#TODO: for Windows only
    set (HAVE_BLAS_ACML 1)
    set (HAVE_BLAS 1)
  elseif (BLA_VENDOR MATCHES "^Intel")
    set (HAVE_MKL 1)#TODO: for Windows only
    set (HAVE_BLAS_MKL 1)
    set (HAVE_BLAS 1)
  elseif (BLA_VENDOR MATCHES "^ATLAS")
    set (HAVE_BLAS_ATLAS 1)
    set (HAVE_BLAS 1)
  elseif( (BLA_VENDOR STREQUAL "Generic") OR (BLA_VENDOR STREQUAL "All") )
    set (HAVE_BLAS 1)
  else()
    message (WARNING "Unknown BLAS vendor: ${BLA_VENDOR}")
    set(BLAS_FOUND FALSE)
  endif()
endif(BLAS_FOUND)

#check for LAPACK library
find_package ( LAPACK )
if (LAPACK_FOUND)
  if ((BLA_VENDOR MATCHES "^ACML") OR (BLA_VENDOR MATCHES "^Intel") OR (BLA_VENDOR MATCHES "^ATLAS") OR ( BLA_VENDOR STREQUAL "Generic" ) OR (BLA_VENDOR STREQUAL "All"))
    set (HAVE_LAPACK 1)
  else()
    message (WARNING "Unknown LAPACK vendor: ${BLA_VENDOR}")
    set(LAPACK_FOUND FALSE)
  endif()
endif(LAPACK_FOUND)

#check for FFT library
find_package ( FFT )
if (FFT_FOUND)
  if (FFT_VENDOR MATCHES "^ACML")
    set (HAVE_FFT_ACML 1)
    set (HAVE_FFT 1)
  elseif (FFT_VENDOR MATCHES "^Intel")
    set (HAVE_FFT_MKL 1)
    set (HAVE_FFT 1)
  elseif (FFT_VENDOR MATCHES "^FFTW3")
    set (HAVE_FFTW3 1)
    set (HAVE_FFT 1)
  else()
    message (WARNING "Unknown FFT vendor: ${FFT_VENDOR}")
    set(FFT_FOUND FALSE)
  endif()
endif(FFT_FOUND)

#check for OpenMP support
find_package ( OpenMP )
if (OPENMP_FOUND)
  set (CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
  if (MINGW)
    set (CMAKE_SHARED_LINKER_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS}")
  endif()
endif()

#check for Doxygen and LaTeX
find_package ( Doxygen )
find_package ( LATEX )

#check for header files
include ( CheckIncludeFile )
include ( CheckIncludeFiles )
CHECK_INCLUDE_FILES ( "sys/time.h;time.h" TIME_WITH_SYS_TIME )
CHECK_INCLUDE_FILE ( "dlfcn.h" HAVE_DLFCN_H )
CHECK_INCLUDE_FILE ( "ieeefp.h" HAVE_IEEEFP_H )
CHECK_INCLUDE_FILE ( "inttypes.h" HAVE_INTTYPES_H )
CHECK_INCLUDE_FILE ( "memory.h" HAVE_MEMORY_H )
CHECK_INCLUDE_FILE ( "stdint.h" HAVE_STDINT_H )
CHECK_INCLUDE_FILE ( "stdlib.h" HAVE_STDLIB_H )
CHECK_INCLUDE_FILE ( "strings.h" HAVE_STRINGS_H )
CHECK_INCLUDE_FILE ( "string.h" HAVE_STRING_H )
CHECK_INCLUDE_FILE ( "sys/stat.h" HAVE_SYS_STAT_H )
CHECK_INCLUDE_FILE ( "sys/types.h" HAVE_SYS_TYPES_H )
CHECK_INCLUDE_FILE ( "unistd.h" HAVE_UNISTD_H )
include ( CheckIncludeFileCXX )
CHECK_INCLUDE_FILE_CXX ( "cmath" HAVE_CMATH )
CHECK_INCLUDE_FILE_CXX ( "complex" HAVE_COMPLEX )

#check for symbols
include (CheckCXXSymbolExists)
CHECK_CXX_SYMBOL_EXISTS ( isfinite "cmath" HAVE_DECL_ISFINITE )
CHECK_CXX_SYMBOL_EXISTS ( isinf "cmath" HAVE_DECL_ISINF )
CHECK_CXX_SYMBOL_EXISTS ( isnan "cmath" HAVE_DECL_ISNAN )
CHECK_CXX_SYMBOL_EXISTS ( signgam "cmath" HAVE_DECL_SIGNGAM )

#check for functions
set(CMAKE_REQUIRED_LIBRARIES m)
include (CheckFunctionExists)
CHECK_FUNCTION_EXISTS ( acosh HAVE_ACOSH )
CHECK_FUNCTION_EXISTS ( asinh HAVE_ASINH )
CHECK_FUNCTION_EXISTS ( atanh HAVE_ATANH )
CHECK_FUNCTION_EXISTS ( cbrt HAVE_CBRT )
CHECK_FUNCTION_EXISTS ( erf HAVE_ERF )
CHECK_FUNCTION_EXISTS ( erfc HAVE_ERFC )
CHECK_FUNCTION_EXISTS ( expm1 HAVE_EXPM1 )
CHECK_FUNCTION_EXISTS ( finite HAVE_FPCLASS )
CHECK_FUNCTION_EXISTS ( fpclass HAVE_FPCLASS )
CHECK_FUNCTION_EXISTS ( isfinite HAVE_ISFINITE )
CHECK_FUNCTION_EXISTS ( isinf HAVE_ISINF )
CHECK_FUNCTION_EXISTS ( isnan HAVE_ISNAN )
CHECK_FUNCTION_EXISTS ( lgamma HAVE_LGAMMA )
CHECK_FUNCTION_EXISTS ( log1p HAVE_LOG1P )
CHECK_FUNCTION_EXISTS ( log2 HAVE_LOG2 )
CHECK_FUNCTION_EXISTS ( rint HAVE_RINT )
CHECK_FUNCTION_EXISTS ( tgamma HAVE_TGAMMA )
include (CheckCXXFunctionExists)
CHECK_CXX_FUNCTION_EXISTS ( std::isfinite HAVE_STD_ISFINITE )
CHECK_CXX_FUNCTION_EXISTS ( std::isinf HAVE_STD_ISINF )
CHECK_CXX_FUNCTION_EXISTS ( std::isnan HAVE_STD_ISNAN )

#package related variables
set (PACKAGE "itpp")
set (PACKAGE_BUGREPORT "http://sourceforge.net/projects/itpp/forums/forum/115656")
set (PACKAGE_NAME "IT++")
set (PACKAGE_VERSION "4.3.1")
set (PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set (PACKAGE_TARNAME ${PACKAGE})
set (PACKAGE_URL "")
if (UNIX)
  EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE PACKAGE_DATE)
  string(REGEX REPLACE "(..)/(..)/..(..).*" "\\1/\\2/\\3" PACKAGE_DATE "${PACKAGE_DATE}")
else()
  set(PACKAGE_DATE "")
endif()

if (WIN32 AND NOT MINGW)
  #TODO: remove this from future versions
  CONFIGURE_FILE ( ${CMAKE_SOURCE_DIR}/itpp/config.h.cmake ${CMAKE_BINARY_DIR}/itpp/config_msvc.h )
else()
  CONFIGURE_FILE ( ${CMAKE_SOURCE_DIR}/itpp/config.h.cmake ${CMAKE_BINARY_DIR}/itpp/config.h )
endif()

#set library type: shared or static
option(ITPP_SHARED_LIB "Building IT++ as shared library" on)
#set library name here to be available for all targets
if (ITPP_SHARED_LIB)
  if ((CMAKE_BUILD_TYPE STREQUAL Release) OR (NOT CMAKE_BUILD_TYPE))
    set (libitpp_target itpp)
    set (CMAKE_CXX_FLAGS "-DNDEBUG ${CMAKE_CXX_FLAGS}")
  else()
    set (libitpp_target itpp_debug)
  endif()
else()
  set (libitpp_target itpp_static)
endif()
CONFIGURE_FILE (${CMAKE_SOURCE_DIR}/itpp/itexports.h.cmake ${CMAKE_BINARY_DIR}/itpp/itexports.h)

#build HTML docs
option(HTML_DOCS "Building HTML documentation with Doxygen if available" on)

#display summary
set(DISP_GENERIC_WARNING FALSE)
if (NOT BLAS_FOUND)
  message(WARNING "BLAS library not found.")
  set(DISP_GENERIC_WARNING TRUE)
endif()
if (NOT LAPACK_FOUND)
  message(WARNING "LAPACK library not found.")
  set(DISP_GENERIC_WARNING TRUE)
endif()
if (NOT FFT_FOUND)
  message(WARNING "FFT library not found.")
  set(DISP_GENERIC_WARNING TRUE)
endif()
if (DISP_GENERIC_WARNING)
  message(WARNING "You can still compile IT++ but the functionality will be reduced.")
else()
  message(STATUS "BLAS vendor: ${BLA_VENDOR}")
endif()

#set external link libraries
set(ITPP_LIBS "")
if (BLAS_LIBRARIES)
  set(ITPP_LIBS ${BLAS_LIBRARIES} ${ITPP_LIBS})
endif()
if (LAPACK_LIBRARIES)
  set(ITPP_LIBS ${LAPACK_LIBRARIES} ${ITPP_LIBS})
endif()
if (FFT_LIBRARIES)
  set(ITPP_LIBS ${FFT_LIBRARIES} ${ITPP_LIBS})
endif()

#add it++ sources
add_subdirectory ( itpp )

#add unit tests
option(OLD_TESTS "Building old unit tests (no longer maintained). Use extras/check_tests.py script to run the tests." off)
if (OLD_TESTS)
  add_subdirectory ( tests )
endif()
if (GTEST_DIR)
  add_subdirectory (gtests)
endif()
