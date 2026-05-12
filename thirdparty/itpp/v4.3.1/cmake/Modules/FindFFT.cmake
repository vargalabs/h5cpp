# File:   FindFFT.cmake
# Brief:  Find the FFT includes and library
# Author: Bogdan Cristea and Jed Brown
#
# Usage: Use the following line in your CMakeLists file: find_package(FFT)
#
# This script tries to find the path to FFT libraries and header files. The
# following variables are set:
#  FFT_VENDOR       - Possible values: ACML, MKL, FFTW3
#  FFT_INCLUDES     - where to find header file
#  FFT_LIBRARIES    - List of libraries 
#  FFT_FOUND        - True if FFT found.
#
# The original version of this script can be found at:
# https://github.com/jedbrown/cmake-modules/blob/master/FindFFTW.cmake
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

#TODO: check for fft function in each library
macro(Check_fft_Libraries LIBRARIES _prefix _name _list)

set(_libraries_work TRUE)
set(${LIBRARIES})
set(_combined_name)
if (NOT _libdir)
  if (WIN32)
    set(_libdir $ENV{LIB})
  elseif (APPLE)
    set(_libdir /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 $ENV{DYLD_LIBRARY_PATH})
  else ()
    set(_libdir /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 $ENV{LD_LIBRARY_PATH})
  endif ()
endif ()

foreach(_library ${_list})
  set(_combined_name ${_combined_name}_${_library})

  if(_libraries_work)
    if (BLA_STATIC)
      if (WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
      endif ( WIN32 )
      if (APPLE)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
      else (APPLE)
        set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
      endif (APPLE)
    else (BLA_STATIC)
      if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
        # for ubuntu's libblas3gf and liblapack3gf packages
        set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .so.3gf)
      endif ()
    endif (BLA_STATIC)
    find_library(${_prefix}_${_library}_LIBRARY
      NAMES ${_library}
      PATHS ${_libdir}
      )
    mark_as_advanced(${_prefix}_${_library}_LIBRARY)
    if (${_prefix}_${_library}_LIBRARY)
      set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARY})
      set(_libraries_work ${${_prefix}_${_library}_LIBRARY})
    endif()
  endif(_libraries_work)
endforeach(_library ${_list})

if(_libraries_work)
  # Test this combination of libraries.
  if(UNIX AND BLA_STATIC)
    set(CMAKE_REQUIRED_LIBRARIES ${_flags} "-Wl,--start-group" ${${LIBRARIES}} ${_blas} "-Wl,--end-group" ${_threads})
  else(UNIX AND BLA_STATIC)
    set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_blas} ${_threads})
  endif(UNIX AND BLA_STATIC)
#  message("DEBUG: CMAKE_REQUIRED_LIBRARIES = ${CMAKE_REQUIRED_LIBRARIES}")
  if (NOT _LANGUAGES_ MATCHES Fortran)
    if (BLA_VENDOR MATCHES "Intel11*")
      #MKL 11
      check_function_exists("${_name}" ${_prefix}${_combined_name}_WORKS)
    else ()
      check_function_exists("${_name}_" ${_prefix}${_combined_name}_WORKS)
    endif()
  else (NOT _LANGUAGES_ MATCHES Fortran)
    check_fortran_function_exists(${_name} ${_prefix}${_combined_name}_WORKS)
  endif (NOT _LANGUAGES_ MATCHES Fortran)
  set(CMAKE_REQUIRED_LIBRARIES)
  mark_as_advanced(${_prefix}${_combined_name}_WORKS)
  set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
  #message("DEBUG: ${LIBRARIES} = ${${LIBRARIES}}")
endif(_libraries_work)

 if(_libraries_work)
   set(${LIBRARIES} ${${LIBRARIES}} ${_blas} ${_threads})
 else(_libraries_work)
    set(${LIBRARIES} FALSE)
 endif(_libraries_work)

endmacro(Check_fft_Libraries)

if (FFT_INCLUDES)
  # Already in cache, be silent
  set (FFT_FIND_QUIETLY TRUE)
else(FFT_INCLUDES)
  set (FFT_FOUND FALSE)
endif (FFT_INCLUDES)

#check for FFTW3
if (NOT FFT_FOUND)
  if (BLA_VENDOR MATCHES "^ACML")
    find_package(BLAS)
    if (BLAS_LIBRARIES)
      set (FFT_VENDOR "ACML")
      set (FFT_LIBRARIES ${BLAS_LIBRARIES})
      get_filename_component (PARENT_PATH ${BLAS_LIBRARIES} PATH )
      find_path (FFT_INCLUDES acml.h PATH "${PARENT_PATH}/../include")
      set (FFT_INCLUDES ${BLAS_INCLUDES})
    endif()
  elseif (BLA_VENDOR MATCHES "^Intel")
    find_package(BLAS)
    if (BLAS_LIBRARIES)
      set (FFT_VENDOR "Intel")
      set (FFT_LIBRARIES ${BLAS_LIBRARIES})
      foreach (lib ${BLAS_LIBRARIES})
        get_filename_component (PARENT_PATH ${lib} PATH )
        find_path (FFT_INCLUDES mkl_dfti.h PATH "${PARENT_PATH}/../../include")
        if (FFT_INCLUDES)
          break()
        endif()
      endforeach()
      if (FFT_INCLUDES)
        find_path(TEMP_INC fftw3.h PATH ${FFT_INCLUDES}/fftw NO_DEFAULT_PATH)
        if (TEMP_INC)
          set(FFT_INCLUDES ${TEMP_INC} ${FFT_INCLUDES})
        else()
          set(FFT_INCLUDES "")
        endif()
      endif()
    endif()
    mark_as_advanced (TEMP_INC)
    #need to find MKL FFT library
    check_fft_libraries(
        FFT_LIBRARIES
        FFT
        DftiCommitDescriptor
        "mkl_rt"
        )
  else()
    #check for generic FFT
    find_path (FFT_INCLUDES fftw3.h)
    find_library (FFT_LIBRARIES NAMES fftw3)
    if (FFT_INCLUDES AND FFT_LIBRARIES)
      set (FFT_VENDOR "FFTW3")
    endif()
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set FFT_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFT DEFAULT_MSG FFT_LIBRARIES FFT_INCLUDES)

if (FFT_FOUND)
  message(STATUS "A library with FFT API found.")
else()
  if(FFT_FIND_REQUIRED)
    message(ERROR_FATAL "A required library with FFT API not found. Please specify library location.")
  else()
    message(STATUS "A library with FFT API not found. Please specify library location.")
  endif()
endif()

if (FFT_LIBRARIES MATCHES "NOTFOUND$")
  set(FFT_LIBRARIES "")
endif()
if (FFT_INCLUDES MATCHES "NOTFOUND$")
  set(FFT_INCLUDES "")
endif()

mark_as_advanced (FFT_LIBRARIES FFT_INCLUDES)
