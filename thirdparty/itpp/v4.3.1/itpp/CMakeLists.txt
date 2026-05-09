#
# \file CMakeLists.txt
# \brief cmake configuration file for generating IT++ library
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

file ( GLOB ITPP_SRCS
       "base/*.cpp"
       "base/algebra/*.cpp"
       "base/bessel/*.cpp"
       "base/math/*.cpp"
       "comm/*.cpp" 
       "fixed/*.cpp"
       "optim/*.cpp"
       "protocol/*.cpp"
       "signal/*.cpp"
       "srccode/*.cpp"
       "stat/*.cpp" )

include_directories ( ${CMAKE_BINARY_DIR}
                      ${CMAKE_SOURCE_DIR}
                      ${BLAS_INCLUDES}
                      ${LAPACK_INCLUDES}
                      ${FFT_INCLUDES}
                    )

if (ITPP_SHARED_LIB)
  add_library (${libitpp_target} SHARED ${ITPP_SRCS})
else()
  add_library (${libitpp_target} STATIC ${ITPP_SRCS})
endif()

#remove duplicates in external libraries list
if (ITPP_LIBS)
  list(REMOVE_DUPLICATES ITPP_LIBS)
endif()

target_link_libraries (${libitpp_target} ${ITPP_LIBS})

set(LIBITPP_VERSION_MAJOR 8)
set(LIBITPP_VERSION_MINOR 2)
set(LIBITPP_VERSION_PATCH 1)
set(LIBITPP_VERSION_STRING ${LIBITPP_VERSION_MAJOR}.${LIBITPP_VERSION_MINOR}.${LIBITPP_VERSION_PATCH})

set_target_properties(${libitpp_target} PROPERTIES VERSION ${LIBITPP_VERSION_STRING}
                                                   SOVERSION ${LIBITPP_VERSION_MAJOR})

#library
install(TARGETS ${libitpp_target}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

#headers
install(DIRECTORY ${CMAKE_SOURCE_DIR}/itpp DESTINATION include FILES_MATCHING PATTERN "*.h")
if (WIN32 AND NOT MINGW)
  install(FILES ${CMAKE_BINARY_DIR}/itpp/config_msvc.h DESTINATION include/itpp)
else()
  install(FILES ${CMAKE_BINARY_DIR}/itpp/config.h DESTINATION include/itpp)
endif()
install(FILES ${CMAKE_BINARY_DIR}/itpp/itexports.h DESTINATION include/itpp)

#extra files (MATLAB and python)
install(FILES ${CMAKE_SOURCE_DIR}/extras/itsave.m
              ${CMAKE_SOURCE_DIR}/extras/itload.m
              ${CMAKE_SOURCE_DIR}/extras/pyitpp.py
        DESTINATION share/itpp)

#html help
if (HTML_DOCS AND DOXYGEN_FOUND)
  configure_file(${CMAKE_SOURCE_DIR}/doc/doxygen_html.cfg.cmake.in ${CMAKE_BINARY_DIR}/Doxyfile @ONLY)
  configure_file(${CMAKE_SOURCE_DIR}/doc/local/index.doc.in ${CMAKE_BINARY_DIR}/doc/local/index.doc @ONLY)
  #make sure that this target is executed last
  add_custom_target(itpp_doc ALL
    ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
  file(COPY ${CMAKE_SOURCE_DIR}/doc/images/itpp_logo.png DESTINATION ${CMAKE_BINARY_DIR}/html)
  install(DIRECTORY ${CMAKE_BINARY_DIR}/html DESTINATION share/doc/itpp)
endif()

#itpp-config script for UNIX-like systems
if (UNIX)
  if (ITPP_LIBS)
    string (REPLACE ";" " " ITPP_LIBS "${ITPP_LIBS}")
  endif()
  configure_file(${CMAKE_SOURCE_DIR}/itpp-config.cmake.in ${CMAKE_BINARY_DIR}/itpp-config @ONLY)
  configure_file(${CMAKE_SOURCE_DIR}/itpp-config.1.cmake.in ${CMAKE_BINARY_DIR}/itpp-config.1 @ONLY)
  configure_file(${CMAKE_SOURCE_DIR}/itpp.pc.cmake.in ${CMAKE_BINARY_DIR}/itpp.pc @ONLY)
  install(FILES ${CMAKE_BINARY_DIR}/itpp-config DESTINATION bin PERMISSIONS OWNER_READ OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  install(FILES ${CMAKE_BINARY_DIR}/itpp-config.1 DESTINATION share/man/man1)
  install(FILES ${CMAKE_BINARY_DIR}/itpp.pc DESTINATION lib/pkgconfig)
endif()
