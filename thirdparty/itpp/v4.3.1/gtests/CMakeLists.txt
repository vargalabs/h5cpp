#
# \file CMakeLists.txt
# \brief cmake configuration file for gtests
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

file ( GLOB ITPP_TEST_SRCS *.cpp )

include_directories ( ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} ${GTEST_DIR} ${GTEST_DIR}/include )

add_executable ( itpp_gtests ${ITPP_TEST_SRCS} ${GTEST_DIR}/src/gtest-all.cc ${GTEST_DIR}/src/gtest_main.cc )
foreach (FILE ${ITPP_TEST_SRCS})
  get_filename_component ( target ${FILE} NAME_WE )
  if (target MATCHES "^fastica")
    set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS -DFASTICA_TEST_FILE=\\""${CMAKE_CURRENT_SOURCE_DIR}"/fastica_test_data.txt\\")
  endif()
  if (target MATCHES "^parser")
    set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS -DPARSER_TEST_FILE=\\""${CMAKE_CURRENT_SOURCE_DIR}"/parser_test_data.txt\\")
  endif()
  if (target MATCHES "^itfile")
    set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS -DITFILE_TEST_FILE=\\""${CMAKE_CURRENT_SOURCE_DIR}"/itfile_test_data.it\\")
  endif()
  if (target MATCHES "^g711")
    set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS -DG711_TEST_FILE=\\""${CMAKE_CURRENT_SOURCE_DIR}"/g711_ref_data.it\\")
  endif()
  if (target MATCHES "^audiofile")
    set_source_files_properties(${target} PROPERTIES COMPILE_FLAGS -DAUDIOFILE_REF_FILE=\\""${CMAKE_CURRENT_SOURCE_DIR}"/audiofile_ref_data.au\\")
  endif()

endforeach (FILE)

if (NOT WIN32)
  target_link_libraries ( itpp_gtests ${ITPP_LIBS} ${libitpp_target} pthread )
else()
  target_link_libraries ( itpp_gtests ${ITPP_LIBS} ${libitpp_target} )
endif()
