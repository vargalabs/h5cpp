#
# \file CMakeLists.txt
# \brief cmake configuration file for old unit tests
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

include_directories ( ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} )

foreach (FILE ${ITPP_TEST_SRCS})
	get_filename_component ( target ${FILE} NAME_WE )
	add_executable ( ${target} ${FILE} )
  if (target MATCHES "^fastica")
    set_target_properties(${target} PROPERTIES COMPILE_FLAGS -DFASTICA_TEST_FILE=\\""${CMAKE_CURRENT_SOURCE_DIR}"/fastica_test_data.txt\\")
  endif()
  if (target MATCHES "^parser")
    set_target_properties(${target} PROPERTIES COMPILE_FLAGS -DPARSER_TEST_FILE=\\""${CMAKE_CURRENT_SOURCE_DIR}"/parser_test_data.txt\\")
  endif()
  if (target MATCHES "^itfile")
    set_target_properties(${target} PROPERTIES COMPILE_FLAGS -DITFILE_TEST_FILE=\\""${CMAKE_CURRENT_SOURCE_DIR}"/itfile_test_data.it\\")
  endif()
  target_link_libraries ( ${target} ${ITPP_LIBS} ${libitpp_target} )
endforeach (FILE)
