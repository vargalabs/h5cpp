#
# \file CMakeLists.txt
# \brief cmake configuration file for tutorial examples
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

file ( GLOB ITPP_EXAMPLES_SRCS *.cpp )

include_directories ( ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} )

foreach (FILE ${ITPP_EXAMPLES_SRCS})
	get_filename_component ( target ${FILE} NAME_WE )
	add_executable ( ${target} ${FILE} )
  target_link_libraries ( ${target} ${ITPP_LIBS} ${libitpp_target} )
endforeach (FILE)
