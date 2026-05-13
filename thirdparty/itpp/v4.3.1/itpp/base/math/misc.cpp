/*!
 * \file
 * \brief Miscellaneous functions - source file
 * \author Tony Ottosson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include <itpp/base/math/misc.h>

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif


namespace itpp
{

std::string itpp_version(void)
{
#ifdef PACKAGE_VERSION
  return std::string(PACKAGE_VERSION);
#else
  return std::string("Warning: Version unknown!");
#endif
}

bool is_bigendian()
{
  int i = 1;
  char *p = reinterpret_cast<char *>(&i);
  if (p[0] == 1) // Lowest address contains the least significant byte
    return false; // LITTLE_ENDIAN
  else
    return true; // BIG_ENDIAN
}

} //namespace itpp
