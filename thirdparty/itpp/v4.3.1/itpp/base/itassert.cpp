/*!
 * \file
 * \brief Error handling functions - source file
 * \author Tobias Ringstrom and Adam Piatyszek
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

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#include <itpp/base/itassert.h>
#include <iostream>
#include <stdexcept>
#include <cstdlib>


namespace itpp
{

//! \cond
static bool warnings_enabled = true;
static bool file_line_info_enabled = true;
static std::ostream *warn = &std::cerr;
//! \endcond

void it_assert_f(std::string ass, std::string msg, std::string file, int line)
{
  std::ostringstream error;
  if (file_line_info_enabled) {
    error << "*** Assertion failed in " << file << " on line " << line
    << ":\n" << msg << " (" << ass << ")\n";
  }
  else {
    error << msg << " (" << ass << ")\n";
  }
  std::cerr << error.str() << std::flush;
#ifdef ITPP_EXCEPTIONS
  throw std::runtime_error(error.str());
#else
  abort();
#endif
}

void it_error_f(std::string msg, std::string file, int line)
{
  std::ostringstream error;
  if (file_line_info_enabled) {
    error << "*** Error in " << file << " on line " << line << ":\n"
    << msg << "\n";
  }
  else {
    error << msg << "\n";
  }
  std::cerr << error.str() << std::flush;
#ifdef ITPP_EXCEPTIONS
  throw std::runtime_error(error.str());
#else
  abort();
#endif
}

void it_info_f(std::string msg)
{
  std::cerr << msg << std::flush;
}

void it_warning_f(std::string msg, std::string file, int line)
{
  if (warnings_enabled) {
    if (file_line_info_enabled) {
      (*warn) << "*** Warning in " << file << " on line " << line << ":\n"
      << msg << std::endl << std::flush;
    }
    else {
      (*warn) << msg << std::endl << std::flush;
    }
  }
}

void it_enable_warnings()
{
  warnings_enabled = true;
}

void it_disable_warnings()
{
  warnings_enabled = false;
}

void it_redirect_warnings(std::ostream *warn_stream)
{
  warn = warn_stream;
}

void it_error_msg_style(error_msg_style style)
{
  switch (style) {
  case Full:
    file_line_info_enabled = true;
    break;
  case Minimum:
    file_line_info_enabled = false;
    break;
  default:
    file_line_info_enabled = true;
  }
}

} //namespace itpp
