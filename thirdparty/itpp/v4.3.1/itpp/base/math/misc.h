/*!
 * \file
 * \brief Miscellaneous functions - header file
 * \author Tony Ottosson, Adam Piatyszek and Conrad Sanderson
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

#ifndef MISC_H
#define MISC_H

#include <complex>
#include <string>
#include <limits>
#include <itpp/itexports.h>

namespace std
{

//! Output stream operator for complex numbers
template <class T>
std::ostream& operator<<(std::ostream &os, const std::complex<T> &x)
{
  os << x.real();
  ios::fmtflags saved_format = os.setf(ios::showpos);
  os << x.imag();
  os.setf(saved_format, ios::showpos);
  return os << 'i';
}

//! Input stream operator for complex numbers
template <class T>
std::istream& operator>>(std::istream &is, std::complex<T> &x)
{
  T re, im;
  char c;
  is >> c;
  if (c == '(') {
    is >> re >> c;
    if (c == ',') {
      is >> im >> c;
      if (c == ')') {
        x = complex<T>(re, im);
      }
      else {
        is.setstate(ios_base::failbit);
      }
    }
    else if (c == ')') {
      x = complex<T>(re, T(0));
    }
    else {
      is.setstate(ios_base::failbit);
    }
  }
  else {
    is.putback(c);
    is >> re;
    if (!is.eof() && ((c = static_cast<char>(is.peek())) == '+' || c == '-')) {
      is >> im >> c;
      if (c == 'i') {
        x = complex<T>(re, im);
      }
      else {
        is.setstate(ios_base::failbit);
      }
    }
    else {
      x = complex<T>(re, T(0));
    }
  }
  return is;
}

} // namespace std


//! itpp namespace
namespace itpp
{

//! Constant Pi
const double pi = 3.14159265358979323846;

//! Constant 2*Pi
const double m_2pi = 2 * pi;

//! Constant eps
const double eps = std::numeric_limits<double>::epsilon();

//! \addtogroup miscfunc
//!@{

//! Return true if x is an integer
inline bool is_int(double x)
{
  double dummy;
  return (modf(x, &dummy) == 0.0);
}

//! Return true if x is an even integer
inline bool is_even(int x) { return ((x&1) == 0); }

//! Returns IT++ library version number, e.g. "3.7.1".
ITPP_EXPORT std::string itpp_version();

//! Returns true if machine endianness is BIG_ENDIAN
ITPP_EXPORT bool is_bigendian();

//! This function is deprecated. Please use is_bigendian() instead.
inline bool check_big_endianness() { return is_bigendian(); }

//!@}

} //namespace itpp


#endif // #ifndef MISC_H
