/*!
 * \file
 * \brief Logarithmic and exponenential functions - header file
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

#ifndef LOG_EXP_H
#define LOG_EXP_H

#include <cmath>
#include <itpp/base/help_functions.h>
#include <itpp/itexports.h>

namespace itpp
{

//!\addtogroup logexpfunc
//!@{

// ----------------------------------------------------------------------
// scalar functions
// ----------------------------------------------------------------------

//! Base-b logarithm
inline double logb(double b, double x)
{
  return (std::log(x) / std::log(b));
}

//! Calculate two to the power of x (2^x); x is integer
inline int pow2i(int x) { return ((x < 0) ? 0 : (1 << x)); }
//! Calculate two to the power of x (2^x)
inline double pow2(double x) { return std::pow(2.0, x); }
//! Calculate two to the power of x (2^x) for large integer x
inline double pow2(int x) 
{
#ifdef _MSC_VER
  //use pow since it seems to be faster on MSVC
  return std::pow(2.0, x); 
#else
  return std::ldexp(1.0,x);
#endif
}

//! Calculate ten to the power of x (10^x)
inline double pow10(double x) { return std::pow(10.0, x); }

//! Decibel of x (10*log10(x))
inline double dB(double x) { return 10.0 * std::log10(x); }
//! Inverse of decibel of x
inline double inv_dB(double x) { return std::pow(10.0, 0.1 * x); }

//! Calculate the number of bits needed to represent an integer \c n
inline int int2bits(int n)
{
  it_assert(n >= 0, "int2bits(): Improper argument value");

  if (n == 0)
    return 1;

  int b = 0;
  while (n) {
    n >>= 1;
    ++b;
  }
  return b;
}

//! Calculate the number of bits needed to represent \c n different values (levels).
inline int levels2bits(int n)
{
  it_assert(n > 0, "levels2bits(): Improper argument value");
  return int2bits(--n);
}

//! Constant definition to speed up trunc_log() and trunc_exp()
const double log_double_max = std::log(std::numeric_limits<double>::max());
//! Constant definition to speed up trunc_log(), trunc_exp() and log_add()
const double log_double_min = std::log(std::numeric_limits<double>::min());

/*!
  \brief Truncated natural logarithm function

  This truncated function provides a solution in the cases when the
  logarithm argument is less or equal to zero or infinity. The function
  checks for such extreme values and use some kind of truncation
  (saturation) before calculating the logarithm.

  The truncated logarithm function can be used for calculation of
  log-likelihood in soft demodulators, when numerical instability problem
  might occur.
*/
inline double trunc_log(double x)
{
  if (std::numeric_limits<double>::is_iec559) {
    if (x == std::numeric_limits<double>::infinity())
      return log_double_max;
    if (x <= 0)
      return log_double_min;
  }
  return std::log(x);
}

/*!
  \brief Truncated exponential function

  This truncated function provides a solution in the case when the
  exponent function results in infinity. The function checks for an
  extreme value and use truncation (saturation) before calculating
  the result.

  The truncated exponential function can be used  when numerical
  instability problem occurs for a standard exp function.
*/
inline double trunc_exp(double x)
{
  if (std::numeric_limits<double>::is_iec559
      && (x >= log_double_max))
    return std::numeric_limits<double>::max();
  return std::exp(x);
}


//! Safe substitute for <tt>log(exp(log_a) + exp(log_b))</tt>
ITPP_EXPORT double log_add(double log_a, double log_b);


// ----------------------------------------------------------------------
// functions on vectors and matrices
// ----------------------------------------------------------------------

//! Exp of the elements of a vector \c x
inline vec exp(const vec &x)
{
  return apply_function<double>(std::exp, x);
}
//! Exp of the elements of a complex vector \c x
inline cvec exp(const cvec &x)
{
  return apply_function<std::complex<double> >(std::exp, x);
}
//! Exp of the elements of a matrix \c m
inline mat exp(const mat &m)
{
  return apply_function<double>(std::exp, m);
}
//! Exp of the elements of a complex matrix \c m
inline cmat exp(const cmat &m)
{
  return apply_function<std::complex<double> >(std::exp, m);
}

//! Calculates x to the power of y (x^y)
inline vec pow(const double x, const vec &y)
{
  return apply_function<double>(std::pow, x, y);
}
//! Calculates x to the power of y (x^y)
inline mat pow(const double x, const mat &y)
{
  return apply_function<double>(std::pow, x, y);
}
//! Calculates x to the power of y (x^y)
inline vec pow(const vec &x, const double y)
{
  return apply_function<double>(std::pow, x, y);
}
//! Calculates x to the power of y (x^y)
inline mat pow(const mat &x, const double y)
{
  return apply_function<double>(std::pow, x, y);
}

//! Calculates two to the power of x (2^x)
inline vec pow2(const vec &x)
{
  return apply_function<double>(pow2, x);
}

//! Calculates two to the power of x (2^x) for integer x
inline vec pow2(const ivec &x)
{
  vec out(x.length());
  for(int i = 0; i < x.length(); i++)
    out(i) = pow2(x(i));
  return out;
}

//! Calculates two to the power of x (2^x)
inline mat pow2(const mat &x)
{
  return apply_function<double>(pow2, x);
}

//! Calculates two to the power of x (2^x) for integer x
inline mat pow2(const imat &x)
{
  mat out(x.rows(), x.cols());
  for(int i = 0; i < x.rows(); i++) {
    for(int j = 0; j < x.cols(); j++) {
      out(i, j) = pow2(x(i, j));
    }
  }
  return out;
}

//! Calculates ten to the power of x (10^x)
inline vec pow10(const vec &x)
{
  return apply_function<double>(pow10, x);
}
//! Calculates ten to the power of x (10^x)
inline mat pow10(const mat &x)
{
  return apply_function<double>(pow10, x);
}

//! The natural logarithm of the elements
inline vec log(const vec &x)
{
  return apply_function<double>(std::log, x);
}
//! The natural logarithm of the elements
inline mat log(const mat &x)
{
  return apply_function<double>(std::log, x);
}
//! The natural logarithm of the elements
inline cvec log(const cvec &x)
{
  return apply_function<std::complex<double> >(std::log, x);
}
//! The natural logarithm of the elements
inline cmat log(const cmat &x)
{
  return apply_function<std::complex<double> >(std::log, x);
}

// Cygwin defines log2 macro conflicting with IT++ functions
#if defined(log2)
#  undef log2
#endif
//! log-2 of the elements
ITPP_EXPORT vec log2(const vec &x);
//! log-2 of the elements
ITPP_EXPORT mat log2(const mat &x);

//! log-10 of the elements
inline vec log10(const vec &x)
{
  return apply_function<double>(std::log10, x);
}
//! log-10 of the elements
inline mat log10(const mat &x)
{
  return apply_function<double>(std::log10, x);
}

//! log-b of \c x
inline vec logb(double b, const vec &x)
{
  return apply_function<double>(itpp::logb, b, x);
}
//! log-b of \c x
inline mat logb(double b, const mat &x)
{
  return apply_function<double>(itpp::logb, b, x);
}

//! Calculates 10*log10(x)
inline vec dB(const vec &x)
{
  return apply_function<double>(dB, x);
}
//! Calculates 10*log10(x)
inline mat dB(const mat &x)
{
  return apply_function<double>(dB, x);
}

//! Calulates the inverse of dB, 10^(x/10)
inline vec inv_dB(const vec &x)
{
  return apply_function<double>(inv_dB, x);
}
//! Calculates the inverse of dB, 10^(x/10)
inline mat inv_dB(const mat &x)
{
  return apply_function<double>(inv_dB, x);
}

//! Calculate the number of bits needed to represent each integer in a vector
inline ivec int2bits(const ivec& v)
{
  return apply_function<int>(int2bits, v);
}

//! Calculate the number of bits needed to represent a numer of levels saved in a vector
inline ivec levels2bits(const ivec& v)
{
  return apply_function<int>(levels2bits, v);
}

//!@}

} // namespace itpp

#endif // #ifndef LOG_EXP_H




