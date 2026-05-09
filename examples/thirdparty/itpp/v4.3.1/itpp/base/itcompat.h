/*!
 * \file
 * \brief IT++ compatibility types and functions
 * \author Adam Piatyszek
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

#ifndef ITCOMPAT_H
#define ITCOMPAT_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

//! \cond

// Microsoft Visual C++ underscore prefixed functions
#if defined(_MSC_VER)
#  include <cfloat>
#  define finite(x) _finite(x)
#  define isfinite(x) _finite(x)
#  define isnan(x) _isnan(x)
#  define fpclass(x) _fpclass(x)
#  define FP_NINF _FPCLASS_NINF
#  define FP_PINF _FPCLASS_PINF
#  define jn(a, b) _jn(a, b)
#  define yn(a, b) _yn(a, b)
#  define j0(a) _j0(a)
#  define j1(a) _j1(a)
#endif // defined(_MSC_VER)


// Solaris uses <ieeefp.h> for declaring isnan() and finite() functions
#if defined(HAVE_IEEEFP_H)
#  include <ieeefp.h>
#endif

// These definitions would collide with IT++ functions
#if defined(min)
#  undef min
#endif
#if defined(max)
#  undef max
#endif
#if defined(log2)
#  undef log2
#endif

namespace std
{

#ifndef HAVE_STD_ISINF
#if (HAVE_DECL_ISINF == 1) || defined(HAVE_ISINF)
inline int isinf(double x) { return ::isinf(x); }
#elif defined(FPCLASS)
inline int isinf(double x)
{
  if (::fpclass(a) == FP_NINF) return -1;
  else if (::fpclass(a) == FP_PINF) return 1;
  else return 0;
}
#else
inline int isinf(double x)
{
  if ((x == x) && ((x - x) != 0.0)) return (x < 0.0 ? -1 : 1);
  else return 0;
}
#endif // #if (HAVE_DECL_ISINF == 1) || defined(HAVE_ISINF)
#endif // #ifndef HAVE_STD_ISINF

#ifndef HAVE_STD_ISNAN
#if (HAVE_DECL_ISNAN == 1) || defined(HAVE_ISNAN)
inline int isnan(double x) { return ::isnan(x); }
#else
inline int isnan(double x) { return ((x != x) ? 1 : 0); }
#endif // #if (HAVE_DECL_ISNAN == 1) || defined(HAVE_ISNAN)
#endif // #ifndef HAVE_STD_ISNAN

#ifndef HAVE_STD_ISFINITE
#if (HAVE_DECL_ISFINITE == 1) || defined(HAVE_ISFINITE)
inline int isfinite(double x) { return ::isfinite(x); }
#elif defined(HAVE_FINITE)
inline int isfinite(double x) { return ::finite(x); }
#else
inline int isfinite(double x)
{
  return ((!std::isnan(x) && !std::isinf(x)) ? 1 : 0);
}
#endif // #if (HAVE_DECL_ISFINITE == 1) || defined(HAVE_ISFINITE)
#endif // #ifndef HAVE_STD_ISFINITE

} // namespace std


#ifndef HAVE_TGAMMA
//! True gamma function
double tgamma(double x);
#endif

#if !defined(HAVE_LGAMMA) || (HAVE_DECL_SIGNGAM != 1)
//! Lograrithm of an absolute gamma function
double lgamma(double x);
//! Global variable needed by \c lgamma function
extern int signgam;
#endif

#ifndef HAVE_CBRT
//! Cubic root
double cbrt(double x);
#endif


#ifndef HAVE_LOG1P
//! Lograrithm of an argument \c x plus one
inline double log1p(double x) { return std::log(1.0 + x); }
#endif

#ifndef HAVE_LOG2
//! Base-2 logarithm
inline double log2(double x)
{
  static const double one_over_log2 = 1.0 / std::log(2.0);
  return std::log(x) * one_over_log2;
}
#endif


#ifndef HAVE_EXPM1
//! C99 exponential minus one (exp(x) - 1.0)
double expm1(double x);
#endif // HAVE_EXPM1


#ifndef HAVE_ERFC
//! Complementary error function
double erfc(double x);
#endif

#ifndef HAVE_ERF
//! Error function
inline double erf(double x) { return (1.0 - ::erfc(x)); }
#endif


#ifndef HAVE_ASINH
//! Arcus sinhyp
double asinh(double x);
#endif

#ifndef HAVE_ACOSH
//! Arcus coshyp
double acosh(double x);
#endif

#ifndef HAVE_ATANH
//! Arcus tanhyp
double atanh(double x);
#endif


#ifndef HAVE_RINT
double rint(double x);
#endif


// Represent GCC version in a concise form
#define GCC_VERSION (__GNUC__ * 10000           \
                     + __GNUC_MINOR__ * 100     \
                     + __GNUC_PATCHLEVEL__)

//! \endcond

#endif // ITCOMPAT_H
