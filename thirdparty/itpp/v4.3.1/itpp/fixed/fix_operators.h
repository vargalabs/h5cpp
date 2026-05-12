/*!
 * \file
 * \brief Definitions of a set of operators for Fix, Fixed, CFix and
 * CFixed classes
 * \author Johan Bergman
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

#ifndef FIX_OPERATORS_H
#define FIX_OPERATORS_H

#include <itpp/fixed/cfix.h>
#include <itpp/fixed/fix_functions.h>
#include <itpp/itexports.h>


namespace itpp
{

//! \addtogroup fixed
//!@{

/////////////////////////////////
// Operators for Fix and Fixed //
/////////////////////////////////

//! Fix + Fix
ITPP_EXPORT Fix operator+(const Fix &x, const Fix &y);
//! Fix - Fix
ITPP_EXPORT Fix operator-(const Fix &x, const Fix &y);
//! Fix * Fix
ITPP_EXPORT Fix operator*(const Fix &x, const Fix &y);
//! Fix / Fix using quantization mode \c TRN
ITPP_EXPORT Fix operator/(const Fix &x, const Fix &y);

//! Fix + int
ITPP_EXPORT Fix operator+(const Fix &x, const int y);
//! Fix - int
ITPP_EXPORT Fix operator-(const Fix &x, const int y);
//! Fix * int
ITPP_EXPORT Fix operator*(const Fix &x, const int y);
//! Fix / int using quantization mode \c TRN
ITPP_EXPORT Fix operator/(const Fix &x, const int y);
//! int + Fix
ITPP_EXPORT Fix operator+(const int x, const Fix &y);
//! int - Fix
ITPP_EXPORT Fix operator-(const int x, const Fix &y);
//! int * Fix
ITPP_EXPORT Fix operator*(const int x, const Fix &y);
//! int / Fix using quantization mode \c TRN
ITPP_EXPORT Fix operator/(const int x, const Fix &y);

//! fixvec + int
inline fixvec operator+(const fixvec &v, const int s) {return v + Fix(s);}
//! int + fixvec
inline fixvec operator+(const int s, const fixvec &v) {return Fix(s) + v;}
//! fixvec - int
inline fixvec operator-(const fixvec &v, const int s) {return v - Fix(s);}
//! int - fixvec
inline fixvec operator-(const int s, const fixvec &v) {return Fix(s) - v;}
//! fixvec * int
inline fixvec operator*(const fixvec &v, const int s) {return v * Fix(s);}
//! int * fixvec
inline fixvec operator*(const int s, const fixvec &v) {return Fix(s) * v;}
//! fixvec / int using quantization mode \c TRN
inline fixvec operator/(const fixvec &v, const int s) {return v / Fix(s);}

//! fixmat + int
inline fixmat operator+(const fixmat &v, const int s) {return v + Fix(s);}
//! int + fixmat
inline fixmat operator+(const int s, const fixmat &v) {return Fix(s) + v;}
//! fixmat - int
inline fixmat operator-(const fixmat &v, const int s) {return v - Fix(s);}
//! int - fixmat
inline fixmat operator-(const int s, const fixmat &v) {return Fix(s) - v;}
//! fixmat * int
inline fixmat operator*(const fixmat &v, const int s) {return v * Fix(s);}
//! int * fixmat
inline fixmat operator*(const int s, const fixmat &v) {return Fix(s) * v;}
//! fixmat / int using quantization mode \c TRN
inline fixmat operator/(const fixmat &v, const int s) {return v / Fix(s);}

//! fixvec + ivec
ITPP_EXPORT fixvec operator+(const fixvec &a, const ivec &b);
//! ivec + fixvec
inline fixvec operator+(const ivec &a, const fixvec &b) {return b + a;}
//! fixvec - ivec
inline fixvec operator-(const fixvec &a, const ivec &b) {return a + (-b);}
//! ivec - fixvec
inline fixvec operator-(const ivec &a, const fixvec &b) {return (-b) + a;}
//! fixvec * ivec
ITPP_EXPORT Fix operator*(const fixvec &a, const ivec &b);
//! ivec * fixvec
inline Fix operator*(const ivec &a, const fixvec &b) {return b*a;}

//! fixmat + imat
ITPP_EXPORT fixmat operator+(const fixmat &a, const imat &b);
//! imat + fixmat
inline fixmat operator+(const imat &a, const fixmat &b) {return b + a;}
//! fixmat - imat
inline fixmat operator-(const fixmat &a, const imat &b) {return a + (-b);}
//! imat - fixmat
inline fixmat operator-(const imat &a, const fixmat &b) {return (-b) + a;}
//! fixmat * imat
ITPP_EXPORT fixmat operator*(const fixmat &a, const imat &b);
//! imat * fixmat
inline fixmat operator*(const imat &a, const fixmat &b) {return b*a;}

///////////////////////////////////
// Operators for CFix and CFixed //
///////////////////////////////////

//! CFix + CFix
ITPP_EXPORT CFix operator+(const CFix &x, const CFix &y);
//! CFix - CFix
ITPP_EXPORT CFix operator-(const CFix &x, const CFix &y);
//! CFix * CFix
ITPP_EXPORT CFix operator*(const CFix &x, const CFix &y);
//! CFix / CFix using quantization mode \c TRN
ITPP_EXPORT CFix operator/(const CFix &x, const CFix &y);

//! CFix + Fix
ITPP_EXPORT CFix operator+(const CFix &x, const Fix &y);
//! CFix - Fix
ITPP_EXPORT CFix operator-(const CFix &x, const Fix &y);
//! CFix * Fix
ITPP_EXPORT CFix operator*(const CFix &x, const Fix &y);
//! CFix / Fix using quantization mode \c TRN
ITPP_EXPORT CFix operator/(const CFix &x, const Fix &y);
//! Fix + CFix
ITPP_EXPORT CFix operator+(const Fix &x, const CFix &y);
//! Fix - CFix
ITPP_EXPORT CFix operator-(const Fix &x, const CFix &y);
//! Fix * CFix
ITPP_EXPORT CFix operator*(const Fix &x, const CFix &y);
//! Fix / CFix using quantization mode \c TRN
ITPP_EXPORT CFix operator/(const Fix &x, const CFix &y);

//! CFix + int
ITPP_EXPORT CFix operator+(const CFix &x, const int y);
//! CFix - int
ITPP_EXPORT CFix operator-(const CFix &x, const int y);
//! CFix * int
ITPP_EXPORT CFix operator*(const CFix &x, const int y);
//! CFix / int using quantization mode \c TRN
ITPP_EXPORT CFix operator/(const CFix &x, const int y);
//! int + CFix
ITPP_EXPORT CFix operator+(const int x, const CFix &y);
//! int - CFix
ITPP_EXPORT CFix operator-(const int x, const CFix &y);
//! int * CFix
ITPP_EXPORT CFix operator*(const int x, const CFix &y);
//! int / CFix using quantization mode \c TRN
ITPP_EXPORT CFix operator/(const int x, const CFix &y);

//! fixvec + CFix
inline cfixvec operator+(const fixvec &v, const CFix &s) {return to<CFix>(v) + s;}
//! CFix + fixvec
inline cfixvec operator+(const CFix &s, const fixvec &v) {return s + to<CFix>(v);}
//! fixvec - CFix
inline cfixvec operator-(const fixvec &v, const CFix &s) {return to<CFix>(v) - s;}
//! CFix - fixvec
inline cfixvec operator-(const CFix &s, const fixvec &v) {return s - to<CFix>(v);}
//! fixvec * CFix
inline cfixvec operator*(const fixvec &v, const CFix &s) {return to<CFix>(v) * s;}
//! CFix * fixvec
inline cfixvec operator*(const CFix &s, const fixvec &v) {return s * to<CFix>(v);}
//! fixvec / CFix using quantization mode \c TRN
inline cfixvec operator/(const fixvec &v, const CFix &s) {return to<CFix>(v) / s;}

//! fixmat + CFix
inline cfixmat operator+(const fixmat &m, const CFix &s) {return to<CFix>(m) + s;}
//! CFix + fixmat
inline cfixmat operator+(const CFix &s, const fixmat &m) {return s + to<CFix>(m);}
//! fixmat - CFix
inline cfixmat operator-(const fixmat &m, const CFix &s) {return to<CFix>(m) - s;}
//! CFix - fixmat
inline cfixmat operator-(const CFix &s, const fixmat &m) {return s - to<CFix>(m);}
//! fixmat * CFix
inline cfixmat operator*(const fixmat &m, const CFix &s) {return to<CFix>(m) * s;}
//! CFix * fixmat
inline cfixmat operator*(const CFix &s, const fixmat &m) {return s * to<CFix>(m);}
//! fixmat / CFix using quantization mode \c TRN
inline cfixmat operator/(const fixmat &m, const CFix &s) {return to<CFix>(m) / s;}

//! ivec + CFix
inline cfixvec operator+(const ivec &v, const CFix &s) {return to<CFix>(to_vec(v)) + s;}
//! CFix + ivec
inline cfixvec operator+(const CFix &s, const ivec &v) {return s + to<CFix>(to_vec(v));}
//! ivec - CFix
inline cfixvec operator-(const ivec &v, const CFix &s) {return to<CFix>(to_vec(v)) - s;}
//! CFix - ivec
inline cfixvec operator-(const CFix &s, const ivec &v) {return s - to<CFix>(to_vec(v));}
//! ivec * CFix
inline cfixvec operator*(const ivec &v, const CFix &s) {return to<CFix>(to_vec(v)) * s;}
//! CFix * ivec
inline cfixvec operator*(const CFix &s, const ivec &v) {return s * to<CFix>(to_vec(v));}
//! ivec / CFix using quantization mode \c TRN
inline cfixvec operator/(const ivec &v, const CFix &s) {return to<CFix>(to_vec(v)) / s;}

//! imat + CFix
inline cfixmat operator+(const imat &m, const CFix &s) {return to<CFix>(to_mat(m)) + s;}
//! CFix + imat
inline cfixmat operator+(const CFix &s, const imat &m) {return s + to<CFix>(to_mat(m));}
//! imat - CFix
inline cfixmat operator-(const imat &m, const CFix &s) {return to<CFix>(to_mat(m)) - s;}
//! CFix - imat
inline cfixmat operator-(const CFix &s, const imat &m) {return s - to<CFix>(to_mat(m));}
//! imat * CFix
inline cfixmat operator*(const imat &m, const CFix &s) {return to<CFix>(to_mat(m)) * s;}
//! CFix * imat
inline cfixmat operator*(const CFix &s, const imat &m) {return s * to<CFix>(to_mat(m));}
//! imat / CFix using quantization mode \c TRN
inline cfixmat operator/(const imat &m, const CFix &s) {return to<CFix>(to_mat(m)) / s;}

//! cfixvec + Fix
inline cfixvec operator+(const cfixvec &v, const Fix &s) {return v + CFix(s);}
//! Fix + cfixvec
inline cfixvec operator+(const Fix &s, const cfixvec &v) {return CFix(s) + v;}
//! cfixvec - Fix
inline cfixvec operator-(const cfixvec &v, const Fix &s) {return v - CFix(s);}
//! Fix - cfixvec
inline cfixvec operator-(const Fix &s, const cfixvec &v) {return CFix(s) - v;}
//! cfixvec * Fix
inline cfixvec operator*(const cfixvec &v, const Fix &s) {return v * CFix(s);}
//! Fix * cfixvec
inline cfixvec operator*(const Fix &s, const cfixvec &v) {return CFix(s) * v;}
//! cfixvec / Fix using quantization mode \c TRN
inline cfixvec operator/(const cfixvec &v, const Fix &s) {return v / CFix(s);}

//! cfixmat + Fix
inline cfixmat operator+(const cfixmat &m, const Fix &s) {return m + CFix(s);}
//! Fix + cfixmat
inline cfixmat operator+(const Fix &s, const cfixmat &m) {return CFix(s) + m;}
//! cfixmat - Fix
inline cfixmat operator-(const cfixmat &m, const Fix &s) {return m - CFix(s);}
//! Fix - cfixmat
inline cfixmat operator-(const Fix &s, const cfixmat &m) {return CFix(s) - m;}
//! cfixmat * Fix
inline cfixmat operator*(const cfixmat &m, const Fix &s) {return m * CFix(s);}
//! Fix * cfixmat
inline cfixmat operator*(const Fix &s, const cfixmat &m) {return CFix(s) * m;}
//! cfixmat / Fix using quantization mode \c TRN
inline cfixmat operator/(const cfixmat &m, const Fix &s) {return m / CFix(s);}

//! cfixvec + int
inline cfixvec operator+(const cfixvec &v, const int s) {return v + CFix(s);}
//! int + cfixvec
inline cfixvec operator+(const int s, const cfixvec &v) {return CFix(s) + v;}
//! cfixvec - int
inline cfixvec operator-(const cfixvec &v, const int s) {return v - CFix(s);}
//! int - cfixvec
inline cfixvec operator-(const int s, const cfixvec &v) {return CFix(s) - v;}
//! cfixvec * int
inline cfixvec operator*(const cfixvec &v, const int s) {return v * CFix(s);}
//! int * cfixvec
inline cfixvec operator*(const int s, const cfixvec &v) {return CFix(s) * v;}
//! cfixvec / int using quantization mode \c TRN
inline cfixvec operator/(const cfixvec &v, const int s) {return v / CFix(s);}

//! cfixmat + int
inline cfixmat operator+(const cfixmat &m, const int s) {return m + CFix(s);}
//! int + cfixmat
inline cfixmat operator+(const int s, const cfixmat &m) {return CFix(s) + m;}
//! cfixmat - int
inline cfixmat operator-(const cfixmat &m, const int s) {return m - CFix(s);}
//! int - cfixmat
inline cfixmat operator-(const int s, const cfixmat &m) {return CFix(s) - m;}
//! cfixmat * int
inline cfixmat operator*(const cfixmat &m, const int s) {return m * CFix(s);}
//! int * cfixmat
inline cfixmat operator*(const int s, const cfixmat &m) {return CFix(s) * m;}
//! cfixmat / int using quantization mode \c TRN
inline cfixmat operator/(const cfixmat &m, const int s) {return m / CFix(s);}

//! cfixvec + fixvec
ITPP_EXPORT cfixvec operator+(const cfixvec &a, const fixvec &b);
//! fixvec + cfixvec
inline cfixvec operator+(const fixvec &a, const cfixvec &b) {return b + a;}
//! cfixvec - fixvec
inline cfixvec operator-(const cfixvec &a, const fixvec &b) {return a + (-b);}
//! fixvec - cfixvec
inline cfixvec operator-(const fixvec &a, const cfixvec &b) {return (-b) + a;}
//! cfixvec * fixvec
ITPP_EXPORT CFix operator*(const cfixvec &a, const fixvec &b);
//! fixvec * cfixvec
inline CFix operator*(const fixvec &a, const cfixvec &b) {return b*a;}

//! cfixmat + fixmat
ITPP_EXPORT cfixmat operator+(const cfixmat &a, const fixmat &b);
//! fixmat + cfixmat
inline cfixmat operator+(const fixmat &a, const cfixmat &b) {return b + a;}
//! cfixmat - fixmat
inline cfixmat operator-(const cfixmat &a, const fixmat &b) {return a + (-b);}
//! fixmat - cfixmat
inline cfixmat operator-(const fixmat &a, const cfixmat &b) {return (-b) + a;}
//! cfixmat * fixmat
ITPP_EXPORT cfixmat operator*(const cfixmat &a, const fixmat &b);
//! fixmat * cfixmat
inline cfixmat operator*(const fixmat &a, const cfixmat &b) {return b*a;}

//! cfixvec + ivec
ITPP_EXPORT cfixvec operator+(const cfixvec &a, const ivec &b);
//! ivec + cfixvec
inline cfixvec operator+(const ivec &a, const cfixvec &b) {return b + a;}
//! cfixvec - ivec
inline cfixvec operator-(const cfixvec &a, const ivec &b) {return a + (-b);}
//! ivec - cfixvec
inline cfixvec operator-(const ivec &a, const cfixvec &b) {return (-b) + a;}
//! cfixvec * ivec
ITPP_EXPORT CFix operator*(const cfixvec &a, const ivec &b);
//! ivec * cfixvec
inline CFix operator*(const ivec &a, const cfixvec &b) {return b*a;}

//! cfixmat + imat
ITPP_EXPORT cfixmat operator+(const cfixmat &a, const imat &b);
//! imat + cfixmat
inline cfixmat operator+(const imat &a, const cfixmat &b) {return b + a;}
//! cfixmat - imat
inline cfixmat operator-(const cfixmat &a, const imat &b) {return a + (-b);}
//! imat - cfixmat
inline cfixmat operator-(const imat &a, const cfixmat &b) {return (-b) + a;}
//! cfixmat * imat
ITPP_EXPORT cfixmat operator*(const cfixmat &a, const imat &b);
//! imat * cfixmat
inline cfixmat operator*(const imat &a, const cfixmat &b) {return b*a;}

//!@}

} // namespace itpp

#endif // #ifndef FIX_OPERATORS_H
