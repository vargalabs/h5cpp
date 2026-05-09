/*!
 * \file
 * \brief Definitions of operators for vectors and matricies of different
 *        types
 * \author Tony Ottosson
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

#ifndef OPERATORS_H
#define OPERATORS_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/converters.h>
#include <itpp/itexports.h>

namespace itpp
{

//---------------------- between scalars and complex<double> -----------------
//! Addition operator for int and complex double
inline std::complex<double> operator+(const int &x, const std::complex<double> &y) {return std::complex<double>(x + y.real(), y.imag());}
//! Addition operator for float and complex double
inline std::complex<double> operator+(const float &x, const std::complex<double> &y) {return std::complex<double>(x + y.real(), y.imag());}
//! Addition operator for int and complex double
inline std::complex<double> operator+(const std::complex<double> &x, const int &y) {return std::complex<double>(x.real() + y, x.imag());}
//! Addition operator for float and complex double
inline std::complex<double> operator+(const std::complex<double> &x, const float &y) {return std::complex<double>(x.real() + y, x.imag());}

//! Subtraction operator for int and complex double
inline std::complex<double> operator-(const int &x, const std::complex<double> &y) {return std::complex<double>(x - y.real(), -y.imag());}
//! Subtraction operator for float and complex double
inline std::complex<double> operator-(const float &x, const std::complex<double> &y) {return std::complex<double>(x - y.real(), -y.imag());}
//! Subtraction operator for int and complex double
inline std::complex<double> operator-(const std::complex<double> &x, const int &y) {return std::complex<double>(x.real() - y, x.imag());}
//! Subtraction operator for float and complex double
inline std::complex<double> operator-(const std::complex<double> &x, const float &y) {return std::complex<double>(x.real() - y, x.imag());}

//! Multiplication operator for int and complex double
inline std::complex<double> operator*(const int &x, const std::complex<double> &y) {return std::complex<double>(x*y.real(), x*y.imag());}
//! Multiplication operator for float and complex double
inline std::complex<double> operator*(const float &x, const std::complex<double> &y) {return std::complex<double>(x*y.real(), x*y.imag());}
//! Multiplication operator for int and complex double
inline std::complex<double> operator*(const std::complex<double> &x, const int &y) {return std::complex<double>(x.real()*y, x.imag()*y);}
//! Multiplication operator for float and complex double
inline std::complex<double> operator*(const std::complex<double> &x, const float &y) {return std::complex<double>(x.real()*y, x.imag()*y);}

//! Division operator for complex double and int
inline std::complex<double> operator/(const std::complex<double> &x, const int &y) {return std::complex<double>(x.real() / y, x.imag() / y);}
//! Division operator for complex double and float
inline std::complex<double> operator/(const std::complex<double> &x, const float &y) {return std::complex<double>(x.real() / y, x.imag() / y);}


//---------------------- between vec and scalar --------------------

/*!
  \relatesalso Vec
  \brief Addition operator for float and vec
*/
inline vec operator+(const float &s, const vec &v) {return static_cast<double>(s) + v;}

/*!
  \relatesalso Vec
  \brief Addition operator for short and vec
*/
inline vec operator+(const short &s, const vec &v) {return static_cast<double>(s) + v;}

/*!
  \relatesalso Vec
  \brief Addition operator for int and vec
*/
inline vec operator+(const int &s, const vec &v) {return static_cast<double>(s) + v;}

/*!
  \relatesalso Vec
  \brief Addition operator for vec and float
*/
inline vec operator+(const vec &v, const float &s) {return static_cast<double>(s) + v;}

/*!
  \relatesalso Vec
  \brief Addition operator for vec and short
*/
inline vec operator+(const vec &v, const short &s) {return static_cast<double>(s) + v;}

/*!
  \relatesalso Vec
  \brief Addition operator for vec and int
*/
inline vec operator+(const vec &v, const int &s) {return static_cast<double>(s) + v;}

/*!
  \relatesalso Vec
  \brief Subtraction operator for float and vec
*/
inline vec operator-(const float &s, const vec &v) {return static_cast<double>(s) - v;}

/*!
  \relatesalso Vec
  \brief Subtraction operator for short and vec
*/
inline vec operator-(const short &s, const vec &v) {return static_cast<double>(s) - v;}

/*!
  \relatesalso Vec
  \brief Subtraction operator for int and vec
*/
inline vec operator-(const int &s, const vec &v) {return static_cast<double>(s) - v;}

/*!
  \relatesalso Vec
  \brief Subtraction operator for vec and float
*/
inline vec operator-(const vec &v, const float &s) {return v -static_cast<double>(s);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for vec and short
*/
inline vec operator-(const vec &v, const short &s) {return v -static_cast<double>(s);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for vec and int
*/
inline vec operator-(const vec &v, const int &s) {return v -static_cast<double>(s);}

/*!
  \relatesalso Vec
  \brief Multiplication operator for float and vec
*/
inline vec operator*(const float &s, const vec &v) {return static_cast<double>(s)*v;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for short and vec
*/
inline vec operator*(const short &s, const vec &v) {return static_cast<double>(s)*v;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for int and vec
*/
inline vec operator*(const int &s, const vec &v) {return static_cast<double>(s)*v;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for complex<double> and vec
*/
cvec operator*(const std::complex<double> &s, const vec &v);

/*!
  \relatesalso Vec
  \brief Multiplication operator for vec and float
*/
inline vec operator*(const vec &v, const float &s) {return static_cast<double>(s)*v;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for vec and short
*/
inline vec operator*(const vec &v, const short &s) {return static_cast<double>(s)*v;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for vec and int
*/
inline vec operator*(const vec &v, const int &s) {return static_cast<double>(s)*v;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for vec and complex<double>
*/
cvec operator*(const vec &v, const std::complex<double> &s);

/*!
  \relatesalso Vec
  \brief Division operator for vec and float
*/
inline vec operator/(const vec &v, const float &s) {return v / static_cast<double>(s);}

/*!
  \relatesalso Vec
  \brief Division operator for vec and short
*/
inline vec operator/(const vec &v, const short &s) {return v / static_cast<double>(s);}

/*!
  \relatesalso Vec
  \brief Division operator for vec and int
*/
inline vec operator/(const vec &v, const int &s) {return v / static_cast<double>(s);}


//---------------------- between ivec and scalar --------------------

/*!
  \relatesalso Vec
  \brief Addition operator for double and ivec
*/
ITPP_EXPORT vec operator+(const double &s, const ivec &v);

/*!
  \relatesalso Vec
  \brief Addition operator for ivec and double
*/
inline vec operator+(const ivec &v, const double &s) { return s + v;}

/*!
  \relatesalso Vec
  \brief Subtraction operator for double and ivec
*/
ITPP_EXPORT vec operator-(const double &s, const ivec &v);

/*!
  \relatesalso Vec
  \brief Subtraction operator for ivec and double
*/
inline vec operator-(const ivec &v, const double &s) { return v + (-s); }

/*!
  \relatesalso Vec
  \brief Multiplication operator for double and ivec
*/
ITPP_EXPORT vec operator*(const double &s, const ivec &v);

/*!
  \relatesalso Vec
  \brief Multiplication operator for ivec and double
*/
inline vec operator*(const ivec &v, const double &s) { return s*v; }

/*!
  \relatesalso Vec
  \brief Division operator for double and ivec
*/
ITPP_EXPORT vec operator/(const double &s, const ivec &v);

/*!
  \relatesalso Vec
  \brief Division operator for ivec and double
*/
ITPP_EXPORT vec operator/(const ivec &v, const double &s);

/*!
  \relatesalso Vec
  \brief Addition operator for complex<double> and ivec
*/
ITPP_EXPORT cvec operator+(const std::complex<double> &s, const ivec &v);

/*!
  \relatesalso Vec
  \brief Addition operator for ivec and complex<double>
*/
inline cvec operator+(const ivec &v, const std::complex<double> &s) { return s + v;}

/*!
  \relatesalso Vec
  \brief Subtraction operator for complex<double> and ivec
*/
ITPP_EXPORT cvec operator-(const std::complex<double> &s, const ivec &v);

/*!
  \relatesalso Vec
  \brief Subtraction operator for ivec and complex<double>
*/
inline cvec operator-(const ivec &v, const std::complex<double> &s) { return v + (-s); }

/*!
  \relatesalso Vec
  \brief Multiplication operator for complex<double> and ivec
*/
ITPP_EXPORT cvec operator*(const std::complex<double> &s, const ivec &v);

/*!
  \relatesalso Vec
  \brief Multiplication operator for ivec and complex<double>
*/
inline cvec operator*(const ivec &v, const std::complex<double> &s) { return s*v; }

/*!
  \relatesalso Vec
  \brief Division operator for complex<double> and ivec
*/
ITPP_EXPORT cvec operator/(const std::complex<double> &s, const ivec &v);

/*!
  \relatesalso Vec
  \brief Division operator for ivec and complex<double>
*/
ITPP_EXPORT cvec operator/(const ivec &v, const std::complex<double> &s);

//---------------------- between cvec and scalar --------------------

/*!
  \relatesalso Vec
  \brief Addition operator for double and cvec
*/
ITPP_EXPORT cvec operator+(const double &s, const cvec &v);

/*!
  \relatesalso Vec
  \brief Addition operator for float and cvec
*/
inline cvec operator+(const float &s, const cvec &v) {return static_cast<double>(s) + v;}

/*!
  \relatesalso Vec
  \brief Addition operator for short and cvec
*/
inline cvec operator+(const short &s, const cvec &v) {return static_cast<double>(s) + v;}

/*!
  \relatesalso Vec
  \brief Addition operator for int and cvec
*/
inline cvec operator+(const int &s, const cvec &v) {return static_cast<double>(s) + v;}

/*!
  \relatesalso Vec
  \brief Addition operator for cvec and float
*/
inline cvec operator+(const cvec &v, const float &s) {return s + v;}

/*!
  \relatesalso Vec
  \brief Addition operator for cvec and double
*/
inline cvec operator+(const cvec &v, const double &s) {return s + v;}

/*!
  \relatesalso Vec
  \brief Addition operator for cvec and short
*/
inline cvec operator+(const cvec &v, const short &s) {return s + v;}

/*!
  \relatesalso Vec
  \brief Addition operator for cvec and int
*/
inline cvec operator+(const cvec &v, const int &s) {return s + v;}

/*!
  \relatesalso Vec
  \brief Subtraction operator for double and cvec
*/
ITPP_EXPORT cvec operator-(const double &s, const cvec &v);

/*!
  \relatesalso Vec
  \brief Subtraction operator for float and cvec
*/
inline cvec operator-(const float &s, const cvec &v) {return static_cast<double>(s) - v;}

/*!
  \relatesalso Vec
  \brief Subtraction operator for short and cvec
*/
inline cvec operator-(const short &s, const cvec &v) {return static_cast<double>(s) - v;}

/*!
  \relatesalso Vec
  \brief Subtraction operator for int and cvec
*/
inline cvec operator-(const int &s, const cvec &v) {return static_cast<double>(s) - v;}

/*!
  \relatesalso Vec
  \brief Subtraction operator for cvec and float
*/
inline cvec operator-(const cvec &v, const float &s) {return v + (-s);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for cvec and double
*/
inline cvec operator-(const cvec &v, const double &s) {return v + (-s);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for cvec and short
*/
inline cvec operator-(const cvec &v, const short &s) {return v + (-s);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for cvec and int
*/
inline cvec operator-(const cvec &v, const int &s) {return v + (-s);}

/*!
  \relatesalso Vec
  \brief Multiplication operator for double and cvec
*/
ITPP_EXPORT cvec operator*(const double &s, const cvec &v);

/*!
  \relatesalso Vec
  \brief Multiplication operator for float and cvec
*/
inline cvec operator*(const float &s, const cvec &v) {return static_cast<double>(s)*v;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for short and cvec
*/
inline cvec operator*(const short &s, const cvec &v) {return static_cast<double>(s)*v;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for int and cvec
*/
inline cvec operator*(const int &s, const cvec &v) {return static_cast<double>(s)*v;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for cvec and float
*/
inline cvec operator*(const cvec &v, const float &s) {return s*v;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for cvec and double
*/
inline cvec operator*(const cvec &v, const double &s) {return s*v;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for cvec and short
*/
inline cvec operator*(const cvec &v, const short &s) {return s*v;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for cvec and int
*/
inline cvec operator*(const cvec &v, const int &s) {return s*v;}

/*!
  \relatesalso Vec
  \brief Division operator for cvec and double
*/
ITPP_EXPORT cvec operator/(const cvec &v, const double &s);

/*!
  \relatesalso Vec
  \brief Division operator for double and cvec
*/
ITPP_EXPORT cvec operator/(const double &s, const cvec &v);

/*!
  \relatesalso Vec
  \brief Division operator for cvec and float
*/
inline cvec operator/(const cvec &v, const float &s) {return v / static_cast<double>(s);}

/*!
  \relatesalso Vec
  \brief Division operator for cvec and short
*/
inline cvec operator/(const cvec &v, const short &s) {return v / static_cast<double>(s);}

/*!
  \relatesalso Vec
  \brief Division operator for cvec and int
*/
inline cvec operator/(const cvec &v, const int &s) {return v / static_cast<double>(s);}

//---------------------- between mat and scalar --------------------

/*!
  \relatesalso Mat
  \brief Addition operator for float and mat
*/
inline mat operator+(const float &s, const mat &m) {return static_cast<double>(s) + m;}

/*!
  \relatesalso Mat
  \brief Addition operator for short and mat
*/
inline mat operator+(const short &s, const mat &m) {return static_cast<double>(s) + m;}

/*!
  \relatesalso Mat
  \brief Addition operator for int and mat
*/
inline mat operator+(const int &s, const mat &m) {return static_cast<double>(s) + m;}

/*!
  \relatesalso Mat
  \brief Addition operator for mat and float
*/
inline mat operator+(const mat &m, const float &s) {return static_cast<double>(s) + m;}

/*!
  \relatesalso Mat
  \brief Addition operator for mat and short
*/
inline mat operator+(const mat &m, const short &s) {return static_cast<double>(s) + m;}

/*!
  \relatesalso Mat
  \brief Addition operator for mat and int
*/
inline mat operator+(const mat &m, const int &s) {return static_cast<double>(s) + m;}

/*!
  \relatesalso Mat
  \brief Subtraction operator for float and mat
*/
inline mat operator-(const float &s, const mat &m) {return static_cast<double>(s) - m;}

/*!
  \relatesalso Mat
  \brief Subtraction operator for short and mat
*/
inline mat operator-(const short &s, const mat &m) {return static_cast<double>(s) - m;}

/*!
  \relatesalso Mat
  \brief Subtraction operator for int and mat
*/
inline mat operator-(const int &s, const mat &m) {return static_cast<double>(s) - m;}

/*!
  \relatesalso Mat
  \brief Subtraction operator for mat and float
*/
inline mat operator-(const mat &m, const float &s) {return m -static_cast<double>(s);}

/*!
  \relatesalso Mat
  \brief Subtraction operator for mat and short
*/
inline mat operator-(const mat &m, const short &s) {return m -static_cast<double>(s);}

/*!
  \relatesalso Mat
  \brief Subtraction operator for mat and int
*/
inline mat operator-(const mat &m, const int &s) {return m -static_cast<double>(s);}

/*!
  \relatesalso Mat
  \brief Multiplication operator for float and mat
*/
inline mat operator*(const float &s, const mat &m) {return static_cast<double>(s)*m;}

/*!
  \relatesalso Mat
  \brief Multiplication operator for short and mat
*/
inline mat operator*(const short &s, const mat &m) {return static_cast<double>(s)*m;}

/*!
  \relatesalso Mat
  \brief Multiplication operator for int and mat
*/
inline mat operator*(const int &s, const mat &m) {return static_cast<double>(s)*m;}

/*!
  \relatesalso Mat
  \brief Multiplication operator for mat and float
*/
inline mat operator*(const mat &m, const float &s) {return static_cast<double>(s)*m;}

/*!
  \relatesalso Mat
  \brief Multiplication operator for mat and short
*/
inline mat operator*(const mat &m, const short &s) {return static_cast<double>(s)*m;}

/*!
  \relatesalso Mat
  \brief Multiplication operator for mat and int
*/
inline mat operator*(const mat &m, const int &s) {return static_cast<double>(s)*m;}

/*!
  \relatesalso Mat
  \brief Division operator for mat and float
*/
inline mat operator/(const mat &m, const float &s) {return m / static_cast<double>(s);}

/*!
  \relatesalso Mat
  \brief Division operator for mat and short
*/
inline mat operator/(const mat &m, const short &s) {return m / static_cast<double>(s);}

/*!
  \relatesalso Mat
  \brief Division operator for mat and int
*/
inline mat operator/(const mat &m, const int &s) {return m / static_cast<double>(s);}

//---------------------- between cmat and scalar --------------------

/*!
  \relatesalso Mat
  \brief Addition operator for double and cmat
*/
ITPP_EXPORT cmat operator+(const double &s, const cmat &m);

/*!
  \relatesalso Mat
  \brief Subtraction operator for double and cmat
*/
ITPP_EXPORT cmat operator-(const double &s, const cmat &m);

/*!
  \relatesalso Mat
  \brief Multiplication operator for double and cmat
*/
ITPP_EXPORT cmat operator*(const double &s, const cmat &m);

/*!
  \relatesalso Mat
  \brief Multiplication operator for complex<double> and mat
*/
ITPP_EXPORT cmat operator*(const std::complex<double> &s, const mat &m);

/*!
  \relatesalso Mat
  \brief Multiplication operator for mat and complex<double>
*/
inline cmat operator*(const mat &m, const std::complex<double> &s) {return s*m;}

/*!
  \relatesalso Mat
  \brief Division operator for cmat and double
*/
ITPP_EXPORT cmat operator/(const cmat &m, const double &s);

//---------------------- between vec and vectors --------------------

/*!
  \relatesalso Vec
  \brief Addition operator for bvec and vec
*/
ITPP_EXPORT vec operator+(const bvec &a, const vec &b);

/*!
  \relatesalso Vec
  \brief Addition operator for svec and vec
*/
ITPP_EXPORT vec operator+(const svec &a, const vec &b);

/*!
  \relatesalso Vec
  \brief Addition operator for ivec and vec
*/
ITPP_EXPORT vec operator+(const ivec &a, const vec &b);

/*!
  \relatesalso Vec
  \brief Addition operator for vec and bvec
*/
inline vec operator+(const vec &a, const bvec &b) {return b + a;}

/*!
  \relatesalso Vec
  \brief Addition operator for vec and svec
*/
inline vec operator+(const vec &a, const svec &b) {return b + a;}

/*!
  \relatesalso Vec
  \brief Addition operator for vec and ivec
*/
inline vec operator+(const vec &a, const ivec &b) {return b + a;}

/*!
  \relatesalso Vec
  \brief Subtraction operator for bvec and vec
*/
inline vec operator-(const bvec &a, const vec &b) {return a + (-b);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for svec and vec
*/
inline vec operator-(const svec &a, const vec &b) {return a + (-b);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for ivec and vec
*/
inline vec operator-(const ivec &a, const vec &b) {return a + (-b);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for vec and bvec
*/
inline vec operator-(const vec &a, const bvec &b) {return a + (-b);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for vec and svec
*/
inline vec operator-(const vec &a, const svec &b) {return a + (-b);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for vec and ivec
*/
inline vec operator-(const vec &a, const ivec &b) {return a + (-b);}

/*!
  \relatesalso Vec
  \brief Multiplication operator for bvec and vec
*/
ITPP_EXPORT double operator*(const bvec &a, const vec &b);

/*!
  \relatesalso Vec
  \brief Multiplication operator for svec and vec
*/
ITPP_EXPORT double operator*(const svec &a, const vec &b);

/*!
  \relatesalso Vec
  \brief Multiplication operator for ivec and vec
*/
ITPP_EXPORT double operator*(const ivec &a, const vec &b);

/*!
  \relatesalso Vec
  \brief Multiplication operator for vec and bvec
*/
inline double operator*(const vec &a, const bvec &b) {return b*a;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for vec and svec
*/
inline double operator*(const vec &a, const svec &b) {return b*a;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for vec and ivec
*/
inline double operator*(const vec &a, const ivec &b) {return b*a;}

//---------------------- between cvec and vectors --------------------

/*!
  \relatesalso Vec
  \brief Addition operator for bvec and cvec
*/
ITPP_EXPORT cvec operator+(const bvec &a, const cvec &b);

/*!
  \relatesalso Vec
  \brief Addition operator for svec and cvec
*/
ITPP_EXPORT cvec operator+(const svec &a, const cvec &b);

/*!
  \relatesalso Vec
  \brief Addition operator for ivec and cvec
*/
ITPP_EXPORT cvec operator+(const ivec &a, const cvec &b);

/*!
  \relatesalso Vec
  \brief Addition operator for cvec and bvec
*/
inline cvec operator+(const cvec &a, const bvec &b) {return b + a;}

/*!
  \relatesalso Vec
  \brief Addition operator for cvec and svec
*/
inline cvec operator+(const cvec &a, const svec &b) {return b + a;}

/*!
  \relatesalso Vec
  \brief Addition operator for cvec and ivec
*/
inline cvec operator+(const cvec &a, const ivec &b) {return b + a;}

/*!
  \relatesalso Vec
  \brief Subtraction operator for bvec and cvec
*/
inline cvec operator-(const bvec &a, const cvec &b) {return a + (-b);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for svec and cvec
*/
inline cvec operator-(const svec &a, const cvec &b) {return a + (-b);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for ivec and cvec
*/
inline cvec operator-(const ivec &a, const cvec &b) {return a + (-b);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for cvec and bvec
*/
inline cvec operator-(const cvec &a, const bvec &b) {return a + (-b);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for cvec and svec
*/
inline cvec operator-(const cvec &a, const svec &b) {return a + (-b);}

/*!
  \relatesalso Vec
  \brief Subtraction operator for cvec and ivec
*/
inline cvec operator-(const cvec &a, const ivec &b) {return a + (-b);}

/*!
  \relatesalso Vec
  \brief Multiplication operator for bvec and cvec
*/
ITPP_EXPORT std::complex<double> operator*(const bvec &a, const cvec &b);

/*!
  \relatesalso Vec
  \brief Multiplication operator for svec and cvec
*/
ITPP_EXPORT std::complex<double> operator*(const svec &a, const cvec &b);

/*!
  \relatesalso Vec
  \brief Multiplication operator for ivec and cvec
*/
ITPP_EXPORT std::complex<double> operator*(const ivec &a, const cvec &b);

/*!
  \relatesalso Vec
  \brief Multiplication operator for cvec and bvec
*/
inline std::complex<double> operator*(const cvec &a, const bvec &b) {return b*a;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for cvec and svec
*/
inline std::complex<double> operator*(const cvec &a, const svec &b) {return b*a;}

/*!
  \relatesalso Vec
  \brief Multiplication operator for cvec and ivec
*/
inline std::complex<double> operator*(const cvec &a, const ivec &b) {return b*a;}

//---------------------- between mat and matricies --------------------

/*!
  \relatesalso Mat
  \brief Addition operator for bmat and mat
*/
ITPP_EXPORT mat operator+(const bmat &a, const mat &b);

/*!
  \relatesalso Mat
  \brief Addition operator for smat and mat
*/
ITPP_EXPORT mat operator+(const smat &a, const mat &b);

/*!
  \relatesalso Mat
  \brief Addition operator for imat and mat
*/
ITPP_EXPORT mat operator+(const imat &a, const mat &b);

/*!
  \relatesalso Mat
  \brief Addition operator for mat and bmat
*/
inline mat operator+(const mat &a, const bmat &b) {return b + a;}

/*!
  \relatesalso Mat
  \brief Addition operator for mat and smat
*/
inline mat operator+(const mat &a, const smat &b) {return b + a;}

/*!
  \relatesalso Mat
  \brief Addition operator for mat and imat
*/
inline mat operator+(const mat &a, const imat &b) {return b + a;}

/*!
  \relatesalso Mat
  \brief Subtraction operator for bmat and mat
*/
inline mat operator-(const bmat &a, const mat &b) {return a + (-b);}

/*!
  \relatesalso Mat
  \brief Subtraction operator for smat and mat
*/
inline mat operator-(const smat &a, const mat &b) {return a + (-b);}

/*!
  \relatesalso Mat
  \brief Subtraction operator for imat and mat
*/
inline mat operator-(const imat &a, const mat &b) {return a + (-b);}

/*!
  \relatesalso Mat
  \brief Subtraction operator for mat and bmat
*/
inline mat operator-(const mat &a, const bmat &b) {return a + (-b);}

/*!
  \relatesalso Mat
  \brief Subtraction operator for mat and smat
*/
inline mat operator-(const mat &a, const smat &b) {return a + (-b);}

/*!
  \relatesalso Mat
  \brief Subtraction operator for mat and imat
*/
inline mat operator-(const mat &a, const imat &b) {return a + (-b);}

//---------------------- between cmat and matricies --------------------

/*!
  \relatesalso Mat
  \brief Addition operator for bmat and cmat
*/
ITPP_EXPORT cmat operator+(const bmat &a, const cmat &b);

/*!
  \relatesalso Mat
  \brief Addition operator for smat and cmat
*/
ITPP_EXPORT cmat operator+(const smat &a, const cmat &b);

/*!
  \relatesalso Mat
  \brief Addition operator for imat and cmat
*/
ITPP_EXPORT cmat operator+(const imat &a, const cmat &b);

/*!
  \relatesalso Mat
  \brief Addition operator for mat and cmat
*/
ITPP_EXPORT cmat operator+(const mat &a, const cmat &b);

/*!
  \relatesalso Mat
  \brief Addition operator for cmat and bmat
*/
inline cmat operator+(const cmat &a, const bmat &b) {return b + a;}

/*!
  \relatesalso Mat
  \brief Addition operator for cmat and smat
*/
inline cmat operator+(const cmat &a, const smat &b) {return b + a;}

/*!
  \relatesalso Mat
  \brief Addition operator for cmat and imat
*/
inline cmat operator+(const cmat &a, const imat &b) {return b + a;}

/*!
  \relatesalso Mat
  \brief Addition operator for cmat and mat
*/
inline cmat operator+(const cmat &a, const mat &b) {return b + a;}

/*!
  \relatesalso Mat
  \brief Subtraction operator for bmat and cmat
*/
inline cmat operator-(const bmat &a, const cmat &b) {return a + (-b);}

/*!
  \relatesalso Mat
  \brief Subtraction operator for smat and cmat
*/
inline cmat operator-(const smat &a, const cmat &b) {return a + (-b);}

/*!
  \relatesalso Mat
  \brief Subtraction operator for imat and cmat
*/
inline cmat operator-(const imat &a, const cmat &b) {return a + (-b);}

/*!
  \relatesalso Mat
  \brief Subtraction operator for mat and cmat
*/
inline cmat operator-(const mat &a, const cmat &b) {return a + (-b);}

/*!
  \relatesalso Mat
  \brief Subtraction operator for cmat and bmat
*/
inline cmat operator-(const cmat &a, const bmat &b) {return a + (-b);}

/*!
  \relatesalso Mat
  \brief Subtraction operator for cmat and smat
*/
inline cmat operator-(const cmat &a, const smat &b) {return a + (-b);}

/*!
  \relatesalso Mat
  \brief Subtraction operator for cmat and imat
*/
inline cmat operator-(const cmat &a, const imat &b) {return a + (-b);}

/*!
  \relatesalso Mat
  \brief Subtraction operator for cmat and mat
*/
inline cmat operator-(const cmat &a, const mat &b) {return a + (-b);}

/*!
  \relatesalso Mat
  \brief Multiplication operator for mat and cmat
*/
inline cmat operator*(const mat &a, const cmat &b) {return to_cmat(a)*b;}

/*!
  \relatesalso Mat
  \brief Multiplication operator for bmat and cmat
*/
inline cmat operator*(const bmat &a, const cmat &b) {return to_cmat(a)*b;}

/*!
  \relatesalso Mat
  \brief Multiplication operator for smat and cmat
*/
inline cmat operator*(const smat &a, const cmat &b) {return to_cmat(a)*b;}

/*!
  \relatesalso Mat
  \brief Multiplication operator for imat and cmat
*/
inline cmat operator*(const imat &a, const cmat &b) {return to_cmat(a)*b;}

/*!
  \relatesalso Mat
  \brief Multiplication operator for cmat and mat
*/
inline cmat operator*(const cmat &a, const mat &b) {return a*to_cmat(b);}

/*!
  \relatesalso Mat
  \brief Multiplication operator for cmat and bmat
*/
inline cmat operator*(const cmat &a, const bmat &b) {return a*to_cmat(b);}

/*!
  \relatesalso Mat
  \brief Multiplication operator for cmat and smat
*/
inline cmat operator*(const cmat &a, const smat &b) {return a*to_cmat(b);}

/*!
  \relatesalso Mat
  \brief Multiplication operator for cmat and imat
*/
inline cmat operator*(const cmat &a, const imat &b) {return a*to_cmat(b);}

} // namespace itpp

#endif // #ifndef OPERATORS_H

