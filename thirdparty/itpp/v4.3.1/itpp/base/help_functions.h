/*!
 * \file
 * \brief Help functions to make functions with vec and mat as arguments
 * \author Tony Ottosson, Adam Piatyszek, Andy Panov
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

#ifndef HELP_FUNCTIONS_H
#define HELP_FUNCTIONS_H

#include <itpp/base/mat.h>
#include <itpp/itexports.h>

namespace itpp
{
//! \cond
namespace details
{
//it would be nice to use binders and functors from <functional> header, but
//bind1st/bind2nd are broken since they can not be used with functions
//getting arguments by reference, so we introduce our own simple functors and
//binders here

//functor made from unary function Ret Ftn(Arg)
template<typename Ret, typename Arg>
class Unary_Function_Wrapper
{
  typedef Ret(*Ftn)(Arg);
  Ftn _f;
public:
  explicit Unary_Function_Wrapper(Ftn f): _f(f) {}
  Ret operator()(Arg v) const { return _f(v);}
};

//functor made from binary function Ret Ftn(Arg, Arg) with 1st arg bound
template<typename Ret, typename Arg>
class Binary_Function_With_Bound_1st
{
  typedef Ret(*Ftn)(Arg, Arg);
  Arg _b;
  Ftn _f;
public:
  Binary_Function_With_Bound_1st(Ftn f, const Arg b): _f(f), _b(b) {}
  Ret operator()(Arg v) const { return _f(_b, v);}
};

//functor made from binary function Ret Ftn(const Arg&, const Arg&) with 1st arg bound
//(partially specialize previous template)
template<typename Ret, typename Arg>
class Binary_Function_With_Bound_1st<Ret, const Arg&>
{
  typedef Ret(*Ftn)(const Arg&, const Arg&);
  const Arg& _b;
  Ftn _f;
public:
  Binary_Function_With_Bound_1st(Ftn f, const Arg& b): _f(f), _b(b) {}
  Ret operator()(const Arg& v) const { return _f(_b, v);}
};

//functor made from binary function Ret Ftn(Arg, Arg) with 2nd arg bound
template<typename Ret, typename Arg>
class Binary_Function_With_Bound_2nd
{
  typedef Ret(*Ftn)(Arg, Arg);
  Arg _b;
  Ftn _f;
public:
  Binary_Function_With_Bound_2nd(Ftn f, Arg b): _f(f), _b(b) {}
  Ret operator()(Arg v) const { return _f(v, _b);}
};

//functor made from binary function Ret Ftn(const Arg&, const Arg&) with 2nd arg bound
//(partially specialize previous template)
template<typename Ret, typename Arg>
class Binary_Function_With_Bound_2nd<Ret, const Arg&>
{
  typedef Ret(*Ftn)(const Arg&, const Arg&);
  const Arg& _b;
  Ftn _f;
public:
  Binary_Function_With_Bound_2nd(Ftn f, const Arg& b): _f(f), _b(b) {}
  Ret operator()(const Arg& v) const { return _f(v, _b);}
};
}
//! \endcond

//! \addtogroup miscfunc
//!@{

//! Help function to apply function object to Vec<T>
template<typename T, typename Ftn>
inline Vec<T> apply_functor(Ftn f, const Vec<T>& v)
{
  Vec<T> out(v.length());
  for(int i = 0; i < v.length(); i++) {
    out(i) = f(v(i));
  }
  return out;
}

//! Help function to call for a function: Vec<T> function(Vec<T>)
template<typename T>
inline Vec<T> apply_function(T(*f)(T), const Vec<T>& v)
{
  return apply_functor(details::Unary_Function_Wrapper<T, T>(f), v);
}

//! Help function to call for a function: Vec<T> function(const Vec<T>&)
template<typename T>
inline Vec<T> apply_function(T(*f)(const T&), const Vec<T>& v)
{
  return apply_functor(details::Unary_Function_Wrapper<T, const T&>(f), v);
}

//! Help function to apply function object to Mat<T>
template<typename T, typename Ftn>
inline Mat<T> apply_functor(Ftn f, const Mat<T>& m)
{
  Mat<T> out(m.rows(), m.cols());
  for(int i = 0; i < m.rows(); i++) {
    for(int j = 0; j < m.cols(); j++) {
      out(i, j) = f(m(i, j));
    }
  }
  return out;
}
//! Help function to call for a function: Mat<T> function(Mat<T>&)
template<typename T>
inline Mat<T> apply_function(T(*f)(T), const Mat<T>& m)
{
  return apply_functor(details::Unary_Function_Wrapper<T, T>(f), m);
}

//! Help function to call for a function: Mat<T> function(const Mat<T>&)
template<typename T>
inline Mat<T> apply_function(T(*f)(const T&), const Mat<T>& m)
{
  return apply_functor(details::Unary_Function_Wrapper<T, const T&>(f), m);
}

//! Help function to call for a function: Vec<T> function(T, Vec<T>)
template<typename T>
inline Vec<T> apply_function(T(*f)(T, T), const T& x, const Vec<T>& v)
{
  return apply_functor(details::Binary_Function_With_Bound_1st<T, T>(f, x), v);
}

//! \brief Help function to call for a function:
//!        Vec<T> function(const T&, const Vec<T>&)
template<typename T>
inline Vec<T> apply_function(T(*f)(const T&, const T&), const T& x,
                             const Vec<T>& v)
{
  return apply_functor(details::Binary_Function_With_Bound_1st<T, const T&>(f, x), v);
}

//! Help function to call for a function: Mat<T> function(T, Mat<T>)
template<typename T>
inline Mat<T> apply_function(T(*f)(T, T), const T& x, const Mat<T>& m)
{
  return apply_functor(details::Binary_Function_With_Bound_1st<T, T>(f, x), m);
}

//! \brief Help function to call for a function:
//!        Mat<T> function(const T&, const Mat<T>&)
template<typename T>
inline Mat<T> apply_function(T(*f)(const T&, const T&), const T& x,
                             const Mat<T>& m)
{
  return apply_functor(details::Binary_Function_With_Bound_1st<T, const T&>(f, x), m);
}

//! Help function to call for a function: Vec<T> function(Vec<T>, T)
template<typename T>
inline Vec<T> apply_function(T(*f)(T, T), const Vec<T>& v, const T& x)
{
  return apply_functor(details::Binary_Function_With_Bound_2nd<T, T>(f, x), v);
}

//! \brief Help function to call for a function:
//!        Vec<T> function(const Vec<T>&, const T&)
template<typename T>
inline Vec<T> apply_function(T(*f)(const T&, const T&), const Vec<T>& v,
                             const T& x)
{
  return apply_functor(details::Binary_Function_With_Bound_2nd<T, const T&>(f, x), v);
}

//! Help function to call for a function: Mat<T> function(Mat<T>, T)
template<typename T>
inline Mat<T> apply_function(T(*f)(T, T), const Mat<T>& m, const T& x)
{
  return apply_functor(details::Binary_Function_With_Bound_2nd<T, T>(f, x), m);
}

//! \brief Help function to call for a function:
//!        Mat<T> function(const Mat<T>&, const T&)
template<typename T>
inline Mat<T> apply_function(T(*f)(const T&, const T&), const Mat<T>& m,
                             const T& x)
{
  return apply_functor(details::Binary_Function_With_Bound_2nd<T, const T&>(f, x), m);
}

//!@}

//! \cond

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec apply_function(double(*f)(double), const vec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec apply_function(std::complex<double> (*f)(const std::complex<double> &),
                                    const cvec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT svec apply_function(short(*f)(short), const svec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec apply_function(int (*f)(int), const ivec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bvec apply_function(bin(*f)(bin), const bvec &v);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat apply_function(double(*f)(double), const mat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat apply_function(std::complex<double> (*f)(const std::complex<double> &),
                                    const cmat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT smat apply_function(short(*f)(short), const smat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT imat apply_function(int (*f)(int), const imat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bmat apply_function(bin(*f)(bin), const bmat &m);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec apply_function(double(*f)(double, double), const double& x, const vec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec apply_function(std::complex<double> (*f)(const std::complex<double> &,
                                    const std::complex<double> &),
                                    const std::complex<double>& x, const cvec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT svec apply_function(short(*f)(short, short), const short& x, const svec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec apply_function(int (*f)(int, int), const int& x, const ivec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bvec apply_function(bin(*f)(bin, bin), const bin& x, const bvec &v);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat apply_function(double(*f)(double, double), const double& x, const mat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat apply_function(std::complex<double> (*f)(const std::complex<double> &,
                                    const std::complex<double> &),
                                    const std::complex<double>& x, const cmat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT smat apply_function(short(*f)(short, short), const short& x, const smat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT imat apply_function(int (*f)(int, int), const int& x, const imat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bmat apply_function(bin(*f)(bin, bin), const bin& x, const bmat &m);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec apply_function(double(*f)(double, double), const vec &v, const double& x);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec apply_function(std::complex<double> (*f)(const std::complex<double> &,
                                    const std::complex<double> &),
                                    const cvec &v, const std::complex<double>& x);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT svec apply_function(short(*f)(short, short), const svec &v, const short& x);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec apply_function(int (*f)(int, int), const ivec &v, const int& x);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bvec apply_function(bin(*f)(bin, bin), const bvec &v, const bin& x);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat apply_function(double(*f)(double, double), const mat &m, const double& x);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat apply_function(std::complex<double> (*f)(const std::complex<double> &,
                                    const std::complex<double> &),
                                    const cmat &m, const std::complex<double>& x);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT smat apply_function(short(*f)(short, short), const smat &m, const short& x);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT imat apply_function(int (*f)(int, int), const imat &m, const int& x);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bmat apply_function(bin(*f)(bin, bin), const bmat &m, const bin& x);

//! \endcond

} // namespace itpp

#endif // #ifndef HELP_FUNCTIONS_H

