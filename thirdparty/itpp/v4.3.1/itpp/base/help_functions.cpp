/*!
 * \file
 * \brief Help functions to make functions with vec and mat as arguments
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

#include <itpp/base/help_functions.h>

//! \cond

namespace itpp
{

template ITPP_EXPORT vec apply_function(double(*f)(double), const vec &v);
template ITPP_EXPORT cvec apply_function(std::complex<double> (*f)(const std::complex<double> &),
                             const cvec &v);
template ITPP_EXPORT svec apply_function(short(*f)(short), const svec &v);
template ITPP_EXPORT ivec apply_function(int (*f)(int), const ivec &v);
template ITPP_EXPORT bvec apply_function(bin(*f)(bin), const bvec &v);

template ITPP_EXPORT mat apply_function(double(*f)(double), const mat &m);
template ITPP_EXPORT cmat apply_function(std::complex<double> (*f)(const std::complex<double> &),
                             const cmat &m);
template ITPP_EXPORT smat apply_function(short(*f)(short), const smat &m);
template ITPP_EXPORT imat apply_function(int (*f)(int), const imat &m);
template ITPP_EXPORT bmat apply_function(bin(*f)(bin), const bmat &m);

template ITPP_EXPORT vec apply_function(double(*f)(double, double), const double& x, const vec &v);
template ITPP_EXPORT cvec apply_function(std::complex<double> (*f)(const std::complex<double> &,
                             const std::complex<double> &),
                             const std::complex<double>& x, const cvec &v);
template ITPP_EXPORT svec apply_function(short(*f)(short, short), const short& x, const svec &v);
template ITPP_EXPORT ivec apply_function(int (*f)(int, int), const int& x, const ivec &v);
template ITPP_EXPORT bvec apply_function(bin(*f)(bin, bin), const bin& x, const bvec &v);

template ITPP_EXPORT mat apply_function(double(*f)(double, double), const double& x, const mat &m);
template ITPP_EXPORT cmat apply_function(std::complex<double> (*f)(const std::complex<double> &,
                             const std::complex<double> &),
                             const std::complex<double>& x, const cmat &m);
template ITPP_EXPORT smat apply_function(short(*f)(short, short), const short& x, const smat &m);
template ITPP_EXPORT imat apply_function(int (*f)(int, int), const int& x, const imat &m);
template ITPP_EXPORT bmat apply_function(bin(*f)(bin, bin), const bin& x, const bmat &m);

template ITPP_EXPORT vec apply_function(double(*f)(double, double), const vec &v, const double& x);
template ITPP_EXPORT cvec apply_function(std::complex<double> (*f)(const std::complex<double> &,
                             const std::complex<double> &),
                             const cvec &v, const std::complex<double>& x);
template ITPP_EXPORT svec apply_function(short(*f)(short, short), const svec &v, const short& x);
template ITPP_EXPORT ivec apply_function(int (*f)(int, int), const ivec &v, const int& x);
template ITPP_EXPORT bvec apply_function(bin(*f)(bin, bin), const bvec &v, const bin& x);

template ITPP_EXPORT mat apply_function(double(*f)(double, double), const mat &m, const double& x);
template ITPP_EXPORT cmat apply_function(std::complex<double> (*f)(const std::complex<double> &,
                             const std::complex<double> &),
                             const cmat &m, const std::complex<double>& x);
template ITPP_EXPORT smat apply_function(short(*f)(short, short), const smat &m, const short& x);
template ITPP_EXPORT imat apply_function(int (*f)(int, int), const imat &m, const int& x);
template ITPP_EXPORT bmat apply_function(bin(*f)(bin, bin), const bmat &m, const bin& x);

} // namespace itpp

//! \endcond
