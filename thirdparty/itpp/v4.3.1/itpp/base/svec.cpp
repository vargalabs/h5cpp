/*!
 * \file
 * \brief Sparse Vector Class implementation
 * \author Tony Ottosson and Tobias Ringstrom
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

#include <itpp/base/svec.h>

//! \cond

namespace itpp
{

// ---------------------------------------------------------------------
// Instantiations
// ---------------------------------------------------------------------

template class ITPP_EXPORT Sparse_Vec<int>;
template class ITPP_EXPORT Sparse_Vec<double>;
template class ITPP_EXPORT Sparse_Vec<std::complex<double> >;

template ITPP_EXPORT sparse_ivec operator+(const sparse_ivec &, const sparse_ivec &);
template ITPP_EXPORT sparse_vec operator+(const sparse_vec &, const sparse_vec &);
template ITPP_EXPORT sparse_cvec operator+(const sparse_cvec &, const sparse_cvec &);

template ITPP_EXPORT int operator*(const sparse_ivec &, const sparse_ivec &);
template ITPP_EXPORT double operator*(const sparse_vec &, const sparse_vec &);
template ITPP_EXPORT std::complex<double> operator*(const sparse_cvec &, const sparse_cvec &);

template ITPP_EXPORT int operator*(const sparse_ivec &, const ivec &);
template ITPP_EXPORT double operator*(const sparse_vec &, const vec &);
template ITPP_EXPORT std::complex<double> operator*(const sparse_cvec &, const cvec &);

template ITPP_EXPORT int operator*(const ivec &, const sparse_ivec &);
template ITPP_EXPORT double operator*(const vec &, const sparse_vec &);
template ITPP_EXPORT std::complex<double> operator*(const cvec &, const sparse_cvec &);

template ITPP_EXPORT sparse_ivec elem_mult(const sparse_ivec &, const sparse_ivec &);
template ITPP_EXPORT sparse_vec elem_mult(const sparse_vec &, const sparse_vec &);
template ITPP_EXPORT sparse_cvec elem_mult(const sparse_cvec &, const sparse_cvec &);

template ITPP_EXPORT ivec elem_mult(const sparse_ivec &, const ivec &);
template ITPP_EXPORT vec elem_mult(const sparse_vec &, const vec &);
template ITPP_EXPORT cvec elem_mult(const sparse_cvec &, const cvec &);

template ITPP_EXPORT sparse_ivec elem_mult_s(const sparse_ivec &, const ivec &);
template ITPP_EXPORT sparse_vec elem_mult_s(const sparse_vec &, const vec &);
template ITPP_EXPORT sparse_cvec elem_mult_s(const sparse_cvec &, const cvec &);

template ITPP_EXPORT ivec elem_mult(const ivec &, const sparse_ivec &);
template ITPP_EXPORT vec elem_mult(const vec &, const sparse_vec &);
template ITPP_EXPORT cvec elem_mult(const cvec &, const sparse_cvec &);

template ITPP_EXPORT sparse_ivec elem_mult_s(const ivec &, const sparse_ivec &);
template ITPP_EXPORT sparse_vec elem_mult_s(const vec &, const sparse_vec &);
template ITPP_EXPORT sparse_cvec elem_mult_s(const cvec &, const sparse_cvec &);

} // namespace itpp

//! \endcond
