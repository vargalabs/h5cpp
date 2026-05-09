/*!
 * \file
 * \brief Various functions on vectors and matrices - source file
 * \author Tony Ottosson, Adam Piatyszek, Conrad Sanderson, Mark Dobossy
 *         and Martin Senst
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

#include <itpp/base/matfunc.h>
#include <itpp/base/algebra/schur.h>
#include <itpp/base/converters.h>
#include <limits>

//! \cond

namespace itpp
{

// Square root of a square matrix (based on Octave sqrtm implementation)
cmat sqrtm(const mat& A)
{
  return sqrtm(to_cmat(A));
}

// Square root of the complex square matrix A
cmat sqrtm(const cmat& A)
{
  cmat U, T;
  schur(A, U, T);

  int n = U.rows();
  cmat R(n, n);

  R.zeros();
  for (int j = 0; j < n; j++)
    R(j, j) = std::sqrt(T(j, j));

  const double fudge = std::sqrt(std::numeric_limits<double>::min());

  for (int p = 0; p < n - 1; p++) {
    for (int i = 0; i < n - (p + 1); i++) {
      const int j = i + p + 1;
      std::complex<double> s = T(i, j);

      for (int k = i + 1; k < j; k++)
        s -= R(i, k) * R(k, j);

      const std::complex<double> d = R(i, i) + R(j, j) + fudge;
      const std::complex<double> conj_d = conj(d);

      R(i, j) = (s * conj_d) / (d * conj_d);
    }
  }

  return U * R * U.H();
}


bool all(const bvec &testvec)
{
  for (int i = 0; i < testvec.length(); i++)
    if (!testvec(i)) return false;
  return true;
}

bool any(const bvec &testvec)
{
  for (int i = 0; i < testvec.length(); i++)
    if (testvec(i)) return true;
  return false;
}

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------

template ITPP_EXPORT int length(const vec &v);
template ITPP_EXPORT int length(const cvec &v);
template ITPP_EXPORT int length(const svec &v);
template ITPP_EXPORT int length(const ivec &v);
template ITPP_EXPORT int length(const bvec &v);

template ITPP_EXPORT double sum(const vec &v);
template ITPP_EXPORT std::complex<double> sum(const cvec &v);
template ITPP_EXPORT short sum(const svec &v);
template ITPP_EXPORT int sum(const ivec &v);
template ITPP_EXPORT bin sum(const bvec &v);

template ITPP_EXPORT double sum_sqr(const vec &v);
template ITPP_EXPORT std::complex<double> sum_sqr(const cvec &v);
template ITPP_EXPORT short sum_sqr(const svec &v);
template ITPP_EXPORT int sum_sqr(const ivec &v);
template ITPP_EXPORT bin sum_sqr(const bvec &v);

template ITPP_EXPORT vec cumsum(const vec &v);
template ITPP_EXPORT cvec cumsum(const cvec &v);
template ITPP_EXPORT svec cumsum(const svec &v);
template ITPP_EXPORT ivec cumsum(const ivec &v);
template ITPP_EXPORT bvec cumsum(const bvec &v);

template ITPP_EXPORT double prod(const vec &v);
template ITPP_EXPORT std::complex<double> prod(const cvec &v);
template ITPP_EXPORT short prod(const svec &v);
template ITPP_EXPORT int prod(const ivec &v);
template ITPP_EXPORT bin prod(const bvec &v);

template ITPP_EXPORT vec cross(const vec &v1, const vec &v2);
template ITPP_EXPORT cvec cross(const cvec &v1, const cvec &v2);
template ITPP_EXPORT ivec cross(const ivec &v1, const ivec &v2);
template ITPP_EXPORT svec cross(const svec &v1, const svec &v2);
template ITPP_EXPORT bvec cross(const bvec &v1, const bvec &v2);

template ITPP_EXPORT vec reverse(const vec &in);
template ITPP_EXPORT cvec reverse(const cvec &in);
template ITPP_EXPORT svec reverse(const svec &in);
template ITPP_EXPORT ivec reverse(const ivec &in);
template ITPP_EXPORT bvec reverse(const bvec &in);

template ITPP_EXPORT vec zero_pad(const vec &v, int n);
template ITPP_EXPORT cvec zero_pad(const cvec &v, int n);
template ITPP_EXPORT ivec zero_pad(const ivec &v, int n);
template ITPP_EXPORT svec zero_pad(const svec &v, int n);
template ITPP_EXPORT bvec zero_pad(const bvec &v, int n);

template ITPP_EXPORT vec zero_pad(const vec &v);
template ITPP_EXPORT cvec zero_pad(const cvec &v);
template ITPP_EXPORT ivec zero_pad(const ivec &v);
template ITPP_EXPORT svec zero_pad(const svec &v);
template ITPP_EXPORT bvec zero_pad(const bvec &v);

template ITPP_EXPORT mat  zero_pad(const mat &, int, int);
template ITPP_EXPORT cmat zero_pad(const cmat &, int, int);
template ITPP_EXPORT imat zero_pad(const imat &, int, int);
template ITPP_EXPORT smat zero_pad(const smat &, int, int);
template ITPP_EXPORT bmat zero_pad(const bmat &, int, int);

template ITPP_EXPORT vec sum(const mat &m, int dim);
template ITPP_EXPORT cvec sum(const cmat &m, int dim);
template ITPP_EXPORT svec sum(const smat &m, int dim);
template ITPP_EXPORT ivec sum(const imat &m, int dim);
template ITPP_EXPORT bvec sum(const bmat &m, int dim);

template ITPP_EXPORT double sumsum(const mat &X);
template ITPP_EXPORT std::complex<double> sumsum(const cmat &X);
template ITPP_EXPORT short sumsum(const smat &X);
template ITPP_EXPORT int sumsum(const imat &X);
template ITPP_EXPORT bin sumsum(const bmat &X);

template ITPP_EXPORT vec sum_sqr(const mat & m, int dim);
template ITPP_EXPORT cvec sum_sqr(const cmat &m, int dim);
template ITPP_EXPORT svec sum_sqr(const smat &m, int dim);
template ITPP_EXPORT ivec sum_sqr(const imat &m, int dim);
template ITPP_EXPORT bvec sum_sqr(const bmat &m, int dim);

template ITPP_EXPORT mat cumsum(const mat &m, int dim);
template ITPP_EXPORT cmat cumsum(const cmat &m, int dim);
template ITPP_EXPORT smat cumsum(const smat &m, int dim);
template ITPP_EXPORT imat cumsum(const imat &m, int dim);
template ITPP_EXPORT bmat cumsum(const bmat &m, int dim);

template ITPP_EXPORT vec prod(const mat &m, int dim);
template ITPP_EXPORT cvec prod(const cmat &v, int dim);
template ITPP_EXPORT svec prod(const smat &m, int dim);
template ITPP_EXPORT ivec prod(const imat &m, int dim);
template ITPP_EXPORT bvec prod(const bmat &m, int dim);

template ITPP_EXPORT vec diag(const mat &in);
template ITPP_EXPORT cvec diag(const cmat &in);

template ITPP_EXPORT void diag(const vec &in, mat &m);
template ITPP_EXPORT void diag(const cvec &in, cmat &m);

template ITPP_EXPORT mat diag(const vec &v, const int K);
template ITPP_EXPORT cmat diag(const cvec &v, const int K);

template ITPP_EXPORT mat bidiag(const vec &, const vec &);
template ITPP_EXPORT cmat bidiag(const cvec &, const cvec &);

template ITPP_EXPORT void bidiag(const vec &, const vec &, mat &);
template ITPP_EXPORT void bidiag(const cvec &, const cvec &, cmat &);

template ITPP_EXPORT void bidiag(const mat &, vec &, vec &);
template ITPP_EXPORT void bidiag(const cmat &, cvec &, cvec &);

template ITPP_EXPORT mat tridiag(const vec &main, const vec &, const vec &);
template ITPP_EXPORT cmat tridiag(const cvec &main, const cvec &, const cvec &);

template ITPP_EXPORT void tridiag(const vec &main, const vec &, const vec &, mat &);
template ITPP_EXPORT void tridiag(const cvec &main, const cvec &, const cvec &, cmat &);

template ITPP_EXPORT void tridiag(const mat &m, vec &, vec &, vec &);
template ITPP_EXPORT void tridiag(const cmat &m, cvec &, cvec &, cvec &);

template ITPP_EXPORT double trace(const mat &in);
template ITPP_EXPORT std::complex<double> trace(const cmat &in);
template ITPP_EXPORT short trace(const smat &in);
template ITPP_EXPORT int trace(const imat &in);
template ITPP_EXPORT bin trace(const bmat &in);

template ITPP_EXPORT void transpose(const mat &m, mat &out);
template ITPP_EXPORT void transpose(const cmat &m, cmat &out);
template ITPP_EXPORT void transpose(const smat &m, smat &out);
template ITPP_EXPORT void transpose(const imat &m, imat &out);
template ITPP_EXPORT void transpose(const bmat &m, bmat &out);

template ITPP_EXPORT mat transpose(const mat &m);
template ITPP_EXPORT cmat transpose(const cmat &m);
template ITPP_EXPORT smat transpose(const smat &m);
template ITPP_EXPORT imat transpose(const imat &m);
template ITPP_EXPORT bmat transpose(const bmat &m);

template ITPP_EXPORT void hermitian_transpose(const mat &m, mat &out);
template ITPP_EXPORT void hermitian_transpose(const cmat &m, cmat &out);
template ITPP_EXPORT void hermitian_transpose(const smat &m, smat &out);
template ITPP_EXPORT void hermitian_transpose(const imat &m, imat &out);
template ITPP_EXPORT void hermitian_transpose(const bmat &m, bmat &out);

template ITPP_EXPORT mat hermitian_transpose(const mat &m);
template ITPP_EXPORT cmat hermitian_transpose(const cmat &m);
template ITPP_EXPORT smat hermitian_transpose(const smat &m);
template ITPP_EXPORT imat hermitian_transpose(const imat &m);
template ITPP_EXPORT bmat hermitian_transpose(const bmat &m);

template ITPP_EXPORT bool is_hermitian(const mat &X);
template ITPP_EXPORT bool is_hermitian(const cmat &X);

template ITPP_EXPORT bool is_unitary(const mat &X);
template ITPP_EXPORT bool is_unitary(const cmat &X);

template ITPP_EXPORT vec rvectorize(const mat &m);
template ITPP_EXPORT cvec rvectorize(const cmat &m);
template ITPP_EXPORT ivec rvectorize(const imat &m);
template ITPP_EXPORT svec rvectorize(const smat &m);
template ITPP_EXPORT bvec rvectorize(const bmat &m);

template ITPP_EXPORT vec cvectorize(const mat &m);
template ITPP_EXPORT cvec cvectorize(const cmat &m);
template ITPP_EXPORT ivec cvectorize(const imat &m);
template ITPP_EXPORT svec cvectorize(const smat &m);
template ITPP_EXPORT bvec cvectorize(const bmat &m);

template ITPP_EXPORT mat reshape(const mat &m, int rows, int cols);
template ITPP_EXPORT cmat reshape(const cmat &m, int rows, int cols);
template ITPP_EXPORT imat reshape(const imat &m, int rows, int cols);
template ITPP_EXPORT smat reshape(const smat &m, int rows, int cols);
template ITPP_EXPORT bmat reshape(const bmat &m, int rows, int cols);

template ITPP_EXPORT mat reshape(const vec &m, int rows, int cols);
template ITPP_EXPORT cmat reshape(const cvec &m, int rows, int cols);
template ITPP_EXPORT imat reshape(const ivec &m, int rows, int cols);
template ITPP_EXPORT smat reshape(const svec &m, int rows, int cols);
template ITPP_EXPORT bmat reshape(const bvec &m, int rows, int cols);

template ITPP_EXPORT mat kron(const mat &X, const mat &Y);
template ITPP_EXPORT cmat kron(const cmat &X, const cmat &Y);
template ITPP_EXPORT imat kron(const imat &X, const imat &Y);
template ITPP_EXPORT smat kron(const smat &X, const smat &Y);
template ITPP_EXPORT bmat kron(const bmat &X, const bmat &Y);

template ITPP_EXPORT vec repmat(const vec &v, int n);
template ITPP_EXPORT cvec repmat(const cvec &v, int n);
template ITPP_EXPORT ivec repmat(const ivec &v, int n);
template ITPP_EXPORT svec repmat(const svec &v, int n);
template ITPP_EXPORT bvec repmat(const bvec &v, int n);

template ITPP_EXPORT mat repmat(const vec &v, int m, int n, bool transpose);
template ITPP_EXPORT cmat repmat(const cvec &v, int m, int n, bool transpose);
template ITPP_EXPORT imat repmat(const ivec &v, int m, int n, bool transpose);
template ITPP_EXPORT smat repmat(const svec &v, int m, int n, bool transpose);
template ITPP_EXPORT bmat repmat(const bvec &v, int m, int n, bool transpose);

template ITPP_EXPORT mat repmat(const mat &data, int m, int n);
template ITPP_EXPORT cmat repmat(const cmat &data, int m, int n);
template ITPP_EXPORT imat repmat(const imat &data, int m, int n);
template ITPP_EXPORT smat repmat(const smat &data, int m, int n);
template ITPP_EXPORT bmat repmat(const bmat &data, int m, int n);

} // namespace itpp

//! \endcond
