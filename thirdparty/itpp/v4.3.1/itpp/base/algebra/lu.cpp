/*!
 * \file
 * \brief Implementation of LU factorisation functions.
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

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#if defined(HAVE_LAPACK)
#  include <itpp/base/algebra/lapack.h>
#endif

#include <itpp/base/algebra/lu.h>
#include <itpp/base/specmat.h>


namespace itpp
{

#if defined(HAVE_LAPACK)

bool lu(const mat &X, mat &L, mat &U, ivec &p)
{
  it_assert_debug(X.rows() == X.cols(), "lu: matrix is not quadratic");
  //int m, n, lda, info;
  //m = n = lda = X.rows();
  int m = X.rows(), info;

  mat A(X);
  L.set_size(m, m, false);
  U.set_size(m, m, false);
  p.set_size(m, false);

  dgetrf_(&m, &m, A._data(), &m, p._data(), &info);

  for (int i = 0; i < m; i++) {
    for (int j = i; j < m; j++) {
      if (i == j) { // diagonal
        L(i, j) = 1;
        U(i, j) = A(i, j);
      }
      else { // upper and lower triangular parts
        L(i, j) = U(j, i) = 0;
        L(j, i) = A(j, i);
        U(i, j) = A(i, j);
      }
    }
  }

  p = p - 1; // Fortran counts from 1

  return (info == 0);
}

// Slower than not using LAPACK when matrix size smaller than approx 20.
bool lu(const cmat &X, cmat &L, cmat &U, ivec &p)
{
  it_assert_debug(X.rows() == X.cols(), "lu: matrix is not quadratic");
  //int m, n, lda, info;
  //m = n = lda = X.rows();
  int m = X.rows(), info;

  cmat A(X);
  L.set_size(m, m, false);
  U.set_size(m, m, false);
  p.set_size(m, false);

  zgetrf_(&m, &m, A._data(), &m, p._data(), &info);

  for (int i = 0; i < m; i++) {
    for (int j = i; j < m; j++) {
      if (i == j) { // diagonal
        L(i, j) = 1;
        U(i, j) = A(i, j);
      }
      else { // upper and lower triangular parts
        L(i, j) = U(j, i) = 0;
        L(j, i) = A(j, i);
        U(i, j) = A(i, j);
      }
    }
  }

  p = p - 1; // Fortran counts from 1

  return (info == 0);
}

#else

bool lu(const mat &X, mat &L, mat &U, ivec &p)
{
  it_error("LAPACK library is needed to use lu() function");
  return false;
}

bool lu(const cmat &X, cmat &L, cmat &U, ivec &p)
{
  it_error("LAPACK library is needed to use lu() function");
  return false;
}

#endif // HAVE_LAPACK


void interchange_permutations(vec &b, const ivec &p)
{
  it_assert(b.size() == p.size(), "interchange_permutations(): dimension mismatch");
  double temp;

  for (int k = 0; k < b.size(); k++) {
    temp = b(k);
    b(k) = b(p(k));
    b(p(k)) = temp;
  }
}

bmat permutation_matrix(const ivec &p)
{
  it_assert(p.size() > 0, "permutation_matrix(): vector must have nonzero size");
  int n = p.size(), k;
  bmat P, identity;
  bvec row_k, row_pk;
  identity = eye_b(n);


  for (k = n - 1; k >= 0; k--) {
    // swap rows k and p(k) in identity
    row_k = identity.get_row(k);
    row_pk = identity.get_row(p(k));
    identity.set_row(k, row_pk);
    identity.set_row(p(k), row_k);

    if (k == n - 1) {
      P = identity;
    }
    else {
      P *= identity;
    }

    // swap back
    identity.set_row(k, row_k);
    identity.set_row(p(k), row_pk);
  }
  return P;
}

} // namespace itpp
