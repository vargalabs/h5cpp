/*!
 * \file
 * \brief Implementation of Cholesky factorisation functions
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

#include <itpp/base/algebra/cholesky.h>


namespace itpp
{

#if defined(HAVE_LAPACK)

bool chol(const mat &X, mat &F)
{

  char uplo = 'U';
  int n, lda, info;
  n = lda = X.rows();

  F = X; // input matrix is overwritten

  dpotrf_(&uplo, &n, F._data(), &lda, &info);

  // Set lower part to zero
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
      F(j, i) = 0;

  return (info == 0);
}

bool chol(const cmat &X, cmat &F)
{
  char uplo = 'U';
  int n, lda, info;
  n = lda = X.rows();

  F = X; // input matrix is overwritten

  zpotrf_(&uplo, &n, F._data(), &lda, &info);

  // Set lower part to zero
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
      F(j, i) = 0;

  return (info == 0);
}

#else // HAVE_LAPACK

bool chol(const mat &X, mat &F)
{
  it_error("LAPACK library is needed to use chol() function");
  return false;
}

bool chol(const cmat &X, cmat &F)
{

  it_error("LAPACK library is needed to use chol() function");
  return false;
}

#endif // HAVE_LAPACK

cmat chol(const cmat &X)
{
  cmat F;
  if (!chol(X, F)) {
    it_warning("cholesky factorization didn't succeed");
  }

  return F;
}

mat chol(const mat &X)
{
  mat F;
  if (!chol(X, F)) {
    it_warning("cholesky factorization didn't succeed");
  }

  return F;
}

} // namespace itpp
