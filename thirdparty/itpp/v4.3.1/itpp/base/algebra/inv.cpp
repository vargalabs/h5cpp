/*!
 * \file
 * \brief Implementation of matrix inversion routines
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

#include <itpp/base/algebra/inv.h>


namespace itpp
{

#if defined(HAVE_LAPACK)

bool inv(const mat &X, mat &Y)
{
  it_assert_debug(X.rows() == X.cols(), "inv: matrix is not square");

  int m = X.rows(), info, lwork;
  lwork = m; // may be choosen better

  ivec p(m);
  Y = X;
  vec work(lwork);

  dgetrf_(&m, &m, Y._data(), &m, p._data(), &info); // LU-factorization
  if (info != 0)
    return false;

  dgetri_(&m, Y._data(), &m, p._data(), work._data(), &lwork, &info);
  return (info == 0);
}

bool inv(const cmat &X, cmat &Y)
{
  it_assert_debug(X.rows() == X.cols(), "inv: matrix is not square");

  int m = X.rows(), info, lwork;
  lwork = m; // may be choosen better

  ivec p(m);
  Y = X;
  cvec work(lwork);

  zgetrf_(&m, &m, Y._data(), &m, p._data(), &info); // LU-factorization
  if (info != 0)
    return false;

  zgetri_(&m, Y._data(), &m, p._data(), work._data(), &lwork, &info);
  return (info == 0);
}

#else

bool inv(const mat &X, mat &Y)
{
  it_error("LAPACK library is needed to use inv() function");
  return false;
}

bool inv(const cmat &X, cmat &Y)
{
  it_error("LAPACK library is needed to use inv() function");
  return false;
}

#endif // HAVE_LAPACK

cmat inv(const cmat &X)
{
  cmat Y;
  inv(X, Y);
  return Y;
}


mat inv(const mat &X)
{
  mat Y;
  inv(X, Y);
  return Y;
}

} // namespace itpp
