/*!
 * \file
 * \brief Implementation of QR factorisation functions
 * \author Tony Ottosson, Simon Wood, Adam Piatyszek and Vasek Smidl
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

#include <itpp/base/algebra/qr.h>
#include <itpp/base/specmat.h>


namespace itpp
{

#if defined(HAVE_LAPACK)

bool qr(const mat &A, mat &Q, mat &R)
{
  int info;
  int m = A.rows();
  int n = A.cols();
  int lwork = n;
  int k = std::min(m, n);
  vec tau(k);
  vec work(lwork);

  R = A;

  // perform workspace query for optimum lwork value
  int lwork_tmp = -1;
  dgeqrf_(&m, &n, R._data(), &m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(work(0));
    work.set_size(lwork, false);
  }
  dgeqrf_(&m, &n, R._data(), &m, tau._data(), work._data(), &lwork, &info);
  Q = R;
  Q.set_size(m, m, true);

  // construct R
  for (int i = 0; i < m; i++)
    for (int j = 0; j < std::min(i, n); j++)
      R(i, j) = 0;

  // perform workspace query for optimum lwork value
  lwork_tmp = -1;
  dorgqr_(&m, &m, &k, Q._data(), &m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(work(0));
    work.set_size(lwork, false);
  }
  dorgqr_(&m, &m, &k, Q._data(), &m, tau._data(), work._data(), &lwork,
          &info);

  return (info == 0);
}

bool qr(const mat &A, mat &R)
{
  int info;
  int m = A.rows();
  int n = A.cols();
  int lwork = n;
  int k = std::min(m, n);
  vec tau(k);
  vec work(lwork);

  R = A;

  // perform workspace query for optimum lwork value
  int lwork_tmp = -1;
  dgeqrf_(&m, &n, R._data(), &m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(work(0));
    work.set_size(lwork, false);
  }
  dgeqrf_(&m, &n, R._data(), &m, tau._data(), work._data(), &lwork, &info);

  // construct R
  for (int i = 0; i < m; i++)
    for (int j = 0; j < std::min(i, n); j++)
      R(i, j) = 0;

  return (info == 0);
}

bool qr(const mat &A, mat &Q, mat &R, bmat &P)
{
  int info;
  int m = A.rows();
  int n = A.cols();
  int lwork = n;
  int k = std::min(m, n);
  vec tau(k);
  vec work(lwork);
  ivec jpvt(n);
  jpvt.zeros();

  R = A;

  // perform workspace query for optimum lwork value
  int lwork_tmp = -1;
  dgeqp3_(&m, &n, R._data(), &m, jpvt._data(), tau._data(), work._data(),
          &lwork_tmp, &info);
  if (info == 0) {
    lwork = static_cast<int>(work(0));
    work.set_size(lwork, false);
  }
  dgeqp3_(&m, &n, R._data(), &m, jpvt._data(), tau._data(), work._data(),
          &lwork, &info);
  Q = R;
  Q.set_size(m, m, true);

  // construct permutation matrix
  P = zeros_b(n, n);
  for (int j = 0; j < n; j++)
    P(jpvt(j) - 1, j) = 1;

  // construct R
  for (int i = 0; i < m; i++)
    for (int j = 0; j < std::min(i, n); j++)
      R(i, j) = 0;

  // perform workspace query for optimum lwork value
  lwork_tmp = -1;
  dorgqr_(&m, &m, &k, Q._data(), &m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(work(0));
    work.set_size(lwork, false);
  }
  dorgqr_(&m, &m, &k, Q._data(), &m, tau._data(), work._data(), &lwork,
          &info);

  return (info == 0);
}



bool qr(const cmat &A, cmat &Q, cmat &R)
{
  int info;
  int m = A.rows();
  int n = A.cols();
  int lwork = n;
  int k = std::min(m, n);
  cvec tau(k);
  cvec work(lwork);

  R = A;

  // perform workspace query for optimum lwork value
  int lwork_tmp = -1;
  zgeqrf_(&m, &n, R._data(), &m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(real(work(0)));
    work.set_size(lwork, false);
  }
  zgeqrf_(&m, &n, R._data(), &m, tau._data(), work._data(), &lwork, &info);

  Q = R;
  Q.set_size(m, m, true);

  // construct R
  for (int i = 0; i < m; i++)
    for (int j = 0; j < std::min(i, n); j++)
      R(i, j) = 0;

  // perform workspace query for optimum lwork value
  lwork_tmp = -1;
  zungqr_(&m, &m, &k, Q._data(), &m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(real(work(0)));
    work.set_size(lwork, false);
  }
  zungqr_(&m, &m, &k, Q._data(), &m, tau._data(), work._data(), &lwork,
          &info);

  return (info == 0);
}

bool qr(const cmat &A, cmat &R)
{
  int info;
  int m = A.rows();
  int n = A.cols();
  int lwork = n;
  int k = std::min(m, n);
  cvec tau(k);
  cvec work(lwork);

  R = A;

  // perform workspace query for optimum lwork value
  int lwork_tmp = -1;
  zgeqrf_(&m, &n, R._data(), &m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(real(work(0)));
    work.set_size(lwork, false);
  }
  zgeqrf_(&m, &n, R._data(), &m, tau._data(), work._data(), &lwork, &info);

  // construct R
  for (int i = 0; i < m; i++)
    for (int j = 0; j < std::min(i, n); j++)
      R(i, j) = 0;

  return (info == 0);
}

bool qr(const cmat &A, cmat &Q, cmat &R, bmat &P)
{
  int info;
  int m = A.rows();
  int n = A.cols();
  int lwork = n;
  int k = std::min(m, n);
  cvec tau(k);
  cvec work(lwork);
  vec rwork(std::max(1, 2*n));
  ivec jpvt(n);
  jpvt.zeros();

  R = A;

  // perform workspace query for optimum lwork value
  int lwork_tmp = -1;
  zgeqp3_(&m, &n, R._data(), &m, jpvt._data(), tau._data(), work._data(),
          &lwork_tmp, rwork._data(), &info);
  if (info == 0) {
    lwork = static_cast<int>(real(work(0)));
    work.set_size(lwork, false);
  }
  zgeqp3_(&m, &n, R._data(), &m, jpvt._data(), tau._data(), work._data(),
          &lwork, rwork._data(), &info);

  Q = R;
  Q.set_size(m, m, true);

  // construct permutation matrix
  P = zeros_b(n, n);
  for (int j = 0; j < n; j++)
    P(jpvt(j) - 1, j) = 1;

  // construct R
  for (int i = 0; i < m; i++)
    for (int j = 0; j < std::min(i, n); j++)
      R(i, j) = 0;

  // perform workspace query for optimum lwork value
  lwork_tmp = -1;
  zungqr_(&m, &m, &k, Q._data(), &m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(real(work(0)));
    work.set_size(lwork, false);
  }
  zungqr_(&m, &m, &k, Q._data(), &m, tau._data(), work._data(), &lwork,
          &info);

  return (info == 0);
}

#else

bool qr(const mat &A, mat &Q, mat &R)
{
  it_error("LAPACK library is needed to use qr() function");
  return false;
}

bool qr(const mat &A, mat &R)
{
  it_error("LAPACK library is needed to use qr() function");
  return false;
}

bool qr(const mat &A, mat &Q, mat &R, bmat &P)
{
  it_error("LAPACK library is needed to use qr() function");
  return false;
}

bool qr(const cmat &A, cmat &Q, cmat &R)
{
  it_error("LAPACK library is needed to use qr() function");
  return false;
}

bool qr(const cmat &A, cmat &R)
{
  it_error("LAPACK library is needed to use qr() function");
  return false;
}

bool qr(const cmat &A, cmat &Q, cmat &R, bmat &P)
{
  it_error("LAPACK library is needed to use qr() function");
  return false;
}

#endif // HAVE_LAPACK

} // namespace itpp
