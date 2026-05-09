/*!
 * \file
 * \brief QR factorisation test program
 * \author Tony Ottosson, Adam Piatyszek and Vasek Smidl
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2013  (see AUTHORS file for a list of contributors)
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

#include <itpp/itstat.h>
#include "gtest/gtest.h"

using namespace itpp;

TEST(QR, All)
{
  // Test of QR factorization routines
  RNG_reset(0);
  static const double tol = 1e-9;

  {
    // QR of Real matrix
    mat Q, R, e;
    mat A = randn(5, 5);
    qr(A, Q, R);
    ASSERT_NEAR(9.14651e-16, norm(A - Q * R), tol);

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn(4, 2);
    qr(A, Q, R);
    ASSERT_NEAR(3.60139e-16, norm(A - Q * R), tol);

    A = randn(2, 4);
    qr(A, Q, R);
    ASSERT_NEAR(5.31202e-16, norm(A - Q * R), tol);
  }

  {
    // QR of Real matrix without Q
    mat R;
    mat A = randn(5, 5);
    qr(A, R);
    ASSERT_NEAR(4.7391e-15, norm(A.T()*A -  R.T()*R), tol);

    A = randn(4, 2);
    qr(A, R);
    ASSERT_NEAR(5.68779e-16, norm(A.T()*A -  R.T()*R), tol);

    A = randn(2, 4);
    qr(A, R);
    ASSERT_NEAR(4.73573e-16, norm(A.T()*A -  R.T()*R), tol);
  }

  {
    // QR of Real matrix with pivoting
    mat Q, R, e;
    bmat P;
    mat A = randn(5, 5);
    qr(A, Q, R, P);
    e = A * to_mat(P) - Q * R;
    ASSERT_NEAR(1.13066e-15, norm(e), tol);

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn(4, 2);
    qr(A, Q, R, P);
    e = A * to_mat(P) - Q * R;
    ASSERT_NEAR(5.05084e-16, norm(e), tol);

    A = randn(2, 4);
    qr(A, Q, R, P);
    e = A * to_mat(P) - Q * R;
    ASSERT_NEAR(2.48662e-16, norm(e), tol);
  }

  {
    // QR of Complex matrix
    cmat A = randn_c(5, 5);
    cmat Q, R, e;

    qr(A, Q, R);
    e = A - Q * R;
    ASSERT_NEAR(1.01366e-15, norm(e), tol);

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn_c(4, 2);
    qr(A, Q, R);
    e = A - Q * R;
    ASSERT_NEAR(2.91757e-16, norm(e), tol);

    A = randn_c(2, 4);
    qr(A, Q, R);
    e = A - Q * R;
    ASSERT_NEAR(1.71413e-16, norm(e), tol);
  }

  {
    // QR of Complex matrix without Q
    cmat A = randn_c(5, 5);
    cmat R, e;

    qr(A, R);
    e = A.H() * A - R.H() * R;
    ASSERT_NEAR(5.20531e-15, norm(e), tol);

    A = randn_c(4, 2);
    qr(A, R);
    e = A.H() * A - R.H() * R;
    ASSERT_NEAR(1.78189e-15, norm(e), tol);

    A = randn_c(2, 4);
    qr(A, R);
    e = A.H() * A - R.H() * R;
    ASSERT_NEAR(5.80022e-16, norm(e), tol);
  }

  {
    // QR of Complex matrix with pivoting
    cmat A = randn_c(5, 5);
    cmat Q, R, e;
    bmat P;

    qr(A, Q, R, P);
    e = A * to_mat(P) - Q * R;
    ASSERT_NEAR(7.59219e-16, norm(e), tol);

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn_c(4, 2);
    qr(A, Q, R, P);
    e = A * to_mat(P) - Q * R;
    ASSERT_NEAR(6.6432e-16, norm(e), tol);

    A = randn_c(2, 4);
    qr(A, Q, R, P);
    e = A * to_mat(P) - Q * R;
    ASSERT_NEAR(1.55268e-16, norm(e), tol);
  }
}
