/*!
 * \file
 * \brief Schur decomposition test program
 * \author Adam Piatyszek
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

TEST(Schur, All)
{
  RNG_reset(0);
  static const double tol = 1e-9;

  int size = 5;

  // Test of Schur decomposition routines
  {
    // Real matrix
    mat A = randn(size, size);
    mat T, U;
    schur(A, U, T);

    ASSERT_NEAR(3.52148e-15, norm(A - (U * T * transpose(U))), tol);
    ASSERT_NEAR(9.92862e-16, norm(eye(size) - (U * transpose(U))), tol);
    double temp_sum = 0;
    for (int i = 2; i < size; i++)
      for (int j = 0; j < i - 1; j++)
        temp_sum += sqr(T(i, j));
    ASSERT_NEAR(0, sqrt(temp_sum), tol);
  }
  {
    // Complex matrix
    cmat A = randn_c(size, size);
    cmat T, U;
    schur(A, U, T);

    ASSERT_NEAR(7.24721e-15, norm(A - (U * T * hermitian_transpose(U))), tol);
    ASSERT_NEAR(1.73911e-15, norm(eye(size) - (U * hermitian_transpose(U))), tol);
    double temp_sum = 0;
    for (int i = 1; i < size; i++)
      for (int j = 0; j < i; j++)
        temp_sum += sqr(T(i, j));
    ASSERT_NEAR(0, sqrt(temp_sum), tol);
  }
}
