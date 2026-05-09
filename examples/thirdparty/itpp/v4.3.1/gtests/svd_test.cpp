/*!
 * \file
 * \brief SVD decomposition test program
 * \author Tony Ottosson and Adam Piatyszek
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

TEST(SVD, All)
{
  RNG_reset(0);
  static const double tol = 1e-4;
  // Test of svd routines

  {
    // Real matrix
    mat A = randn(5, 5);
    mat U, V;
    vec S;
    svd(A, U, S, V);
    ASSERT_NEAR(2.38835e-15, norm(A - U * diag(S) * transpose(V)), tol);

  }
  {
    // Complex matrix
    cmat A = randn_c(5, 5);
    cmat U, V;
    vec S;
    svd(A, U, S, V);
    ASSERT_NEAR(6.20568e-15, norm(A - U * diag(S) * hermitian_transpose(V)), tol);
  }
}
