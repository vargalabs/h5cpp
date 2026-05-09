/*!
 * \file
 * \brief Eigenvalue decomposition test program
 * \author Tony Ottosson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2012  (see AUTHORS file for a list of contributors)
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
using namespace std;


TEST(Eigen, All)
{
  double actual_norm = 0;
  const double eps = 1e-12;
  {
    //Real symmetric matrix
    mat A = randn(5, 5);
    A = transpose(A) * A; // make it symmetic
    mat V;
    vec d;
    eig_sym(A, d, V);
    actual_norm = norm(A * V - V * diag(d));
    ASSERT_NEAR(0, actual_norm, eps);
  }

  {
    //Real non-symmetric matrix
    mat A = randn(5, 5);
    cmat V;
    cvec d;
    eig(A, d, V);
    actual_norm = norm(A * V - V * diag(d));
    ASSERT_NEAR(0, actual_norm, eps);
  }

  {
    //Complex hermitian matrix
    cmat A = randn_c(5, 5);
    A = transpose(conj(A)) * A; // make it hermitian
    cmat V;
    vec d;
    eig_sym(A, d, V);
    actual_norm = norm(A * V - V * diag(d));
    ASSERT_NEAR(0, actual_norm, eps);
  }

  {
    //Complex non-hermitian matrix
    cmat A = randn_c(5, 5);
    cmat V;
    cvec d;
    eig(A, d, V);
    actual_norm = norm(A * V - V * diag(d));
    ASSERT_NEAR(0, actual_norm, eps);
  }
}
