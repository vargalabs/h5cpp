/*!
 * \file
 * \brief Cholesky factorisation test program
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


TEST(Cholesky, All)
{
  ostringstream ss(ostringstream::out);
  const string ref[] = {"X = [[3.1227 -2.07284 4.95845 -0.526723 0.837649]\n\
 [-2.07284 4.04154 -3.41073 2.58253 -2.42765]\n\
 [4.95845 -3.41073 15.1747 -0.518676 0.767615]\n\
 [-0.526723 2.58253 -0.518676 2.69632 -1.01246]\n\
 [0.837649 -2.42765 0.767615 -1.01246 2.36107]]",
 "X = [[3.84646+0i 0.263101+2.41269i 1.41977-0.156761i -4.66083+1.4819i -1.61119+0.105067i]\n\
 [0.263101-2.41269i 6.52894+0i 1.83059+1.97006i -1.05794+4.01142i -1.20091+2.68544i]\n\
 [1.41977+0.156761i 1.83059-1.97006i 5.8672+0i -0.433947+1.42927i -1.39383+0.665716i]\n\
 [-4.66083-1.4819i -1.05794-4.01142i -0.433947-1.42927i 9.4098+0i 1.13971-1.30921i]\n\
 [-1.61119-0.105067i -1.20091-2.68544i -1.39383-0.665716i 1.13971+1.30921i 3.75041+0i]]"};

  int k = 0;

  RNG_reset(0);

  {
    mat X, F;
    bool ok = false;

    while (!ok) {
      X = randn(5, 5);
      X = transpose(X) * X; // create a symmetric matrix
      // Make diagonal real and positive
      for (int i = 0; i < X.cols(); i++)
        X(i, i) = std::abs(X(i, i));

      ok = chol(X, F);
      ASSERT_TRUE(ok);
      ss << "X = " << round_to_zero(X);
      ASSERT_TRUE(ss.str() == ref[k++]);
      ss.str("");
      ASSERT_EQ(0, round_to_zero(norm(X - transpose(F) * F)));
    }
  }

  {
    cmat X, F;
    bool ok = false;

    while (!ok) {
      X = randn_c(5, 5);
      X = hermitian_transpose(X) * X; // create a symmetric matrix
      // Make diagonal real and positive
      for (int i = 0; i < X.cols(); i++)
        X(i, i) = std::abs(real(X(i, i)))
                  ;
      ok = chol(X, F);
      ASSERT_TRUE(ok);
      ss << "X = " << round_to_zero(X);
      ASSERT_TRUE(ss.str() == ref[k++]);
      ASSERT_NEAR(12.683, round_to_zero(norm(X - transpose(F) * F)), 1e-3);
    }
  }
}
