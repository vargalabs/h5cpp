/*!
 * \file
 * \brief Miscellaneous statistical routines test program
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
using namespace std;

static
const double tol = 1e-4;
static
void assert_vec_p(const vec &ref, const vec &act, int line)
{
  ASSERT_EQ(ref.length(), act.length()) << line;
  for (int n = 0; n < ref.length(); ++n) {
    ASSERT_NEAR(ref(n), act(n), tol) << line;
  }
}
#define assert_vec(ref, act) assert_vec_p(ref, act, __LINE__)

TEST(Stat, All)
{
  RNG_reset(0);

  // Test of statistical routines
  vec a = randn(5);
  ASSERT_NEAR(1.02823, max(a), tol);
  ASSERT_NEAR(-1.27491, min(a), tol);
  ASSERT_NEAR(-0.255701, mean(a), tol);
  ASSERT_NEAR(0.505336, geometric_mean(abs(a)), tol);
  ASSERT_NEAR(1.76712, norm(a), tol);
  ASSERT_NEAR(1.76712, norm(a, 2), tol);
  ASSERT_NEAR(3.33497, norm(a, 1), tol);
  ASSERT_NEAR(1.76712, norm(a, "fro"), tol);
  ASSERT_NEAR(3.1227, energy(a), tol);
  ASSERT_NEAR(0.698947, variance(a), tol);
  ASSERT_NEAR(0, moment(a, 1), tol);
  ASSERT_NEAR(0.559158, moment(a, 2), tol);
  ASSERT_NEAR(0.207182, moment(a, 3), tol);
  ASSERT_NEAR(0.528541, skewness(a), tol);
  ASSERT_NEAR(-1.7716, kurtosisexcess(a), tol);
  ASSERT_NEAR(1.2284, kurtosis(a), tol);

  mat A = randn(5, 5);
  vec ref = "1.29011 0.249247 0.875791 1.29677 1.44749";
  assert_vec(ref, max(A));
  ref = "1.29011 0.249247 0.875791 1.29677 1.44749";
  assert_vec(ref, max(A, 1));
  ref = "0.382527 1.29011 1.44749 1.29677 0.875791";
  assert_vec(ref, max(A, 2));
  ref = "-0.878434 -3.7327 -1.27628 -0.610395 -0.466618";
  assert_vec(ref, min(A));
  ref = "-0.878434 -3.7327 -1.27628 -0.610395 -0.466618";
  assert_vec(ref, min(A, 1));
  ref = "-1.27628 -3.7327 -0.675383 -0.940147 -0.19548";
  assert_vec(ref, min(A, 2));
  ASSERT_NEAR(-0.0769627, mean(A), tol);
  ASSERT_NEAR(0.518696, geometric_mean(abs(A)), tol);
  ASSERT_NEAR(4.11346, norm(A), tol);
  ASSERT_NEAR(4.11346, norm(A, 2), tol);
  ASSERT_NEAR(5.56682, norm(A, 1), tol);
  ASSERT_NEAR(5.23492, norm(A, "fro"), tol);
}
