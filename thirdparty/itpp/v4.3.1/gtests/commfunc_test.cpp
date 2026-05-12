/*!
 * \file
 * \brief Commfunc test program
 * \author Erik G. Larsson
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

#include <itpp/itcomm.h>
#include "gtest/gtest.h"

using namespace std;
using namespace itpp;

TEST(CommFunc, All)
{
  const double tol = 1e-9;

  bmat expected_gray = "0 0 0 0 0;\
                        0 0 0 0 1;\
                        0 0 0 1 1;\
                        0 0 0 1 0;\
                        0 0 1 1 0;\
                        0 0 1 1 1;\
                        0 0 1 0 1;\
                        0 0 1 0 0;\
                        0 1 1 0 0;\
                        0 1 1 0 1;\
                        0 1 1 1 1;\
                        0 1 1 1 0;\
                        0 1 0 1 0;\
                        0 1 0 1 1;\
                        0 1 0 0 1;\
                        0 1 0 0 0;\
                        1 1 0 0 0;\
                        1 1 0 0 1;\
                        1 1 0 1 1;\
                        1 1 0 1 0;\
                        1 1 1 1 0;\
                        1 1 1 1 1;\
                        1 1 1 0 1;\
                        1 1 1 0 0;\
                        1 0 1 0 0;\
                        1 0 1 0 1;\
                        1 0 1 1 1;\
                        1 0 1 1 0;\
                        1 0 0 1 0;\
                        1 0 0 1 1;\
                        1 0 0 0 1;\
                        1 0 0 0 0";
  ASSERT_TRUE(expected_gray == graycode(5));

  RNG_reset(0);
  bvec b1 = randb(10);
  bvec b2 = randb(10);
  ASSERT_EQ(8, hamming_distance(b1, b2));
  ASSERT_EQ(7, weight(b1));

  vec expected_water = "0 0 0.001";
  vec actual_water = waterfilling("1 2 3", 0.001);
  int i = 0;
  for (i = 0; i < actual_water.length(); ++i) {
    ASSERT_NEAR(expected_water[i], actual_water[i], tol);
  }
  expected_water = "0 0 0.010000000000000009";
  actual_water = waterfilling("1 2 3", 0.01);
  for (i = 0; i < actual_water.length(); ++i) {
    ASSERT_NEAR(expected_water[i], actual_water[i], tol);
  }
  expected_water = "0 0 0.1";
  actual_water = waterfilling("1 2 3", 0.1);
  for (i = 0; i < actual_water.length(); ++i) {
    ASSERT_NEAR(expected_water[i], actual_water[i], tol);
  }
  expected_water = "0 0.16666666666666663 0.33333333333333331";
  actual_water = waterfilling("1 2 3", 0.5);
  for (i = 0; i < actual_water.length(); ++i) {
    ASSERT_NEAR(expected_water[i], actual_water[i], tol);
  }
  expected_water = "0 0.41666666666666652 0.58333333333333326";
  actual_water = waterfilling("1 2 3", 1);
  for (i = 0; i < actual_water.length(); ++i) {
    ASSERT_NEAR(expected_water[i], actual_water[i], tol);
  }
  expected_water = "1.2777777777777777 1.7777777777777777 1.9444444444444444";
  actual_water = waterfilling("1 2 3", 5);
  for (i = 0; i < actual_water.length(); ++i) {
    ASSERT_NEAR(expected_water[i], actual_water[i], tol);
  }
  expected_water = "32.944444444444443 33.444444444444443 33.611111111111107";
  actual_water = waterfilling("1 2 3", 100);
  for (i = 0; i < actual_water.length(); ++i) {
    ASSERT_NEAR(expected_water[i], actual_water[i], tol);
  }
}
