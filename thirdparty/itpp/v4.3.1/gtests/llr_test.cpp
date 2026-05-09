/*!
 * \file
 * \brief Test program for the LLR class
 * \author Erik G. Larsson
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

#include <itpp/itcomm.h>
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;

TEST(LLR, All)
{
  const double tol = 1e-6;

  LLR_calc_unit lcu1; // standard table resolution
  LLR_calc_unit lcu2(10, 7, 9);  // low table resolution
  LLR_calc_unit lcu3(2, 15, 0);  // low table resolution and low LLR granuality
  LLR_calc_unit lcu4(10, 0, 0);  // this gives logexp=logmax
  ivec ref = "12 300 7";
  ASSERT_TRUE(ref == lcu1.get_Dint());
  ref = "10 7 9";
  ASSERT_TRUE(ref == lcu2.get_Dint());
  ref = "2 15 0";
  ASSERT_TRUE(ref == lcu3.get_Dint());
  ref = "10 0 0";
  ASSERT_TRUE(ref == lcu4.get_Dint());

  // Testing Jacobian logarithm with four different resolutions
  mat ref_m = "0.693115 0.693359 0.75 0;"
            "0.647461 0.693359 0.75 0;"
            "0.60376 0.693359 0.5 0;"
            "0.562256 0.693359 0.5 0;"
            "0.523193 0.693359 0.5 0;"
            "0.474121 0.473633 0.5 0;"
            "0.439697 0.473633 0.5 0;"
            "0.407471 0.473633 0.5 0;"
            "0.376953 0.473633 0.5 0;"
            "0.348389 0.473633 0.25 0;"
            "0.313232 0.313477 0.25 0;"
            "0.288818 0.313477 0.25 0;"
            "0.266113 0.313477 0.25 0;"
            "0.245117 0.313477 0.25 0;"
            "0.225342 0.313477 0.25 0;"
            "0.201416 0.201172 0.25 0;"
            "0.185059 0.201172 0.25 0;"
            "0.169678 0.201172 0.25 0;"
            "0.155762 0.201172 0.25 0;"
            "0.142578 0.201172 0.25 0;"
            "0.126953 0.126953 0.25 0;"
            "0.116211 0.126953 0.25 0;"
            "0.106445 0.126953 0 0;"
            "0.097168 0.126953 0 0;"
            "0.0888672 0.126953 0 0;"
            "0.0788574 0.0791016 0 0;"
            "0.0720215 0.0791016 0 0;"
            "0.065918 0.0791016 0 0;"
            "0.0600586 0.0791016 0 0;"
            "0.0549316 0.0791016 0 0;"
            "0.048584 0.0488281 0 0;"
            "0.0444336 0.0488281 0 0;"
            "0.0405273 0.0488281 0 0;"
            "0.0368652 0.0488281 0 0;"
            "0.0336914 0.0488281 0 0;"
            "0.0297852 0 0 0;"
            "0.0270996 0 0 0;"
            "0.0246582 0 0 0;"
            "0.0224609 0 0 0;"
            "0.0205078 0 0 0;"
            "0.0180664 0 0 0;"
            "0.0166016 0 0 0;"
            "0.0151367 0 0 0;"
            "0.0136719 0 0 0;"
            "0.0124512 0 0 0;"
            "0.0109863 0 0 0;"
            "0.0100098 0 0 0;"
            "0.00927734 0 0 0;"
            "0.00830078 0 0 0;"
            "0.00756836 0 0 0;"
            "0.00683594 0 0 0;"
            "0.00610352 0 0 0;"
            "0.00561523 0 0 0;"
            "0.00512695 0 0 0;"
            "0.00463867 0 0 0;"
            "0.00415039 0 0 0;"
            "0.00366211 0 0 0;"
            "0.00341797 0 0 0;"
            "0.00317383 0 0 0;"
            "0.00268555 0 0 0;"
            "0.00244141 0 0 0;"
            "0.00219727 0 0 0;"
            "0.00195312 0 0 0;"
            "0.00195312 0 0 0;"
            "0.00170898 0 0 0;"
            "0.00146484 0 0 0;"
            "0.00146484 0 0 0;"
            "0.0012207 0 0 0;"
            "0.0012207 0 0 0;"
            "0.000976562 0 0 0;"
            "0.000976562 0 0 0;"
            "0.000732422 0 0 0;"
            "0.000732422 0 0 0;"
            "0.000732422 0 0 0;"
            "0.000732422 0 0 0;"
            "0.000488281 0 0 0;"
            "0.000488281 0 0 0;"
            "0.000488281 0 0 0;"
            "0.000488281 0 0 0;"
            "0.000488281 0 0 0;"
            "0.000244141 0 0 0;"
            "0.000244141 0 0 0;"
            "0.000244141 0 0 0;"
            "0.000244141 0 0 0;"
            "0.000244141 0 0 0;"
            "0.000244141 0 0 0;"
            "0.000244141 0 0 0;"
            "0.000244141 0 0 0;"
            "0.000244141 0 0 0;"
            "0.000244141 0 0 0;"
            "0.000244141 0 0 0;"
            "0 0 0 0;"
            "0 0 0 0;"
            "0 0 0 0;"
            "0 0 0 0;"
            "0 0 0 0;"
            "0 0 0 0;"
            "0 0 0 0;"
            "0 0 0 0;"
            "0 0 0 0;"
            "0 0 0 0";
  int n = 0;
  for (double x = 0.0; x < 10; x += 0.1, ++n) {
    ASSERT_NEAR(ref_m(n,0), lcu1.to_double(lcu1.logexp(lcu1.to_qllr(x))), tol);
    ASSERT_NEAR(ref_m(n,1), lcu2.to_double(lcu2.logexp(lcu2.to_qllr(x))), tol);
    ASSERT_NEAR(ref_m(n,2), lcu3.to_double(lcu3.logexp(lcu3.to_qllr(x))), tol);
    ASSERT_NEAR(ref_m(n,3), lcu4.to_double(lcu4.logexp(lcu4.to_qllr(x))), tol);
  }

  ASSERT_NEAR(0.75, lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(100.0), lcu1.to_qllr(0.75))), tol);
  ASSERT_NEAR(-0.75, lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(100.0), lcu1.to_qllr(-0.75))), tol);
  ASSERT_NEAR(-0.75, lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(-100.0), lcu1.to_qllr(0.75))), tol);
  ASSERT_NEAR(0.75, lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(-100.0), lcu1.to_qllr(-0.75))), tol);
  ASSERT_NEAR(0, lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(0.0), lcu1.to_qllr(0.75))), tol);
  ASSERT_NEAR(0, lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(0.0), lcu1.to_qllr(-0.75))), tol);
  ASSERT_NEAR(-1.177978515625, lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(3.75), lcu1.to_qllr(-1.25))), tol);
  ASSERT_NEAR(-1.177978515625, lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(-1.25), lcu1.to_qllr(3.75))), tol);
  ASSERT_NEAR(-1.177978515625, lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(-3.75), lcu1.to_qllr(1.25))), tol);
  ASSERT_NEAR(-1.177978515625, lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(1.25), lcu1.to_qllr(-3.75))), tol);
  ASSERT_NEAR(1.177978515625, lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(3.75), lcu1.to_qllr(1.25))), tol);
  ASSERT_NEAR(1.177978515625, lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(1.25), lcu1.to_qllr(3.75))), tol);
}

