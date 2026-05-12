/*!
 * \file
 * \brief LDPC class test program
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

TEST(LDPC, All)
{
  RNG_reset(0);
  int n;

  LDPC_Parity_Regular H;
  H.generate(200, 3, 6, "rand", "100 6");
  int girth = H.cycle_removal_MGW(6);
  ASSERT_EQ(6, girth);
  LDPC_Generator_Systematic G;
  G.construct(&H);
  LDPC_Code C(&H, &G);
  C.save_code("ldpc_test.codec");
  LDPC_Code C1("ldpc_test.codec", &G);
  ASSERT_EQ(200, C.get_nvar());
  ASSERT_EQ(100, C.get_ncheck());
  ASSERT_DOUBLE_EQ(0.5, C.get_rate());
  bvec bitsin = randb(C.get_nvar() - C.get_ncheck());
  bvec bitsout;
  C.encode(bitsin, bitsout);
  ASSERT_TRUE(C.syndrome_check(bitsout)) << "syndrome check failed";

  double EbN0db = 1.5;
  double N0 = pow(10.0, -EbN0db / 10.0) / C.get_rate();
  double sigma = sqrt(N0 / 2.0);
  vec x = 1.0 + sigma * randn(C.get_nvar());
  QLLRvec LLRin = C.get_llrcalc().to_qllr(2.0 * x / (N0 / 2.0));
  QLLRvec LLRout(C.get_nvar());
  C.bp_decode(LLRin, LLRout);
  QLLRvec LLRout_ref = "27277 29922 2454 37509 36082 1560 25684 39929 29353 29133 26578 42497 39556 22777 "
  "47940 32903 35425 17665 35097 33419 18795 31183 30935 27616 20303";
  QLLRvec LLRout_act = LLRout.left(25);
  for (n = 0; n < LLRout_act.length(); ++n) {
    ASSERT_EQ(LLRout_ref(n), LLRout_act(n));
  }

  // BLDPC code
  {
    imat A = "0 -1 -1 0; -1 1 4 -1; -1 2 -1 6";
    imat B = "1; -1; 2";
    imat T = "0 -1 -1; 0 0 -1; -1 0 0";
    imat C = "3 -1 5 -1";
    imat D = "3";
    imat E = "-1 -1 0";
    // base matrix
    imat H_b = concat_vertical(concat_horizontal(concat_horizontal(A, B), T),
                               concat_horizontal(concat_horizontal(C, D), E));
    int Z = 4; // expansion factor

    BLDPC_Parity H(H_b, Z);

    BLDPC_Generator G(&H);
    bvec in_bits = randb(H.get_nvar() - H.get_ncheck());
    bvec codeword;
    G.encode(in_bits, codeword);
    bvec cw_ref = "0 0 0 1 0 0 1 0 0 1 0 0 0 1 0 1 1 0 0 1 0 1 1 1 0 1 1 1 1 1 0 0";
    for (n = 0; n < cw_ref.length(); ++n) {
      ASSERT_EQ(cw_ref(n), codeword(n));
    }
  }
}
