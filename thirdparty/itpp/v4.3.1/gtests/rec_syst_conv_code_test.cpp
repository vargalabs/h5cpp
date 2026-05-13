/*!
 * \file
 * \brief Recursive systematic convolutional codes class test program
 * \author Pal Frenger
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

#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include "gtest/gtest.h"

using namespace itpp;

TEST(RSCC, All)
{
  RNG_reset(0);

  bvec uncoded_bits = "0 1 1 0 1 0 1", tail_bits(4), decoded_bits(11);
  bmat parity_bits(11, 1);
  vec received_systematic(11), extrinsic_output(11), L(11), symbols;
  mat received_parity(11, 1);

  double Ec = 1.0, N0 = 1e-10, Lc = 4.0 * sqrt(Ec) / N0;
  Normal_RNG noise(0.0, N0 / 2.0);

  BPSK bpsk;
  Rec_Syst_Conv_Code rscc;
  // Fixed by ediap: ivec("string") constructor uses only decimal notation
  // rscc.set_generator_polynomials("031 027",5);
  rscc.set_generator_polynomials("25 23", 5);
  rscc.set_awgn_channel_parameters(Ec, N0);

  rscc.encode_tail(uncoded_bits, tail_bits, parity_bits);
  bpsk.modulate_bits(concat(uncoded_bits, tail_bits), symbols);
  received_systematic = symbols + noise(11);

  bpsk.modulate_bits(parity_bits.get_col(0), symbols);
  received_parity.set_col(0, symbols + noise(11));

  vec extrinsic_input = zeros(11);
  rscc.map_decode(received_systematic, received_parity, extrinsic_input,
                  extrinsic_output);

  L = Lc * received_systematic + extrinsic_output;
  for (int k = 0; k < 11; k++) {
    (L(k) > 0) ? (decoded_bits(k) = bin(0)) : (decoded_bits(k) = bin(1));
  }
  bvec ref = concat(uncoded_bits, tail_bits);
  ASSERT_TRUE(ref == decoded_bits);
}
