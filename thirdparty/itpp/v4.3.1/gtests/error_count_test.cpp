/*!
 * \file
 * \brief BERC and BLER error counters test program
 * \author Adam Piatyszek
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

using namespace itpp;
using namespace std;


TEST(ErrorCount, All)
{
  //BERC and BLER test

  const double eps = 1e-12;
  RNG_reset(0);

  const int block_size = 20;
  const ivec error_pos = "3 5 9 17";
  BERC berc;
  BLERC blerc(block_size);

  for (int i = 0; i < 100; ++i) {
    bvec input = randb(block_size);
    bvec output = input;
    // introduce some errors
    if (i < 80)
      for (int j = 0; j < error_pos.size(); ++j)
        output(error_pos(j)) = !output(error_pos(j));
    // extend the output vector by one bit
    output = concat(output, bin(1));
    // count errors
    berc.count(input, output);
    blerc.count(input, output);
  }

  ASSERT_NEAR(0.8, blerc.get_errorrate(), eps);
  ASSERT_NEAR(80, blerc.get_errors(), eps);
  ASSERT_NEAR(20, blerc.get_corrects(), eps);
  ASSERT_EQ(100, blerc.get_total_blocks());

  ASSERT_NEAR(0.16, berc.get_errorrate(), eps);
  ASSERT_NEAR(320, berc.get_errors(), eps);
  ASSERT_NEAR(1680, berc.get_corrects(), eps);
  ASSERT_EQ(2000, berc.get_total_bits());
}
