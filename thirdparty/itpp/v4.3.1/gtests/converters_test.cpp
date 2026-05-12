/*!
 * \file
 * \brief Tests of miscellaneous conversion functions
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

#include <itpp/itbase.h>
#include "gtest/gtest.h"

using namespace itpp;
using namespace std;


TEST (Converters, All)
{
  vec v1 = "-4:0.5:5.5";
  vec v2 = randn(10) * 20 - 10;
  vec v = concat(v1, v2);

  ASSERT_EQ(v1.length()+v2.length(), v.length());
  ASSERT_TRUE(v1 == v.left(v1.length()));
  ASSERT_TRUE(v2 == v.right(v2.length()));

  vec expect = round_to_zero(round(v));
  vec actual = to_vec(round_i(v));
  ASSERT_TRUE(expect == actual);
  expect = round_to_zero(ceil(v));
  actual = to_vec(ceil_i(v));
  ASSERT_TRUE(expect == actual);
  expect = round_to_zero(floor(v));
  actual = to_vec(floor_i(v));
  ASSERT_TRUE(expect == actual);

  RNG_reset(0);
  bvec b = randb(15);

  int i = bin2dec(b);
  ASSERT_EQ(31523, i);
  ASSERT_TRUE(b == dec2bin(i));
  i = bin2dec(b, false);
  ASSERT_EQ(25199, i);
  ASSERT_TRUE(b == dec2bin(i, false));

  ivec iv = bin2oct(b);
  ivec expect_i = "7 5 4 4 3";
  ASSERT_TRUE(expect_i == iv);
  ASSERT_TRUE(b == oct2bin(iv));
  ASSERT_TRUE(b == oct2bin(iv, 1));

  iv = bin2pol(b);
  expect_i = 1-2*to_ivec(b);
  ASSERT_TRUE(expect_i == iv);
  ASSERT_TRUE(b == pol2bin(iv));
}
