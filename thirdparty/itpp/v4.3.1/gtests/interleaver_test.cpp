/*!
 * \file
 * \brief Interleaver classes test program
 * \author Pal Frenger and Adam Piatyszek
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
using namespace std;

static
void assert_ivec(const ivec &exp, const ivec &act)
{
  ASSERT_EQ(exp.length(), act.length());
  for (int n = 0; n < exp.length(); ++n)
  {
    ASSERT_EQ(exp[n], act[n]);
  }
}

TEST (Interleaver, All)
{
  //Declare scalars and vectors:
  int rows, cols, order, depth;
  ivec input, output, deinterleaved;
  ivec ref;

  RNG_reset(0);

  //Declare the interleavers.
  Block_Interleaver<int> block_interleaver;
  Cross_Interleaver<int> cross_interleaver;
  Sequence_Interleaver<int> sequence_interleaver;

  //Testing Block_Interleaver
  rows = 4;
  cols = 5;
  block_interleaver.set_rows(rows);
  block_interleaver.set_cols(cols);
  input = "1:20";
  output = block_interleaver.interleave(input);
  deinterleaved = block_interleaver.deinterleave(output);
  ref = "1 5 9 13 17 2 6 10 14 18 3 7 11 15 19 4 8 12 16 20";
  assert_ivec(ref, output);
  assert_ivec(input, deinterleaved);

  //Testing Cross_Interleaver
  order = 5;
  cross_interleaver.set_order(order);
  input = "1:25";
  output = cross_interleaver.interleave(input);
  deinterleaved = cross_interleaver.deinterleave(output);
  ref = "1 0 0 0 0 6 2 0 0 0 11 7 3 0 0 16 12 8 4 0 21 17 13 9 5 0 22 18 14 10 "
  "0 0 23 19 15 0 0 0 24 20 0 0 0 0 25 0 0 0 0 0";
  assert_ivec(ref, output);
  assert_ivec(input, deinterleaved);

  //Testing Sequence_Interleaver
  depth = 25;
  sequence_interleaver.set_interleaver_depth(depth);
  sequence_interleaver.randomize_interleaver_sequence();
  input = "1:25";
  output = sequence_interleaver.interleave(input);
  deinterleaved = sequence_interleaver.deinterleave(output);
  ref = "19 1 7 6 10 18 2 14 3 23 4 15 16 21 11 24 13 9 12 25 17 22 5 8 20";
  assert_ivec(ref, output);
  assert_ivec(input, deinterleaved);
}
