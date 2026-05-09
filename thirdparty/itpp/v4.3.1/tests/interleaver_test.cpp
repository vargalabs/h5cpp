/*!
 * \file
 * \brief Interleaver classes test program
 * \author Pal Frenger and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
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

using namespace itpp;
using namespace std;

int main()
{
  //Declare scalars and vectors:
  int rows, cols, order, depth;
  ivec input, output, deinterleaved;

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
  cout << "Testing Block_Interleaver:" << endl;
  cout << "input = " << input << endl;
  cout << "output = " << output << endl;
  cout << "deinterleaved = " << deinterleaved << endl;
  cout << "===============================================================" << endl;

  //Testing Cross_Interleaver
  order = 5;
  cross_interleaver.set_order(order);
  input = "1:25";
  output = cross_interleaver.interleave(input);
  deinterleaved = cross_interleaver.deinterleave(output);
  cout << "Testing Cross_Interleaver:" << endl;
  cout << "input = " << input << endl;
  cout << "output = " << output << endl;
  cout << "deinterleaved = " << deinterleaved << endl;
  cout << "===============================================================" << endl;

  //Testing Sequence_Interleaver
  depth = 25;
  sequence_interleaver.set_interleaver_depth(depth);
  sequence_interleaver.randomize_interleaver_sequence();
  input = "1:25";
  output = sequence_interleaver.interleave(input);
  deinterleaved = sequence_interleaver.deinterleave(output);
  cout << "Testing Sequence_Interleaver:" << endl;
  cout << "input = " << input << endl;
  cout << "output = " << output << endl;
  cout << "deinterleaved = " << deinterleaved << endl;
  cout << "===============================================================" << endl;

  //Exit program:
  return 0;
}
