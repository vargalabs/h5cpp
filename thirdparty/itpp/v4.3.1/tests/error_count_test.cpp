/*!
 * \file
 * \brief BERC and BLER error counters test program
 * \author Adam Piatyszek
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

#include <itpp/itcomm.h>

using namespace itpp;
using namespace std;


int main()
{
  cout << "==========================" << endl;
  cout << "    BERC and BLER test    " << endl;
  cout << "==========================" << endl << endl;

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
  cout << "BLER = " << blerc.get_errorrate() << endl
       << "Block errors = " << blerc.get_errors() << endl
       << "Correct blocks = " << blerc.get_corrects() << endl
       << "Total blocks = " << blerc.get_total_blocks() << endl << endl;
  blerc.clear();

  cout << "BER = " << berc.get_errorrate() << endl
       << "Bit errors = " << berc.get_errors() << endl
       << "Correct bits = " << berc.get_corrects() << endl
       << "Total bits = " << berc.get_total_bits() << endl;

  return 0;
}
