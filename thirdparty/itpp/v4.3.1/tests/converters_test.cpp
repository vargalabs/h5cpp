/*!
 * \file
 * \brief Tests of miscellaneous conversion functions
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

#include <itpp/itbase.h>

using namespace itpp;
using namespace std;


int main()
{
  cout << "==================================" << endl
       << "   Test of conversion functions" << endl
       << "==================================" << endl;

  vec v = "-4:0.5:5.5";
  v = concat(v, randn(10) * 20 - 10);

  cout << "v = " << v << endl << endl;

  cout << "round(v) = " << round_to_zero(round(v)) << endl;
  cout << "round_i(v) = " << round_i(v) << endl << endl;

  cout << "ceil(v) = " << round_to_zero(ceil(v)) << endl;
  cout << "ceil_i(v) = " << ceil_i(v) << endl << endl;

  cout << "floor(v) = " << round_to_zero(floor(v)) << endl;
  cout << "floor_i(v) = " << floor_i(v) << endl << endl;

  bvec b = randb(15);

  cout << "b = " << b << endl << endl;

  int i = bin2dec(b);
  cout << "bin2dec(b) = " << i << endl;
  cout << "dec2bin(bin2dec(b)) = " << dec2bin(i) << endl;
  i = bin2dec(b, false);
  cout << "bin2dec(b, false) = " << i << endl;
  cout << "dec2bin(bin2dec(b, false), false) = " << dec2bin(i, false)
       << endl << endl;

  ivec iv = bin2oct(b);
  cout << "bin2oct(b) = " << iv << endl;
  cout << "oct2bin(bin2oct(b)) = " << oct2bin(iv) << endl;
  cout << "oct2bin(bin2oct(b), 1) = " << oct2bin(iv, 1) << endl << endl;

  iv = bin2pol(b);
  cout << "bin2pol(b) = " << iv << endl;
  cout << "pol2bin(bin2pol(b)) = " << pol2bin(iv) << endl << endl;

  return 0;
}
