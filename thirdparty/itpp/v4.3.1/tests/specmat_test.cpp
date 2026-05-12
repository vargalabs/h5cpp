/*!
 * \file
 * \brief Test program of special vectors and matrices
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
#include <iomanip>

using namespace std;
using namespace itpp;

int main()
{
  cout << "=================================" << endl;
  cout << "    Test of specmat routines" << endl;
  cout << "=================================" << endl;

  cout.setf(ios::fixed);
  cout.precision(3);

  bvec b1 = randb(7);
  bvec b2 = randb(13);
  cout << "b1 = " << b1 << endl;
  cout << "b2 = " << b2 << endl;
  cout << "toeplitz(b1, b2) =" << endl << toeplitz(b1, b2) << endl;
  cout << "toeplitz(b1) =" << endl << toeplitz(b1) << endl << endl;

  vec v1 = randn(5);
  vec v2 = randn(7);
  cout << "v1 = " << v1 << endl;
  cout << "v2 = " << v2 << endl;
  cout << "toeplitz(v1, v2) =" << endl << toeplitz(v1, v2) << endl;
  cout << "toeplitz(v1) =" << endl << toeplitz(v1) << endl << endl;

  cvec c1 = randn_c(4);
  cvec c2 = randn_c(3);
  cout << "c1 = " << c1 << endl;
  cout << "c2 = " << c2 << endl;
  cout << "toeplitz(c1, c2) =" << endl << toeplitz(c1, c2) << endl;
  cout << "toeplitz(c1, c1) =" << endl << toeplitz(c1, c1) << endl;
  cout << "toeplitz(c1) =" << endl << toeplitz(c1) << endl << endl;

  return 0;
}
