/*!
 * \file
 * \brief Windowing functions test program
 * \author Tony Ottosson, Tobias Ringstrom, Pal Frenger, Adam Piatyszek
 *         and Kumar Appaiah
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

#include <itpp/itsignal.h>
#include <iomanip>

using namespace itpp;
using namespace std;

int main()
{
  // This is a hack for improper rounding under MinGW
  cout.setf(ios::fixed);
  cout.precision(8);

  cout << "================================" << endl;
  cout << "    Test of window functions    " << endl;
  cout << "================================" << endl;

  cout << "hamming(32) = " << hamming(32) << endl;
  cout << "hamming(128) = " << hamming(128) << endl;

  cout << "hanning(32) = " << hanning(32) << endl;
  cout << "hanning(128) = " << hanning(128) << endl;

  cout << "hann(32) = " << hann(32) << endl;
  cout << "hann(128) = " << hann(128) << endl;

  cout << "blackman(32) = " << blackman(32) << endl;
  cout << "blackman(128) = " << blackman(128) << endl;

  cout << "triang(32) = " << triang(32) << endl;
  cout << "triang(128) = " << triang(128) << endl;

  cout << "chebwin(32, 50) = " << chebwin(32, 50) << endl;
  cout << "chebwin(33, 20) = " << chebwin(33, 20) << endl;
  cout << "chebwin(127, 25) = " << chebwin(127, 25) << endl;
  cout << "chebwin(128, 25) = " << chebwin(128, 25) << endl;

  return 0;
}
