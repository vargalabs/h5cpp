/*!
 * \file
 * \brief Frequency filter test program
 * \author Simon Wood and Adam Piatyszek
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

using namespace itpp;
using namespace std;


int main()
{
  vec b = "1 2 3 4";
  vec x(20);
  x.zeros();
  x(0) = 1;

  // Define a filter object for doubles
  Freq_Filt<double> FF(b, x.length());

  // Filter the data
  vec y = FF.filter(x);

  // Check the FFT and block sizes that were used
  int fftsize = FF.get_fft_size();
  int blksize = FF.get_blk_size();

  cout << fftsize << endl;
  cout << blksize << endl;

  cout << round_to_zero(y) << endl;

  // Test streaming mode
  x = linspace(0, 10, 100);
  Freq_Filt<double> FFS(b, x.length());
  vec y1 = FFS.filter(x(0, 49), 1);
  vec y2 = FFS.filter(x(50, 99), 1);

  cout << round_to_zero(concat(y1, y2)) << endl;

  return 0;
}
