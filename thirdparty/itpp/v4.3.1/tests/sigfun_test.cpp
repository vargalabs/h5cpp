/*!
 * \file
 * \brief Signal processing routines test program
 * \author Tony Ottosson, Thomas Eriksson, Pal Frenger, Tobias Ringstrom
 *         and Adam Piatyszek
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
  RNG_reset(12345);

  vec x = randn(10);
  vec y = randn(10);
  int max_lag;

  cout << "x = " << x << endl;
  cout << "y = " << y << endl;

  cout << "Testing cross correlation" << endl;
  cout << "--------------------------------------------------------------------------" << endl;

  cout << "xcorr(x,y) = " << xcorr(x, y) << endl;
  max_lag = 0;
  cout << "xcorr(x,y," << max_lag << ") = " << xcorr(x, y, max_lag) << endl;
  max_lag = 5;
  cout << "xcorr(x,y," << max_lag << ") = " << xcorr(x, y, max_lag) << endl;
  max_lag = -1;
  cout << "xcorr(x,y," << max_lag << ") = " << xcorr(x, y, max_lag) << endl;

  max_lag = 0;
  cout << "xcorr(x,y," << max_lag << ",biased) = " << xcorr(x, y, max_lag, "biased") << endl;
  max_lag = 5;
  cout << "xcorr(x,y," << max_lag << ",biased) = " << xcorr(x, y, max_lag, "biased") << endl;
  max_lag = -1;
  cout << "xcorr(x,y," << max_lag << ",biased) = " << xcorr(x, y, max_lag, "biased") << endl;

  max_lag = 0;
  cout << "xcorr(x,y," << max_lag << ",unbiased) = " << xcorr(x, y, max_lag, "unbiased") << endl;
  max_lag = 5;
  cout << "xcorr(x,y," << max_lag << ",unbiased) = " << xcorr(x, y, max_lag, "unbiased") << endl;
  max_lag = -1;
  cout << "xcorr(x,y," << max_lag << ",unbiased) = " << xcorr(x, y, max_lag, "unbiased") << endl;

  max_lag = 0;
  cout << "xcorr(x,y," << max_lag << ",coeff) = " << xcorr(x, y, max_lag, "coeff") << endl;
  max_lag = 5;
  cout << "xcorr(x,y," << max_lag << ",coeff) = " << xcorr(x, y, max_lag, "coeff") << endl;
  max_lag = -1;
  cout << "xcorr(x,y," << max_lag << ",coeff) = " << xcorr(x, y, max_lag, "coeff") << endl;

  cout << "Testing auto correlation" << endl;
  cout << "--------------------------------------------------------------------------" << endl;

  cout << "xcorr(x) = " << xcorr(x) << endl;
  max_lag = 0;
  cout << "xcorr(x," << max_lag << ") = " << xcorr(x, max_lag) << endl;
  max_lag = 5;
  cout << "xcorr(x," << max_lag << ") = " << xcorr(x, max_lag) << endl;
  max_lag = -1;
  cout << "xcorr(x," << max_lag << ") = " << xcorr(x, max_lag) << endl;

  max_lag = 0;
  cout << "xcorr(x," << max_lag << ",biased) = " << xcorr(x, max_lag, "biased") << endl;
  max_lag = 5;
  cout << "xcorr(x," << max_lag << ",biased) = " << xcorr(x, max_lag, "biased") << endl;
  max_lag = -1;
  cout << "xcorr(x," << max_lag << ",biased) = " << xcorr(x, max_lag, "biased") << endl;

  max_lag = 0;
  cout << "xcorr(x," << max_lag << ",unbiased) = " << xcorr(x, max_lag, "unbiased") << endl;
  max_lag = 5;
  cout << "xcorr(x," << max_lag << ",unbiased) = " << xcorr(x, max_lag, "unbiased") << endl;
  max_lag = -1;
  cout << "xcorr(x," << max_lag << ",unbiased) = " << xcorr(x, max_lag, "unbiased") << endl;

  max_lag = 0;
  cout << "xcorr(x," << max_lag << ",coeff) = " << xcorr(x, max_lag, "coeff") << endl;
  max_lag = 5;
  cout << "xcorr(x," << max_lag << ",coeff) = " << xcorr(x, max_lag, "coeff") << endl;
  max_lag = -1;
  cout << "xcorr(x," << max_lag << ",coeff) = " << xcorr(x, max_lag, "coeff") << endl;

  return 0;
}
