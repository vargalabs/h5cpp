/*!
 * \file
 * \brief Transforms test program
 * \author Tony Ottosson, Thomas Eriksson, Simon Wood and Adam Piatyszek
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
  cout << "========================" << endl;
  cout << "   Test of Transforms   " << endl;
  cout << "========================" << endl << endl;

  int N, s;
  {
    vec x, z;
    cvec y;

    N = 16;
    x = randn(N);
    cout << "Test 1: FFT/IFFT; Real input vector x; N = " << N << endl
         << "        fft_real(x, y), ifft_real(y, z):" << endl << endl;

    cout << "x = " << round_to_zero(x) << endl;
    fft_real(x, y);
    cout << "y = " << round_to_zero(y) << endl;
    ifft_real(y, z);
    cout << "z = " << round_to_zero(z) << endl << endl;

    N = 15;
    s = N - 4;
    x = randn(s);
    cout << "Test 2: FFT/IFFT; Real input vector x of size s = " << s
         << "; N = " << N << endl
         << "        y = fft_real(x, N), z = ifft_real(y, N):" << endl << endl;

    cout << "x = " << round_to_zero(x) << endl;
    y = fft_real(x, N);
    cout << "y = " << round_to_zero(y) << endl;
    z = ifft_real(y, N);
    cout << "z = " << round_to_zero(z.left(s)) << endl << endl;
  }
  {
    cvec x, y, z;

    N = 16;
    x = randn_c(N);
    cout << "Test 3: FFT/IFFT; Complex input vector x; N = " << N << endl
         << "        fft(x, y), ifft(y, z):" << endl << endl;

    cout << "x = " << round_to_zero(x) << endl;
    fft(x, y);
    cout << "y = " << round_to_zero(y) << endl;
    ifft(y, z);
    cout << "z = " << round_to_zero(z) << endl << endl;

    N = 16;
    s = N - 7;
    x = randn_c(s);
    cout << "Test 4: FFT/IFFT; Complex input vector x of size s = " << s
         << "; N = " << N << endl
         << "        y = fft(x, N), z = ifft(y, N):" << endl << endl;

    cout << "x = " << round_to_zero(x) << endl;
    y = fft(x, N);
    cout << "y = " << round_to_zero(y) << endl;
    z = ifft(y, N);
    cout << "z = " << round_to_zero(z.left(s)) << endl << endl;
  }
  {
    vec x, y, z;

    N = 8;
    x = randn(N);
    cout << "Test 5: DCT/IDCT; Real input vector; N = " << N << endl
         << "        dct(x, y), idct(y, z):" << endl << endl;

    cout << "x = " << round_to_zero(x) << endl;
    dct(x, y);
    cout << "y = " << round_to_zero(y) << endl;
    idct(y, z);
    cout << "z = " << round_to_zero(z) << endl << endl;

    N = 11;
    x = randn(N);
    cout << "Test 6: DCT/IDCT; Real input vector; N = " << N << endl
         << "        dct(x, y), idct(y, z):" << endl << endl;

    cout << "x = " << round_to_zero(x) << endl;
    dct(x, y);
    cout << "y = " << round_to_zero(y) << endl;
    idct(y, z);
    cout << "z = " << round_to_zero(z) << endl;
  }

  return 0;
}
