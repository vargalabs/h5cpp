/*!
 * \file
 * \brief Filter design test program
 * \author Tony Ottosson and Adam Piatyszek
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
  cout << "====================================" << endl;
  cout << "   Test of filter design routines   " << endl;
  cout << "====================================" << endl;

  {
    cout << "Stabilisation of real filters" << endl;
    vec p = "0.7 3.0 -0.4";
    vec p2 = polystab(p);
    cout << "Polynomial, p = " << p << endl;
    cout << "p2 = polystab(p) = " << round_to_zero(p2) << endl;

    p = randn(7);
    p2 = polystab(p);
    cout << "Polynomial, p = " << p << endl;
    cout << "p2 = polystab(p) = " << round_to_zero(p2) << endl;

    p = poly(vec("1.1 0.7 0.2"));
    p2 = polystab(p);
    cout << "Polynomial, p = " << p << endl;
    cout << "p2 = polystab(p) = " << round_to_zero(p2) << endl;
  }

  {
    cout << "Stabilisation of complex filters" << endl;
    cvec p = "0.7 3.0 -0.4";
    cvec p2 = polystab(p);
    cout << "Polynomial, p = " << p << endl;
    cout << "p2 = polystab(p) = " << round_to_zero(p2) << endl;

    p = randn_c(7);
    p2 = polystab(p);
    cout << "Polynomial, p = " << p << endl;
    cout << "p2 = polystab(p) = " << round_to_zero(p2) << endl;

    p = poly(cvec("1.1 0.7 0.2"));
    p2 = polystab(p);
    cout << "Polynomial, p = " << p << endl;
    cout << "p2 = polystab(p) = " << round_to_zero(p2) << endl;

    cvec a = randn_c(4);
    cvec b = randn_c(6);
    cvec h = freqz(b, a, 32);

    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "freqz(b,a,32) = " << round_to_zero(h) << endl;
  }
  {
    cout << "Yulewalk filter design" << endl;
    vec f = "0 0.5 0.6 1";
    vec m = "1 1 0 0";
    vec a, b, R;

    cout << "f = " << f << endl;
    cout << "m = " << m << endl;
    cout << "filter_design_autocorrelation(32, f, m, R): " << endl;
    filter_design_autocorrelation(256, f, m, R);

    cout << "R = " << fixed << R << endl;

    cout << "arma_estimator(8, 8, R, b, a): " << endl;
    arma_estimator(8, 8, R, b, a);
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;

    vec n = "0:1:256";
    double fd = 0.1;
    R = besselj(0, 2 * pi * fd * n);
    cout << "R = " << R << endl;
    arma_estimator(8, 8, R, b, a);
    cout << "arma_estimator(8, 8, R, b, a): " << endl;
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;

  }

  return 0;
}
