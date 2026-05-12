/*!
 * \file
 * \brief Polynomial routines test program
 * \author Tony Ottosson, Adam Piatyszek and Kumar Appaiah
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
  cout << "===================================" << endl;
  cout << "    Test of polynomial routines    " << endl;
  cout << "===================================" << endl;

  cout.setf(ios::fixed);
  cout.precision(6);

  {
    cout << "Real polynomials" << endl;
    vec r = randn(3);
    vec p = poly(r);
    cvec r2 = roots(p);
    cout << "Roots, r = " << r << endl;
    cout << "Polynomial, p = " << p << endl;
    cout << "r = roots(p) = " << r2 << endl;

    r = randn(7);
    p = poly(r);
    r2 = roots(p);
    cout << "Roots, r = " << r << endl;
    cout << "Polynomial, p = " << p << endl;
    cout << "r = roots(p) = " << r2 << endl;

    vec x = randn(9);
    vec y = polyval(p, x);
    cout << "x = " << x << endl;
    cout << "polyval(p, x) = " << y << endl;
  }

  {
    cout << "Complex polynomials" << endl;
    cvec r = randn_c(3);
    cvec p = poly(r);
    cvec r2 = roots(p);
    cout << "Roots, r = " << r << endl;
    cout << "Polynomial, p = " << p << endl;
    cout << "r = roots(p) = " << r2 << endl;

    r = randn_c(7);
    p = poly(r);
    r2 = roots(p);
    cout << "Roots, r = " << r << endl;
    cout << "Polynomial, p = " << p << endl;
    cout << "r = roots(p) = " << r2 << endl;

    cvec x = randn_c(9);
    cvec y = polyval(p, x);
    cout << "x = " << x << endl;
    cout << "polyval(p, x) = " << y << endl;
  }

  {
    cout << "Chebyshev polynomial" << endl;
    vec x = randn(8);
    cout << "x = " << x << endl;
    cout << "cheb(10, x) = " << cheb(10, x) << endl;
    cout << "cheb(15, x) = " << cheb(15, x) << endl;
  }

  return 0;
}
