/*!
 * \file
 * \brief BLAS aided routines test program
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

using namespace itpp;
using namespace std;

int main()
{
  cout.setf(ios::fixed);
  cout.precision(4);

  // dot() tests
  {
    vec a = randn(10);
    vec b = randn(10);
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "dot(a,b) = " << dot(a, b) << endl << endl;

    cvec c = randn_c(6);
    cvec d = randn_c(6);
    cout << "c = " << c << endl;
    cout << "d = " << d << endl;
    cout << "dot(c,d) = " << dot(c, d) << endl << endl;
  }

  // outer_product() tests
  {
    vec a = randn(4);
    vec b = randn(7);
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "outer_product(a,b) = " << outer_product(a, b) << endl;

    cvec c = randn_c(4);
    cvec d = randn_c(7);
    cout << "c = " << c << endl;
    cout << "d = " << d << endl;
    cout << "outer_product(c,d) = " << outer_product(c, d) << endl;
    cout << "outer_product(c,d,true) = " << outer_product(c, d, true) << endl << endl;
  }

  // Mat *= Mat operator test
  {
    mat M1 = randn(3, 5);
    mat N1 = randn(5, 2);
    cout << "M = " << M1 << endl;
    cout << "N = " << N1 << endl;
    M1 *= N1;
    cout << "M *= N;\nM = " << M1 << endl << endl;

    cmat M2 = randn_c(4, 4);
    cmat N2 = randn_c(4, 2);
    cout << "M = " << M2 << endl;
    cout << "N = " << N2 << endl;
    M2 *= N2;
    cout << "M *= N;\nM = " << M2 << endl << endl;
  }

  // Vec = Mat * Vec operator test
  {
    mat M1 = randn(3, 4);
    vec v1 = randn(4);
    cout << "M = " << M1 << endl;
    cout << "v = " << v1 << endl;
    cout << "out = M * v = " << M1 * v1 << endl << endl;

    cmat M2 = randn_c(3, 2);
    cvec v2 = randn_c(2);
    cout << "M = " << M2 << endl;
    cout << "v = " << v2 << endl;
    cout << "out = M * v = " << M2 * v2 << endl << endl;
  }

  return 0;
}
