/*!
 * \file
 * \brief Determinant routines test program
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

#include <itpp/itstat.h>

using namespace itpp;
using namespace std;


int main()
{
  cout << "===================================" << endl;
  cout << "    Test of Determinant routines   " << endl;
  cout << "===================================" << endl;

  {
    cout << "Real matrix" << endl;
    mat X = randn(5, 5);
    double d;
    d = det(X);
    cout << "X = " << round_to_zero(X) << endl;
    cout << "det(X) = " << round_to_zero(d) << endl;

    X = randn(5, 5);
    d = det(X);
    cout << "X = " << round_to_zero(X) << endl;
    cout << "det(X) = " << round_to_zero(d) << endl;
  }

  {
    cout << endl << "Complex matrix" << endl;
    cmat X = randn_c(5, 5);
    complex<double> d;
    d = det(X);
    cout << "X = " << round_to_zero(X) << endl;
    cout << "det(X) = " << round_to_zero(d) << endl;

    X = randn_c(5, 5);
    d = det(X);
    cout << "X = " << round_to_zero(X) << endl;
    cout << "det(X) = " << round_to_zero(d) << endl;
  }

  return 0;
}
