/*!
 * \file
 * \brief SVD decomposition test program
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
  cout << "================================" << endl;
  cout << "      Test of svd routines      " << endl;
  cout << "================================" << endl;

  {
    cout << "Real matrix" << endl;
    mat A = randn(5, 5);
    mat U, V;
    vec S;
    svd(A, U, S, V);

    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A - U*diag(S)*V^T) = "
         << round_to_zero(norm(A - U * diag(S) * transpose(V))) << endl;

  }
  {
    cout << endl << "Complex matrix" << endl;
    cmat A = randn_c(5, 5);
    cmat U, V;
    vec S;
    svd(A, U, S, V);

    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A - U*diag(S)*V^H) = "
         << round_to_zero(norm(A - U * diag(S) * hermitian_transpose(V)))
         << endl;
  }

  return 0;
}
