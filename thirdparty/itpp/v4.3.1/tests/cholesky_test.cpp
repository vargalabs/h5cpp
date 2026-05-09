/*!
 * \file
 * \brief Cholesky factorisation test program
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
  cout << "    Test of Cholesky routines   " << endl;
  cout << "================================" << endl;

  {
    cout << "Real matrix" << endl;
    mat X, F;
    bool ok = false;

    while (!ok) {
      X = randn(5, 5);
      X = transpose(X) * X; // create a symmetric matrix
      // Make diagonal real and positive
      for (int i = 0; i < X.cols(); i++)
        X(i, i) = std::abs(X(i, i));

      ok = chol(X, F);
      cout << "X = " << round_to_zero(X) << endl;
      if (!ok)
        cout << "matrix is not positive definite" << endl;
      else {
        cout << "norm(X - F^T*F) = "
             << round_to_zero(norm(X - transpose(F) * F)) << endl;
      }
    }
  }

  {
    cout << endl << "Complex matrix" << endl;
    cmat X, F;
    bool ok = false;

    while (!ok) {
      X = randn_c(5, 5);
      X = hermitian_transpose(X) * X; // create a symmetric matrix
      // Make diagonal real and positive
      for (int i = 0; i < X.cols(); i++)
        X(i, i) = std::abs(real(X(i, i)))
                  ;
      ok = chol(X, F);
      cout << "X = " << round_to_zero(X) << endl;

      if (!ok)
        cout << "matrix is not positive definite" << endl;
      else {
        cout << "norm(X - F^H*F) = "
             << round_to_zero(norm(X - hermitian_transpose(F) * F)) << endl;
      }
    }
  }

  return 0;
}
