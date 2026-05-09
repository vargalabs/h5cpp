/*!
 * \file
 * \brief Schur decomposition test program
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

#include <itpp/itstat.h>

using namespace itpp;
using namespace std;


int main()
{
  int size = 5;
  double thres = 1e-13;

  cout << "==========================================" << endl;
  cout << "   Test of Schur decomposition routines   " << endl;
  cout << "==========================================" << endl;
  {
    cout << "Real matrix" << endl;
    mat A = randn(size, size);
    mat T, U;
    schur(A, U, T);

    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A - U*T*U^T) = "
         << round_to_zero(norm(A - (U * T * transpose(U))), thres) << endl;
    cout << "norm(I - U*U^T) = "
         << round_to_zero(norm(eye(size) - (U * transpose(U))), thres) << endl;
    double temp_sum = 0;
    for (int i = 2; i < size; i++)
      for (int j = 0; j < i - 1; j++)
        temp_sum += sqr(T(i, j));
    cout << "norm(lower triangular part of T) = "
         << round_to_zero(sqrt(temp_sum), thres) << endl;
  }
  {
    cout << endl << "Complex matrix" << endl;
    cmat A = randn_c(size, size);
    cmat T, U;
    schur(A, U, T);

    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A - U*T*U^H) = "
         << round_to_zero(norm(A - (U * T * hermitian_transpose(U))), thres)
         << endl;
    cout << "norm(I - U*U^H) = "
         << round_to_zero(norm(eye(size) - (U * hermitian_transpose(U))), thres)
         << endl;
    double temp_sum = 0;
    for (int i = 1; i < size; i++)
      for (int j = 0; j < i; j++)
        temp_sum += sqr(T(i, j));
    cout << "norm(lower triangular part of T) = "
         << round_to_zero(sqrt(temp_sum), thres) << endl;
  }

  return 0;
}
