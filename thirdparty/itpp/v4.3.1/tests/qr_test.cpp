/*!
 * \file
 * \brief QR factorisation test program
 * \author Tony Ottosson, Adam Piatyszek and Vasek Smidl
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
  cout << "=================================" << endl;
  cout << "Test of QR factorization routines" << endl;
  cout << "=================================" << endl;

  {
    cout << "QR of Real matrix" << endl;
    mat Q, R, e;
    mat A = randn(5, 5);
    qr(A, Q, R);
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A - Q*R) = " << round_to_zero(norm(A - Q * R)) << endl
         << endl;

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn(4, 2);
    qr(A, Q, R);
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A - Q*R) = " << round_to_zero(norm(A - Q * R)) << endl
         << endl;

    A = randn(2, 4);
    qr(A, Q, R);
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A - Q*R) = " << round_to_zero(norm(A - Q * R)) << endl
         << endl;
  }

  {
    cout << "QR of Real matrix without Q" << endl;
    mat R;
    mat A = randn(5, 5);
    qr(A, R);
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A.T()*A - R.T()*R) = " << round_to_zero(norm(A.T()*A -  R.T()*R)) << endl
         << endl;

    A = randn(4, 2);
    qr(A, R);
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A.T()*A - R.T()*R) = " << round_to_zero(norm(A.T()*A -  R.T()*R)) << endl
         << endl;

    A = randn(2, 4);
    qr(A, R);
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A.T()*A - R.T()*R) = " << round_to_zero(norm(A.T()*A -  R.T()*R)) << endl
         << endl;
  }

  {
    cout << "QR of Real matrix with pivoting" << endl;
    mat Q, R, e;
    bmat P;
    mat A = randn(5, 5);
    qr(A, Q, R, P);
    e = A * to_mat(P) - Q * R;
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A*mat(P) - Q*R) = " << round_to_zero(norm(e)) << endl << endl;

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn(4, 2);
    qr(A, Q, R, P);
    e = A * to_mat(P) - Q * R;
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A*mat(P) - Q*R) = " << round_to_zero(norm(e)) << endl << endl;

    A = randn(2, 4);
    qr(A, Q, R, P);
    e = A * to_mat(P) - Q * R;
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A*mat(P) - Q*R) = " << round_to_zero(norm(e)) << endl << endl;
  }

  {
    cout << "QR of Complex matrix" << endl;
    cmat A = randn_c(5, 5);
    cmat Q, R, e;

    qr(A, Q, R);
    e = A - Q * R;
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A - Q*R) = " << round_to_zero(norm(e)) << endl << endl;

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn_c(4, 2);
    qr(A, Q, R);
    e = A - Q * R;
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A - Q*R) = " << round_to_zero(norm(e)) << endl << endl;

    A = randn_c(2, 4);
    qr(A, Q, R);
    e = A - Q * R;
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A - Q*R) = " << round_to_zero(norm(e)) << endl << endl;
  }

  {
    cout << "QR of Complex matrix without Q" << endl;
    cmat A = randn_c(5, 5);
    cmat R, e;

    qr(A, R);
    e = A.H() * A - R.H() * R;
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A.H()*A - R.H()*R) = " << round_to_zero(norm(e)) << endl << endl;

    A = randn_c(4, 2);
    qr(A, R);
    e = A.H() * A - R.H() * R;
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A.H()*A - R.H()*R) = " << round_to_zero(norm(e)) << endl << endl;

    A = randn_c(2, 4);
    qr(A, R);
    e = A.H() * A - R.H() * R;
    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A.H()*A - R.H()*R) = " << round_to_zero(norm(e)) << endl << endl;
  }

  {
    cout << "QR of Complex matrix with pivoting" << endl;
    cmat A = randn_c(5, 5);
    cmat Q, R, e;
    bmat P;

    qr(A, Q, R, P);
    e = A * to_mat(P) - Q * R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A*mat(P) - Q*R) = " << round_to_zero(norm(e)) << endl << endl;

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn_c(4, 2);
    qr(A, Q, R, P);
    e = A * to_mat(P) - Q * R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A*mat(P) - Q*R) = " << round_to_zero(norm(e)) << endl << endl;

    A = randn_c(2, 4);
    qr(A, Q, R, P);
    e = A * to_mat(P) - Q * R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A*mat(P) - Q*R) = " << round_to_zero(norm(e)) << endl << endl;
  }

  return 0;
}
