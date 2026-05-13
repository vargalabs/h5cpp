/*!
 * \file
 * \brief Linear systems of equations solving test program
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

#include <itpp/itbase.h>

using namespace itpp;
using namespace std;


int main()
{
  cout << "==========================================" << endl;
  cout << "   Solving linear systems of equations" << endl;
  cout << "==========================================" << endl;

  {
    cout << "Real systems:" << endl << endl;
    mat A, B, X;
    vec b, x;

    A = randn(4, 4);
    b = randn(4);
    x = ls_solve(A, b);

    cout << "Square system: Ax=b" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "b=" << b << endl
         << "x=" << x << endl << endl;

    A = randn(4, 4);
    B = randn(4, 2);
    X = ls_solve(A, B);

    cout << "Square system: AX=B" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "B=" << B << endl
         << "X=" << X << endl << endl;



    A = randn(4, 4);
    A = A.transpose() * A;
    b = randn(4);
    x = ls_solve(A, b);

    cout << "Square system (chol): Ax=b" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "b=" << b << endl
         << "x=" << x << endl << endl;


    A = randn(4, 4);
    A = A.transpose() * A;
    B = randn(4, 2);
    X = ls_solve(A, B);

    cout << "Square system (Chol): AX=B" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "B=" << B << endl
         << "X=" << X << endl << endl;



    A = randn(4, 2);
    b = randn(4);
    x = ls_solve_od(A, b);

    cout << "Overdetermined system: Ax=b" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "b=" << b << endl
         << "x=" << x << endl << endl;

    A = randn(4, 2);
    B = randn(4, 3);
    X = ls_solve_od(A, B);

    cout << "Overdetermined system: AX=B" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "B=" << B << endl
         << "X=" << X << endl << endl;

    A = randn(2, 4);
    b = randn(2);
    x = ls_solve_ud(A, b);

    cout << "Underdetermined system: Ax=b" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "b=" << b << endl
         << "x=" << x << endl << endl;

    A = randn(2, 4);
    B = randn(2, 3);
    X = ls_solve_ud(A, B);

    cout << "Underdetermined system: AX=B" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "B=" << B << endl
         << "X=" << X << endl << endl;

  }


  {
    cout << "Complex systems:" << endl << endl;
    cmat A, B, X;
    cvec b, x;

    A = randn_c(4, 4);
    b = randn_c(4);
    x = ls_solve(A, b);

    cout << "Square system: Ax=b" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "b=" << b << endl
         << "x=" << x << endl << endl;

    A = randn_c(4, 4);
    B = randn_c(4, 2);
    X = ls_solve(A, B);

    cout << "Square system: AX=B" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "B=" << B << endl
         << "X=" << X << endl << endl;



    A = randn_c(4, 4);
    A = A.transpose() * A;
    b = randn_c(4);
    x = ls_solve(A, b);

    cout << "Square system (chol): Ax=b" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "b=" << b << endl
         << "x=" << x << endl << endl;


    A = randn_c(4, 4);
    A = A.transpose() * A;
    B = randn_c(4, 2);
    X = ls_solve(A, B);

    cout << "Square system (Chol): AX=B" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "B=" << B << endl
         << "X=" << X << endl << endl;



    A = randn_c(4, 2);
    b = randn_c(4);
    x = ls_solve_od(A, b);

    cout << "Overdetermined system: Ax=b" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "b=" << b << endl
         << "x=" << x << endl << endl;

    A = randn_c(4, 2);
    B = randn_c(4, 3);
    X = ls_solve_od(A, B);

    cout << "Overdetermined system: AX=B" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "B=" << B << endl
         << "X=" << X << endl << endl;

    A = randn_c(2, 4);
    b = randn_c(2);
    x = ls_solve_ud(A, b);

    cout << "Underdetermined system: Ax=b" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "b=" << b << endl
         << "x=" << x << endl << endl;

    A = randn_c(2, 4);
    B = randn_c(2, 3);
    X = ls_solve_ud(A, B);

    cout << "Underdetermined system: AX=B" << endl
         << "============================" << endl
         << "A=" << A << endl
         << "B=" << B << endl
         << "X=" << X << endl << endl;

  }

  return 0;
}
