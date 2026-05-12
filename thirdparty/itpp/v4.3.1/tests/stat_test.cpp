/*!
 * \file
 * \brief Miscellaneous statistical routines test program
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
#include <iomanip>

using namespace itpp;
using namespace std;


int main()
{
  cout << "=================================" << endl;
  cout << "  Test of statistical routines   " << endl;
  cout << "=================================" << endl;

  vec a = randn(5);

  cout << "a = " << a << endl << endl;

  cout << "max(a) = " << max(a) << endl;
  cout << "min(a) = " << min(a) << endl;
  cout << "mean(a) = " << mean(a) << endl;
  cout << "geometric_mean(abs(a)) = " << geometric_mean(abs(a)) << endl;
  cout << "norm(a) = " << norm(a) << endl;
  cout << "norm(a, 2) = " << norm(a, 2) << endl;
  cout << "norm(a, 1) = " << norm(a, 1) << endl;
  cout << "norm(a, \"fro\") = " << norm(a, "fro") << endl;
  cout << "energy(a) = " << energy(a) << endl;
  cout << "variance(a) = " << variance(a) << endl;
  cout << "moment(a, 1) = " << round_to_zero(moment(a, 1)) << endl;
  cout << "moment(a, 2) = " << moment(a, 2) << endl;
  cout << "moment(a, 3) = " << moment(a, 3) << endl;
  cout << "skewness(a) = " << skewness(a) << endl;
  cout << "kurtosisexcess(a) = " << kurtosisexcess(a) << endl;
  cout << "kurtosis(a) = " << kurtosis(a) << endl << endl;

  mat A = randn(5, 5);
  cout << "A = " << A << endl << endl;

  cout << "max(A) = " << max(A) << endl;
  cout << "max(A, 1) = " << max(A, 1) << endl;
  cout << "max(A, 2) = " << max(A, 2) << endl;
  cout << "min(A) = " << min(A) << endl;
  cout << "min(A, 1) = " << min(A, 1) << endl;
  cout << "min(A, 2) = " << min(A, 2) << endl;
  cout << "mean(A) = " << mean(A) << endl;
  cout << "geometric_mean(abs(A)) = " << geometric_mean(abs(A)) << endl;
  cout << "norm(A) = " << norm(A) << endl;
  cout << "norm(A, 2) = " << norm(A, 2) << endl;
  cout << "norm(A, 1) = " << norm(A, 1) << endl;
  cout << "norm(A, \"fro\") = " << norm(A, "fro") << endl << endl;

  return 0;
}
