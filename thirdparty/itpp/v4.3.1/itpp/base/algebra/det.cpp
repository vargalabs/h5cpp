/*!
 * \file
 * \brief Implementation of determinant calculations
 * \author Tony Ottosson
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

#include <itpp/base/algebra/det.h>
#include <itpp/base/algebra/lu.h>


namespace itpp
{


/* Determinant of square matrix.
   Calculate determinant of inmatrix (Uses LU-factorisation)
   (See Theorem 3.2.1 p. 97 in Golub & van Loan, "Matrix Computations").

   det(X) = det(P')*det(L)*det(U) = det(P')*1*prod(diag(U))
*/
double det(const mat &X)
{
  it_assert_debug(X.rows() == X.cols(), "det : Only square matrices");

  mat L, U;
  ivec p;
  double s = 1.0;
  int i;

  lu(X, L, U, p); // calculate LU-factorisation

  double temp = U(0, 0);
  for (i = 1;i < X.rows();i++) {
    temp *= U(i, i);
  }

  // Calculate det(P'). Equal to (-1)^(no row changes)
  for (i = 0; i < p.size(); i++)
    if (i != p(i))
      s *= -1.0;

  return temp*s;
}


/* Determinant of complex square matrix.
   Calculate determinant of inmatrix (Uses LU-factorisation)
   (See Theorem 3.2.1 p. 97 in Golub & van Loan, "Matrix Computations").

   det(X) = det(P')*det(L)*det(U) = det(P')*1*prod(diag(U))

   Needs LU-factorization of complex matrices (LAPACK)
*/
std::complex<double> det(const cmat &X)
{
  it_assert_debug(X.rows() == X.cols(), "det : Only square matrices");

  int i;
  cmat L, U;
  ivec p;
  double s = 1.0;

  lu(X, L, U, p); // calculate LU-factorisation

  std::complex<double> temp = U(0, 0);
  for (i = 1;i < X.rows();i++) {
    temp *= U(i, i);
  }

  // Calculate det(P'). Equal to (-1)^(no row changes)
  for (i = 0; i < p.size(); i++)
    if (i != p(i))
      s *= -1.0;

  return temp*s;
}


} // namespace itpp
