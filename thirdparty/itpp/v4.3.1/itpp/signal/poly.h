/*!
 * \file
 * \brief Polynomial functions
 * \author Tony Ottosson, Kumar Appaiah and Adam Piatyszek
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

#ifndef POLY_H
#define POLY_H

#include <itpp/base/vec.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \brief Create a polynomial of the given roots
  \ingroup poly

  Create a polynomial \c p with roots \c r

  @{
*/
ITPP_EXPORT void poly(const vec &r, vec &p);
inline vec poly(const vec &r) { vec temp; poly(r, temp); return temp; }
ITPP_EXPORT void poly(const cvec &r, cvec &p);
inline cvec poly(const cvec &r) { cvec temp; poly(r, temp); return temp; }
/*! @} */


/*!
  \brief Calculate the roots of the polynomial
  \ingroup poly

  Calculate the roots \c r of the polynomial \c p

  @{
*/
ITPP_EXPORT void roots(const vec &p, cvec &r);
inline cvec roots(const vec &p) { cvec temp; roots(p, temp); return temp; }
ITPP_EXPORT void roots(const cvec &p, cvec &r);
inline cvec roots(const cvec &p) { cvec temp; roots(p, temp); return temp; }
/*! @} */


/*!
  \brief Evaluate polynomial
  \ingroup poly

  Evaluate the polynomial \c p (of length \f$N+1\f$ at the points \c x
  The output is given by
  \f[
  p_0 x^N + p_1 x^{N-1} + \ldots + p_{N-1} x + p_N
  \f]

  @{
*/
ITPP_EXPORT vec polyval(const vec &p, const vec &x);
ITPP_EXPORT cvec polyval(const vec &p, const cvec &x);
ITPP_EXPORT cvec polyval(const cvec &p, const vec &x);
ITPP_EXPORT cvec polyval(const cvec &p, const cvec &x);
/*! @} */

/*!
  \brief Chebyshev polynomial of the first kind
  \ingroup poly

  Chebyshev polynomials of the first kind can be defined as follows:
  \f[
  T(x) = \left\{
  \begin{array}{ll}
  \cos(n\arccos(x)),& |x| \leq 0 \\
  \cosh(n\mathrm{arccosh}(x)),& x > 1 \\
  (-1)^n \cosh(n\mathrm{arccosh}(-x)),& x < -1
  \end{array}
  \right.
  \f]

  \param n order of the Chebyshev polynomial
  \param x value at which the Chebyshev polynomial is to be evaluated

  \author Kumar Appaiah, Adam Piatyszek (code review)
*/
ITPP_EXPORT double cheb(int n, double x);

/*!
  \brief Chebyshev polynomial of the first kind
  \ingroup poly

  Chebyshev polynomials of the first kind can be defined as follows:
  \f[
  T(x) = \left\{
  \begin{array}{ll}
  \cos(n\arccos(x)),& |x| \leq 0 \\
  \cosh(n\mathrm{arccosh}(x)),& x > 1 \\
  (-1)^n \cosh(n\mathrm{arccosh}(-x)),& x < -1
  \end{array}
  \right.
  \f]

  \param n order of the Chebyshev polynomial
  \param x vector of values at which the Chebyshev polynomial is to
  be evaluated
  \return values of the Chebyshev polynomial evaluated for each
  element of \c x

  \author Kumar Appaiah, Adam Piatyszek (code review)
*/
ITPP_EXPORT vec cheb(int n, const vec &x);

/*!
  \brief Chebyshev polynomial of the first kind
  \ingroup poly

  Chebyshev polynomials of the first kind can be defined as follows:
  \f[
  T(x) = \left\{
  \begin{array}{ll}
  \cos(n\arccos(x)),& |x| \leq 0 \\
  \cosh(n\mathrm{arccosh}(x)),& x > 1 \\
  (-1)^n \cosh(n\mathrm{arccosh}(-x)),& x < -1
  \end{array}
  \right.
  \f]

  \param n order of the Chebyshev polynomial
  \param x matrix of values at which the Chebyshev polynomial is to
  be evaluated
  \return values of the Chebyshev polynomial evaluated for each
  element in \c x.

  \author Kumar Appaiah, Adam Piatyszek (code review)
*/
ITPP_EXPORT mat cheb(int n, const mat &x);
} // namespace itpp

#endif // #ifndef POLY_H
