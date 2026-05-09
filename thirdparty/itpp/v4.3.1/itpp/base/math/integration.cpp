/*!
 * \file
 * \brief Implementation of numerical integration
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

#include <itpp/base/math/integration.h>

//! \cond
namespace itpp
{

//wrapper to convert double(*f)(double) to function object
class Integrand_Wrapper
{
  typedef double(*Ftn)(double);
  Ftn _f;
public:
  explicit Integrand_Wrapper(Ftn f): _f(f) {}
  double operator()(double x) const {return _f(x);}
};


template double quad(Integrand_Wrapper, double, double, double);
template double quadl(Integrand_Wrapper, double, double, double);

//--------------------- quad() ----------------------------------------
double quad(double(*f)(double), double a, double b,
            double tol)
{
  return quad(Integrand_Wrapper(f), a, b, tol);
}

//--------------------- quadl() ----------------------------------------
double quadl(double(*f)(double), double a, double b, double tol)
{
  return quadl(Integrand_Wrapper(f), a, b, tol);
}


} // namespace itpp

//! \endcond


