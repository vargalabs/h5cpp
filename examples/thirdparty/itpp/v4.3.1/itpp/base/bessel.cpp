/*!
 * \file
 * \brief Implementation of Bessel functions
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

#include <itpp/base/bessel.h>
#include <itpp/base/bessel/bessel_internal.h>
#include <itpp/base/itcompat.h>

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif


namespace itpp
{

// Bessel function of order nu
double besselj(int nu, double x) { return jn(nu, x); }

vec besselj(int nu, const vec &x)
{
  vec out(x.size());
  for (int i = 0; i < x.size(); i++)
    out(i) = jn(nu, x(i));

  return out;
}

// Bessel function of order nu. nu is real.
double besselj(double nu, double x) { return jv(nu, x); }

vec besselj(double nu, const vec &x)
{
  vec out(x.size());
  for (int i = 0; i < x.size(); i++)
    out(i) = jv(nu, x(i));

  return out;
}

// Bessel function of second kind of order nu
double bessely(int nu, double x) { return yn(nu, x); }

vec bessely(int nu, const vec &x)
{
  vec out(x.size());
  for (int i = 0; i < x.size(); i++)
    out(i) = yn(nu, x(i));

  return out;
}
// Bessel function of second kind of order nu
double bessely(double nu, double x) { return yv(nu, x); }

vec bessely(double nu, const vec &x)
{
  vec out(x.size());
  for (int i = 0; i < x.size(); i++)
    out(i) = yv(nu, x(i));

  return out;
}

// Modified Bessel function of order nu
double besseli(double nu, double x) { return iv(nu, x); }

vec besseli(double nu, const vec &x)
{
  vec out(x.size());
  for (int i = 0; i < x.size(); i++)
    out(i) = iv(nu, x(i));

  return out;
}

// Modified Bessel function of second kind of order n
double besselk(int n, double x) { return kn(n, x); }

vec besselk(int nu, const vec &x)
{
  vec out(x.size());
  for (int i = 0; i < x.size(); i++)
    out(i) = kn(nu, x(i));

  return out;
}

} // namespace itpp
