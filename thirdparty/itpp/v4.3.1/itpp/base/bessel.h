/*!
 * \file
 * \brief Definitions of Bessel functions
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

#ifndef BESSEL_H
#define BESSEL_H

#include <itpp/base/vec.h>
#include <itpp/itexports.h>

namespace itpp
{

/*! \addtogroup besselfunctions
 */

/*!
  \ingroup besselfunctions
  \brief Bessel function of first kind of order \a nu for \a nu integer

  The bessel function of first kind is defined as:
  \f[
  J_{\nu}(x) = \sum_{k=0}^{\infty} \frac{ (-1)^{k} }{k! \Gamma(\nu+k+1) } \left(\frac{x}{2}\right)^{\nu+2k}
  \f]
  where \f$\nu\f$ is the order and \f$ 0 < x < \infty \f$.
*/
ITPP_EXPORT double besselj(int nu, double x);

/*!
  \ingroup besselfunctions
  \brief Bessel function of first kind of order \a nu for \a nu integer
*/
ITPP_EXPORT vec besselj(int nu, const vec &x);

/*!
  \ingroup besselfunctions
  \brief Bessel function of first kind of order \a nu. \a nu is real.
*/
ITPP_EXPORT double besselj(double nu, double x);

/*!
  \ingroup besselfunctions
  \brief Bessel function of first kind of order \a nu. \a nu is real.
*/
ITPP_EXPORT vec besselj(double nu, const vec &x);

/*!
  \ingroup besselfunctions
  \brief Bessel function of second kind of order \a nu. \a nu is integer.

  The Bessel function of second kind is defined as:
  \f[
  Y_{\nu}(x) = \frac{J_{\nu}(x) \cos(\nu\pi) - J_{-\nu}(x)}{\sin(\nu\pi)}
  \f]
  where \f$\nu\f$ is the order and \f$ 0 < x < \infty \f$.
*/
ITPP_EXPORT double bessely(int nu, double x);

/*!
  \ingroup besselfunctions
  \brief Bessel function of second kind of order \a nu. \a nu is integer.
*/
ITPP_EXPORT vec bessely(int nu, const vec &x);

/*!
  \ingroup besselfunctions
  \brief Bessel function of second kind of order \a nu. \a nu is real.
*/
ITPP_EXPORT double bessely(double nu, double x);

/*!
  \ingroup besselfunctions
  \brief Bessel function of second kind of order \a nu. \a nu is real.
*/
ITPP_EXPORT vec bessely(double nu, const vec &x);

/*!
  \ingroup besselfunctions
  \brief Modified Bessel function of first kind of order \a nu. \a nu is \a double. \a x is \a double.

  The Modified Bessel function of first kind is defined as:
  \f[
  I_{\nu}(x) = i^{-\nu} J_{\nu}(ix)
  \f]
  where \f$\nu\f$ is the order and \f$ 0 < x < \infty \f$.
*/
ITPP_EXPORT double besseli(double nu, double x);

/*!
  \ingroup besselfunctions
  \brief Modified Bessel function of first kind of order \a nu. \a nu is \a double. \a x is \a double.
*/
ITPP_EXPORT vec besseli(double nu, const vec &x);

/*!
  \ingroup besselfunctions
  \brief Modified Bessel function of second kind of order \a nu. \a nu is double. \a x is double.

  The Modified Bessel function of second kind is defined as:
  \f[
  K_{\nu}(x) = \frac{\pi}{2} i^{\nu+1} [J_{\nu}(ix) + i Y_{\nu}(ix)]
  \f]
  where \f$\nu\f$ is the order and \f$ 0 < x < \infty \f$.
*/
ITPP_EXPORT double besselk(int nu, double x);

/*!
  \ingroup besselfunctions
  \brief Modified Bessel function of second kind of order \a nu. \a nu is double. \a x is double.
*/
ITPP_EXPORT vec besselk(int nu, const vec &x);

} //namespace itpp

#endif // #ifndef BESSEL_H
