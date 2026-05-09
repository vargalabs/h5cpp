/*!
 * \file
 * \brief Definitions of window functions
 * \author Tony Ottosson, Tobias Ringstrom, Pal Frenger, Adam Piatyszek
 *         and Kumar Appaiah
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

#ifndef WINDOW_H
#define WINDOW_H

#include <itpp/base/vec.h>
#include <itpp/itexports.h>


namespace itpp
{

/*!
  \addtogroup windfunc
*/

/*!\addtogroup windfunc
  \brief Windowing functions

*/
//!@{

/*! \brief Hamming window

The \c n size Hamming window is a vector \f$w\f$ where the \f$i\f$th component is
\f[
w_i = 0.54 - 0.46 \cos(2\pi i/(n-1))
\f]
*/
ITPP_EXPORT vec hamming(int size);


/*! \brief Hanning window

The \c n size Hanning window is a vector \f$w\f$ where the \f$i\f$th component is
\f[
w_i = 0.5(1 - \cos(2\pi (i+1)/(n+1))
\f]

Observe that this function is not the same as the hann() function which is defined
as in matlab.
*/
ITPP_EXPORT vec hanning(int n);

/*! \brief Hanning window compatible with matlab

The \c n size Hanning window is a vector \f$w\f$ where the \f$i\f$th component is
\f[
w_i = 0.5(1 - \cos(2\pi i/(n-1))
\f]
*/
ITPP_EXPORT vec hann(int n);


/*! \brief Blackman window

The \c n size Blackman window is a vector \f$w\f$ where the \f$i\f$th component is
\f[
w_i = 0.42 - 0.5\cos(2\pi i/(n-1)) + 0.08\cos(4\pi i/(n-1))
\f]
*/
ITPP_EXPORT vec blackman(int n);

/*! \brief Triangular window


The \c n size triangle window is a vector \f$w\f$ where the \f$i\f$th component is
\f[
w_i = w_{n-i-1} = \frac{2(i+1)}{n+1}
\f]
for \c n odd and for \c n even
\f[
w_i = w_{n-i-1} = \frac{2i+1}{n}
\f]
*/
ITPP_EXPORT vec triang(int n);

/*! \brief Square root window

The square-root of the Triangle window.
sqrt_win(n) = sqrt(triang(n))
*/
ITPP_EXPORT vec sqrt_win(int n);


/*!
  \brief Dolph-Chebyshev window

  The length \c n Dolph-Chebyshev window is a vector \f$w\f$ whose \f$i\f$th
  transform component is given by
  \f[
  W[k] = \frac{T_M\left(\beta \cos\left(\frac{\pi k}{M}\right)
  \right)}{T_M(\beta)},k = 0, 1, 2, \ldots, M - 1
  \f]
  where \c T_n(x) is the order \c n Chebyshev polynomial of the first kind.

  \param n length of the Doplh-Chebyshev window
  \param at attenutation of side lobe (in dB)
  \return symmetric length \c n Doplh-Chebyshev window

  \author Kumar Appaiah and Adam Piatyszek (code review)
*/
ITPP_EXPORT vec chebwin(int n, double at);
//!@}


} //namespace itpp

#endif // #ifndef WINDOW_H
