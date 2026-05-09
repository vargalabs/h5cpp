/*!
 * \file
 * \brief Filter design functions
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

#ifndef FILTER_DESIGN_H
#define FILTER_DESIGN_H

#include <itpp/base/vec.h>
#include <itpp/itexports.h>


namespace itpp
{

/*!
  \addtogroup filters
*/



/*!
  \brief Polynomial Stabilization
  \ingroup filters
  \author Tony Ottosson

  Stabilizes the polynomial transfer function by replacing all roots outside
  the unit cirlce with their reflection inside the unit circle.

  @{
*/
ITPP_EXPORT void polystab(const vec &a, vec &out);
inline vec polystab(const vec &a) { vec temp; polystab(a, temp); return temp; }
ITPP_EXPORT void polystab(const cvec &a, cvec &out);
inline cvec polystab(const cvec &a) { cvec temp; polystab(a, temp); return temp; }
/*! @} */

/*!
  \brief Frequency response of filter
  \ingroup filters
  \author Tony Ottosson

  Calculates the N-point frequency response of the supplied digital filter over the frequencies w.
  If w is not given the response is evaluated over the range 0 to \f$\pi\f$ with N values.
  The default value of N is 512.

  If \c w is supplied polyval() is used. Otherwise the calculation is based on the fft.

  @{
*/
ITPP_EXPORT void freqz(const cvec &b, const cvec& a, const int N, cvec &h, vec &w);
ITPP_EXPORT cvec freqz(const cvec &b, const cvec& a, const int N = 512);
ITPP_EXPORT cvec freqz(const cvec &b, const cvec& a, const vec &w);

ITPP_EXPORT void freqz(const vec &b, const vec& a, const int N, cvec &h, vec &w);
ITPP_EXPORT cvec freqz(const vec &b, const vec& a, const int N = 512);
ITPP_EXPORT cvec freqz(const vec &b, const vec& a, const vec &w);
/*! @} */



/*!
  \brief Calculate autocorrelation from the specified frequency-response (suitable for filter design)
  \ingroup filters
  \author Tony Ottosson

  Calculates the autocorrelation function of size \c N corresponding to the specified frequency response.
  Useful as a first step in designing filters.

  The vectors \c f and \c m is the frequency response. The frequencies should be between 0 and 1.0
  (equal to half the sample rate) in increasing order. Both 0.0 and 1.0 must be included. The frequency
  response is upsampled to 512 points and the autocorrelation is ifft of the power
  magnitude response of the upsampled frequency response.
*/
ITPP_EXPORT void filter_design_autocorrelation(const int N, const vec &f, const vec &m, vec &R);


/*!
  \brief Estimation of AR-part in an ARMA model given the autocorrelation
  \ingroup filters
  \author Tony Ottosson

  Estimates the AR-part of an ARMA model from the given autocorrelation. The AR part is of order \c n.
  The overdetermined modified Yule-Walker equations are used.

  If \f$N>n\f$ then the system is overdetermined and a least squares solution is used.
  As a rule of thumb use \f$N = 4 n\f$

  The parameter \c m is the order of the MA-part such that \f$R(k) = 0, \forall \|k\| > m\f$.

  The supplied autocorrelation should at least be of size \c N.

  References:
  Stoica and Moses, Introduction to spectral analysis, Prentice Hall, 1997.
*/
ITPP_EXPORT void modified_yule_walker(const int m, const int n, const int N, const vec &R, vec &a);



/*!
  \brief Estimation of ARMA model given the autocorrelation
  \ingroup filters
  \author Tony Ottosson

  Estimates an ARMA model from the given autocorrelation. The AR part is of order \c n and the MA part
  is of order \c m.

  The AR part (the denominator) is calcuated using the modified Yule-Walker equations. The
  the MA part (the nominator) is calculated by calculating the inverse magnitude spectrum using FFTs of
  size 512 which is an AR-system. This AR-system is then solved using the Levinson-Durbin algorithm.

  The supplied autocorrelation is windowed using a Hamming window of size \f$2(m+n)\f$ and hence should at
  least be of that size.

  References:
  [1] Stoica and Moses, Introduction to spectral analysis, Prentice Hall, 1997.
  [2] B. Friedlander and B. Porat, The modified Yule-Walker method of ARMA spectral estimation,
  IEEE Trans. Aerospace and Electronic Systems, Vol. AES-20, No. 2, pp. 158--173, March 1984.

*/
ITPP_EXPORT void arma_estimator(const int m, const int n, const vec &R, vec &b, vec &a);


/*!
  \brief ARMA filter design using a least-squares fit to the specified frequency-response
  \ingroup filters
  \author Tony Ottosson

  The arma_estimator() function is used to calculate the a and b coefficients.

  The vectors \c f and \c m is the frequency response. The frequencies should be between 0 and 1.0
  (equal to half the sample rate) in increasing order. Both 0.0 and 1.0 must be included. The
  filter_design_autocorrelation() fucnction is used to interpolate the frequency response and calculate
  the corresponding autocorrelation.

  Observe: this function will not always give exactly the same result as the matlab yulewalk function.
*/
ITPP_EXPORT void yulewalk(const int N, const vec &f, const vec &m, vec &b, vec &a);


} // namespace itpp

#endif // #ifndef FILTER_DESIGN_H
