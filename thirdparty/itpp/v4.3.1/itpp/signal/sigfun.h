/*!
 * \file
 * \brief Definitions of signal processing functions
 * \author Tony Ottosson, Thomas Eriksson, Pal Frenger, and Tobias Ringstrom
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

#ifndef SIGFUN_H
#define SIGFUN_H

#include <itpp/base/vec.h>
#include <itpp/itexports.h>


namespace itpp
{

/*!
 * \addtogroup sigproc
 * @{
 */

/*!
  \brief Cross-correlation calculation

  \c z=xcorr(x,y,max_lag), where \a x and \a y are length \a M vectors \a (M>1), returns the length \c 2*max_lag+1
  cross-correlation sequence \a z. (lags: \c -max_lag,...,0,...,max_lag)

  For \c max_lag=-1 the cross-correlation sequence is of length \c 2*M-1, i.e., the cross-correlation for all possible lags.

  Scaling options \a scaleopt:
  \arg \a none \c => No scaling of the cross-correlation vector
  \arg \a biased \c => Scales the cross-correlation vector by \c 1/M
  \arg \a unbiased \c => Scales the cross-correlation vector by \c 1/(M-abs(lag))
  \arg \a coeff \c => Normalises the cross-correlation to 1 for zero lag.

  \note \c max_lag \c <= \c M-1
  \note If \c x and \c y are of different length, the shortest one is zero-padded

  \param x (Input) Vector of samples
  \param y (Input) Vector of samples
  \param out (Output) The cross correlation between \a x and \a y.
  \param max_lag (Input) Maximum lag for which the cross-correlation is calculated. The output vector is of size \c 2*maxlag+1.
  Default value: \c max_lag=-1 calculates the cross-correlations for all possible lags
  \param scaleopt (Input) Indicates how the cross-correlation function should be scaled. Default value: \c "none" indicates that no scaling is done

  @{
*/
ITPP_EXPORT void xcorr_old(const vec &x, const vec &y, vec &out, const int max_lag = -1, const std::string scaleopt = "none");
ITPP_EXPORT void xcorr(const vec &x, const vec &y, vec &out, const int max_lag = -1, const std::string scaleopt = "none");
/*! @} */

/*!
  \brief Cross-correlation calculation

  \c z=xcorr(x,y,max_lag), where \a x and \a y are length \a M vectors \a (M>1), returns the length \c 2*max_lag+1
  cross-correlation sequence \a z. (lags: \c -max_lag,...,0,...,max_lag)

  For \c max_lag=-1 the cross-correlation sequence is of length \c 2*M-1, i.e., the cross-correlation for all possible lags.

  Scaling options \a scaleopt:
  \arg \a none \c => No scaling of the cross-correlation vector
  \arg \a biased \c => Scales the cross-correlation vector by \c 1/M
  \arg \a unbiased \c => Scales the cross-correlation vector by \c 1/(M-abs(lag))
  \arg \a coeff \c => Normalises the cross-correlation to 1 for zero lag.

  \note \c max_lag \c <= \c M-1
  \note If \c x and \c y are of different length, the shortest one is zero-padded

  \param x (Input) Vector of samples
  \param y (Input) Vector of samples
  \param max_lag (Input) Maximum lag for which the cross-correlation is calculated. The output vector is of size \c 2*maxlag+1.
  Default value: \c max_lag=-1 calculates the cross-correlations for all possible lags
  \param scaleopt (Input) Indicates how the cross-correlation function should be scaled. Default
  value: \c "none" indicates that no scaling is done
  \returns The cross correlation between \a x and \a y.

  @{
*/
ITPP_EXPORT vec xcorr_old(const vec &x, const vec &y, const int max_lag = -1, const std::string scaleopt = "none");
ITPP_EXPORT vec xcorr(const vec &x, const vec &y, const int max_lag = -1, const std::string scaleopt = "none");
/*! @} */

/*!
  \brief Cross Correlation

  \code r = xcorr(x,y) \endcode returns the cross-correlation vector \b r.
*/
ITPP_EXPORT cvec xcorr(const cvec &x, const cvec &y, const int max_lag = -1, const std::string scaleopt = "none");


/*!
  \brief Auto-correlation calculation

  \c z=xcorr(x,max_lag), where \a x and is a length \a M vector \a (M>1), returns the length \c 2*max_lag+1 auto-correlation
  sequence \a z. (lags: \c -max_lag,...,0,...,max_lag)

  For \c max_lag=-1 the auto-correlation sequence is of length \c 2*M-1, i.e., the cross correlation for all possible lags.

  Scaling options \a scaleopt:
  \arg \a none \c => No scaling of the auto-correlation vector
  \arg \a biased \c => Scales the auto-correlation vector by \c 1/M
  \arg \a unbiased \c => Scales the auto-correlation vector by \c 1/(M-abs(lag))
  \arg \a coeff \c => Normalises the auto-correlation so that \c acf(x)=1 for zero lag.

  \note \c max_lag \c <= \c M-1

  \param x (Input) Vector of samples
  \param max_lag (Input) Maximum lag for which the auto-correlation is calculated. The output vector is of size \c 2*maxlag+1.
  Default value \c max_lag=-1 calculates the auto-correlations for all possible lags.
  \param scaleopt (Input) Indicates how the auto-correlation function should be scaled.
  Default value: \c "none" indicates that no scaling is done.
  \returns The auto-correlation of \a x.

  @{
*/
ITPP_EXPORT vec xcorr_old(const vec &x, const int max_lag = -1, const std::string scaleopt = "none");
ITPP_EXPORT vec xcorr(const vec &x, const int max_lag = -1, const std::string scaleopt = "none");
/*! @} */

/*!
  \brief Cross Correlation

  \code r = xcorr(x) \endcode returns the auto-correlation vecotr \b r.
*/
ITPP_EXPORT cvec xcorr(const cvec &x, const int max_lag = -1, const std::string scaleopt = "none");

/*!
  \brief Cross Correlation

  \code xcorr(x,y,out) \endcode Computes the cross-correlatin and returns in vector \b out
*/
ITPP_EXPORT void xcorr(const cvec &x, const cvec &y, cvec &out, const int max_lag = -1, const std::string scaleopt = "none",
           bool autoflag = true);

/*!
  \brief Covariance matrix calculation

  Calculates the covariance matrix of the observations in the matrix \f$X\f$. Each
  row is an observation and each column represents a variable.

  The covariance is normalized with the number of observations \f$N\f$.
  The mean value is removed before calculation.

  Set is_zero_mean if X already has zero mean.
*/
ITPP_EXPORT mat cov(const mat &X, bool is_zero_mean = false);

//vec cov(const vec &x, short order);

/*!
  \brief Power spectrum calculation

  Calculates the power spectrum using the Welch method and a Hanning window.
*/
ITPP_EXPORT vec spectrum(const vec &v, int nfft = 256, int noverlap = 0);

/*!
  \brief Power spectrum calculation

  Calculates the power spectrum using using the Welch method and the supplied window w.
*/
ITPP_EXPORT vec spectrum(const vec &v, const vec &w, int noverlap = 0);

/*!
  \brief Power spectrum calculation of a filter

  Calculates the power spectrum of a filter with transfer function a(z)
*/
ITPP_EXPORT vec filter_spectrum(const vec &a, int nfft = 256);

/*!
  \brief Power spectrum calculation of a filter

  Calculates the power spectrum of a filter with transfer function a(z)/b(z)
*/
ITPP_EXPORT vec filter_spectrum(const vec &a, const vec &b, int nfft = 256);

/*! @} */

} // namespace itpp

#endif // #ifndef SIGFUN_H
