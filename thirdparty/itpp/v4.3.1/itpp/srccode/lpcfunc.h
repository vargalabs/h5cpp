/*!
 * \file
 * \brief Implementation of linear prediction functions, and conversion
 * between common representations of linear predictive parameters
 * \author Thomas Eriksson
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
 * Here we define functions for linear predictive analysis and
 * coding. Several functions for computation of lpc parameters are
 * given, together with functions to convert between various common
 * representations of the lpc polynomial coefficients. The
 * functionality is the same as in the MATLAB functions with the same
 * names.
 *
 * The short term for the various parameter types are:
 * poly - LPC polynomial coefficients
 * ac - autocorrelation coefficients
 * rc - reflection coefficients
 * lar - log area ratios
 * is - inverse sine parameters
 * lsf - line spectral frequencies
 */

#ifndef LPCFUNC_H
#define LPCFUNC_H

#include <itpp/base/vec.h>
#include <itpp/itexports.h>


namespace itpp
{

/*! \addtogroup lpc
 */
//!@{

//! Returns a chirped version of the input vector
ITPP_EXPORT vec chirp(const vec &a, double factor);
//! Spectral distortion between two vectors, in dB.
ITPP_EXPORT double sd(const vec &In1, const vec &In2);
//! Spectral distortion between two vectors, in dB, up to highest frequency highestfreq.
ITPP_EXPORT double sd(const vec &In1, const vec &In2, double highestfreq);
//! Computes reflection coefficients from autocorrelation, using the Le-Roux-Guegen algorithm.
ITPP_EXPORT vec lerouxguegenrc(const vec &R, int order);
//! Levinson    - Levinson-Durbin recursion.
ITPP_EXPORT vec levinson(const vec &R2, int order);
//! Computes the autocorrelation function
ITPP_EXPORT vec autocorr(const vec &x, int order);
//!    lpc         - Linear Predictive Coefficients using autocorrelation method.
ITPP_EXPORT vec lpc(const vec &x, int order);
//!    schurrc     - Schur algorithm.
ITPP_EXPORT vec schurrc(const vec &R, int order);
//!    ac2rc       - Autocorrelation sequence to reflection coefficients conversion.
ITPP_EXPORT vec ac2rc(const vec &ac);
//!    ac2poly     - Autocorrelation sequence to prediction polynomial conversion.
ITPP_EXPORT vec ac2poly(const vec &ac);
//!    is2rc       - Inverse sine parameters to reflection coefficients conversion.
ITPP_EXPORT vec is2rc(const vec &is);
//!    lar2rc      - Log area ratios to reflection coefficients conversion.
ITPP_EXPORT vec lar2rc(const vec &lar);
//!    lsf2poly    - Line spectral frequencies to prediction polynomial conversion.
ITPP_EXPORT vec lsf2poly(const vec &lsf);
//!    poly2ac     - Prediction polynomial to autocorrelation sequence conversion.
ITPP_EXPORT vec poly2ac(const vec &poly);
//!    poly2lsf    - Prediction polynomial to line spectral frequencies conversion.
ITPP_EXPORT vec poly2lsf(const vec &poly);
//!    poly2rc     - Prediction polynomial to reflection coefficients conversion.
ITPP_EXPORT vec poly2rc(const vec &poly);
//!    poly2cepstrum     - Prediction polynomial to cepstrum conversion.
ITPP_EXPORT vec poly2cepstrum(const vec &a);
//!    poly2cepstrum     - Prediction polynomial to cepstrum conversion, to the specified order.
ITPP_EXPORT vec poly2cepstrum(const vec &a, int num);
//!    cepstrum2poly     - Cepstrum to prediction polynomial conversion.
ITPP_EXPORT vec cepstrum2poly(const vec &c);
//!    rc2ac       - Reflection coefficients to autocorrelation sequence conversion.
ITPP_EXPORT vec rc2ac(const vec &rc);
//!    rc2is       - Reflection coefficients to inverse sine parameters conversion.
ITPP_EXPORT vec rc2is(const vec &rc);
//!    rc2lar      - Reflection coefficients to log area ratios conversion.
ITPP_EXPORT vec rc2lar(const vec &rc);
//!    rc2poly     - Reflection coefficients to prediction polynomial conversion.
ITPP_EXPORT vec rc2poly(const vec &rc);

//!@}

} // namespace itpp

#endif // #ifndef LPCFUNC_H
