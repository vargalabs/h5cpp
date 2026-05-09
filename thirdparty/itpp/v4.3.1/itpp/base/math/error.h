/*!
 * \file
 * \brief Error functions - header file
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

#ifndef ERROR_H
#define ERROR_H

#include <itpp/base/help_functions.h>
#include <itpp/itexports.h>

namespace itpp
{

//!\addtogroup errorfunc
//!@{

/*!
 * \brief Error function for complex argument
 * \author Adam Piatyszek
 *
 * This function calculates a well known error function \c erf(z)
 * for complex \c z. The implementation is based on unofficial
 * implementation for Octave. Here is a part of the author's note
 * from original sources:
 *
 * Put together by John Smith john at arrows dot demon dot co dot uk,
 * using ideas by others.
 *
 * Calculate \c erf(z) for complex \c z.
 * Three methods are implemented; which one is used depends on z.
 *
 * The code includes some hard coded constants that are intended to
 * give about 14 decimal places of accuracy. This is appropriate for
 * 64-bit floating point numbers.
 */
ITPP_EXPORT std::complex<double> erf(const std::complex<double>& z);

//! Inverse of error function
ITPP_EXPORT double erfinv(double x);

//! Q-function
ITPP_EXPORT double Qfunc(double x);


// ----------------------------------------------------------------------
// functions for matrices and vectors
// ----------------------------------------------------------------------

//! Error function
ITPP_EXPORT vec erf(const vec &x);
//! Error function
ITPP_EXPORT mat erf(const mat &x);
//! Error function
ITPP_EXPORT cvec erf(const cvec &x);
//! Error function
ITPP_EXPORT cmat erf(const cmat &x);

//! Inverse of error function
ITPP_EXPORT vec erfinv(const vec &x);
//! Inverse of error function
ITPP_EXPORT mat erfinv(const mat &x);

//! Complementary error function
ITPP_EXPORT vec erfc(const vec &x);
//! Complementary error function
ITPP_EXPORT mat erfc(const mat &x);

//! Q-function
ITPP_EXPORT vec Qfunc(const vec &x);
//! Q-function
ITPP_EXPORT mat Qfunc(const mat &x);
//!@}

} // namespace itpp

#endif // #ifndef ERROR_H




