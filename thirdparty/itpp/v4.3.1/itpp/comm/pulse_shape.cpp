/*!
 * \file
 * \brief Pulse shaping classes - source file
 * \author Tony Ottosson, Hakan Eriksson and Adam Piatyszek
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

#include <itpp/comm/pulse_shape.h>


namespace itpp
{

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------

template class ITPP_EXPORT Pulse_Shape<double, double, double>;
template class ITPP_EXPORT Pulse_Shape < std::complex<double>, double,
std::complex<double> >;
template class ITPP_EXPORT Pulse_Shape < std::complex<double>, std::complex<double>,
std::complex<double> >;

template class ITPP_EXPORT Root_Raised_Cosine<double>;
template class ITPP_EXPORT Root_Raised_Cosine<std::complex<double> >;

template class ITPP_EXPORT Raised_Cosine<double>;
template class ITPP_EXPORT Raised_Cosine<std::complex<double> >;

//! \endcond

} // namespace itpp
