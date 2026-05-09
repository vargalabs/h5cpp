/*!
 * \file
 * \brief Implementation of interleaver classes
 * \author Pal Frenger
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

#include <itpp/comm/interleave.h>


namespace itpp
{

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------

template class ITPP_EXPORT Block_Interleaver<double>;
template class ITPP_EXPORT Block_Interleaver<short>;
template class ITPP_EXPORT Block_Interleaver<int>;
template class ITPP_EXPORT Block_Interleaver<std::complex<double> >;
template class ITPP_EXPORT Block_Interleaver<bin>;

template class ITPP_EXPORT Cross_Interleaver<double>;
template class ITPP_EXPORT Cross_Interleaver<short>;
template class ITPP_EXPORT Cross_Interleaver<int>;
template class ITPP_EXPORT Cross_Interleaver<std::complex<double> >;
template class ITPP_EXPORT Cross_Interleaver<bin>;

template class ITPP_EXPORT Sequence_Interleaver<double>;
template class ITPP_EXPORT Sequence_Interleaver<short>;
template class ITPP_EXPORT Sequence_Interleaver<int>;
template class ITPP_EXPORT Sequence_Interleaver<std::complex<double> >;
template class ITPP_EXPORT Sequence_Interleaver<bin>;

//! \endcond

} // namespace itpp
