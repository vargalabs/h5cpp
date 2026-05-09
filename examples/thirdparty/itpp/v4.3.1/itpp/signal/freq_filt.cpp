/*!
 * \file
 * \brief Implementation of Freq_Filt Class
 * \author Simon Wood
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

#include <itpp/signal/freq_filt.h>


//! \cond

namespace itpp
{

template  class ITPP_EXPORT Vec<double>;
template  class ITPP_EXPORT Vec<short>;
template  class ITPP_EXPORT Vec<int>;
template  class ITPP_EXPORT Vec<std::complex<double> >;

template class ITPP_EXPORT Freq_Filt<double>;
template class ITPP_EXPORT Freq_Filt<std::complex<double> >;
template class ITPP_EXPORT Freq_Filt<short>;
template class ITPP_EXPORT Freq_Filt<int>;

} // namespace itpp

//! \endcond
