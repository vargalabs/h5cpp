/*!
 * \file
 * \brief Import/Export definitions for some templates defined in base folder.
 *
 * \author Andy Panov
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2013  (see AUTHORS file for a list of contributors)
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

#ifndef BASEEXPORTS_H
#define BASEEXPORTS_H

#include <itpp/itexports.h>

namespace itpp
{

//! \cond

#if (defined(_MSC_VER) && defined(ITPP_SHARED_LIB))
//MSVC needs to explicitely instantiate required templates while building the
//shared library. Also, these definitions are marked as imported when library is
//linked with user's code.
template class ITPP_EXPORT Array<bool>;
template class ITPP_EXPORT Array<std::string>;
template class ITPP_EXPORT Array<bmat>;
template class ITPP_EXPORT Array<mat>;
template class ITPP_EXPORT Array<vec>;
template class ITPP_EXPORT Array<ivec>;
template class ITPP_EXPORT Array<cvec>;
template class ITPP_EXPORT Array<Array<int> >;
template class ITPP_EXPORT Array<Vec<unsigned int> >;
template class ITPP_EXPORT Array<Array<vec> >;
template class ITPP_EXPORT Array<Array<cvec> >;
#endif

//! \endcond

}

#endif
