/*!
 * \file
 * \brief Definitions of special operations on vectors and matricies optimized
 * for speed
 * \author Tony Ottosson and Tobias Ringstrom
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

#ifndef FASTMATH_H
#define FASTMATH_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \relatesalso Mat
  \brief Calculates m=m-v*v'*m
*/
ITPP_EXPORT void sub_v_vT_m(mat &m, const vec &v);

/*!
  \relatesalso Mat
  \brief Calculates m=m-m*v*v'
*/
ITPP_EXPORT void sub_m_v_vT(mat &m, const vec &v);

} // namespace itpp

#endif // #ifndef FASTMATH_H
