/*!
 * \file
 * \brief Trigonometric and hyperbolic functions - source file
 * \author Tony Ottosson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2008  (see AUTHORS file for a list of contributors)
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

#include <itpp/base/math/trig_hyp.h>
#include <itpp/base/itcompat.h>

namespace itpp
{

//! Inverse sine hyperbolic function
vec asinh(const vec &x) { return apply_function<double>(::asinh, x); }
//! Inverse sine hyperbolic function
mat asinh(const mat &x) { return apply_function<double>(::asinh, x); }
//! Inverse cosine hyperbolic function
vec acosh(const vec &x) { return apply_function<double>(::acosh, x); }
//! Inverse cosine hyperbolic function
mat acosh(const mat &x) { return apply_function<double>(::acosh, x); }
//! Inverse tan hyperbolic function
vec atanh(const vec &x) { return apply_function<double>(::atanh, x); }
//! Inverse tan hyperbolic function
mat atanh(const mat &x) { return apply_function<double>(::atanh, x); }

} // namespace itpp
