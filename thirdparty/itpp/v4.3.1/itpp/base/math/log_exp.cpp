/*!
 * \file
 * \brief Logarithmic and exponenential functions - source file
 * \author Tony Ottosson, Adam Piatyszek and Conrad Sanderson
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

#include <itpp/base/math/log_exp.h>
#include <itpp/base/itcompat.h>

namespace itpp
{

// log-2 of the elements
vec log2(const vec &x) { return apply_function<double>(::log2, x); }
mat log2(const mat &x) { return apply_function<double>(::log2, x); }

// Safe substitute for <tt>log(exp(log_a) + exp(log_b))</tt>
double log_add(double log_a, double log_b)
{
  if (log_a < log_b) {
    double tmp = log_a;
    log_a = log_b;
    log_b = tmp;
  }
  double negdelta = log_b - log_a;
  if ((negdelta < log_double_min) || std::isnan(negdelta))
    return log_a;
  else
    return (log_a + log1p(std::exp(negdelta)));
}

}
