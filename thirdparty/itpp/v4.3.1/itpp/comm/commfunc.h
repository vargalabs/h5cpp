/*!
 * \file
 * \brief Definitions of some specific functions useful in communications
 * \author Tony Ottosson and Erik G. Larsson
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

#ifndef COMMFUNC_H
#define COMMFUNC_H

#include <itpp/base/mat.h>
#include <itpp/base/vec.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \brief Generate Gray code of blocklength m.
  \ingroup misccommfunc

  The codes are contained as binary codewords {0,1} in the rows of the
  returned matrix.
  See also the \c gray() function in \c math/scalfunc.h.
*/
ITPP_EXPORT bmat graycode(int m);

/*!
  \brief Calculate the Hamming distance between \a a and \a b
  \ingroup misccommfunc
*/
ITPP_EXPORT int hamming_distance(const bvec &a, const bvec &b);

/*!
  \brief Calculate the Hamming weight of \a a
  \ingroup misccommfunc
*/
ITPP_EXPORT int weight(const bvec &a);

/*!
 * \brief Compute the water-filling solution
 * \ingroup misccommfunc
 *
 * This function computes the solution of the water-filling problem
 * \f[
 * \max_{p_0,...,p_{n-1}} \sum_{i=0}^{n-1} \log\left(1+p_i\alpha_i\right)
 * \f]
 * subject to
 * \f[
 * \sum_{i=0}^{n-1} p_i \le P
 * \f]
 *
 * \param alpha vector of \f$\alpha_0,...,\alpha_{n-1}\f$ gains (must have
 * strictly positive elements)
 * \param P power constraint
 * \return vector of power allocations \f$p_0,...,p_{n-1}\f$
 *
 * The computational complexity of the method is \f$O(n^2)\f$ at most
 */
ITPP_EXPORT vec waterfilling(const vec& alpha, double P);

} // namespace itpp

#endif // #ifndef COMMFUNC_H
