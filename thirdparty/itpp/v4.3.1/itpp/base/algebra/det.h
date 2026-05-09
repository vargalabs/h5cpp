/*!
 * \file
 * \brief Definitions of determinant calculations
 * \author Tony Ottosson
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

#ifndef DET_H
#define DET_H

#include <itpp/base/mat.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \brief Determinant of real square matrix.
  \ingroup determinant

  Calculate determinant of the real matrix \f$\mathbf{X}\f$

  Uses LU-factorisation.
  \f[
  \det(\mathbf{X}) = \det(\mathbf{P}^T \mathbf{L}) \det(\mathbf{U}) = \det(\mathbf{P}^T) \prod(\mathrm{diag}(\mathbf{U}))
  \f]
  and the determinant of the permuation matrix is \f$ \pm 1\f$ depending on the number of row permutations
*/
ITPP_EXPORT double det(const mat &X);


/*!
  \brief Determinant of complex square matrix.
  \ingroup determinant

  Calculate determinant of the complex matrix \f$\mathbf{X}\f$

  Uses LU-factorisation.
  \f[
  \det(\mathbf{X}) = \det(\mathbf{P}^T \mathbf{L}) \det(\mathbf{U}) = \det(\mathbf{P}^T) \prod(\mathrm{diag}(\mathbf{U}))
  \f]
  and the determinant of the permuation matrix is \f$ \pm 1\f$ depending on the number of row permutations
*/
ITPP_EXPORT std::complex<double> det(const cmat &X);


} // namespace itpp

#endif // #ifndef DET_H
