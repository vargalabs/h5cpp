/*!
 * \file
 * \brief Definitions of Cholesky factorisation functions
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

#ifndef CHOLESKY_H
#define CHOLESKY_H

#include <itpp/base/mat.h>
#include <itpp/itexports.h>

namespace itpp
{

/*! \addtogroup matrixdecomp
 */
//!@{

/*!
  \brief Cholesky factorisation of real symmetric and positive definite matrix

  The Cholesky factorisation of a real symmetric positive-definite matrix \f$\mathbf{X}\f$
  of size \f$n \times n\f$ is given by
  \f[
  \mathbf{X} = \mathbf{F}^T \mathbf{F}
  \f]
  where \f$\mathbf{F}\f$ is an upper triangular \f$n \times n\f$ matrix.

  Returns true if calculation succeeded. False otherwise.
*/
ITPP_EXPORT bool chol(const mat &X, mat &F);

/*!
  \brief Cholesky factorisation of real symmetric and positive definite matrix

  The Cholesky factorisation of a real symmetric positive-definite matrix \f$\mathbf{X}\f$
  of size \f$n \times n\f$ is given by
  \f[
  \mathbf{X} = \mathbf{F}^T \mathbf{F}
  \f]
  where \f$\mathbf{F}\f$ is an upper triangular \f$n \times n\f$ matrix.
*/
ITPP_EXPORT mat chol(const mat &X);


/*!
  \brief Cholesky factorisation of complex hermitian and positive-definite matrix

  The Cholesky factorisation of a hermitian positive-definite matrix \f$\mathbf{X}\f$
  of size \f$n \times n\f$ is given by
  \f[
  \mathbf{X} = \mathbf{F}^H \mathbf{F}
  \f]
  where \f$\mathbf{F}\f$ is an upper triangular \f$n \times n\f$ matrix.

  Returns true if calculation succeeded. False otherwise.

  If \c X is positive definite, true is returned and \c F=chol(X)
  produces an upper triangular \c F. If also \c X is symmetric then \c F'*F = X.
  If \c X is not positive definite, false is returned.
*/
ITPP_EXPORT bool chol(const cmat &X, cmat &F);

/*!
  \brief Cholesky factorisation of complex hermitian and positive-definite matrix

  The Cholesky factorisation of a hermitian positive-definite matrix \f$\mathbf{X}\f$
  of size \f$n \times n\f$ is given by
  \f[
  \mathbf{X} = \mathbf{F}^H \mathbf{F}
  \f]
  where \f$\mathbf{F}\f$ is an upper triangular \f$n \times n\f$ matrix.
*/
ITPP_EXPORT cmat chol(const cmat &X);

//!@}

} // namespace itpp

#endif // #ifndef CHOLESKY_H
