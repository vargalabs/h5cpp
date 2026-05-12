/*!
 * \file
 * \brief Definitions of Schur decomposition functions
 * \author Adam Piatyszek
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

#ifndef SCHUR_H
#define SCHUR_H

#include <itpp/base/mat.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
 * \ingroup matrixdecomp
 * \brief Schur decomposition of a real matrix
 *
 * This function computes the Schur form of a square real matrix
 * \f$ \mathbf{A} \f$. The Schur decomposition satisfies the
 * following equation:
 * \f[ \mathbf{U} \mathbf{T} \mathbf{U}^{T} = \mathbf{A} \f]
 * where: \f$ \mathbf{U} \f$ is a unitary, \f$ \mathbf{T} \f$ is upper
 * quasi-triangular, and \f$ \mathbf{U}^{T} \f$ is the transposed
 * \f$ \mathbf{U} \f$ matrix.
 *
 * The upper quasi-triangular matrix may have \f$ 2 \times 2 \f$ blocks on
 * its diagonal.
 *
 * Uses the LAPACK routine DGEES.
 */
ITPP_EXPORT bool schur(const mat &A, mat &U, mat &T);

/*!
 * \ingroup matrixdecomp
 * \brief Schur decomposition of a real matrix
 *
 * This function computes the Schur form of a square real matrix
 * \f$ \mathbf{A} \f$. The Schur decomposition satisfies the
 * following equation:
 * \f[ \mathbf{U} \mathbf{T} \mathbf{U}^{T} = \mathbf{A} \f]
 * where: \f$ \mathbf{U} \f$ is a unitary, \f$ \mathbf{T} \f$ is upper
 * quasi-triangular, and \f$ \mathbf{U}^{T} \f$ is the transposed
 * \f$ \mathbf{U} \f$ matrix.
 *
 * The upper quasi-triangular matrix may have \f$ 2 \times 2 \f$ blocks on
 * its diagonal.
 *
 * \return  Real Schur matrix \f$ \mathbf{T} \f$
 *
 * uses the LAPACK routine DGEES.
 */
ITPP_EXPORT mat schur(const mat &A);


/*!
 * \ingroup matrixdecomp
 * \brief Schur decomposition of a complex matrix
 *
 * This function computes the Schur form of a square complex matrix
 * \f$ \mathbf{A} \f$. The Schur decomposition satisfies
 * the following equation:
 * \f[ \mathbf{U} \mathbf{T} \mathbf{U}^{H} = \mathbf{A} \f]
 * where: \f$ \mathbf{U} \f$ is a unitary, \f$ \mathbf{T} \f$ is upper
 * triangular, and \f$ \mathbf{U}^{H} \f$ is the Hermitian
 * transposition of the \f$ \mathbf{U} \f$ matrix.
 *
 * Uses the LAPACK routine ZGEES.
 */
ITPP_EXPORT bool schur(const cmat &A, cmat &U, cmat &T);

/*!
 * \ingroup matrixdecomp
 * \brief Schur decomposition of a complex matrix
 *
 * This function computes the Schur form of a square complex matrix
 * \f$ \mathbf{A} \f$. The Schur decomposition satisfies
 * the following equation:
 * \f[ \mathbf{U} \mathbf{T} \mathbf{U}^{H} = \mathbf{A} \f]
 * where: \f$ \mathbf{U} \f$ is a unitary, \f$ \mathbf{T} \f$ is upper
 * triangular, and \f$ \mathbf{U}^{H} \f$ is the Hermitian
 * transposition of the \f$ \mathbf{U} \f$ matrix.
 *
 * \return  Complex Schur matrix \f$ \mathbf{T} \f$
 *
 * Uses the LAPACK routine ZGEES.
 */
ITPP_EXPORT cmat schur(const cmat &A);


} // namespace itpp

#endif // #ifndef SCHUR_H
