/*!
 * \file
 * \brief Definitions of Singular Value Decompositions
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

#ifndef SVD_H
#define SVD_H

#include <itpp/base/mat.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
 * \ingroup matrixdecomp
 * \brief Get singular values \c s of a real matrix \c A using SVD
 *
 * This function calculates singular values \f$s\f$ from the SVD
 * decomposition of a real matrix \f$A\f$. The SVD algorithm computes the
 * decomposition of a real \f$m \times n\f$ matrix \f$\mathbf{A}\f$ so
 * that
 * \f[
 * \mathrm{diag}(\mathbf{U}^T \mathbf{A} \mathbf{V}) = \mathbf{s}
 * = \sigma_1, \ldots, \sigma_p
 * \f]
 * where \f$\sigma_1 \geq \sigma_2 \geq \ldots \sigma_p \geq 0\f$ are the
 * singular values of \f$\mathbf{A}\f$.
 * Or put differently:
 * \f[
 * \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^T
 * \f]
 * where \f$ \mathrm{diag}(\mathbf{S}) = \mathbf{s} \f$
 *
 * \note An external LAPACK library is required by this function.
 */
ITPP_EXPORT bool svd(const mat &A, vec &s);

/*!
 * \ingroup matrixdecomp
 * \brief Get singular values \c s of a complex matrix \c A using SVD
 *
 * This function calculates singular values \f$s\f$ from the SVD
 * decomposition of a complex matrix \f$A\f$. The SVD algorithm computes
 * the decomposition of a complex \f$m \times n\f$ matrix \f$\mathbf{A}\f$
 * so that
 * \f[
 * \mathrm{diag}(\mathbf{U}^H \mathbf{A} \mathbf{V}) = \mathbf{s}
 * = \sigma_1, \ldots, \sigma_p
 * \f]
 * where \f$\sigma_1 \geq \sigma_2 \geq \ldots \sigma_p \geq 0\f$
 * are the singular values of \f$\mathbf{A}\f$.
 * Or put differently:
 * \f[
 * \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^H
 * \f]
 * where \f$ \mathrm{diag}(\mathbf{S}) = \mathbf{s} \f$
 *
 * \note An external LAPACK library is required by this function.
 */
ITPP_EXPORT bool svd(const cmat &A, vec &s);

/*!
   * \ingroup matrixdecomp
   * \brief Return singular values of a real matrix \c A using SVD
   *
   * This function returns singular values from the SVD decomposition
   * of a real matrix \f$A\f$. The SVD algorithm computes the decomposition
   * of a real \f$m \times n\f$ matrix \f$\mathbf{A}\f$ so that
   * \f[
   * \mathrm{diag}(\mathbf{U}^T \mathbf{A} \mathbf{V}) = \mathbf{s}
   * = \sigma_1, \ldots, \sigma_p
   * \f]
   * where \f$\sigma_1 \geq \sigma_2 \geq \ldots \sigma_p \geq 0\f$ are the
   * singular values of \f$\mathbf{A}\f$.
   * Or put differently:
   * \f[
   * \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^T
   * \f]
   * where \f$ \mathrm{diag}(\mathbf{S}) = \mathbf{s} \f$
   *
   * \note An external LAPACK library is required by this function.
   */
ITPP_EXPORT vec svd(const mat &A);

/*!
 * \ingroup matrixdecomp
 * \brief Return singular values of a complex matrix \c A using SVD
 *
 * This function returns singular values from the SVD
 * decomposition of a complex matrix \f$A\f$. The SVD algorithm computes
 * the decomposition of a complex \f$m \times n\f$ matrix \f$\mathbf{A}\f$
 * so that
 * \f[
 * \mathrm{diag}(\mathbf{U}^H \mathbf{A} \mathbf{V}) = \mathbf{s}
 * = \sigma_1, \ldots, \sigma_p
 * \f]
 * where \f$\sigma_1 \geq \sigma_2 \geq \ldots \sigma_p \geq 0\f$
 * are the singular values of \f$\mathbf{A}\f$.
 * Or put differently:
 * \f[
 * \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^H
 * \f]
 * where \f$ \mathrm{diag}(\mathbf{S}) = \mathbf{s} \f$
 *
 * \note An external LAPACK library is required by this function.
 */
ITPP_EXPORT vec svd(const cmat &A);

/*!
 * \ingroup matrixdecomp
 * \brief Perform Singular Value Decomposition (SVD) of a real matrix \c A
 *
 * This function returns two orthonormal matrices \f$U\f$ and \f$V\f$
 * and a vector of singular values \f$s\f$.
 * The SVD algorithm computes the decomposition of a real \f$m \times n\f$
 * matrix \f$\mathbf{A}\f$ so that
 * \f[
 * \mathrm{diag}(\mathbf{U}^T \mathbf{A} \mathbf{V}) = \mathbf{s}
 * = \sigma_1, \ldots, \sigma_p
 * \f]
 * where the elements of \f$\mathbf{s}\f$, \f$\sigma_1 \geq \sigma_2 \geq
 * \ldots \sigma_p \geq 0\f$ are the singular values of \f$\mathbf{A}\f$.
 * Or put differently:
 * \f[
 * \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^T
 * \f]
 * where \f$ \mathrm{diag}(\mathbf{S}) = \mathbf{s} \f$
 *
 * \note An external LAPACK library is required by this function.
 */
ITPP_EXPORT bool svd(const mat &A, mat &U, vec &s, mat &V);

/*!
 * \ingroup matrixdecomp
 * \brief Perform Singular Value Decomposition (SVD) of a complex matrix \c A
 *
 * This function returns two orthonormal matrices \f$U\f$ and \f$V\f$
 * and a vector of singular values \f$s\f$.
 * The SVD algorithm computes the decomposition of a complex \f$m \times n\f$
 * matrix \f$\mathbf{A}\f$ so that
 * \f[
 * \mathrm{diag}(\mathbf{U}^H \mathbf{A} \mathbf{V}) = \mathbf{s}
 * = \sigma_1, \ldots, \sigma_p
 * \f]
 * where the elements of \f$\mathbf{s}\f$, \f$\sigma_1 \geq \sigma_2 \geq
 * \ldots \sigma_p \geq 0\f$ are the singular values of \f$\mathbf{A}\f$.
 * Or put differently:
 * \f[
 * \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^H
 * \f]
 * where \f$ \mathrm{diag}(\mathbf{S}) = \mathbf{s} \f$
 *
 * \note An external LAPACK library is required by this function.
 */
ITPP_EXPORT bool svd(const cmat &A, cmat &U, vec &s, cmat &V);


} // namespace itpp

#endif // #ifndef SVD_H
