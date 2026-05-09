/*!
 * \file
 * \brief Definitions of eigenvalue decomposition functions
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

#ifndef EIGEN_H
#define EIGEN_H

#include <itpp/base/mat.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \ingroup matrixdecomp
  \brief Calculates the eigenvalues and eigenvectors of a symmetric real matrix

  The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
  \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the real and symmetric \f$n \times n\f$
  matrix \f$\mathbf{A}\f$ satisfies
  \f[
  \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
  \f]
  The eigenvectors are the columns of the matrix V.
  True is returned if the calculation was successful. Otherwise false.

  Uses the LAPACK routine DSYEV.
*/
ITPP_EXPORT bool eig_sym(const mat &A, vec &d, mat &V);

/*!
  \ingroup matrixdecomp
  \brief Calculates the eigenvalues of a symmetric real matrix

  The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
  \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the real and symmetric \f$n \times n\f$
  matrix \f$\mathbf{A}\f$ satisfies
  \f[
  \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
  \f]
  True is returned if the calculation was successful. Otherwise false.

  Uses the LAPACK routine DSYEV.
*/
ITPP_EXPORT bool eig_sym(const mat &A, vec &d);

/*!
  \ingroup matrixdecomp
  \brief Calculates the eigenvalues of a symmetric real matrix

  The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
  \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the real and symmetric \f$n \times n\f$
  matrix \f$\mathbf{A}\f$ satisfies
  \f[
  \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
  \f]

  Uses the LAPACK routine DSYEV.
*/
ITPP_EXPORT vec eig_sym(const mat &A);

/*!
  \ingroup matrixdecomp
  \brief Calculates the eigenvalues and eigenvectors of a hermitian complex matrix

  The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
  \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the complex and hermitian \f$n \times n\f$
  matrix \f$\mathbf{A}\f$ satisfies
  \f[
  \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
  \f]
  The eigenvectors are the columns of the matrix V.
  True is returned if the calculation was successful. Otherwise false.

  Uses the LAPACK routine ZHEEV.
*/
ITPP_EXPORT bool eig_sym(const cmat &A, vec &d, cmat &V);

/*!
  \ingroup matrixdecomp
  \brief Calculates the eigenvalues of a hermitian complex matrix

  The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
  \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the complex and hermitian \f$n \times n\f$
  matrix \f$\mathbf{A}\f$ satisfies
  \f[
  \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
  \f]
  True is returned if the calculation was successful. Otherwise false.

  Uses the LAPACK routine ZHEEV.
*/
ITPP_EXPORT bool eig_sym(const cmat &A, vec &d);

/*!
  \ingroup matrixdecomp
  \brief Calculates the eigenvalues of a hermitian complex matrix

  The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
  \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the complex and hermitian \f$n \times n\f$
  matrix \f$\mathbf{A}\f$ satisfies
  \f[
  \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
  \f]

  Uses the LAPACK routine ZHEEV.
*/
ITPP_EXPORT vec eig_sym(const cmat &A);

/*!
  \ingroup matrixdecomp
  \brief Calculates the eigenvalues and eigenvectors of a real non-symmetric matrix

  The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
  \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the real \f$n \times n\f$
  matrix \f$\mathbf{A}\f$ satisfies
  \f[
  \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
  \f]
  The eigenvectors are the columns of the matrix V.
  True is returned if the calculation was successful. Otherwise false.

  Uses the LAPACK routine DGEEV.
*/
ITPP_EXPORT bool eig(const mat &A, cvec &d, cmat &V);

/*!
  \ingroup matrixdecomp
  \brief Calculates the eigenvalues of a real non-symmetric matrix

  The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
  \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the real \f$n \times n\f$
  matrix \f$\mathbf{A}\f$ satisfies
  \f[
  \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
  \f]
  True is returned if the calculation was successful. Otherwise false.

  Uses the LAPACK routine DGEEV.
*/
ITPP_EXPORT bool eig(const mat &A, cvec &d);

/*!
  \ingroup matrixdecomp
  \brief Calculates the eigenvalues of a real non-symmetric matrix

  The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
  \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the real \f$n \times n\f$
  matrix \f$\mathbf{A}\f$ satisfies
  \f[
  \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
  \f]

  Uses the LAPACK routine DGEEV.
*/
ITPP_EXPORT cvec eig(const mat &A);

/*!
  \ingroup matrixdecomp
  \brief Calculates the eigenvalues and eigenvectors of a complex non-hermitian matrix

  The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
  \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the complex \f$n \times n\f$
  matrix \f$\mathbf{A}\f$ satisfies
  \f[
  \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
  \f]
  The eigenvectors are the columns of the matrix V.
  True is returned if the calculation was successful. Otherwise false.

  Uses the LAPACK routine ZGEEV.
*/
ITPP_EXPORT bool eig(const cmat &A, cvec &d, cmat &V);

/*!
  \ingroup matrixdecomp
  \brief Calculates the eigenvalues of a complex non-hermitian matrix

  The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
  \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the complex \f$n \times n\f$
  matrix \f$\mathbf{A}\f$ satisfies
  \f[
  \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
  \f]
  True is returned if the calculation was successful. Otherwise false.

  Uses the LAPACK routine ZGEEV.
*/
ITPP_EXPORT bool eig(const cmat &A, cvec &d);

/*!
  \ingroup matrixdecomp
  \brief Calculates the eigenvalues of a complex non-hermitian matrix

  The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
  \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the complex \f$n \times n\f$
  matrix \f$\mathbf{A}\f$ satisfies
  \f[
  \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
  \f]

  Uses the LAPACK routine ZGEEV.
*/
ITPP_EXPORT cvec eig(const cmat &A);

} // namespace itpp

#endif // #ifndef EIGEN_H
