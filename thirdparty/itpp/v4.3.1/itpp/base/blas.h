/*!
 * \file
 * \brief BLAS header functions. For internal use only.
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

#ifndef BLAS_H
#define BLAS_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#include <complex>

//! \cond

#if defined(_MSC_VER) && (defined(HAVE_ACML) || defined(HAVE_MKL))
#  define dswap_ DSWAP
#  define zswap_ ZSWAP
#  define dscal_ DSCAL
#  define zscal_ ZSCAL
#  define dcopy_ DCOPY
#  define zcopy_ ZCOPY
#  define daxpy_ DAXPY
#  define zaxpy_ ZAXPY
#  define ddot_  DDOT
#  define dgemv_ DGEMV
#  define zgemv_ ZGEMV
#  define dger_  DGER
#  define zgeru_ ZGERU
#  define zgerc_ ZGERC
#  define dgemm_ DGEMM
#  define zgemm_ ZGEMM
#endif // #if defined(_MSC_VER) && (defined(HAVE_ACML) || defined(HAVE_MKL))


namespace blas
{

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

  // ----------------------------------------------------------------------
  // BLAS 1 functions
  // ----------------------------------------------------------------------

  void dswap_(const int *n,
              double *x, const int *incx,
              double *y, const int *incy);

  void zswap_(const int *n,
              std::complex<double> *x, const int *incx,
              std::complex<double> *y, const int *incy);

  void dscal_(const int *n,
              const double *alpha,
              double *x, const int *incx);

  void zscal_(const int *n,
              const std::complex<double> *alpha,
              std::complex<double> *x, const int *incx);

  void dcopy_(const int *n,
              const double *x, const int *incx,
              double *y, const int *incy);

  void zcopy_(const int *n,
              const std::complex<double> *x, const int *incx,
              std::complex<double> *y, const int *incy);

  void daxpy_(const int *n,
              const double *alpha,
              const double *x, const int *incx,
              double *y, const int *incy);

  void zaxpy_(const int *n,
              const std::complex<double> *alpha,
              const std::complex<double> *x, const int *incx,
              std::complex<double> *y, const int *incy);

  double ddot_(const int *n,
               const double *x, const int *incx,
               const double *y, const int *incy);

  // ----------------------------------------------------------------------
  // BLAS 2 functions
  // ----------------------------------------------------------------------

  void dgemv_(const char *transA, const int *m, const int *n,
              const double *alpha,
              const double *A, const int *ldA,
              const double *x, const int *incx,
              const double *beta,
              double *y, const int *incy);

  void zgemv_(const char *transA, const int *m, const int *n,
              const std::complex<double> *alpha,
              const std::complex<double> *A, const int *ldA,
              const std::complex<double> *x, const int *incx,
              const std::complex<double> *beta,
              std::complex<double> *y, const int *incy);

  void dger_(const int *m, const int *n,
             const double *alpha,
             const double *x, const int *incx,
             const double *y, const int *incy,
             double *A, const int *ldA);

  void zgeru_(const int *m, const int *n,
              const std::complex<double> *alpha,
              const std::complex<double> *x, const int *inxx,
              const std::complex<double> *y, const int *incy,
              std::complex<double> *A, const int *ldA);

  void zgerc_(const int *m, const int *n,
              const std::complex<double> *alpha,
              const std::complex<double> *x, const int *inxx,
              const std::complex<double> *y, const int *incy,
              std::complex<double> *A, const int *ldA);

  // ----------------------------------------------------------------------
  // BLAS 3 functions
  // ----------------------------------------------------------------------

  void dgemm_(const char *transA, const char *transB,
              const int *m, const int *n, const int *k,
              const double *alpha,
              const double *A, const int *ldA,
              const double *B, const int *ldB,
              const double *beta,
              double *C, const int *ldC);

  void zgemm_(const char *transA, const char *transB,
              const int *m, const int *n, const int *k,
              const std::complex<double> *alpha,
              const std::complex<double> *A, const int *ldA,
              const std::complex<double> *B, const int *ldB,
              const std::complex<double> *beta,
              std::complex<double> *C, const int *ldC);

#ifdef __cplusplus
}
#endif /* __cplusplus */

} // namespace blas

//! \endcond

#endif // #ifndef BLAS_H
