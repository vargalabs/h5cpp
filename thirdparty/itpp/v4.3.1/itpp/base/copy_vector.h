/*!
 * \file
 * \brief Vector copy functions for internal use
 * \author Tony Ottosson and Adam Piatyszek
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

#ifndef COPY_VECTOR_H
#define COPY_VECTOR_H

#include <itpp/base/binary.h>
#include <complex>
#include <cstring>
#include <itpp/itexports.h>

//! \cond

namespace itpp
{


/*
  Copy vector x to vector y. Both vectors are of size n
*/
inline void copy_vector(int n, const int *x, int *y)
{
  memcpy(y, x, n * sizeof(int));
}
inline void copy_vector(int n, const short *x, short *y)
{
  memcpy(y, x, n * sizeof(short));
}
inline void copy_vector(int n, const bin *x, bin *y)
{
  memcpy(y, x, n * sizeof(bin));
}

ITPP_EXPORT void copy_vector(int n, const double *x, double *y);
ITPP_EXPORT void copy_vector(int n, const std::complex<double> *x,
                 std::complex<double> *y);

template<class T> inline
void copy_vector(int n, const T *x, T *y)
{
  for (int i = 0; i < n; i++)
    y[i] = x[i];
}


/*
  Copy vector x to vector y. Both vectors are of size n
  vector x elements are stored linearly with element increament incx
  vector y elements are stored linearly with element increament incx
*/
ITPP_EXPORT void copy_vector(int n, const double *x, int incx, double *y, int incy);
ITPP_EXPORT void copy_vector(int n, const std::complex<double> *x, int incx,
                 std::complex<double> *y, int incy);

template<class T> inline
void copy_vector(int n, const T *x, int incx, T *y, int incy)
{
  for (int i = 0; i < n; i++)
    y[i*incy] = x[i*incx];
}


/*
  Swap vector x and vector y. Both vectors are of size n
*/
inline void swap_vector(int n, int *x, int *y)
{
  for (int i = 0; i < n; i++)
    std::swap(x[i], y[i]);
}
inline void swap_vector(int n, short *x, short *y)
{
  for (int i = 0; i < n; i++)
    std::swap(x[i], y[i]);
}
inline void swap_vector(int n, bin *x, bin *y)
{
  for (int i = 0; i < n; i++)
    std::swap(x[i], y[i]);
}

ITPP_EXPORT void swap_vector(int n, double *x, double *y);
ITPP_EXPORT void swap_vector(int n, std::complex<double> *x, std::complex<double> *y);

template<class T> inline
void swap_vector(int n, T *x, T *y)
{
  T tmp;
  for (int i = 0; i < n; i++) {
    tmp = y[i];
    y[i] = x[i];
    x[i] = tmp;
  }
}


/*
  Swap vector x and vector y. Both vectors are of size n
  vector x elements are stored linearly with element increament incx
  vector y elements are stored linearly with element increament incx
*/
inline void swap_vector(int n, int *x, int incx, int *y, int incy)
{
  for (int i = 0; i < n; i++)
    std::swap(x[i*incx], y[i*incy]);
}
inline void swap_vector(int n, short *x, int incx, short *y, int incy)
{
  for (int i = 0; i < n; i++)
    std::swap(x[i*incx], y[i*incy]);
}
inline void swap_vector(int n, bin *x, int incx, bin *y, int incy)
{
  for (int i = 0; i < n; i++)
    std::swap(x[i*incx], y[i*incy]);
}

ITPP_EXPORT void swap_vector(int n, double *x, int incx, double *y, int incy);
ITPP_EXPORT void swap_vector(int n, std::complex<double> *x, int incx,
                 std::complex<double> *y, int incy);

template<class T> inline
void swap_vector(int n, T *x, int incx, T *y, int incy)
{
  T tmp;
  for (int i = 0; i < n; i++) {
    tmp = y[i*incy];
    y[i*incy] = x[i*incx];
    x[i*incx] = tmp;
  }
}


/*
 * Realise scaling operation: x = alpha*x
 */
ITPP_EXPORT void scal_vector(int n, double alpha, double *x);
ITPP_EXPORT void scal_vector(int n, std::complex<double> alpha, std::complex<double> *x);

template<typename T> inline
void scal_vector(int n, T alpha, T *x)
{
  if (alpha != T(1)) {
    for (int i = 0; i < n; ++i) {
      x[i] *= alpha;
    }
  }
}


/*
 * Realise scaling operation: x = alpha*x
 * Elements of x are stored linearly with increament incx
 */
ITPP_EXPORT void scal_vector(int n, double alpha, double *x, int incx);
ITPP_EXPORT void scal_vector(int n, std::complex<double> alpha, std::complex<double> *x,
                 int incx);

template<typename T> inline
void scal_vector(int n, T alpha, T *x, int incx)
{
  if (alpha != T(1)) {
    for (int i = 0; i < n; ++i) {
      x[i*incx] *= alpha;
    }
  }
}


/*
 * Realise the following equation on vectors: y = alpha*x + y
 */
ITPP_EXPORT void axpy_vector(int n, double alpha, const double *x, double *y);

ITPP_EXPORT void axpy_vector(int n, std::complex<double> alpha,
                 const std::complex<double> *x, std::complex<double> *y);

template<typename T> inline
void axpy_vector(int n, T alpha, const T *x, T *y)
{
  if (alpha != T(1)) {
    for (int i = 0; i < n; ++i) {
      y[i] += alpha * x[i];
    }
  }
  else {
    for (int i = 0; i < n; ++i) {
      y[i] += x[i];
    }
  }
}


/*
 * Realise the following equation on vectors: y = alpha*x + y
 * Elements of x are stored linearly with increment incx
 * and elements of y are stored linearly with increment incx
 */
ITPP_EXPORT void axpy_vector(int n, double alpha, const double *x, int incx, double *y,
                 int incy);
ITPP_EXPORT void axpy_vector(int n, std::complex<double> alpha,
                 const std::complex<double> *x, int incx,
                 std::complex<double> *y, int incy);

template<typename T> inline
void axpy_vector(int n, T alpha, const T *x, int incx, T *y, int incy)
{
  if (alpha != T(1)) {
    for (int i = 0; i < n; ++i) {
      y[i*incy] += alpha * x[i*incx];
    }
  }
  else {
    for (int i = 0; i < n; ++i) {
      y[i*incy] += x[i*incx];
    }
  }
}


} // namespace itpp

//! \endcond

#endif // #ifndef COPY_VECTOR_H
