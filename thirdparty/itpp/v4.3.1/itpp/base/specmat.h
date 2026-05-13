/*!
 * \file
 * \brief Definitions of special vectors and matrices
 * \author Tony Ottosson, Tobias Ringstrom, Pal Frenger, Adam Piatyszek
 *         and Erik G. Larsson
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

#ifndef SPECMAT_H
#define SPECMAT_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/converters.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \brief Return a integer vector with indicies where bvec == 1
  \ingroup miscfunc
*/
ITPP_EXPORT ivec find(const bvec &invector);

/*!
  \addtogroup specmat
*/

//!\addtogroup specmat
//!@{

//! A float vector of ones
ITPP_EXPORT vec ones(int size);
//! A Binary vector of ones
ITPP_EXPORT bvec ones_b(int size);
//! A Int vector of ones
ITPP_EXPORT ivec ones_i(int size);
//! A float Complex vector of ones
ITPP_EXPORT cvec ones_c(int size);

//! A float (rows,cols)-matrix of ones
ITPP_EXPORT mat ones(int rows, int cols);
//! A Binary (rows,cols)-matrix of ones
ITPP_EXPORT bmat ones_b(int rows, int cols);
//! A Int (rows,cols)-matrix of ones
ITPP_EXPORT imat ones_i(int rows, int cols);
//! A Double Complex (rows,cols)-matrix of ones
ITPP_EXPORT cmat ones_c(int rows, int cols);

//! A Double vector of zeros
ITPP_EXPORT vec zeros(int size);
//! A Binary vector of zeros
ITPP_EXPORT bvec zeros_b(int size);
//! A Int vector of zeros
ITPP_EXPORT ivec zeros_i(int size);
//! A Double Complex vector of zeros
ITPP_EXPORT cvec zeros_c(int size);

//! A Double (rows,cols)-matrix of zeros
ITPP_EXPORT mat zeros(int rows, int cols);
//! A Binary (rows,cols)-matrix of zeros
ITPP_EXPORT bmat zeros_b(int rows, int cols);
//! A Int (rows,cols)-matrix of zeros
ITPP_EXPORT imat zeros_i(int rows, int cols);
//! A Double Complex (rows,cols)-matrix of zeros
ITPP_EXPORT cmat zeros_c(int rows, int cols);

//! A Double (size,size) unit matrix
ITPP_EXPORT mat eye(int size);
//! A Binary (size,size) unit matrix
ITPP_EXPORT bmat eye_b(int size);
//! A Int (size,size) unit matrix
ITPP_EXPORT imat eye_i(int size);
//! A Double Complex (size,size) unit matrix
ITPP_EXPORT cmat eye_c(int size);
//! A non-copying version of the eye function.
template <class T>
void eye(int size, Mat<T> &m)
{
    m.set_size(size, size, false);
    m = T(0);
    for (int i = size - 1; i >= 0; i--)
        m(i, i) = T(1);
}

//! Impulse vector
ITPP_EXPORT vec impulse(int size);

//! linspace (works in the same way as the MATLAB version)
ITPP_EXPORT vec linspace(double from, double to, int length = 100);

//! linspace_fixed_step (works in the same way as "from:step:to" in MATLAB)
template<class T>
Vec<T> linspace_fixed_step(T from, T to, T step = 1)
{
    int points = 0;
    if (0 != step) {
        points = itpp::floor_i(double(to-from)/step)+1;
    }
    if (0 >= points) {
        return Vec<T>(0);
    }

    Vec<T> output(points);
    output(0) = from;
    for (int n = 1; n < points; ++n) {
    	output(n) = output(n-1)+step;
    }
    return output;
}

/*! \brief Zig-zag space function (variation on linspace)

This function is a variation on linspace().  It traverses the points
in different order. For example
\code
zigzag_space(-5,5,3)
\endcode
gives the vector
\code
[-5 5 0 -2.5 2.5 -3.75 -1.25 1.25 3.75]
\endcode
and
\code
zigzag_space(-5,5,4)
\endcode
gives
the vector
\code
[-5 5 0 -2.5 2.5 -3.75 -1.25 1.25 3.75 -4.375 -3.125 -1.875 -0.625 0.625 1.875 3.125 4.375]
\endcode
and so on.

I.e. the function samples the interval [t0,t1] with finer and finer
density and with points uniformly distributed over the interval,
rather than from left to right (as does linspace).

The result is a vector of length 1+2^K.
*/
ITPP_EXPORT vec zigzag_space(double t0, double t1, int K = 5);

/*!
 * \brief Hadamard matrix
 *
 * This function constructs a \a size by \a size Hadammard matrix, where
 * \a size is a power of 2.
 */
ITPP_EXPORT imat hadamard(int size);

/*!
  \brief Jacobsthal matrix.

  Constructs an p by p matrix Q where p is a prime (not checked).
  The elements in Q {qij} is given by qij=X(j-i), where X(x) is the
  Legendre symbol given as:

  <ul>
  <li> X(x)=0 if x is a multiple of p, </li>
  <li> X(x)=1 if x is a quadratic residue modulo p, </li>
  <li> X(x)=-1 if x is a quadratic nonresidue modulo p. </li>
  </ul>

  See Wicker "Error Control Systems for digital communication and storage", p. 134
  for more information on these topics. Do not check that p is a prime.
*/
ITPP_EXPORT imat jacobsthal(int p);

/*!
  \brief Conference matrix.

  Constructs an n by n matrix C, where n=p^m+1=2 (mod 4) and p is a odd prime (not checked).
  This code only work with m=1, that is n=p+1 and p odd prime. The valid sizes
  of n is then n=6, 14, 18, 30, 38, ... (and not 10, 26, ...).
  C has the property that C*C'=(n-1)I, that is it has orthogonal rows and columns
  in the same way as Hadamard matricies. However, one element in each row (on the
  diagonal) is zeros. The others are {-1,+1}.

  For more details see pp. 55-58 in MacWilliams & Sloane "The theory of error correcting codes",
  North-Holland, 1977.
*/
ITPP_EXPORT imat conference(int n);

/*!
 * \brief Generate Toeplitz matrix from two vectors \c c and \c r.
 *
 * Returns the Toeplitz matrix constructed given the first column C, and
 * (optionally) the first row R. If the first element of C is not the same
 * as the first element of R, the first element of C is used. If the second
 * argument is omitted, the first row is taken to be the same as the first
 * column and a symmetric (Hermitian) Toeplitz matrix is created.
 *
 * An example square Toeplitz matrix has the form:
 * \verbatim
 *       c(0)    r(1)     r(2)   ...   r(n)
 *       c(1)    c(0)     r(1)        r(n-1)
 *       c(2)    c(1)     c(0)        r(n-2)
 *        .                             .
 *        .                             .
 *        .                             .
 *
 *       c(n)   c(n-1)   c(n-2)  ...   c(0)
 * \endverbatim
 *
 * \author Adam Piatyszek
 */
template <typename Num_T>
const Mat<Num_T> toeplitz(const Vec<Num_T> &c, const Vec<Num_T> &r)
{
    int n_rows = c.size();
    int n_cols = r.size();
    Mat<Num_T> output(n_rows, n_cols);
    for (int i = 0; i < n_rows; ++i) {
        int j_limit = std::min(n_cols, n_rows - i);
        for (int j = 0; j < j_limit; ++j) {
            output(i + j, j) = c(i);
        }
    }
    for (int j = 1; j < n_cols; ++j) {
        int i_limit = std::min(n_rows, n_cols - j);
        for (int i = 0; i < i_limit; ++i) {
            output(i, i + j) = r(j);
        }
    }
    return output;
}

//! Generate symmetric Toeplitz matrix from vector \c c.
template <typename Num_T>
const Mat<Num_T> toeplitz(const Vec<Num_T> &c)
{
    int s = c.size();
    Mat<Num_T> output(s, s);
    for (int i = 0; i < s; ++i) {
        for (int j = 0; j < s - i; ++j) {
            output(i + j, j) = c(i);
        }
    }
    for (int j = 1; j < s; ++j) {
        for (int i = 0; i < s - j; ++i) {
            output(i, i + j) = c(j);
        }
    }
    return output;
}

//! Generate symmetric Toeplitz matrix from vector \c c (complex valued)
ITPP_EXPORT const cmat toeplitz(const cvec &c);


//!@}


/*!
  \brief Create a rotation matrix that rotates the given plane \c angle radians. Note that the order of the planes are important!
  \ingroup miscfunc
*/
ITPP_EXPORT mat rotation_matrix(int dim, int plane1, int plane2, double angle);

/*!
  \brief Calcualte the Householder vector
  \ingroup miscfunc
*/
ITPP_EXPORT void house(const vec &x, vec &v, double &beta);

/*!
  \brief Calculate the Givens rotation values
  \ingroup miscfunc
*/
ITPP_EXPORT void givens(double a, double b, double &c, double &s);

/*!
  \brief Calculate the Givens rotation matrix
  \ingroup miscfunc
*/
ITPP_EXPORT void givens(double a, double b, mat &m);

/*!
  \brief Calculate the Givens rotation matrix
  \ingroup miscfunc
*/
ITPP_EXPORT mat givens(double a, double b);

/*!
  \brief Calculate the transposed Givens rotation matrix
  \ingroup miscfunc
*/
ITPP_EXPORT void givens_t(double a, double b, mat &m);

/*!
  \brief Calculate the transposed Givens rotation matrix
  \ingroup miscfunc
*/
ITPP_EXPORT mat givens_t(double a, double b);

/*!
  \relates Vec
  \brief Vector of length 1
*/
template <class T>
Vec<T> vec_1(T v0)
{
    Vec<T> v(1);
    v(0) = v0;
    return v;
}

/*!
  \relates Vec
  \brief Vector of length 2
*/
template <class T>
Vec<T> vec_2(T v0, T v1)
{
    Vec<T> v(2);
    v(0) = v0;
    v(1) = v1;
    return v;
}

/*!
  \relates Vec
  \brief Vector of length 3
*/
template <class T>
Vec<T> vec_3(T v0, T v1, T v2)
{
    Vec<T> v(3);
    v(0) = v0;
    v(1) = v1;
    v(2) = v2;
    return v;
}

/*!
  \relates Mat
  \brief Matrix of size 1 by 1
*/
template <class T>
Mat<T> mat_1x1(T m00)
{
    Mat<T> m(1, 1);
    m(0, 0) = m00;
    return m;
}

/*!
  \relates Mat
  \brief Matrix of size 1 by 2
*/
template <class T>
Mat<T> mat_1x2(T m00, T m01)
{
    Mat<T> m(1, 2);
    m(0, 0) = m00;
    m(0, 1) = m01;
    return m;
}

/*!
  \relates Mat
  \brief Matrix of size 2 by 1
*/
template <class T>
Mat<T> mat_2x1(T m00,
               T m10)
{
    Mat<T> m(2, 1);
    m(0, 0) = m00;
    m(1, 0) = m10;
    return m;
}

/*!
  \relates Mat
  \brief Matrix of size 2 by 2
*/
template <class T>
Mat<T> mat_2x2(T m00, T m01,
               T m10, T m11)
{
    Mat<T> m(2, 2);
    m(0, 0) = m00;
    m(0, 1) = m01;
    m(1, 0) = m10;
    m(1, 1) = m11;
    return m;
}

/*!
  \relates Mat
  \brief Matrix of size 1 by 3
*/
template <class T>
Mat<T> mat_1x3(T m00, T m01, T m02)
{
    Mat<T> m(1, 3);
    m(0, 0) = m00;
    m(0, 1) = m01;
    m(0, 2) = m02;
    return m;
}

/*!
  \relates Mat
  \brief Matrix of size 3 by 1
*/
template <class T>
Mat<T> mat_3x1(T m00,
               T m10,
               T m20)
{
    Mat<T> m(3, 1);
    m(0, 0) = m00;
    m(1, 0) = m10;
    m(2, 0) = m20;
    return m;
}

/*!
  \relates Mat
  \brief Matrix of size 2 by 3
*/
template <class T>
Mat<T> mat_2x3(T m00, T m01, T m02,
               T m10, T m11, T m12)
{
    Mat<T> m(2, 3);
    m(0, 0) = m00;
    m(0, 1) = m01;
    m(0, 2) = m02;
    m(1, 0) = m10;
    m(1, 1) = m11;
    m(1, 2) = m12;
    return m;
}

/*!
  \relates Mat
  \brief Matrix of size 3 by 2
*/
template <class T>
Mat<T> mat_3x2(T m00, T m01,
               T m10, T m11,
               T m20, T m21)
{
    Mat<T> m(3, 2);
    m(0, 0) = m00;
    m(0, 1) = m01;
    m(1, 0) = m10;
    m(1, 1) = m11;
    m(2, 0) = m20;
    m(2, 1) = m21;
    return m;
}

/*!
  \relates Mat
  \brief Matrix of size 3 by 3
*/
template <class T>
Mat<T> mat_3x3(T m00, T m01, T m02,
               T m10, T m11, T m12,
               T m20, T m21, T m22)
{
    Mat<T> m(3, 3);
    m(0, 0) = m00;
    m(0, 1) = m01;
    m(0, 2) = m02;
    m(1, 0) = m10;
    m(1, 1) = m11;
    m(1, 2) = m12;
    m(2, 0) = m20;
    m(2, 1) = m21;
    m(2, 2) = m22;
    return m;
}

} //namespace itpp

#endif // #ifndef SPECMAT_H
