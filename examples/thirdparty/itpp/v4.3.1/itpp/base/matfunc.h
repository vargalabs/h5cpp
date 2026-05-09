/*!
 * \file
 * \brief Various functions on vectors and matrices - header file
 * \author Tony Ottosson, Adam Piatyszek, Conrad Sanderson, Mark Dobossy
 *         and Martin Senst
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

#ifndef MATFUNC_H
#define MATFUNC_H

#include <itpp/base/mat.h>
#include <itpp/base/math/log_exp.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/algebra/inv.h>
#include <itpp/base/algebra/svd.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \addtogroup matrix_functions
  \brief Functions on vectors and matrices
*/
//!@{

//! Length of vector
template<class T>
int length(const Vec<T> &v) { return v.length(); }

//! Length of vector
template<class T>
int size(const Vec<T> &v) { return v.length(); }

//! Sum of all elements in the vector
template<class T>
T sum(const Vec<T> &v)
{
  T M = 0;

  for (int i = 0;i < v.length();i++)
    M += v[i];

  return M;
}

/*!
 * \brief Sum of elements in the matrix \c m, either along columns or rows
 *
 * <tt>sum(m) = sum(m, 1)</tt> returns a vector where the elements are sum
 * over each column, whereas <tt>sum(m, 2)</tt> returns a vector where the
 * elements are sum over each row.
 */
template<class T>
Vec<T> sum(const Mat<T> &m, int dim = 1)
{
  it_assert((dim == 1) || (dim == 2), "sum: dimension need to be 1 or 2");
  Vec<T> out;

  if (dim == 1) {
    out.set_size(m.cols(), false);

    for (int i = 0; i < m.cols(); i++)
      out(i) = sum(m.get_col(i));
  }
  else {
    out.set_size(m.rows(), false);

    for (int i = 0; i < m.rows(); i++)
      out(i) = sum(m.get_row(i));
  }

  return out;
}


//! Sum of all elements in the given matrix. Fast version of sum(sum(X))
template<class T>
T sumsum(const Mat<T> &X)
{
  const T * X_data = X._data();
  const int X_datasize = X._datasize();
  T acc = 0;

  for (int i = 0;i < X_datasize;i++)
    acc += X_data[i];

  return acc;
}


//! Sum of square of the elements in a vector
template<class T>
T sum_sqr(const Vec<T> &v)
{
  T M = 0;

  for (int i = 0; i < v.length(); i++)
    M += v[i] * v[i];

  return M;
}

/*!
 * \brief Sum of the square of elements in the matrix \c m
 *
 * <tt>sum(m) = sum(m, 1)</tt> returns a vector where the elements are sum
 * squared over each column, whereas <tt>sum(m, 2)</tt> returns a vector
 * where the elements are sum squared over each row
 */
template<class T>
Vec<T> sum_sqr(const Mat<T> &m, int dim = 1)
{
  it_assert((dim == 1) || (dim == 2), "sum_sqr: dimension need to be 1 or 2");
  Vec<T> out;

  if (dim == 1) {
    out.set_size(m.cols(), false);

    for (int i = 0; i < m.cols(); i++)
      out(i) = sum_sqr(m.get_col(i));
  }
  else {
    out.set_size(m.rows(), false);

    for (int i = 0; i < m.rows(); i++)
      out(i) = sum_sqr(m.get_row(i));
  }

  return out;
}

//! Cumulative sum of all elements in the vector
template<class T>
Vec<T> cumsum(const Vec<T> &v)
{
  Vec<T> out(v.size());

  out(0) = v(0);
  for (int i = 1; i < v.size(); i++)
    out(i) = out(i - 1) + v(i);

  return out;
}

/*!
 * \brief Cumulative sum of elements in the matrix \c m
 *
 * <tt>cumsum(m) = cumsum(m, 1)</tt> returns a matrix where the elements
 * are sums over each column, whereas <tt>cumsum(m, 2)</tt> returns a
 * matrix where the elements are sums over each row
 */
template<class T>
Mat<T> cumsum(const Mat<T> &m, int dim = 1)
{
  it_assert((dim == 1) || (dim == 2), "cumsum: dimension need to be 1 or 2");
  Mat<T> out(m.rows(), m.cols());

  if (dim == 1) {
    for (int i = 0; i < m.cols(); i++)
      out.set_col(i, cumsum(m.get_col(i)));
  }
  else {
    for (int i = 0; i < m.rows(); i++)
      out.set_row(i, cumsum(m.get_row(i)));
  }

  return out;
}

//! The product of all elements in the vector
template<class T>
T prod(const Vec<T> &v)
{
  it_assert(v.size() >= 1, "prod: size of vector should be at least 1");
  T out = v(0);

  for (int i = 1; i < v.size(); i++)
    out *= v(i);

  return out;
}

/*!
 * \brief Product of elements in the matrix \c m
 *
 * <tt>prod(m) = prod(m, 1)</tt> returns a vector where the elements are
 * products over each column, whereas <tt>prod(m, 2)</tt> returns a vector
 * where the elements are products over each row
 */
template<class T>
Vec<T> prod(const Mat<T> &m, int dim = 1)
{
  it_assert((dim == 1) || (dim == 2), "prod: dimension need to be 1 or 2");
  Vec<T> out(m.cols());

  if (dim == 1) {
    it_assert((m.cols() >= 1) && (m.rows() >= 1),
              "prod: number of columns should be at least 1");
    out.set_size(m.cols(), false);

    for (int i = 0; i < m.cols(); i++)
      out(i) = prod(m.get_col(i));
  }
  else {
    it_assert((m.cols() >= 1) && (m.rows() >= 1),
              "prod: number of rows should be at least 1");
    out.set_size(m.rows(), false);

    for (int i = 0; i < m.rows(); i++)
      out(i) = prod(m.get_row(i));
  }
  return out;
}

//! Vector cross product. Vectors need to be of size 3
template<class T>
Vec<T> cross(const Vec<T> &v1, const Vec<T> &v2)
{
  it_assert((v1.size() == 3) && (v2.size() == 3),
            "cross: vectors should be of size 3");

  Vec<T> r(3);

  r(0) = v1(1) * v2(2) - v1(2) * v2(1);
  r(1) = v1(2) * v2(0) - v1(0) * v2(2);
  r(2) = v1(0) * v2(1) - v1(1) * v2(0);

  return r;
}


//! Zero-pad a vector to size n
template<class T>
Vec<T> zero_pad(const Vec<T> &v, int n)
{
  it_assert(n >= v.size(), "zero_pad() cannot shrink the vector!");
  Vec<T> v2(n);
  v2.set_subvector(0, v);
  if (n > v.size())
    v2.set_subvector(v.size(), n - 1, T(0));

  return v2;
}

//! Zero-pad a vector to the nearest greater power of two
template<class T>
Vec<T> zero_pad(const Vec<T> &v)
{
  int n = pow2i(levels2bits(v.size()));

  return (n == v.size()) ? v : zero_pad(v, n);
}

//! Zero-pad a matrix to size rows x cols
template<class T>
Mat<T> zero_pad(const Mat<T> &m, int rows, int cols)
{
  it_assert((rows >= m.rows()) && (cols >= m.cols()),
            "zero_pad() cannot shrink the matrix!");
  Mat<T> m2(rows, cols);
  m2.set_submatrix(0, 0, m);
  if (cols > m.cols()) // Zero
    m2.set_submatrix(0, m.rows() - 1, m.cols(), cols - 1, T(0));
  if (rows > m.rows()) // Zero
    m2.set_submatrix(m.rows(), rows - 1, 0, cols - 1, T(0));

  return m2;
}


//! Return zero if indexing outside the vector \c v otherwise return the
//! element \c index
template<class T>
T index_zero_pad(const Vec<T> &v, const int index)
{
  if (index >= 0 && index < v.size())
    return v(index);
  else
    return T(0);
}


//! Transposition of the matrix \c m returning the transposed matrix in \c out
template<class T>
void transpose(const Mat<T> &m, Mat<T> &out) { out = m.T(); }

//! Transposition of the matrix \c m
template<class T>
Mat<T> transpose(const Mat<T> &m) { return m.T(); }


//! Hermitian transpose (complex conjugate transpose) of the matrix \c m
//! returning the transposed matrix in \c out
template<class T>
void hermitian_transpose(const Mat<T> &m, Mat<T> &out) { out = m.H(); }

//! Hermitian transpose (complex conjugate transpose) of the matrix \c m
template<class T>
Mat<T> hermitian_transpose(const Mat<T> &m) { return m.H(); }



/*!
 * \brief Returns true if matrix \c X is hermitian, false otherwise
 * \author M. Szalay
 *
 * A square matrix \f$\mathbf{X}\f$ is hermitian if
 * \f[
 * \mathbf{X} = \mathbf{X}^H
 * \f]
 */
template<class Num_T>
bool is_hermitian(const Mat<Num_T>& X)
{

  if (X == X.H())
    return true;
  else
    return false;
}

/*!
 * \brief Returns true if matrix \c X is unitary, false otherwise
 * \author M. Szalay
 *
 * A square matrix \f$\mathbf{X}\f$ is unitary if
 * \f[
 * \mathbf{X}^H = \mathbf{X}^{-1}
 * \f]
 */
template<class Num_T>
bool is_unitary(const Mat<Num_T>& X)
{

  if (inv(X) == X.H())
    return true;
  else
    return false;
}


/*!
 * \relates Vec
 * \brief Creates a vector with \c n copies of the vector \c v
 * \author Martin Senst
 *
 * \param v Vector to be repeated
 * \param n Number of times to repeat \c v
 */
template<class T>
Vec<T> repmat(const Vec<T> &v, int n)
{
  it_assert(n > 0, "repmat(): Wrong repetition parameter");
  int data_length = v.length();
  it_assert(data_length > 0, "repmat(): Input vector can not be empty");
  Vec<T> assembly(data_length * n);
  for (int j = 0; j < n; ++j) {
    assembly.set_subvector(j * data_length, v);
  }
  return assembly;
}


/*!
 * \relates Mat
 * \brief Creates a matrix with \c m by \c n copies of the matrix \c data
 * \author Mark Dobossy
 *
 * \param data Matrix to be repeated
 * \param m Number of times to repeat data vertically
 * \param n Number of times to repeat data horizontally
 */
template<class T>
Mat<T> repmat(const Mat<T> &data, int m, int n)
{
  it_assert((m > 0) && (n > 0), "repmat(): Wrong repetition parameters");
  int data_rows = data.rows();
  int data_cols = data.cols();
  it_assert((data_rows > 0) && (data_cols > 0), "repmat(): Input matrix can "
            "not be empty");
  Mat<T> assembly(data_rows*m, data_cols*n);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      assembly.set_submatrix(i*data_rows, j*data_cols, data);
    }
  }
  return assembly;
}

/*!
 * \relates Mat
 * \brief Returns a matrix with \c m by \c n copies of the vector \c data
 * \author Adam Piatyszek
 *
 * \param v Vector to be repeated
 * \param m Number of times to repeat data vertically
 * \param n Number of times to repeat data horizontally
 * \param transpose Specifies the input vector orientation (column vector
 * by default)
 */
template<class T> inline
Mat<T> repmat(const Vec<T> &v, int m, int n, bool transpose = false)
{
  return repmat((transpose ? v.T() : Mat<T>(v)), m, n);
}


/*!
 * \brief Computes the Kronecker product of two matrices
 *
 * <tt>K = kron(X, Y)</tt> returns the Kronecker tensor product of \c X
 * and \c Y. The result is a large array formed by taking all possible
 * products between the elements of \c X and those of \c Y. If \c X is
 * <tt>(m x n)</tt> and \c Y is <tt>(p x q)</tt>, then <tt>kron(X, Y)</tt>
 * is <tt>(m*p x n*q)</tt>.
 *
 * \author Adam Piatyszek
 */
template<class Num_T>
Mat<Num_T> kron(const Mat<Num_T>& X, const Mat<Num_T>& Y)
{
  Mat<Num_T> result(X.rows() * Y.rows(), X.cols() * Y.cols());

  for (int i = 0; i < X.rows(); i++)
    for (int j = 0; j < X.cols(); j++)
      result.set_submatrix(i * Y.rows(), j * Y.cols(), X(i, j) * Y);

  return result;
}


/*!
 * \brief Square root of the complex square matrix \c A
 *
 * This function computes the matrix square root of the complex square
 * matrix \c A. The implementation is based on the Matlab/Octave \c
 * sqrtm() function.
 *
 * Ref: N. J. Higham, "Numerical Analysis Report No. 336", Manchester
 * Centre for Computational Mathematics, Manchester, England, January 1999
 *
 * \author Adam Piatyszek
 */
ITPP_EXPORT cmat sqrtm(const cmat& A);

/*!
 * \brief Square root of the real square matrix \c A
 *
 * This function computes the matrix square root of the real square matrix
 * \c A. Please note that the returned matrix is complex. The
 * implementation is based on the Matlab/Octave \c sqrtm() function.
 *
 * Ref: N. J. Higham, "Numerical Analysis Report No. 336", Manchester
 * Centre for Computational Mathematics, Manchester, England, January 1999
 *
 * \author Adam Piatyszek
 */
ITPP_EXPORT cmat sqrtm(const mat& A);


/*!
 * \brief Calculate the rank of matrix \c m
 * \author Martin Senst
 *
 * \param m Input matrix
 * \param tol Tolerance used for comparing the singular values with zero.
 *            If negative, it is automatically determined.
 */
template<class T>
int rank(const Mat<T> &m, double tol = -1.0)
{
  int rows = m.rows();
  int cols = m.cols();
  if ((rows == 0) || (cols == 0))
    return 0;

  vec sing_val = svd(m);

  if (tol < 0.0) { // Calculate default tolerance
    tol = eps * sing_val(0) * (rows > cols ? rows : cols);
  }

  // Count number of nonzero singular values
  int r = 0;
  while ((r < sing_val.length()) && (sing_val(r) > tol)) {
    r++;
  }

  return r;
}

//! Specialisation of rank() function
template<> inline
int rank(const imat &m, double tol)
{
  return rank(to_mat(m), tol);
}

//! Specialisation of rank() function
template<> inline
int rank(const smat &m, double tol)
{
  return rank(to_mat(m), tol);
}

//! Specialisation of rank() function
template<> inline
int rank(const bmat &, double)
{
  it_error("rank(bmat): Function not implemented for GF(2) algebra");
  return 0;
}

//!@}



// -------------------- Diagonal matrix functions -------------------------

//! \addtogroup diag
//!@{

/*!
 * \brief Create a diagonal matrix using vector \c v as its diagonal
 *
 * All other matrix elements except the ones on its diagonal are set to
 * zero. An optional parameter \c K can be used to shift the diagonal in
 * the resulting matrix. By default \c K is equal to zero.
 *
 * The size of the diagonal matrix will be \f$n+|K| \times n+|K|\f$, where
 * \f$n\f$ is the length of the input vector \c v.
 */
template<class T>
Mat<T> diag(const Vec<T> &v, const int K = 0)
{
  Mat<T> m(v.size() + std::abs(K), v.size() + std::abs(K));
  m = T(0);
  if (K > 0)
    for (int i = v.size() - 1; i >= 0; i--)
      m(i, i + K) = v(i);
  else
    for (int i = v.size() - 1; i >= 0; i--)
      m(i - K, i) = v(i);

  return m;
}

/*!
 * \brief Create a diagonal matrix using vector \c v as its diagonal
 *
 * All other matrix elements except the ones on its diagonal are set to
 * zero.
 *
 * The size of the diagonal matrix will be \f$n \times n\f$, where \f$n\f$
 * is the length of the input vector \c v.
 */
template<class T>
void diag(const Vec<T> &v, Mat<T> &m)
{
  m.set_size(v.size(), v.size(), false);
  m = T(0);
  for (int i = v.size() - 1; i >= 0; i--)
    m(i, i) = v(i);
}

/*!
 * \brief Get the diagonal elements of the input matrix \c m
 *
 * The size of the output vector with diagonal elements will be
 * \f$n = min(r, c)\f$, where \f$r \times c\f$ are the dimensions of
 * matrix \c m.
*/
template<class T>
Vec<T> diag(const Mat<T> &m)
{
  Vec<T> t(std::min(m.rows(), m.cols()));

  for (int i = 0; i < t.size(); i++)
    t(i) = m(i, i);

  return t;
}

/*!
  \brief Returns a matrix with the elements of the input vector \c main on
  the diagonal and the elements of the input vector \c sup on the diagonal
  row above.

  If the number of elements in the vector \c main is \f$n\f$, then the
  number of elements in the input vector \c sup must be \f$n-1\f$. The
  size of the return matrix will be \f$n \times n\f$.
*/
template<class T>
Mat<T> bidiag(const Vec<T> &main, const Vec<T> &sup)
{
  it_assert(main.size() == sup.size() + 1, "bidiag()");

  int n = main.size();
  Mat<T> m(n, n);
  m = T(0);
  for (int i = 0; i < n - 1; i++) {
    m(i, i) = main(i);
    m(i, i + 1) = sup(i);
  }
  m(n - 1, n - 1) = main(n - 1);

  return m;
}

/*!
  \brief Returns in the output variable \c m a matrix with the elements of
  the input vector \c main on the diagonal and the elements of the input
  vector \c sup on the diagonal row above.

  If the number of elements in the vector \c main is \f$n\f$, then the
  number of elements in the input vector \c sup must be \f$n-1\f$. The
  size of the output matrix \c m will be \f$n \times n\f$.
*/
template<class T>
void bidiag(const Vec<T> &main, const Vec<T> &sup, Mat<T> &m)
{
  it_assert(main.size() == sup.size() + 1, "bidiag()");

  int n = main.size();
  m.set_size(n, n);
  m = T(0);
  for (int i = 0; i < n - 1; i++) {
    m(i, i) = main(i);
    m(i, i + 1) = sup(i);
  }
  m(n - 1, n - 1) = main(n - 1);
}

/*!
  \brief Returns the main diagonal and the diagonal row above in the two
  output vectors \c main and \c sup.

  The input matrix \c in must be a square \f$n \times n\f$ matrix. The
  length of the output vector \c main will be \f$n\f$ and the length of
  the output vector \c sup will be \f$n-1\f$.
*/
template<class T>
void bidiag(const Mat<T> &m, Vec<T> &main, Vec<T> &sup)
{
  it_assert(m.rows() == m.cols(), "bidiag(): Matrix must be square!");

  int n = m.cols();
  main.set_size(n);
  sup.set_size(n - 1);
  for (int i = 0; i < n - 1; i++) {
    main(i) = m(i, i);
    sup(i) = m(i, i + 1);
  }
  main(n - 1) = m(n - 1, n - 1);
}

/*!
  \brief Returns a matrix with the elements of \c main on the diagonal,
  the elements of \c sup on the diagonal row above, and the elements of \c
  sub on the diagonal row below.

  If the length of the input vector \c main is \f$n\f$ then the lengths of
  the vectors \c sup and \c sub must equal \f$n-1\f$. The size of the
  return matrix will be \f$n \times n\f$.
*/
template<class T>
Mat<T> tridiag(const Vec<T> &main, const Vec<T> &sup, const Vec<T> &sub)
{
  it_assert(main.size() == sup.size() + 1 && main.size() == sub.size() + 1, "bidiag()");

  int n = main.size();
  Mat<T> m(n, n);
  m = T(0);
  for (int i = 0; i < n - 1; i++) {
    m(i, i) = main(i);
    m(i, i + 1) = sup(i);
    m(i + 1, i) = sub(i);
  }
  m(n - 1, n - 1) = main(n - 1);

  return m;
}

/*!
  \brief Returns in the output matrix \c m a matrix with the elements of
  \c main on the diagonal, the elements of \c sup on the diagonal row
  above, and the elements of \c sub on the diagonal row below.

  If the length of the input vector \c main is \f$n\f$ then the lengths of
  the vectors \c sup and \c sub must equal \f$n-1\f$. The size of the
  output matrix \c m will be \f$n \times n\f$.
*/
template<class T>
void tridiag(const Vec<T> &main, const Vec<T> &sup, const Vec<T> &sub, Mat<T> &m)
{
  it_assert(main.size() == sup.size() + 1 && main.size() == sub.size() + 1, "bidiag()");

  int n = main.size();
  m.set_size(n, n);
  m = T(0);
  for (int i = 0; i < n - 1; i++) {
    m(i, i) = main(i);
    m(i, i + 1) = sup(i);
    m(i + 1, i) = sub(i);
  }
  m(n - 1, n - 1) = main(n - 1);
}

/*!
  \brief Returns the main diagonal, the diagonal row above, and the
  diagonal row below int the output vectors \c main, \c sup, and \c sub.

  The input matrix \c m must be a square \f$n \times n\f$ matrix. The
  length of the output vector \c main will be \f$n\f$ and the length of
  the output vectors \c sup and \c sup will be \f$n-1\f$.
*/
template<class T>
void tridiag(const Mat<T> &m, Vec<T> &main, Vec<T> &sup, Vec<T> &sub)
{
  it_assert(m.rows() == m.cols(), "tridiag(): Matrix must be square!");

  int n = m.cols();
  main.set_size(n);
  sup.set_size(n - 1);
  sub.set_size(n - 1);
  for (int i = 0; i < n - 1; i++) {
    main(i) = m(i, i);
    sup(i) = m(i, i + 1);
    sub(i) = m(i + 1, i);
  }
  main(n - 1) = m(n - 1, n - 1);
}


/*!
  \brief The trace of the matrix \c m, i.e. the sum of the diagonal elements.
*/
template<class T>
T trace(const Mat<T> &m)
{
  return sum(diag(m));
}

//!@}


// ----------------- reshaping vectors and matrices ------------------------

//! \addtogroup reshaping
//!@{

//! Reverse the input vector
template<class T>
Vec<T> reverse(const Vec<T> &in)
{
  int i, s = in.length();

  Vec<T> out(s);
  for (i = 0;i < s;i++)
    out[i] = in[s-1-i];
  return out;
}

//! Row vectorize the matrix [(0,0) (0,1) ... (N-1,N-2) (N-1,N-1)]
template<class T>
Vec<T> rvectorize(const Mat<T> &m)
{
  int i, j, n = 0, r = m.rows(), c = m.cols();
  Vec<T> v(r * c);

  for (i = 0; i < r; i++)
    for (j = 0; j < c; j++)
      v(n++) = m(i, j);

  return v;
}

//! Column vectorize the matrix [(0,0) (1,0) ... (N-2,N-1) (N-1,N-1)]
template<class T>
Vec<T> cvectorize(const Mat<T> &m)
{
  int i, j, n = 0, r = m.rows(), c = m.cols();
  Vec<T> v(r * c);

  for (j = 0; j < c; j++)
    for (i = 0; i < r; i++)
      v(n++) = m(i, j);

  return v;
}

/*!
  \brief Reshape the matrix into an rows*cols matrix

  The data is taken columnwise from the original matrix and written
  columnwise into the new matrix.
*/
template<class T>
Mat<T> reshape(const Mat<T> &m, int rows, int cols)
{
  it_assert_debug(m.rows()*m.cols() == rows*cols, "Mat<T>::reshape: Sizes must match");
  Mat<T> temp(rows, cols);
  int i, j, ii = 0, jj = 0;
  for (j = 0; j < m.cols(); j++) {
    for (i = 0; i < m.rows(); i++) {
      temp(ii++, jj) = m(i, j);
      if (ii == rows) {
        jj++;
        ii = 0;
      }
    }
  }
  return temp;
}

/*!
  \brief Reshape the vector into an rows*cols matrix

  The data is element by element from the vector and written columnwise
  into the new matrix.
*/
template<class T>
Mat<T> reshape(const Vec<T> &v, int rows, int cols)
{
  it_assert_debug(v.size() == rows*cols, "Mat<T>::reshape: Sizes must match");
  Mat<T> temp(rows, cols);
  int i, j, ii = 0;
  for (j = 0; j < cols; j++) {
    for (i = 0; i < rows; i++) {
      temp(i, j) = v(ii++);
    }
  }
  return temp;
}

//!@}


//! Returns \a true if all elements are ones and \a false otherwise
ITPP_EXPORT bool all(const bvec &testvec);
//! Returns \a true if any element is one and \a false otherwise
ITPP_EXPORT bool any(const bvec &testvec);

//! \cond

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT int length(const vec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT int length(const cvec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT int length(const svec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT int length(const ivec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT int length(const bvec &v);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT double sum(const vec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT std::complex<double> sum(const cvec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT short sum(const svec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT int sum(const ivec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bin sum(const bvec &v);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT double sum_sqr(const vec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT std::complex<double> sum_sqr(const cvec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT short sum_sqr(const svec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT int sum_sqr(const ivec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bin sum_sqr(const bvec &v);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec cumsum(const vec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec cumsum(const cvec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT svec cumsum(const svec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec cumsum(const ivec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bvec cumsum(const bvec &v);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT double prod(const vec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT std::complex<double> prod(const cvec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT short prod(const svec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT int prod(const ivec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bin prod(const bvec &v);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec cross(const vec &v1, const vec &v2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec cross(const cvec &v1, const cvec &v2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec cross(const ivec &v1, const ivec &v2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT svec cross(const svec &v1, const svec &v2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bvec cross(const bvec &v1, const bvec &v2);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec reverse(const vec &in);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec reverse(const cvec &in);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT svec reverse(const svec &in);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec reverse(const ivec &in);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bvec reverse(const bvec &in);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec zero_pad(const vec &v, int n);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec zero_pad(const cvec &v, int n);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec zero_pad(const ivec &v, int n);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT svec zero_pad(const svec &v, int n);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bvec zero_pad(const bvec &v, int n);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec zero_pad(const vec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec zero_pad(const cvec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec zero_pad(const ivec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT svec zero_pad(const svec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bvec zero_pad(const bvec &v);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat zero_pad(const mat &, int, int);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat zero_pad(const cmat &, int, int);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT imat zero_pad(const imat &, int, int);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT smat zero_pad(const smat &, int, int);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bmat zero_pad(const bmat &, int, int);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec sum(const mat &m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec sum(const cmat &m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT svec sum(const smat &m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec sum(const imat &m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bvec sum(const bmat &m, int dim);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT double sumsum(const mat &X);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT std::complex<double> sumsum(const cmat &X);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT short sumsum(const smat &X);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT int sumsum(const imat &X);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bin sumsum(const bmat &X);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec sum_sqr(const mat & m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec sum_sqr(const cmat &m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT svec sum_sqr(const smat &m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec sum_sqr(const imat &m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bvec sum_sqr(const bmat &m, int dim);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat cumsum(const mat &m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat cumsum(const cmat &m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT smat cumsum(const smat &m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT imat cumsum(const imat &m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bmat cumsum(const bmat &m, int dim);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec prod(const mat &m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec prod(const cmat &v, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT svec prod(const smat &m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec prod(const imat &m, int dim);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bvec prod(const bmat &m, int dim);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec diag(const mat &in);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec diag(const cmat &in);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void diag(const vec &in, mat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void diag(const cvec &in, cmat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat diag(const vec &v, const int K);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat diag(const cvec &v, const int K);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat bidiag(const vec &, const vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat bidiag(const cvec &, const cvec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void bidiag(const vec &, const vec &, mat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void bidiag(const cvec &, const cvec &, cmat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void bidiag(const mat &, vec &, vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void bidiag(const cmat &, cvec &, cvec &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat tridiag(const vec &main, const vec &, const vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat tridiag(const cvec &main, const cvec &, const cvec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void tridiag(const vec &main, const vec &, const vec &, mat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void tridiag(const cvec &main, const cvec &, const cvec &, cmat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void tridiag(const mat &m, vec &, vec &, vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void tridiag(const cmat &m, cvec &, cvec &, cvec &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT double trace(const mat &in);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT std::complex<double> trace(const cmat &in);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT short trace(const smat &in);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT int trace(const imat &in);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bin trace(const bmat &in);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void transpose(const mat &m, mat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void transpose(const cmat &m, cmat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void transpose(const smat &m, smat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void transpose(const imat &m, imat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void transpose(const bmat &m, bmat &out);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat transpose(const mat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat transpose(const cmat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT smat transpose(const smat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT imat transpose(const imat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bmat transpose(const bmat &m);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void hermitian_transpose(const mat &m, mat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void hermitian_transpose(const cmat &m, cmat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void hermitian_transpose(const smat &m, smat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void hermitian_transpose(const imat &m, imat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT void hermitian_transpose(const bmat &m, bmat &out);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat hermitian_transpose(const mat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat hermitian_transpose(const cmat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT smat hermitian_transpose(const smat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT imat hermitian_transpose(const imat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bmat hermitian_transpose(const bmat &m);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bool is_hermitian(const mat &X);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bool is_hermitian(const cmat &X);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bool is_unitary(const mat &X);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bool is_unitary(const cmat &X);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec rvectorize(const mat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec rvectorize(const cmat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec rvectorize(const imat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT svec rvectorize(const smat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bvec rvectorize(const bmat &m);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec cvectorize(const mat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec cvectorize(const cmat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec cvectorize(const imat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT svec cvectorize(const smat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bvec cvectorize(const bmat &m);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat reshape(const mat &m, int rows, int cols);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat reshape(const cmat &m, int rows, int cols);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT imat reshape(const imat &m, int rows, int cols);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT smat reshape(const smat &m, int rows, int cols);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bmat reshape(const bmat &m, int rows, int cols);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat reshape(const vec &m, int rows, int cols);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat reshape(const cvec &m, int rows, int cols);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT imat reshape(const ivec &m, int rows, int cols);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT smat reshape(const svec &m, int rows, int cols);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bmat reshape(const bvec &m, int rows, int cols);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat kron(const mat &X, const mat &Y);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat kron(const cmat &X, const cmat &Y);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT imat kron(const imat &X, const imat &Y);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT smat kron(const smat &X, const smat &Y);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bmat kron(const bmat &X, const bmat &Y);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec repmat(const vec &v, int n);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec repmat(const cvec &v, int n);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec repmat(const ivec &v, int n);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT svec repmat(const svec &v, int n);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bvec repmat(const bvec &v, int n);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat repmat(const vec &v, int m, int n, bool transpose);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat repmat(const cvec &v, int m, int n, bool transpose);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT imat repmat(const ivec &v, int m, int n, bool transpose);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT smat repmat(const svec &v, int m, int n, bool transpose);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bmat repmat(const bvec &v, int m, int n, bool transpose);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat repmat(const mat &data, int m, int n);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat repmat(const cmat &data, int m, int n);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT imat repmat(const imat &data, int m, int n);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT smat repmat(const smat &data, int m, int n);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT bmat repmat(const bmat &data, int m, int n);

//! \endcond

} // namespace itpp

#endif // #ifndef MATFUNC_H
