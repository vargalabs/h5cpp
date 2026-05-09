/*!
 * \file
 * \brief Sparse Matrix Class Definitions
 * \author Tony Ottosson and Tobias Ringstrom
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

#ifndef SMAT_H
#define SMAT_H

#include <itpp/base/svec.h>
#include <itpp/itexports.h>

namespace itpp
{

// Declaration of class Vec
template <class T> class Vec;
// Declaration of class Mat
template <class T> class Mat;
// Declaration of class Sparse_Vec
template <class T> class Sparse_Vec;
// Declaration of class Sparse_Mat
template <class T> class Sparse_Mat;

// ------------------------ Sparse_Mat Friends -------------------------------------

//! m1+m2 where m1 and m2 are sparse matrices
template <class T>
Sparse_Mat<T> operator+(const Sparse_Mat<T> &m1, const Sparse_Mat<T> &m2);

//! c*m where c is a scalar and m is a sparse matrix
template <class T>
Sparse_Mat<T> operator*(const T &c, const Sparse_Mat<T> &m);

//! m1*m2 where m1 and m2 are sparse matrices
template <class T>
Sparse_Mat<T> operator*(const Sparse_Mat<T> &m1, const Sparse_Mat<T> &m2);

//! m*v where m is a sparse matrix and v is a sparse vector
template <class T>
Sparse_Vec<T> operator*(const Sparse_Mat<T> &m, const Sparse_Vec<T> &v);

//! m*v where m is a sparse matrix and v is a full column vector
template <class T>
Vec<T> operator*(const Sparse_Mat<T> &m, const Vec<T> &v);

//! v'*m where m is a sparse matrix and v is a full column vector
template <class T>
Vec<T> operator*(const Vec<T> &v, const Sparse_Mat<T> &m);

//! m'*m where m is a sparse matrix
template <class T>
Mat<T> trans_mult(const Sparse_Mat<T> &m);

//! m'*m where m is a sparse matrix
template <class T>
Sparse_Mat<T> trans_mult_s(const Sparse_Mat<T> &m);

//! m1'*m2 where m1 and m2 are sparse matrices
template <class T>
Sparse_Mat<T> trans_mult(const Sparse_Mat<T> &m1, const Sparse_Mat<T> &m2);

//! m'*v where m is a sparse matrix and v is a full column vector
template <class T>
Vec<T> trans_mult(const Sparse_Mat<T> &m, const Vec<T> &v);

//! m1*m2' where m1 and m2 are sparse matrices
template <class T>
Sparse_Mat<T> mult_trans(const Sparse_Mat<T> &m1, const Sparse_Mat<T> &m2);

/*!
  \brief Templated Sparse Matrix Class
  \author Tony Ottosson and Tobias Ringstrom

  A sparse matrix is a matrix where most elements are zero. The
  maximum number of non-zero elements in each column is a parameter
  to the constructor.

  The implementation is based on representing all columns as sparse
  vectors. Thus, column access generally is much faster than row
  access. The elements in each vector are stored in random order,
  i.e. they are not sorted.
*/
template <class T>
class Sparse_Mat
{
public:

  //! Default constructor
  Sparse_Mat();

  /*!
    \brief Initiate an empty sparse matrix

    A Sparse_Mat consists of colums that have the type Sparse_Vec. The maximum number of non-zero elements is each column
    is denoted \c row_data_init.

    \param rows Number of rows in the matrix
    \param cols Number of columns in the matrix
    \param row_data_init The maximum number of non-zero elements in each column (default value is 200)
  */
  Sparse_Mat(int rows, int cols, int row_data_init = 200);

  //! Initiate a new sparse matrix. The elements of \c m are copied into the new sparse matrix
  Sparse_Mat(const Sparse_Mat<T> &m);

  //! Initiate a new sparse matrix from a dense matrix. The elements of \c m are copied into the new sparse matrix
  Sparse_Mat(const Mat<T> &m);

  /*!
    \brief Initiate a new sparse matrix from a dense matrix. Elements of \c m larger than \c epsilon are copied into the new sparse matrix.

    \note If the type T is double complex, then the elements of \c m larger than \c abs(epsilon) are copied into the new sparse matrix.
  */
  Sparse_Mat(const Mat<T> &m, T epsilon);

  //! Destructor
  ~Sparse_Mat();

  /*!
    \brief Set the size of the sparse matrix

    A Sparse_Mat consists of colums that have the type Sparse_Vec. The maximum number of non-zero elements is each column
    is denoted \c row_data_init, with default value =-1 indicating that the number of data elements is not changed.

    \param rows Number of rows in the matrix
    \param cols Number of columns in the matrix
    \param row_data_init The maximum number of non-zero elements in each column (default value -1 \c => allocated size for the data is not changed)
  */
  void set_size(int rows, int cols, int row_data_init = -1);

  //! Returns the number of rows of the sparse matrix
  int rows() const { return n_rows; }

  //! Returns the number of columns of the sparse matrix
  int cols() const { return n_cols; }

  //! The number of non-zero elements in the sparse matrix
  int nnz();

  //! Returns the density of the sparse matrix: (number of non-zero elements)/(total number of elements)
  double density();

  //! Set the maximum number of non-zero elements in each column equal to the actual number of non-zero elements in each column
  void compact();

  //! Returns a full, dense matrix in \c m
  void full(Mat<T> &m) const;

  //! Returns a full, dense matrix
  Mat<T> full() const;

  //! Returns element of row \c r and column \c c
  T operator()(int r, int c) const;

  //! Set element (\c r, \c c ) equal to \c v
  void set(int r, int c, T v);

  //! Set a new element with index (\c r, \c c ) equal to \c v
  void set_new(int r, int c, T v);

  //! Add the element in row \c r and column \c c with \c v
  void add_elem(const int r, const int c, const T v);

  //! Set the sparse matrix to the all zero matrix (removes all non-zero elements)
  void zeros();

  //! Set the element in row \c r and column \c c to zero (i.e. clear that element if it contains a non-zero value)
  void zero_elem(const int r, const int c);

  //! Clear all non-zero elements of the sparse matrix
  void clear();

  //! Clear the element in row \c r and column \c c (if it contains a non-zero value)
  void clear_elem(const int r, const int c);

  //! Set submatrix defined by rows r1,r2 and columns c1,c2 to matrix m
  void set_submatrix(int r1, int r2, int c1, int c2, const Mat<T> &m);

  //! Set submatrix defined by upper-left element (\c r,\c c) and the size of matrix \c m to \c m
  void set_submatrix(int r, int c, const Mat<T>& m);

  //! Returns the sub-matrix from rows \c r1 to \c r2 and columns \c c1 to \c c2
  Sparse_Mat<T> get_submatrix(int r1, int r2, int c1, int c2) const;

  //! Returns the sub-matrix from columns \c c1 to \c c2 (all rows)
  Sparse_Mat<T> get_submatrix_cols(int c1, int c2) const;

  //! Returns column \c c of the Sparse_Mat in the Sparse_Vec \c v
  void get_col(int c, Sparse_Vec<T> &v) const;

  //! Returns column \c c of the Sparse_Mat
  Sparse_Vec<T> get_col(int c) const;

  //! Set column \c c of the Sparse_Mat
  void set_col(int c, const Sparse_Vec<T> &v);

  /*! Transpose the sparse matrix, return the result in \c m

  Note: this function can be slow for large matrices.
   */
  void transpose(Sparse_Mat<T> &m) const;

  /*! Returns the transpose of the sparse matrix

  Note: this function can be slow for large matrices.
  */
  Sparse_Mat<T> transpose() const;

  /*! Returns the transpose of the sparse matrix

  Note: this function can be slow for large matrices.
  */
  // Sparse_Mat<T> T() const { return this->transpose(); };

  //! Assign sparse matrix the value and dimensions of the sparse matrix \c m
  void operator=(const Sparse_Mat<T> &m);

  //! Assign sparse matrix the value and dimensions of the dense matrix \c m
  void operator=(const Mat<T> &m);

  //! Returns the sign inverse of all elements in the sparse matrix
  Sparse_Mat<T> operator-() const;

  //! Compare two sparse matricies. False if wrong sizes or different values
  bool operator==(const Sparse_Mat<T> &m) const;

  //! Add sparse matrix \c v to all non-zero elements of the sparse matrix
  void operator+=(const Sparse_Mat<T> &v);

  //! Add matrix \c v to all non-zero elements of the sparse matrix
  void operator+=(const Mat<T> &v);

  //! Subtract sparse matrix \c v from all non-zero elements of the sparse matrix
  void operator-=(const Sparse_Mat<T> &v);

  //! Subtract matrix \c v from all non-zero elements of the sparse matrix
  void operator-=(const Mat<T> &v);

  //! Multiply all non-zero elements of the sparse matrix with the scalar \c v
  void operator*=(const T &v);

  //! Divide all non-zero elements of the sparse matrix with the scalar \c v
  void operator/=(const T &v);

  //! Addition m1+m2 where m1 and m2 are sparse matrices
  friend Sparse_Mat<T> operator+<>(const Sparse_Mat<T> &m1, const Sparse_Mat<T> &m2);

  //! Multiplication c*m where c is a scalar and m is a sparse matrix
  friend Sparse_Mat<T> operator*<>(const T &c, const Sparse_Mat<T> &m);

  //! Multiplication m1*m2 where m1 and m2 are sparse matrices
  friend Sparse_Mat<T> operator*<>(const Sparse_Mat<T> &m1, const Sparse_Mat<T> &m2);

  //! Multiplication m*v where m is a sparse matrix and v is a sparse vector
  friend Sparse_Vec<T> operator*<>(const Sparse_Mat<T> &m, const Sparse_Vec<T> &v);

  //! Multiplication m*v where m is a sparse matrix and v is a full column vector
  friend Vec<T> operator*<>(const Sparse_Mat<T> &m, const Vec<T> &v);

  //! Multiplication v'*m where m is a sparse matrix and v is a full column vector
  friend Vec<T> operator*<>(const Vec<T> &v, const Sparse_Mat<T> &m);

  //! Multiplication m'*m where m is a sparse matrix. Returns a full, dense matrix
  friend Mat<T> trans_mult <>(const Sparse_Mat<T> &m);

  //! Multiplication m'*m where m is a sparse matrix, Returns a sparse matrix
  friend Sparse_Mat<T> trans_mult_s <>(const Sparse_Mat<T> &m);

  //! Multiplication m1'*m2 where m1 and m2 are sparse matrices
  friend Sparse_Mat<T> trans_mult <>(const Sparse_Mat<T> &m1, const Sparse_Mat<T> &m2);

  //! Multiplication m'*v where m is a sparse matrix and v is a full column vector
  friend Vec<T> trans_mult <>(const Sparse_Mat<T> &m, const Vec<T> &v);

  //! Multiplication m1*m2' where m1 and m2 are sparse matrices
  friend Sparse_Mat<T> mult_trans <>(const Sparse_Mat<T> &m1, const Sparse_Mat<T> &m2);

private:
  void init();
  void alloc_empty();
  void alloc(int row_data_size = 200);
  void free();

  int n_rows, n_cols;
  Sparse_Vec<T> *col;
};

/*!
  \relates Sparse_Mat
  \brief Sparse integer matrix
*/
typedef Sparse_Mat<int> sparse_imat;

/*!
  \relates Sparse_Mat
  \brief Sparse double matrix
*/
typedef Sparse_Mat<double> sparse_mat;

/*!
  \relates Sparse_Mat
  \brief Sparse complex<double> matrix
*/
typedef Sparse_Mat<std::complex<double> > sparse_cmat;

//---------------------- Implementation starts here --------------------------------

template <class T>
void Sparse_Mat<T>::init()
{
  n_rows = 0;
  n_cols = 0;
  col = 0;
}

template <class T>
void Sparse_Mat<T>::alloc_empty()
{
  if (n_cols == 0)
    col = 0;
  else
    col = new Sparse_Vec<T>[n_cols];
}

template <class T>
void Sparse_Mat<T>::alloc(int row_data_init)
{
  if (n_cols == 0)
    col = 0;
  else
    col = new Sparse_Vec<T>[n_cols];
  for (int c = 0; c < n_cols; c++)
    col[c].set_size(n_rows, row_data_init);
}

template <class T>
void Sparse_Mat<T>::free()
{
  delete [] col;
  col = 0;
}

template <class T>
Sparse_Mat<T>::Sparse_Mat()
{
  init();
}

template <class T>
Sparse_Mat<T>::Sparse_Mat(int rows, int cols, int row_data_init)
{
  init();
  n_rows = rows;
  n_cols = cols;
  alloc(row_data_init);
}

template <class T>
Sparse_Mat<T>::Sparse_Mat(const Sparse_Mat<T> &m)
{
  init();
  n_rows = m.n_rows;
  n_cols = m.n_cols;
  alloc_empty();

  for (int c = 0; c < n_cols; c++)
    col[c] = m.col[c];
}

template <class T>
Sparse_Mat<T>::Sparse_Mat(const Mat<T> &m)
{
  init();
  n_rows = m.rows();
  n_cols = m.cols();
  alloc();

  for (int c = 0; c < n_cols; c++) {
    for (int r = 0; r < n_rows; r++) {
      //if (abs(m(r,c)) != T(0))
      if (m(r, c) != T(0))
        col[c].set_new(r, m(r, c));
    }
    col[c].compact();
  }
}

template <class T>
Sparse_Mat<T>::Sparse_Mat(const Mat<T> &m, T epsilon)
{
  init();
  n_rows = m.rows();
  n_cols = m.cols();
  alloc();

  for (int c = 0; c < n_cols; c++) {
    for (int r = 0; r < n_rows; r++) {
      if (std::abs(m(r, c)) > std::abs(epsilon))
        col[c].set_new(r, m(r, c));
    }
    col[c].compact();
  }
}

template <class T>
Sparse_Mat<T>::~Sparse_Mat()
{
  free();
}

template <class T>
void Sparse_Mat<T>::set_size(int rows, int cols, int row_data_init)
{
  n_rows = rows;

  //Allocate new memory for data if the number of columns has changed or if row_data_init != -1
  if (cols != n_cols || row_data_init != -1) {
    n_cols = cols;
    free();
    alloc(row_data_init);
  }
}

template <class T>
int Sparse_Mat<T>::nnz()
{
  int n = 0;
  for (int c = 0; c < n_cols; c++)
    n += col[c].nnz();

  return n;
}

template <class T>
double Sparse_Mat<T>::density()
{
  //return static_cast<double>(nnz())/(n_rows*n_cols);
  return double(nnz()) / (n_rows*n_cols);
}

template <class T>
void Sparse_Mat<T>::compact()
{
  for (int c = 0; c < n_cols; c++)
    col[c].compact();
}

template <class T>
void Sparse_Mat<T>::full(Mat<T> &m) const
{
  m.set_size(n_rows, n_cols);
  m = T(0);
  for (int c = 0; c < n_cols; c++) {
    for (int p = 0; p < col[c].nnz(); p++)
      m(col[c].get_nz_index(p), c) = col[c].get_nz_data(p);
  }
}

template <class T>
Mat<T> Sparse_Mat<T>::full() const
{
  Mat<T> r(n_rows, n_cols);
  full(r);
  return r;
}

template <class T>
T Sparse_Mat<T>::operator()(int r, int c) const
{
  it_assert_debug(r >= 0 && r<n_rows && c >= 0 && c < n_cols, "Incorrect input indexes given");
  return col[c](r);
}

template <class T>
void Sparse_Mat<T>::set(int r, int c, T v)
{
  it_assert_debug(r >= 0 && r<n_rows && c >= 0 && c < n_cols, "Incorrect input indexes given");
  col[c].set(r, v);
}

template <class T>
void Sparse_Mat<T>::set_new(int r, int c, T v)
{
  it_assert_debug(r >= 0 && r<n_rows && c >= 0 && c < n_cols, "Incorrect input indexes given");
  col[c].set_new(r, v);
}

template <class T>
void Sparse_Mat<T>::add_elem(int r, int c, T v)
{
  it_assert_debug(r >= 0 && r<n_rows && c >= 0 && c < n_cols, "Incorrect input indexes given");
  col[c].add_elem(r, v);
}

template <class T>
void Sparse_Mat<T>::zeros()
{
  for (int c = 0; c < n_cols; c++)
    col[c].zeros();
}

template <class T>
void Sparse_Mat<T>::zero_elem(const int r, const int c)
{
  it_assert_debug(r >= 0 && r<n_rows && c >= 0 && c < n_cols, "Incorrect input indexes given");
  col[c].zero_elem(r);
}

template <class T>
void Sparse_Mat<T>::clear()
{
  for (int c = 0; c < n_cols; c++)
    col[c].clear();
}

template <class T>
void Sparse_Mat<T>::clear_elem(const int r, const int c)
{
  it_assert_debug(r >= 0 && r<n_rows && c >= 0 && c < n_cols, "Incorrect input indexes given");
  col[c].clear_elem(r);
}

template <class T>
void Sparse_Mat<T>::set_submatrix(int r1, int r2, int c1, int c2, const Mat<T>& m)
{
  if (r1 == -1) r1 = n_rows - 1;
  if (r2 == -1) r2 = n_rows - 1;
  if (c1 == -1) c1 = n_cols - 1;
  if (c2 == -1) c2 = n_cols - 1;

  it_assert_debug(r1 >= 0 && r2 >= 0 && r1 < n_rows && r2 < n_rows &&
                  c1 >= 0 && c2 >= 0 && c1 < n_cols && c2 < n_cols, "Sparse_Mat<Num_T>::set_submatrix(): index out of range");

  it_assert_debug(r2 >= r1 && c2 >= c1, "Sparse_Mat<Num_T>::set_submatrix: r2<r1 or c2<c1");
  it_assert_debug(m.rows() == r2 - r1 + 1 && m.cols() == c2 - c1 + 1, "Mat<Num_T>::set_submatrix(): sizes don't match");

  for (int i = 0 ; i < m.rows() ; i++) {
    for (int j = 0 ; j < m.cols() ; j++) {
      set(r1 + i, c1 + j, m(i, j));
    }
  }
}

template <class T>
void Sparse_Mat<T>::set_submatrix(int r, int c, const Mat<T>& m)
{
  it_assert_debug(r >= 0 && r + m.rows() <= n_rows &&
                  c >= 0 && c + m.cols() <= n_cols, "Sparse_Mat<Num_T>::set_submatrix(): index out of range");

  for (int i = 0 ; i < m.rows() ; i++) {
    for (int j = 0 ; j < m.cols() ; j++) {
      set(r + i, c + j, m(i, j));
    }
  }
}

template <class T>
Sparse_Mat<T> Sparse_Mat<T>::get_submatrix(int r1, int r2, int c1, int c2) const
{
  it_assert_debug(r1 <= r2 && r1 >= 0 && r1 < n_rows && c1 <= c2 && c1 >= 0 && c1 < n_cols,
                  "Sparse_Mat<T>::get_submatrix(): illegal input variables");

  Sparse_Mat<T> r(r2 - r1 + 1, c2 - c1 + 1);

  for (int c = c1; c <= c2; c++)
    r.col[c-c1] = col[c].get_subvector(r1, r2);
  r.compact();

  return r;
}

template <class T>
Sparse_Mat<T> Sparse_Mat<T>::get_submatrix_cols(int c1, int c2) const
{
  it_assert_debug(c1 <= c2 && c1 >= 0 && c1 < n_cols, "Sparse_Mat<T>::get_submatrix_cols()");
  Sparse_Mat<T> r(n_rows, c2 - c1 + 1, 0);

  for (int c = c1; c <= c2; c++)
    r.col[c-c1] = col[c];
  r.compact();

  return r;
}

template <class T>
void Sparse_Mat<T>::get_col(int c, Sparse_Vec<T> &v) const
{
  it_assert(c >= 0 && c < n_cols, "Sparse_Mat<T>::get_col()");
  v = col[c];
}

template <class T>
Sparse_Vec<T> Sparse_Mat<T>::get_col(int c) const
{
  it_assert(c >= 0 && c < n_cols, "Sparse_Mat<T>::get_col()");
  return col[c];
}

template <class T>
void Sparse_Mat<T>::set_col(int c, const Sparse_Vec<T> &v)
{
  it_assert(c >= 0 && c < n_cols, "Sparse_Mat<T>::set_col()");
  col[c] = v;
}

template <class T>
void Sparse_Mat<T>::transpose(Sparse_Mat<T> &m) const
{
  m.set_size(n_cols, n_rows);
  for (int c = 0; c < n_cols; c++) {
    for (int p = 0; p < col[c].nnz(); p++)
      m.col[col[c].get_nz_index(p)].set_new(c, col[c].get_nz_data(p));
  }
}

template <class T>
Sparse_Mat<T> Sparse_Mat<T>::transpose() const
{
  Sparse_Mat<T> m;
  transpose(m);
  return m;
}

template <class T>
void Sparse_Mat<T>::operator=(const Sparse_Mat<T> &m)
{
  free();
  n_rows = m.n_rows;
  n_cols = m.n_cols;
  alloc_empty();

  for (int c = 0; c < n_cols; c++)
    col[c] = m.col[c];
}

template <class T>
void Sparse_Mat<T>::operator=(const Mat<T> &m)
{
  free();
  n_rows = m.rows();
  n_cols = m.cols();
  alloc();

  for (int c = 0; c < n_cols; c++) {
    for (int r = 0; r < n_rows; r++) {
      if (m(r, c) != T(0))
        col[c].set_new(r, m(r, c));
    }
    col[c].compact();
  }
}

template <class T>
Sparse_Mat<T> Sparse_Mat<T>::operator-() const
{
  Sparse_Mat r(n_rows, n_cols, 0);

  for (int c = 0; c < n_cols; c++) {
    r.col[c].resize_data(col[c].nnz());
    for (int p = 0; p < col[c].nnz(); p++)
      r.col[c].set_new(col[c].get_nz_index(p), -col[c].get_nz_data(p));
  }

  return r;
}

template <class T>
bool Sparse_Mat<T>::operator==(const Sparse_Mat<T> &m) const
{
  if (n_rows != m.n_rows || n_cols != m.n_cols)
    return false;
  for (int c = 0; c < n_cols; c++) {
    if (!(col[c] == m.col[c]))
      return false;
  }
  // If they passed all tests, they must be equal
  return true;
}

template <class T>
void Sparse_Mat<T>::operator+=(const Sparse_Mat<T> &m)
{
  it_assert_debug(m.rows() == n_rows && m.cols() == n_cols, "Addition of unequal sized matrices is not allowed");

  Sparse_Vec<T> v;
  for (int c = 0; c < n_cols; c++) {
    m.get_col(c, v);
    col[c] += v;
  }
}

template <class T>
void Sparse_Mat<T>::operator+=(const Mat<T> &m)
{
  it_assert_debug(m.rows() == n_rows && m.cols() == n_cols, "Addition of unequal sized matrices is not allowed");

  for (int c = 0; c < n_cols; c++)
    col[c] += (m.get_col(c));
}

template <class T>
void Sparse_Mat<T>::operator-=(const Sparse_Mat<T> &m)
{
  it_assert_debug(m.rows() == n_rows && m.cols() == n_cols, "Subtraction of unequal sized matrices is not allowed");

  Sparse_Vec<T> v;
  for (int c = 0; c < n_cols; c++) {
    m.get_col(c, v);
    col[c] -= v;
  }
}

template <class T>
void Sparse_Mat<T>::operator-=(const Mat<T> &m)
{
  it_assert_debug(m.rows() == n_rows && m.cols() == n_cols, "Subtraction of unequal sized matrices is not allowed");

  for (int c = 0; c < n_cols; c++)
    col[c] -= (m.get_col(c));
}

template <class T>
void Sparse_Mat<T>::operator*=(const T &m)
{
  for (int c = 0; c < n_cols; c++)
    col[c] *= m;
}

template <class T>
void Sparse_Mat<T>::operator/=(const T &m)
{
  for (int c = 0; c < n_cols; c++)
    col[c] /= m;
}

template <class T>
Sparse_Mat<T> operator+(const Sparse_Mat<T> &m1, const Sparse_Mat<T> &m2)
{
  it_assert_debug(m1.n_cols == m2.n_cols && m1.n_rows == m2.n_rows , "Sparse_Mat<T> + Sparse_Mat<T>");

  Sparse_Mat<T> m(m1.n_rows, m1.n_cols, 0);

  for (int c = 0; c < m.n_cols; c++)
    m.col[c] = m1.col[c] + m2.col[c];

  return m;
}

// This function added by EGL, May'05
template <class T>
Sparse_Mat<T> operator*(const T &c, const Sparse_Mat<T> &m)
{
  int i, j;
  Sparse_Mat<T> ret(m.n_rows, m.n_cols);
  for (j = 0; j < m.n_cols; j++) {
    for (i = 0; i < m.col[j].nnz(); i++) {
      T x = c * m.col[j].get_nz_data(i);
      int k = m.col[j].get_nz_index(i);
      ret.set_new(k, j, x);
    }
  }
  return ret;
}

template <class T>
Sparse_Mat<T> operator*(const Sparse_Mat<T> &m1, const Sparse_Mat<T> &m2)
{
  it_assert_debug(m1.n_cols == m2.n_rows, "Sparse_Mat<T> * Sparse_Mat<T>");

  Sparse_Mat<T> ret(m1.n_rows, m2.n_cols);

  for (int c = 0; c < m2.n_cols; c++) {
    Sparse_Vec<T> &m2colc = m2.col[c];
    for (int p2 = 0; p2 < m2colc.nnz(); p2++) {
      Sparse_Vec<T> &mcol = m1.col[m2colc.get_nz_index(p2)];
      T x = m2colc.get_nz_data(p2);
      for (int p1 = 0; p1 < mcol.nnz(); p1++) {
        int r = mcol.get_nz_index(p1);
        T inc = mcol.get_nz_data(p1) * x;
        ret.col[c].add_elem(r, inc);
      }
    }
  }
  // old code
  /*       for (int c=0; c<m2.n_cols; c++) { */
  /*  for (int p2=0; p2<m2.col[c].nnz(); p2++) { */
  /*    Sparse_Vec<T> &mcol = m1.col[m2.col[c].get_nz_index(p2)]; */
  /*    for (int p1=0; p1<mcol.nnz(); p1++) { */
  /*      int r = mcol.get_nz_index(p1); */
  /*      T inc = mcol.get_nz_data(p1) * m2.col[c].get_nz_data(p2); */
  /*      ret.col[c].add_elem(r,inc); */
  /*    } */
  /*  } */
  /*       } */
  ret.compact();
  return ret;
}


// This is apparently buggy.
/*   template <class T> */
/*     Sparse_Mat<T> operator*(const Sparse_Mat<T> &m1, const Sparse_Mat<T> &m2) */
/*     { */
/*       it_assert_debug(m1.n_cols == m2.n_rows, "Sparse_Mat<T> * Sparse_Mat<T>"); */

/*       Sparse_Mat<T> ret(m1.n_rows, m2.n_cols); */
/*       ivec occupied_by(ret.n_rows), pos(ret.n_rows); */
/*       for (int rp=0; rp<m1.n_rows; rp++) */
/*  occupied_by[rp] = -1; */
/*       for (int c=0; c<ret.n_cols; c++) { */
/*  Sparse_Vec<T> &m2col = m2.col[c]; */
/*  for (int p2=0; p2<m2col.nnz(); p2++) { */
/*    Sparse_Vec<T> &m1col = m1.col[m2col.get_nz_index(p2)]; */
/*    for (int p1=0; p1<m1col.nnz(); p1++) { */
/*      int r = m1col.get_nz_index(p1); */
/*      T inc = m1col.get_nz_data(p1) * m2col.get_nz_data(p2); */
/*      if (occupied_by[r] == c) { */
/*        int index=ret.col[c].get_nz_index(pos[r]); */
/*        ret.col[c].add_elem(index,inc); */
/*      } */
/*      else { */
/*        occupied_by[r] = c; */
/*        pos[r] = ret.col[c].nnz(); */
/*        ret.col[c].set_new(r, inc); */
/*      } */
/*    } */
/*  } */
/*       } */
/*       ret.compact(); */

/*       return ret; */
/*     } */


// This function added by EGL, May'05
template <class T>
Sparse_Vec<T> operator*(const Sparse_Mat<T> &m, const Sparse_Vec<T> &v)
{
  it_assert_debug(m.n_cols == v.size(), "Sparse_Mat<T> * Sparse_Vec<T>");

  Sparse_Vec<T> ret(m.n_rows);

  /* The two lines below added because the input parameter "v" is
  declared const, but the some functions (e.g., nnz()) change
  the vector... Is there a better workaround? */
  Sparse_Vec<T> vv = v;

  for (int p2 = 0; p2 < vv.nnz(); p2++) {
    Sparse_Vec<T> &mcol = m.col[vv.get_nz_index(p2)];
    T x = vv.get_nz_data(p2);
    for (int p1 = 0; p1 < mcol.nnz(); p1++) {
      int r = mcol.get_nz_index(p1);
      T inc = mcol.get_nz_data(p1) * x;
      ret.add_elem(r, inc);
    }
  }
  ret.compact();
  return ret;
}


template <class T>
Vec<T> operator*(const Sparse_Mat<T> &m, const Vec<T> &v)
{
  it_assert_debug(m.n_cols == v.size(), "Sparse_Mat<T> * Vec<T>");

  Vec<T> r(m.n_rows);
  r.clear();

  for (int c = 0; c < m.n_cols; c++) {
    for (int p = 0; p < m.col[c].nnz(); p++)
      r(m.col[c].get_nz_index(p)) += m.col[c].get_nz_data(p) * v(c);
  }

  return r;
}

template <class T>
Vec<T> operator*(const Vec<T> &v, const Sparse_Mat<T> &m)
{
  it_assert_debug(v.size() == m.n_rows, "Vec<T> * Sparse_Mat<T>");

  Vec<T> r(m.n_cols);
  r.clear();

  for (int c = 0; c < m.n_cols; c++)
    r[c] = v * m.col[c];

  return r;
}

template <class T>
Mat<T> trans_mult(const Sparse_Mat<T> &m)
{
  Mat<T> ret(m.n_cols, m.n_cols);
  Vec<T> col;
  for (int c = 0; c < ret.cols(); c++) {
    m.col[c].full(col);
    for (int r = 0; r < c; r++) {
      T tmp = m.col[r] * col;
      ret(r, c) = tmp;
      ret(c, r) = tmp;
    }
    ret(c, c) = m.col[c].sqr();
  }

  return ret;
}

template <class T>
Sparse_Mat<T> trans_mult_s(const Sparse_Mat<T> &m)
{
  Sparse_Mat<T> ret(m.n_cols, m.n_cols);
  Vec<T> col;
  T tmp;
  for (int c = 0; c < ret.n_cols; c++) {
    m.col[c].full(col);
    for (int r = 0; r < c; r++) {
      tmp = m.col[r] * col;
      if (tmp != T(0)) {
        ret.col[c].set_new(r, tmp);
        ret.col[r].set_new(c, tmp);
      }
    }
    tmp = m.col[c].sqr();
    if (tmp != T(0))
      ret.col[c].set_new(c, tmp);
  }

  return ret;
}

template <class T>
Sparse_Mat<T> trans_mult(const Sparse_Mat<T> &m1, const Sparse_Mat<T> &m2)
{
  it_assert_debug(m1.n_rows == m2.n_rows, "trans_mult()");

  Sparse_Mat<T> ret(m1.n_cols, m2.n_cols);
  Vec<T> col;
  for (int c = 0; c < ret.n_cols; c++) {
    m2.col[c].full(col);
    for (int r = 0; r < ret.n_rows; r++)
      ret.col[c].set_new(r, m1.col[r] * col);
  }

  return ret;
}

template <class T>
Vec<T> trans_mult(const Sparse_Mat<T> &m, const Vec<T> &v)
{
  Vec<T> r(m.n_cols);
  for (int c = 0; c < m.n_cols; c++)
    r(c) = m.col[c] * v;

  return r;
}

template <class T>
Sparse_Mat<T> mult_trans(const Sparse_Mat<T> &m1, const Sparse_Mat<T> &m2)
{
  return trans_mult(m1.transpose(), m2.transpose());
}

//! Convert a dense matrix \c m into its sparse representation
template <class T>
inline Sparse_Mat<T> sparse(const Mat<T> &m, T epsilon)
{
  Sparse_Mat<T> s(m, epsilon);
  return s;
}

//! Convert a sparse matrix \c s into its dense representation
template <class T>
inline Mat<T> full(const Sparse_Mat<T> &s)
{
  Mat<T> m;
  s.full(m);
  return m;
}

//! Transpose a sparse matrix \c s
template <class T>
inline Sparse_Mat<T> transpose(const Sparse_Mat<T> &s)
{
  Sparse_Mat<T> m;
  s.transpose(m);
  return m;
}

//! \cond

// ---------------------------------------------------------------------
// Instantiations
// ---------------------------------------------------------------------

ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Sparse_Mat<int>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Sparse_Mat<double>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Sparse_Mat<std::complex<double> >;

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_imat operator+(const sparse_imat &, const sparse_imat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_mat operator+(const sparse_mat &, const sparse_mat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_cmat operator+(const sparse_cmat &, const sparse_cmat &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_imat operator*(const sparse_imat &, const sparse_imat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_mat operator*(const sparse_mat &, const sparse_mat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_cmat operator*(const sparse_cmat &, const sparse_cmat &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec operator*(const ivec &, const sparse_imat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec operator*(const vec &, const sparse_mat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec operator*(const cvec &, const sparse_cmat &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec operator*(const sparse_imat &, const ivec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec operator*(const sparse_mat &, const vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec operator*(const sparse_cmat &, const cvec &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT imat trans_mult(const sparse_imat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT mat trans_mult(const sparse_mat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cmat trans_mult(const sparse_cmat &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_imat trans_mult_s(const sparse_imat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_mat trans_mult_s(const sparse_mat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_cmat trans_mult_s(const sparse_cmat &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_imat trans_mult(const sparse_imat &, const sparse_imat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_mat trans_mult(const sparse_mat &, const sparse_mat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_cmat trans_mult(const sparse_cmat &, const sparse_cmat &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec trans_mult(const sparse_imat &, const ivec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec trans_mult(const sparse_mat &, const vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec trans_mult(const sparse_cmat &, const cvec &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_imat mult_trans(const sparse_imat &, const sparse_imat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_mat mult_trans(const sparse_mat &, const sparse_mat &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_cmat mult_trans(const sparse_cmat &, const sparse_cmat &);

//! \endcond

} // namespace itpp

#endif // #ifndef SMAT_H

