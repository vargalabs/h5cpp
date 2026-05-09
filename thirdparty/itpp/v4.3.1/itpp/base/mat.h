/*!
 * \file
 * \brief Matrix Class Definitions
 * \author Tony Ottosson, Tobias Ringstrom, Adam Piatyszek and Conrad Sanderson
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

#ifndef MAT_H
#define MAT_H

#include <itpp/base/itassert.h>
#include <itpp/base/math/misc.h>
#include <itpp/base/factory.h>
#include <itpp/itexports.h>

namespace itpp
{

// Declaration of Vec
template<class Num_T> class Vec;
// Declaration of Mat
template<class Num_T> class Mat;
// Declaration of bin
class bin;

//! Horizontal concatenation of two matrices
template<class Num_T>
Mat<Num_T> concat_horizontal(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
//! Vertical concatenation of two matrices
template<class Num_T>
Mat<Num_T> concat_vertical(const Mat<Num_T> &m1, const Mat<Num_T> &m2);

//! Addition of two matrices
template<class Num_T>
Mat<Num_T> operator+(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
//! Addition of a matrix and a scalar
template<class Num_T>
Mat<Num_T> operator+(const Mat<Num_T> &m, Num_T t);
//! Addition of a scalar and a matrix
template<class Num_T>
Mat<Num_T> operator+(Num_T t, const Mat<Num_T> &m);

//! Subtraction of two matrices
template<class Num_T>
Mat<Num_T> operator-(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
//! Subtraction of matrix and scalar
template<class Num_T>
Mat<Num_T> operator-(const Mat<Num_T> &m, Num_T t);
//! Subtraction of scalar and matrix
template<class Num_T>
Mat<Num_T> operator-(Num_T t, const Mat<Num_T> &m);
//! Negation of matrix
template<class Num_T>
Mat<Num_T> operator-(const Mat<Num_T> &m);

//! Multiplication of two matrices
template<class Num_T>
Mat<Num_T> operator*(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
//! Multiplication of matrix and vector
template<class Num_T>
Vec<Num_T> operator*(const Mat<Num_T> &m, const Vec<Num_T> &v);
//! Multiplication of matrix and scalar
template<class Num_T>
Mat<Num_T> operator*(const Mat<Num_T> &m, Num_T t);
//! Multiplication of scalar and matrix
template<class Num_T>
Mat<Num_T> operator*(Num_T t, const Mat<Num_T> &m);

//! Element wise multiplication of two matrices
template<class Num_T>
Mat<Num_T> elem_mult(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
//! Element wise multiplication of two matrices, storing the result in matrix \c out
template<class Num_T>
void elem_mult_out(const Mat<Num_T> &m1, const Mat<Num_T> &m2,
                   Mat<Num_T> &out);
//! Element wise multiplication of three matrices, storing the result in matrix \c out
template<class Num_T>
void elem_mult_out(const Mat<Num_T> &m1, const Mat<Num_T> &m2,
                   const Mat<Num_T> &m3, Mat<Num_T> &out);
//! Element wise multiplication of four matrices, storing the result in matrix \c out
template<class Num_T>
void elem_mult_out(const Mat<Num_T> &m1, const Mat<Num_T> &m2,
                   const Mat<Num_T> &m3, const Mat<Num_T> &m4,
                   Mat<Num_T> &out);
//! In-place element wise multiplication of two matrices. Fast version of B = elem_mult(A, B).
template<class Num_T>
void elem_mult_inplace(const Mat<Num_T> &m1, Mat<Num_T> &m2);
//! Element wise multiplication of two matrices, followed by summation of the resultant elements. Fast version of sumsum(elem_mult(A, B)).
template<class Num_T>
Num_T elem_mult_sum(const Mat<Num_T> &m1, const Mat<Num_T> &m2);

//! Element-wise division by a scalar
template<class Num_T>
Mat<Num_T> operator/(const Mat<Num_T> &m, Num_T t);
//! Element-wise division (\c t is the dividend, elements of \c m are divisors)
template<class Num_T>
Mat<Num_T> operator/(Num_T t, const Mat<Num_T> &m);

//! Element wise division of two matrices
template<class Num_T>
Mat<Num_T> elem_div(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
//! Element wise division of two matrices, storing the result in matrix \c out
template<class Num_T>
void elem_div_out(const Mat<Num_T> &m1, const Mat<Num_T> &m2,
                  Mat<Num_T> &out);
//! Element wise division of two matrices, followed by summation of the resultant elements. Fast version of sumsum(elem_div(A, B)).
template<class Num_T>
Num_T elem_div_sum(const Mat<Num_T> &m1, const Mat<Num_T> &m2);

// -------------------------------------------------------------------------------------
// Declaration of Mat
// -------------------------------------------------------------------------------------

/*!
  \ingroup arr_vec_mat
  \brief Matrix Class (Templated)
  \author Tony Ottosson, Tobias Ringstrom, Adam Piatyszek and Conrad Sanderson

  Matrices can be of arbitrarily types, but conversions and functions are
  prepared for \c bin, \c short, \c int, \c double, and \c complex<double>
  vectors and these are predefined as: \c bmat, \c smat, \c imat, \c mat,
  and \c cmat. \c double and \c complex<double> are usually \c double and
  \c complex<double> respectively. However, this can be changed when
  compiling the it++ (see installation notes for more details). (Note: for
  binary matrices, an alternative to the bmat class is \c GF2mat and
  \c GF2mat_dense, which offer a more memory efficient representation and
  additional functions for linear algebra.)

  Examples:

  Matrix Constructors:
  When constructing a matrix without dimensions (memory) use
  \code mat temp; \endcode
  For construction of a matrix of a given size use
  \code mat temp(rows, cols); \endcode
  It is also possible to assign the constructed matrix the value and dimension
  of another matrix by
  \code vec temp(inmatrix); \endcode
  If you have explicit values you would like to assign to the matrix it is
  possible to do this using strings as:
  \code
  mat a("0 0.7;5 9.3"); // that is a = [0, 0.7; 5, 9.3]
  mat a="0 0.7;5 9.3";  // the constructor are called implicitly
  \endcode
  It is also possible to change dimension by
  \code temp.set_size(new_rows, new_cols, false); \endcode
  where \c false is used to indicate that the old values in \c temp
  is not copied. If you like to preserve the values use \c true.

  There are a number of methods to access parts of a matrix. Examples are
  \code
  a(5,3);     // Element number (5,3)
  a(5,9,3,5);  // Sub-matrix from rows 5, 6, 7, 8, 9 the columns 3, 4, and 5
  a.get_row(10);  // Row 10
  a.get_col(10); // Column 10
  \endcode

  It is also possible to modify parts of a vector as e.g. in
  \code
  a.set_row(5, invector);    // Set row 5 to \c invector
  a.set_col(3, invector); // Set column 3 to \c invector
  a.copy_col(1, 5); // Copy column 5 to column 1
  a.swap_cols(1, 5); // Swap the contents of columns 1 and 5
  \endcode

  It is of course also possible to perform the common linear algebra
  methods such as addition, subtraction, and matrix multiplication. Observe
  though, that vectors are assumed to be column-vectors in operations with
  matrices.

  Most elementary functions such as sin(), cosh(), log(), abs(), ..., are
  available as operations on the individual elements of the matrices. Please
  see the individual functions for more details.

  By default, the Mat elements are created using the default constructor for
  the element type. This can be changed by specifying a suitable Factory in
  the Mat constructor call; see Detailed Description for Factory.
*/
template<class Num_T>
class Mat
{
public:
  //! The type of the matrix values
  typedef Num_T value_type;

  //! Default constructor. An element factory \c f can be specified
  explicit Mat(const Factory &f = DEFAULT_FACTORY);
  //! Create a matrix of size (rows, cols). An element factory \c f can be specified.
  Mat(int rows, int cols, const Factory &f = DEFAULT_FACTORY);
  //! Copy constructor
  Mat(const Mat<Num_T> &m);
  //! Constructor, similar to the copy constructor, but also takes an element factory \c f as argument
  Mat(const Mat<Num_T> &m, const Factory &f);
  //! Construct a matrix from a column vector \c v. An element factory \c f can be specified.
  Mat(const Vec<Num_T> &v, const Factory &f = DEFAULT_FACTORY);
  //! Set matrix equal to values in string \c str. An element factory \c f can be specified.
  Mat(const std::string &str, const Factory &f = DEFAULT_FACTORY);
  //! Set matrix equal to values in string \c str. An element factory \c f can be specified.
  Mat(const char *str, const Factory &f = DEFAULT_FACTORY);
  /*!
   * \brief Constructor taking a C-array as input. An element factory \c f
   * can be specified.
   *
   * By default the matrix is stored as a row-major matrix (i.e. listing
   * elements in sequence beginning with the first column).
   */
  Mat(const Num_T *c_array, int rows, int cols, bool row_major = true,
      const Factory &f = DEFAULT_FACTORY);

  //! Destructor
  ~Mat();

  //! The number of columns
  int cols() const { return no_cols; }
  //! The number of rows
  int rows() const { return no_rows; }
  //! The number of elements
  int size() const { return datasize; }
  //! Set size of matrix. If copy = true then keep the data before resizing.
  void set_size(int rows, int cols, bool copy = false);
  //! Set matrix equal to the all zero matrix
  void zeros();
  //! Set matrix equal to the all zero matrix
  void clear() { zeros(); }
  //! Set matrix equal to the all one matrix
  void ones();
  //! Set matrix equal to values in the string \c str
  void set(const std::string &str);
  //! Set matrix equal to values in the string \c str
  void set(const char *str);

  //! Get element (r,c) from matrix
  const Num_T &operator()(int r, int c) const;
  //! Get element (r,c) from matrix
  Num_T &operator()(int r, int c);
  //! Get element \c i using linear addressing (by rows)
  const Num_T &operator()(int i) const;
  //! Get element \c i using linear addressing (by rows)
  Num_T &operator()(int i);
  //! Get element (r,c) from matrix
  const Num_T &get(int r, int c) const;
  //! Get element \c i using linear addressing (by rows)
  const Num_T &get(int i) const;
  //! Set element (r,c) of matrix
  void set(int r, int c, Num_T t);

  /*!
    \brief Sub-matrix from row \c r1 to row \c r2 and columns \c c1 to \c c2.

    Value -1 indicates the last row and column, respectively.
  */
  Mat<Num_T> operator()(int r1, int r2, int c1, int c2) const;
  /*!
    \brief Sub-matrix from row \c r1 to row \c r2 and columns \c c1 to \c c2.

    Value -1 indicates the last row and column, respectively.
  */
  Mat<Num_T> get(int r1, int r2, int c1, int c2) const;

  //! Get row \c r
  Vec<Num_T> get_row(int r) const;
  //! Get rows \c r1 through \c r2
  Mat<Num_T> get_rows(int r1, int r2) const;
  //! Get the rows specified by \c indexlist
  Mat<Num_T> get_rows(const Vec<int> &indexlist) const;
  //! Get column \c c
  Vec<Num_T> get_col(int c) const;
  //! Get columns \c c1 through \c c2
  Mat<Num_T> get_cols(int c1, int c2) const;
  //! Get the columns specified by \c indexlist
  Mat<Num_T> get_cols(const Vec<int> &indexlist) const;
  //! Set row \c r to vector \c v
  void set_row(int r, const Vec<Num_T> &v);
  //! Set column \c c to vector \c v
  void set_col(int c, const Vec<Num_T> &v);
  //! Set rows to matrix \c m, staring from row \c r
  void set_rows(int r, const Mat<Num_T> &m);
  //! Set columns to matrix \c m, starting from column \c c
  void set_cols(int c, const Mat<Num_T> &m);
  //! Copy row \c from onto row \c to
  void copy_row(int to, int from);
  //! Copy column \c from onto column \c to
  void copy_col(int to, int from);
  //! Swap the rows \c r1 and \c r2
  void swap_rows(int r1, int r2);
  //! Swap the columns \c c1 and \c c2
  void swap_cols(int c1, int c2);

  //! This function is deprecated. Please use set_submatrix(int r, int c, const Mat<> &m) instead.
  void set_submatrix(int r1, int r2, int c1, int c2, const Mat<Num_T> &m);
  //! Set submatrix defined by upper-left element (r,c) and the size of matrix m to m
  void set_submatrix(int r, int c, const Mat<Num_T> &m);
  //! Set all elements of submatrix defined by rows r1,r2 and columns c1,c2 to value t
  void set_submatrix(int r1, int r2, int c1, int c2, Num_T t);

  //! Delete row number \c r
  void del_row(int r);
  //! Delete rows from \c r1 to \c r2
  void del_rows(int r1, int r2);
  //! Delete column number \c c
  void del_col(int c);
  //! Delete columns from \c c1 to \c c2
  void del_cols(int c1, int c2);
  //! Insert vector \c v at row number \c r. The matrix can be empty.
  void ins_row(int r, const Vec<Num_T> &v);
  //! Insert vector \c v at column number \c c. The matrix can be empty.
  void ins_col(int c, const Vec<Num_T> &v);
  //! Append vector \c v to the bottom of the matrix. The matrix can be empty.
  void append_row(const Vec<Num_T> &v);
  //! Append vector \c v to the right side of the matrix. The matrix can be empty.
  void append_col(const Vec<Num_T> &v);

  //! Matrix transpose
  Mat<Num_T> transpose() const;
  //! Matrix transpose
  Mat<Num_T> T() const { return this->transpose(); }
  //! Hermitian matrix transpose (conjugate transpose)
  Mat<Num_T> hermitian_transpose() const;
  //! Hermitian matrix transpose (conjugate transpose)
  Mat<Num_T> H() const { return this->hermitian_transpose(); }

  //! Concatenate the matrices \c m1 and \c m2 horizontally
  friend Mat<Num_T> concat_horizontal<>(const Mat<Num_T> &m1,
                                        const Mat<Num_T> &m2);
  //! Concatenate the matrices \c m1 and \c m2 vertically
  friend Mat<Num_T> concat_vertical<>(const Mat<Num_T> &m1,
                                      const Mat<Num_T> &m2);

  //! Set all elements of the matrix equal to \c t
  Mat<Num_T>& operator=(Num_T t);
  //! Set matrix equal to \c m
  Mat<Num_T>& operator=(const Mat<Num_T> &m);
  //! Set matrix equal to the vector \c v, assuming column vector
  Mat<Num_T>& operator=(const Vec<Num_T> &v);
  //! Set matrix equal to values in the string \c str
  Mat<Num_T>& operator=(const std::string &str);
  //! Set matrix equal to values in the string \c str
  Mat<Num_T>& operator=(const char *str);

  //! Addition of matrices
  Mat<Num_T>& operator+=(const Mat<Num_T> &m);
  //! Addition of scalar to matrix
  Mat<Num_T>& operator+=(Num_T t);
  //! Addition of two matrices
  friend Mat<Num_T> operator+<>(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
  //! Addition of matrix and scalar
  friend Mat<Num_T> operator+<>(const Mat<Num_T> &m, Num_T t);
  //! Addition of scalar and matrix
  friend Mat<Num_T> operator+<>(Num_T t, const Mat<Num_T> &m);

  //! Subtraction of matrix
  Mat<Num_T>& operator-=(const Mat<Num_T> &m);
  //! Subtraction of scalar from matrix
  Mat<Num_T>& operator-=(Num_T t);
  //! Subtraction of \c m2 from \c m1
  friend Mat<Num_T> operator-<>(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
  //! Subtraction of scalar from matrix
  friend Mat<Num_T> operator-<>(const Mat<Num_T> &m, Num_T t);
  //! Subtract matrix from scalar
  friend Mat<Num_T> operator-<>(Num_T t, const Mat<Num_T> &m);
  //! Subtraction of matrix
  friend Mat<Num_T> operator-<>(const Mat<Num_T> &m);

  //! Matrix multiplication
  Mat<Num_T>& operator*=(const Mat<Num_T> &m);
  //! Multiplication by a scalar
  Mat<Num_T>& operator*=(Num_T t);

  //! Element wise multiplication of two matrices
  friend Mat<Num_T> elem_mult<>(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
  //! Element wise multiplication of two matrices, storing the result in matrix \c out
  friend void elem_mult_out<>(const Mat<Num_T> &m1, const Mat<Num_T> &m2,
                              Mat<Num_T> &out);
  //! Element wise multiplication of three matrices, storing the result in matrix \c out
  friend void elem_mult_out<>(const Mat<Num_T> &m1, const Mat<Num_T> &m2,
                              const Mat<Num_T> &m3, Mat<Num_T> &out);
  //! Element wise multiplication of four matrices, storing the result in matrix \c out
  friend void elem_mult_out<>(const Mat<Num_T> &m1, const Mat<Num_T> &m2,
                              const Mat<Num_T> &m3, const Mat<Num_T> &m4,
                              Mat<Num_T> &out);
  //! In-place element wise multiplication of two matrices. Fast version of B = elem_mult(A, B).
  friend void elem_mult_inplace<>(const Mat<Num_T> &m1, Mat<Num_T> &m2);
  //! Element wise multiplication of two matrices, followed by summation of the resultant elements. Fast version of sumsum(elem_mult(A, B)).
  friend Num_T elem_mult_sum<>(const Mat<Num_T> &m1, const Mat<Num_T> &m2);

  //! Division by a scalar
  Mat<Num_T>& operator/=(Num_T t);
  //! Element-wise division with the current matrix
  Mat<Num_T>& operator/=(const Mat<Num_T> &m);

  //! Element-wise division by a scalar
  friend Mat<Num_T> operator/<>(const Mat<Num_T> &m, Num_T t);
  //! Element-wise division (\c t is the dividend, elements of \c m are divisors)
  friend Mat<Num_T> operator/<>(Num_T t, const Mat<Num_T> &m);

  //! Element wise division of two matrices
  friend Mat<Num_T> elem_div<>(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
  //! Element wise division of two matrices, storing the result in matrix \c out
  friend void elem_div_out<>(const Mat<Num_T> &m1, const Mat<Num_T> &m2,
                             Mat<Num_T> &out);
  //! Element wise division of two matrices, followed by summation of the resultant elements. Fast version of sumsum(elem_div(A, B)).
  friend Num_T elem_div_sum<>(const Mat<Num_T> &m1, const Mat<Num_T> &m2);

  //! Compare two matrices. False if wrong sizes or different values
  bool operator==(const Mat<Num_T> &m) const;
  //! Compare two matrices. True if different
  bool operator!=(const Mat<Num_T> &m) const;

  //! Get element (r,c) from matrix without boundary check (not recommended to use)
  Num_T &_elem(int r, int c) { return data[r+c*no_rows]; }
  //! Get element (r,c) from matrix without boundary check (not recommended to use)
  const Num_T &_elem(int r, int c) const { return data[r+c*no_rows]; }
  //! Get element \c i using linear addressing (by rows) without boundary check (not recommended to use)
  Num_T &_elem(int i) { return data[i]; }
  //! Get element \c i using linear addressing (by rows) without boundary check (not recommended to use)
  const Num_T &_elem(int i) const { return data[i]; }

  //! Access of the internal data structure (not recommended to use)
  Num_T *_data() { return data; }
  //! Access to the internal data structure (not recommended to use)
  const Num_T *_data() const { return data; }
  //! Access to the internal data structure (not recommended to use)
  int _datasize() const { return datasize; }

protected:
  //! Allocate memory for the matrix
  void alloc(int rows, int cols);
  //! Free the memory space of the matrix
  void free();

  /*! Protected integer variables
   * @{ */
  int datasize, no_rows, no_cols;
  /*! @} */
  //! Protected data pointer
  Num_T *data;
  //! Element factory (set to DEFAULT_FACTORY to use Num_T default constructors only)
  const Factory &factory;

private:
  //! Check whether element (r,c) is within the matrix
  bool in_range(int r, int c) const {
    return ((r >= 0) && (r < no_rows) && (c >= 0) && (c < no_cols));
  }
  //! Check whether row \c r is in the allowed range
  bool row_in_range(int r) const { return ((r >= 0) && (r < no_rows)); }
  //! Check whether column \c c is in the allowed range
  bool col_in_range(int c) const { return ((c >= 0) && (c < no_cols)); }
  //! Check whether element \c i is in the allowed range
  bool in_range(int i) const { return ((i >= 0) && (i < datasize)); }
};

// -------------------------------------------------------------------------------------
// Type definitions of mat, cmat, imat, smat, and bmat
// -------------------------------------------------------------------------------------

/*!
  \relates Mat
  \brief Default Matrix Type
*/
typedef Mat<double> mat;

/*!
  \relates Mat
  \brief Default Complex Matrix Type
*/
typedef Mat<std::complex<double> > cmat;

/*!
  \relates Mat
  \brief Integer matrix
*/
typedef Mat<int> imat;

/*!
  \relates Mat
  \brief short int matrix
*/
typedef Mat<short int> smat;

/*!
  \relates Mat
  \relates GF2mat
  \relates GF2mat_sparse
  \brief bin matrix
*/
typedef Mat<bin> bmat;

} //namespace itpp


#include <itpp/base/vec.h>

namespace itpp
{

// ----------------------------------------------------------------------
// Declaration of input and output streams for Mat
// ----------------------------------------------------------------------

/*!
  \relatesalso Mat
  \brief Output stream for matrices
*/
template <class Num_T>
std::ostream &operator<<(std::ostream &os, const Mat<Num_T> &m);

/*!
  \relatesalso Mat
  \brief Input stream for matrices

  The input can be on the form "1 2 3; 4 5 6" or "[[1 2 3][4 5 6]]", i.e. with
  brackets or semicolons as row delimiters. The first form is compatible with
  the set method, while the second form is compatible with the ostream
  operator. The elements on a row can be separated by blank space or commas.
  Rows that are shorter than the longest row are padded with zero elements.
  "[]" means an empty matrix.
*/
template <class Num_T>
std::istream &operator>>(std::istream &is, Mat<Num_T> &m);

// ----------------------------------------------------------------------
// Implementation of templated Mat members and friends
// ----------------------------------------------------------------------

template<class Num_T> inline
void Mat<Num_T>::alloc(int rows, int cols)
{
  if ((rows > 0) && (cols > 0)) {
    datasize = rows * cols;
    no_rows = rows;
    no_cols = cols;
    create_elements(data, datasize, factory);
  }
  else {
    data = 0;
    datasize = 0;
    no_rows = 0;
    no_cols = 0;
  }
}

template<class Num_T> inline
void Mat<Num_T>::free()
{
  destroy_elements(data, datasize);
  datasize = 0;
  no_rows = 0;
  no_cols = 0;
}


template<class Num_T> inline
Mat<Num_T>::Mat(const Factory &f) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(f) {}

template<class Num_T> inline
Mat<Num_T>::Mat(int rows, int cols, const Factory &f) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(f)
{
  it_assert_debug((rows >= 0) && (cols >= 0), "Mat<>::Mat(): Wrong size");
  alloc(rows, cols);
}

template<class Num_T> inline
Mat<Num_T>::Mat(const Mat<Num_T> &m) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(m.factory)
{
  alloc(m.no_rows, m.no_cols);
  copy_vector(m.datasize, m.data, data);
}

template<class Num_T> inline
Mat<Num_T>::Mat(const Mat<Num_T> &m, const Factory &f) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(f)
{
  alloc(m.no_rows, m.no_cols);
  copy_vector(m.datasize, m.data, data);
}

template<class Num_T> inline
Mat<Num_T>::Mat(const Vec<Num_T> &v, const Factory &f) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(f)
{
  int size = v.size();
  alloc(size, 1);
  copy_vector(size, v._data(), data);
}

template<class Num_T> inline
Mat<Num_T>::Mat(const std::string &str, const Factory &f) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(f)
{
  set(str);
}

template<class Num_T> inline
Mat<Num_T>::Mat(const char *str, const Factory &f) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(f)
{
  set(std::string(str));
}

template<class Num_T>
Mat<Num_T>::Mat(const Num_T *c_array, int rows, int cols, bool row_major,
                const Factory &f):
    datasize(0), no_rows(0), no_cols(0), data(0), factory(f)
{
  alloc(rows, cols);
  if (!row_major)
    copy_vector(datasize, c_array, data);
  else
    for (int i = 0; i < rows; i++)
      for (int j = 0; j < cols; j++)
        data[i+j*no_rows] = c_array[i*no_cols+j];
}

template<class Num_T> inline
Mat<Num_T>::~Mat()
{
  free();
}


template<class Num_T>
void Mat<Num_T>::set_size(int rows, int cols, bool copy)
{
  it_assert_debug((rows >= 0) && (cols >= 0),
                  "Mat<>::set_size(): Wrong size");
  // check if we have to resize the current matrix
  if ((no_rows == rows) && (no_cols == cols))
    return;
  // check if one of dimensions is zero
  if ((rows == 0) || (cols == 0)) {
    free();
    return;
  }
  // conditionally copy previous matrix content
  if (copy) {
    // create a temporary pointer to the allocated data
    Num_T* tmp = data;
    // store the current number of elements and number of rows
    int old_datasize = datasize;
    int old_rows = no_rows;
    // check the boundaries of the copied data
    int min_r = (no_rows < rows) ? no_rows : rows;
    int min_c = (no_cols < cols) ? no_cols : cols;
    // allocate new memory
    alloc(rows, cols);
    // copy the previous data into the allocated memory
    for (int i = 0; i < min_c; ++i) {
      copy_vector(min_r, &tmp[i*old_rows], &data[i*no_rows]);
    }
    // fill-in the rest of matrix with zeros
    for (int i = min_r; i < rows; ++i)
      for (int j = 0; j < cols; ++j)
        data[i+j*rows] = Num_T(0);
    for (int j = min_c; j < cols; ++j)
      for (int i = 0; i < min_r; ++i)
        data[i+j*rows] = Num_T(0);
    // delete old elements
    destroy_elements(tmp, old_datasize);
  }
  // if possible, reuse the allocated memory
  else if (datasize == rows * cols) {
    no_rows = rows;
    no_cols = cols;
  }
  // finally release old memory and allocate a new one
  else {
    free();
    alloc(rows, cols);
  }
}

template<class Num_T> inline
void Mat<Num_T>::zeros()
{
  for (int i = 0; i < datasize; i++)
    data[i] = Num_T(0);
}

template<class Num_T> inline
void Mat<Num_T>::ones()
{
  for (int i = 0; i < datasize; i++)
    data[i] = Num_T(1);
}

template<class Num_T> inline
const Num_T& Mat<Num_T>::operator()(int r, int c) const
{
  it_assert_debug(in_range(r, c),
                  "Mat<>::operator(): Indexing out of range");
  return data[r+c*no_rows];
}

template<class Num_T> inline
Num_T& Mat<Num_T>::operator()(int r, int c)
{
  it_assert_debug(in_range(r, c),
                  "Mat<>::operator(): Indexing out of range");
  return data[r+c*no_rows];
}

template<class Num_T> inline
Num_T& Mat<Num_T>::operator()(int i)
{
  it_assert_debug(in_range(i), "Mat<>::operator(): Index out of range");
  return data[i];
}

template<class Num_T> inline
const Num_T& Mat<Num_T>::operator()(int i) const
{
  it_assert_debug(in_range(i), "Mat<>::operator(): Index out of range");
  return data[i];
}

template<class Num_T> inline
const Num_T& Mat<Num_T>::get(int r, int c) const
{
  return (*this)(r, c);
}

template<class Num_T> inline
const Num_T& Mat<Num_T>::get(int i) const
{
  return (*this)(i);
}

template<class Num_T> inline
void Mat<Num_T>::set(int r, int c, Num_T t)
{
  it_assert_debug(in_range(r, c), "Mat<>::set(): Indexing out of range");
  data[r+c*no_rows] = t;
}


template<class Num_T>
void Mat<Num_T>::set(const std::string &str)
{
  // actual row counter
  int rows = 0;
  // number of rows to allocate next time (8, 16, 32, 64, etc.)
  int maxrows = 8;

  // clean the current matrix content
  free();

  // variable to store the start of a current vector
  std::string::size_type beg = 0;
  std::string::size_type end = 0;
  while (end != std::string::npos) {
    // find next occurrence of a semicolon in string str
    end = str.find(';', beg);
    // parse first row into a vector v
    Vec<Num_T> v(str.substr(beg, end - beg));
    int v_size = v.size();

    // this check is necessary to parse the following two strings as the
    // same matrix: "1 0 1; ; 1 1; " and "1 0 1; 0 0 0; 1 1 0"
    if ((end != std::string::npos) || (v_size > 0)) {
      // matrix empty -> insert v as a first row and allocate maxrows
      if (rows == 0) {
        set_size(maxrows, v_size, true);
        set_row(rows++, v);
      }
      else {
        // check if we need to resize the matrix
        if ((rows == maxrows) || (v_size != no_cols)) {
          // we need to add new rows
          if (rows == maxrows) {
            maxrows *= 2;
          }
          // check if we need to add new columns
          if (v_size > no_cols) {
            set_size(maxrows, v_size, true);
          }
          else {
            set_size(maxrows, no_cols, true);
            // set the size of the parsed vector to the number of columns
            v.set_size(no_cols, true);
          }
        }
        // set the parsed vector as the next row
        set_row(rows++, v);
      }
    }
    // update the starting position of the next vector in the parsed
    // string
    beg = end + 1;
  } // if ((end != std::string::npos) || (v.size > 0))

  set_size(rows, no_cols, true);
}

template<class Num_T> inline
void Mat<Num_T>::set(const char *str)
{
  set(std::string(str));
}

template<class Num_T> inline
Mat<Num_T> Mat<Num_T>::operator()(int r1, int r2, int c1, int c2) const
{
  if (r1 == -1) r1 = no_rows - 1;
  if (r2 == -1) r2 = no_rows - 1;
  if (c1 == -1) c1 = no_cols - 1;
  if (c2 == -1) c2 = no_cols - 1;

  it_assert_debug((r1 >= 0) && (r1 <= r2) && (r2 < no_rows) &&
                  (c1 >= 0) && (c1 <= c2) && (c2 < no_cols),
                  "Mat<>::operator()(r1, r2, c1, c2): Wrong indexing");

  Mat<Num_T> s(r2 - r1 + 1, c2 - c1 + 1);

  for (int i = 0;i < s.no_cols;i++)
    copy_vector(s.no_rows, data + r1 + (c1 + i)*no_rows, s.data + i*s.no_rows);

  return s;
}

template<class Num_T> inline
Mat<Num_T> Mat<Num_T>::get(int r1, int r2, int c1, int c2) const
{
  return (*this)(r1, r2, c1, c2);
}

template<class Num_T> inline
Vec<Num_T> Mat<Num_T>::get_row(int r) const
{
  it_assert_debug(row_in_range(r), "Mat<>::get_row(): Index out of range");
  Vec<Num_T> a(no_cols);

  copy_vector(no_cols, data + r, no_rows, a._data(), 1);
  return a;
}

template<class Num_T>
Mat<Num_T> Mat<Num_T>::get_rows(int r1, int r2) const
{
  it_assert_debug((r1 >= 0) && (r1 <= r2) && (r2 < no_rows),
                  "Mat<>::get_rows(): Wrong indexing");
  Mat<Num_T> m(r2 - r1 + 1, no_cols);

  for (int i = 0; i < m.rows(); i++)
    copy_vector(no_cols, data + i + r1, no_rows, m.data + i, m.no_rows);

  return m;
}

template<class Num_T>
Mat<Num_T> Mat<Num_T>::get_rows(const Vec<int> &indexlist) const
{
  Mat<Num_T> m(indexlist.size(), no_cols);

  for (int i = 0;i < indexlist.size();i++) {
    it_assert_debug(row_in_range(indexlist(i)),
                    "Mat<>::get_rows(indexlist): Indexing out of range");
    copy_vector(no_cols, data + indexlist(i), no_rows, m.data + i, m.no_rows);
  }

  return m;
}

template<class Num_T> inline
Vec<Num_T> Mat<Num_T>::get_col(int c) const
{
  it_assert_debug(col_in_range(c), "Mat<>::get_col(): Index out of range");
  Vec<Num_T> a(no_rows);

  copy_vector(no_rows, data + c*no_rows, a._data());

  return a;
}

template<class Num_T>
Mat<Num_T> Mat<Num_T>::get_cols(int c1, int c2) const
{
  it_assert_debug((c1 >= 0) && (c1 <= c2) && (c2 < no_cols),
                  "Mat<>::get_cols(): Wrong indexing");
  Mat<Num_T> m(no_rows, c2 - c1 + 1);

  for (int i = 0; i < m.cols(); i++)
    copy_vector(no_rows, data + (i + c1)*no_rows, m.data + i*m.no_rows);

  return m;
}

template<class Num_T>
Mat<Num_T> Mat<Num_T>::get_cols(const Vec<int> &indexlist) const
{
  Mat<Num_T> m(no_rows, indexlist.size());

  for (int i = 0; i < indexlist.size(); i++) {
    it_assert_debug(col_in_range(indexlist(i)),
                    "Mat<>::get_cols(indexlist): Indexing out of range");
    copy_vector(no_rows, data + indexlist(i)*no_rows, m.data + i*m.no_rows);
  }

  return m;
}

template<class Num_T> inline
void Mat<Num_T>::set_row(int r, const Vec<Num_T> &v)
{
  it_assert_debug(row_in_range(r), "Mat<>::set_row(): Index out of range");
  it_assert_debug(v.size() == no_cols,
                  "Mat<>::set_row(): Wrong size of input vector");
  copy_vector(v.size(), v._data(), 1, data + r, no_rows);
}

template<class Num_T> inline
void Mat<Num_T>::set_col(int c, const Vec<Num_T> &v)
{
  it_assert_debug(col_in_range(c), "Mat<>::set_col(): Index out of range");
  it_assert_debug(v.size() == no_rows,
                  "Mat<>::set_col(): Wrong size of input vector");
  copy_vector(v.size(), v._data(), data + c*no_rows);
}


template<class Num_T>
void Mat<Num_T>::set_rows(int r, const Mat<Num_T> &m)
{
  it_assert_debug(row_in_range(r), "Mat<>::set_rows(): Index out of range");
  it_assert_debug(no_cols == m.cols(),
                  "Mat<>::set_rows(): Column sizes do not match");
  it_assert_debug(m.rows() + r <= no_rows,
                  "Mat<>::set_rows(): Not enough rows");

  for (int i = 0; i < m.rows(); ++i) {
    copy_vector(no_cols, m.data + i, m.no_rows, data + i + r, no_rows);
  }
}

template<class Num_T>
void Mat<Num_T>::set_cols(int c, const Mat<Num_T> &m)
{
  it_assert_debug(col_in_range(c), "Mat<>::set_cols(): Index out of range");
  it_assert_debug(no_rows == m.rows(),
                  "Mat<>::set_cols(): Row sizes do not match");
  it_assert_debug(m.cols() + c <= no_cols,
                  "Mat<>::set_cols(): Not enough colums");

  for (int i = 0; i < m.cols(); ++i) {
    copy_vector(no_rows, m.data + i*no_rows, data + (i + c)*no_rows);
  }
}


template<class Num_T> inline
void Mat<Num_T>::copy_row(int to, int from)
{
  it_assert_debug(row_in_range(to) && row_in_range(from),
                  "Mat<>::copy_row(): Indexing out of range");
  if (from == to)
    return;

  copy_vector(no_cols, data + from, no_rows, data + to, no_rows);
}

template<class Num_T> inline
void Mat<Num_T>::copy_col(int to, int from)
{
  it_assert_debug(col_in_range(to) && col_in_range(from),
                  "Mat<>::copy_col(): Indexing out of range");
  if (from == to)
    return;

  copy_vector(no_rows, data + from*no_rows, data + to*no_rows);
}

template<class Num_T> inline
void Mat<Num_T>::swap_rows(int r1, int r2)
{
  it_assert_debug(row_in_range(r1) && row_in_range(r2),
                  "Mat<>::swap_rows(): Indexing out of range");
  if (r1 == r2)
    return;

  swap_vector(no_cols, data + r1, no_rows, data + r2, no_rows);
}

template<class Num_T> inline
void Mat<Num_T>::swap_cols(int c1, int c2)
{
  it_assert_debug(col_in_range(c1) && col_in_range(c2),
                  "Mat<>::swap_cols(): Indexing out of range");
  if (c1 == c2)
    return;

  swap_vector(no_rows, data + c1*no_rows, data + c2*no_rows);
}

template<class Num_T>
void Mat<Num_T>::set_submatrix(int r1, int, int c1, int, const Mat<Num_T> &m)
{
  it_warning("Mat<>::set_submatrix(r1, r2, r3, r4, m): This function is "
             "deprecated and might be removed from future IT++ releases. "
             "Please use Mat<>::set_submatrix(r, c, m) function instead.");
  set_submatrix(r1, c1, m);
}

template<class Num_T> inline
void Mat<Num_T>::set_submatrix(int r, int c, const Mat<Num_T> &m)
{
  it_assert_debug((r >= 0) && (r + m.no_rows <= no_rows) &&
                  (c >= 0) && (c + m.no_cols <= no_cols),
                  "Mat<>::set_submatrix(): Indexing out of range "
                  "or wrong input matrix");
  for (int i = 0; i < m.no_cols; i++)
    copy_vector(m.no_rows, m.data + i*m.no_rows, data + (c + i)*no_rows + r);
}



template<class Num_T> inline
void Mat<Num_T>::set_submatrix(int r1, int r2, int c1, int c2, Num_T t)
{
  if (r1 == -1) r1 = no_rows - 1;
  if (r2 == -1) r2 = no_rows - 1;
  if (c1 == -1) c1 = no_cols - 1;
  if (c2 == -1) c2 = no_cols - 1;
  it_assert_debug((r1 >= 0) && (r1 <= r2) && (r2 < no_rows) &&
                  (c1 >= 0) && (c1 <= c2) && (c2 < no_cols),
                  "Mat<>::set_submatrix(): Wrong indexing");
  for (int i = c1; i <= c2; i++) {
    int pos = i * no_rows + r1;
    for (int j = r1; j <= r2; j++)
      data[pos++] = t;
  }
}

template<class Num_T>
void Mat<Num_T>::del_row(int r)
{
  it_assert_debug(row_in_range(r), "Mat<>::del_row(): Index out of range");
  Mat<Num_T> Temp(*this);
  set_size(no_rows - 1, no_cols, false);
  for (int i = 0 ; i < r ; i++) {
    copy_vector(no_cols, &Temp.data[i], no_rows + 1, &data[i], no_rows);
  }
  for (int i = r ; i < no_rows ; i++) {
    copy_vector(no_cols, &Temp.data[i+1], no_rows + 1, &data[i], no_rows);
  }

}

template<class Num_T>
void Mat<Num_T>::del_rows(int r1, int r2)
{
  it_assert_debug((r1 >= 0) && (r1 <= r2) && (r2 < no_rows),
                  "Mat<>::del_rows(): Indexing out of range");
  Mat<Num_T> Temp(*this);
  int no_del_rows = r2 - r1 + 1;
  set_size(no_rows - no_del_rows, no_cols, false);
  for (int i = 0; i < r1 ; ++i) {
    copy_vector(no_cols, &Temp.data[i], Temp.no_rows, &data[i], no_rows);
  }
  for (int i = r2 + 1; i < Temp.no_rows; ++i) {
    copy_vector(no_cols, &Temp.data[i], Temp.no_rows, &data[i-no_del_rows],
                no_rows);
  }
}

template<class Num_T>
void Mat<Num_T>::del_col(int c)
{
  it_assert_debug(col_in_range(c), "Mat<>::del_col(): Index out of range");
  Mat<Num_T> Temp(*this);

  set_size(no_rows, no_cols - 1, false);
  copy_vector(c*no_rows, Temp.data, data);
  copy_vector((no_cols - c)*no_rows, &Temp.data[(c+1)*no_rows], &data[c*no_rows]);
}

template<class Num_T>
void Mat<Num_T>::del_cols(int c1, int c2)
{
  it_assert_debug((c1 >= 0) && (c1 <= c2) && (c2 < no_cols),
                  "Mat<>::del_cols(): Indexing out of range");
  Mat<Num_T> Temp(*this);
  int n_deleted_cols = c2 - c1 + 1;
  set_size(no_rows, no_cols - n_deleted_cols, false);
  copy_vector(c1*no_rows, Temp.data, data);
  copy_vector((no_cols - c1)*no_rows, &Temp.data[(c2+1)*no_rows], &data[c1*no_rows]);
}

template<class Num_T>
void Mat<Num_T>::ins_row(int r, const Vec<Num_T> &v)
{
  it_assert_debug((r >= 0) && (r <= no_rows),
                  "Mat<>::ins_row(): Index out of range");
  it_assert_debug((v.size() == no_cols) || (no_rows == 0),
                  "Mat<>::ins_row(): Wrong size of the input vector");

  if (no_cols == 0) {
    no_cols = v.size();
  }

  Mat<Num_T> Temp(*this);
  set_size(no_rows + 1, no_cols, false);

  for (int i = 0 ; i < r ; i++) {
    copy_vector(no_cols, &Temp.data[i], no_rows - 1, &data[i], no_rows);
  }
  copy_vector(no_cols, v._data(), 1, &data[r], no_rows);
  for (int i = r + 1 ; i < no_rows ; i++) {
    copy_vector(no_cols, &Temp.data[i-1], no_rows - 1, &data[i], no_rows);
  }
}

template<class Num_T>
void Mat<Num_T>::ins_col(int c, const Vec<Num_T> &v)
{
  it_assert_debug((c >= 0) && (c <= no_cols),
                  "Mat<>::ins_col(): Index out of range");
  it_assert_debug((v.size() == no_rows) || (no_cols == 0),
                  "Mat<>::ins_col(): Wrong size of the input vector");

  if (no_rows == 0) {
    no_rows = v.size();
  }

  Mat<Num_T> Temp(*this);
  set_size(no_rows, no_cols + 1, false);

  copy_vector(c*no_rows, Temp.data, data);
  copy_vector(no_rows, v._data(), &data[c*no_rows]);
  copy_vector((no_cols - c - 1)*no_rows, &Temp.data[c*no_rows], &data[(c+1)*no_rows]);
}

template<class Num_T> inline
void Mat<Num_T>::append_row(const Vec<Num_T> &v)
{
  ins_row(no_rows, v);
}

template<class Num_T> inline
void Mat<Num_T>::append_col(const Vec<Num_T> &v)
{
  ins_col(no_cols, v);
}

template<class Num_T>
Mat<Num_T> Mat<Num_T>::transpose() const
{
  Mat<Num_T> temp(no_cols, no_rows);
  for (int i = 0; i < no_rows; ++i) {
    copy_vector(no_cols, &data[i], no_rows, &temp.data[i * no_cols], 1);
  }
  return temp;
}

template<class Num_T>
Mat<Num_T> Mat<Num_T>::hermitian_transpose() const
{
  Mat<Num_T> temp(no_cols, no_rows);
  for (int i = 0; i < no_rows; ++i) {
    copy_vector(no_cols, &data[i], no_rows, &temp.data[i * no_cols], 1);
  }
  return temp;
}

//! \cond
template<>
ITPP_EXPORT cmat cmat::hermitian_transpose() const;
//! \endcond

template<class Num_T>
Mat<Num_T> concat_horizontal(const Mat<Num_T> &m1, const Mat<Num_T> &m2)
{
  // if one of the input matrix is empty just copy the other one as a result
  if (m1.no_cols == 0)
    return m2;
  if (m2.no_cols == 0)
    return m1;
  it_assert_debug(m1.no_rows == m2.no_rows,
                  "Mat<>::concat_horizontal(): Wrong sizes");
  int no_rows = m1.no_rows;
  Mat<Num_T> temp(no_rows, m1.no_cols + m2.no_cols);
  for (int i = 0; i < m1.no_cols; ++i) {
    copy_vector(no_rows, &m1.data[i * no_rows], &temp.data[i * no_rows]);
  }
  for (int i = 0; i < m2.no_cols; ++i) {
    copy_vector(no_rows, &m2.data[i * no_rows], &temp.data[(m1.no_cols + i)
                * no_rows]);
  }
  return temp;
}

template<class Num_T>
Mat<Num_T> concat_vertical(const Mat<Num_T> &m1, const Mat<Num_T> &m2)
{
  // if one of the input matrix is empty just copy the other one as a result
  if (m1.no_rows == 0)
    return m2;
  if (m2.no_rows == 0)
    return m1;
  it_assert_debug(m1.no_cols == m2.no_cols,
                  "Mat<>::concat_vertical(): Wrong sizes");
  int no_cols = m1.no_cols;
  Mat<Num_T> temp(m1.no_rows + m2.no_rows, no_cols);
  for (int i = 0; i < no_cols; ++i) {
    copy_vector(m1.no_rows, &m1.data[i * m1.no_rows],
                &temp.data[i * temp.no_rows]);
    copy_vector(m2.no_rows, &m2.data[i * m2.no_rows],
                &temp.data[i * temp.no_rows + m1.no_rows]);
  }
  return temp;
}

template<class Num_T> inline
Mat<Num_T>& Mat<Num_T>::operator=(Num_T t)
{
  for (int i = 0; i < datasize; i++)
    data[i] = t;
  return *this;
}

template<class Num_T> inline
Mat<Num_T>& Mat<Num_T>::operator=(const Mat<Num_T> &m)
{
  if (this != &m) {
    set_size(m.no_rows, m.no_cols, false);
    if (m.datasize != 0)
      copy_vector(m.datasize, m.data, data);
  }
  return *this;
}

template<class Num_T> inline
Mat<Num_T>& Mat<Num_T>::operator=(const Vec<Num_T> &v)
{
  it_assert_debug(((no_rows == 1) && (no_cols == v.size()))
                  || ((no_cols == 1) && (no_rows == v.size())),
                  "Mat<>::operator=(): Wrong size of the input vector");
  set_size(v.size(), 1, false);
  copy_vector(v.size(), v._data(), data);
  return *this;
}

template<class Num_T> inline
Mat<Num_T>& Mat<Num_T>::operator=(const std::string &str)
{
  set(str);
  return *this;
}

template<class Num_T> inline
Mat<Num_T>& Mat<Num_T>::operator=(const char *str)
{
  set(std::string(str));
  return *this;
}

//-------------------- Templated friend functions --------------------------

template<class Num_T>
Mat<Num_T>& Mat<Num_T>::operator+=(const Mat<Num_T> &m)
{
  if (datasize == 0)
    operator=(m);
  else {
    int i, j, m_pos = 0, pos = 0;
    it_assert_debug(m.no_rows == no_rows && m.no_cols == no_cols, "Mat<Num_T>::operator+=: wrong sizes");
    for (i = 0; i < no_cols; i++) {
      for (j = 0; j < no_rows; j++)
        data[pos+j] += m.data[m_pos+j];
      pos += no_rows;
      m_pos += m.no_rows;
    }
  }
  return *this;
}

template<class Num_T> inline
Mat<Num_T>& Mat<Num_T>::operator+=(Num_T t)
{
  for (int i = 0; i < datasize; i++)
    data[i] += t;
  return *this;
}

template<class Num_T>
Mat<Num_T> operator+(const Mat<Num_T> &m1, const Mat<Num_T> &m2)
{
  Mat<Num_T> r(m1.no_rows, m1.no_cols);
  int i, j, m1_pos = 0, m2_pos = 0, r_pos = 0;

  it_assert_debug((m1.no_rows == m2.no_rows) && (m1.no_cols == m2.no_cols),
                  "Mat<>::operator+(): Wrong sizes");

  for (i = 0; i < r.no_cols; i++) {
    for (j = 0; j < r.no_rows; j++)
      r.data[r_pos+j] = m1.data[m1_pos+j] + m2.data[m2_pos+j];
    // next column
    m1_pos += m1.no_rows;
    m2_pos += m2.no_rows;
    r_pos += r.no_rows;
  }

  return r;
}


template<class Num_T>
Mat<Num_T> operator+(const Mat<Num_T> &m, Num_T t)
{
  Mat<Num_T> r(m.no_rows, m.no_cols);

  for (int i = 0; i < r.datasize; i++)
    r.data[i] = m.data[i] + t;

  return r;
}

template<class Num_T>
Mat<Num_T> operator+(Num_T t, const Mat<Num_T> &m)
{
  Mat<Num_T> r(m.no_rows, m.no_cols);

  for (int i = 0; i < r.datasize; i++)
    r.data[i] = t + m.data[i];

  return r;
}

template<class Num_T>
Mat<Num_T>& Mat<Num_T>::operator-=(const Mat<Num_T> &m)
{
  int i, j, m_pos = 0, pos = 0;

  if (datasize == 0) {
    set_size(m.no_rows, m.no_cols, false);
    for (i = 0; i < no_cols; i++) {
      for (j = 0; j < no_rows; j++)
        data[pos+j] = -m.data[m_pos+j];
      // next column
      m_pos += m.no_rows;
      pos += no_rows;
    }
  }
  else {
    it_assert_debug((m.no_rows == no_rows) && (m.no_cols == no_cols),
                    "Mat<>::operator-=(): Wrong sizes");
    for (i = 0; i < no_cols; i++) {
      for (j = 0; j < no_rows; j++)
        data[pos+j] -= m.data[m_pos+j];
      // next column
      m_pos += m.no_rows;
      pos += no_rows;
    }
  }
  return *this;
}

template<class Num_T>
Mat<Num_T> operator-(const Mat<Num_T> &m1, const Mat<Num_T> &m2)
{
  Mat<Num_T> r(m1.no_rows, m1.no_cols);
  int i, j, m1_pos = 0, m2_pos = 0, r_pos = 0;
  it_assert_debug((m1.no_rows == m2.no_rows) && (m1.no_cols == m2.no_cols),
                  "Mat<>::operator-(): Wrong sizes");

  for (i = 0; i < r.no_cols; i++) {
    for (j = 0; j < r.no_rows; j++)
      r.data[r_pos+j] = m1.data[m1_pos+j] - m2.data[m2_pos+j];
    // next column
    m1_pos += m1.no_rows;
    m2_pos += m2.no_rows;
    r_pos += r.no_rows;
  }

  return r;
}

template<class Num_T> inline
Mat<Num_T>& Mat<Num_T>::operator-=(Num_T t)
{
  for (int i = 0; i < datasize; i++)
    data[i] -= t;
  return *this;
}

template<class Num_T>
Mat<Num_T> operator-(const Mat<Num_T> &m, Num_T t)
{
  Mat<Num_T> r(m.no_rows, m.no_cols);
  int i, j, m_pos = 0, r_pos = 0;

  for (i = 0; i < r.no_cols; i++) {
    for (j = 0; j < r.no_rows; j++)
      r.data[r_pos+j] = m.data[m_pos+j] - t;
    // next column
    m_pos += m.no_rows;
    r_pos += r.no_rows;
  }

  return r;
}

template<class Num_T>
Mat<Num_T> operator-(Num_T t, const Mat<Num_T> &m)
{
  Mat<Num_T> r(m.no_rows, m.no_cols);
  int i, j, m_pos = 0, r_pos = 0;

  for (i = 0; i < r.no_cols; i++) {
    for (j = 0; j < r.no_rows; j++)
      r.data[r_pos+j] = t - m.data[m_pos+j];
    // next column
    m_pos += m.no_rows;
    r_pos += r.no_rows;
  }

  return r;
}

template<class Num_T>
Mat<Num_T> operator-(const Mat<Num_T> &m)
{
  Mat<Num_T> r(m.no_rows, m.no_cols);
  int i, j, m_pos = 0, r_pos = 0;

  for (i = 0; i < r.no_cols; i++) {
    for (j = 0; j < r.no_rows; j++)
      r.data[r_pos+j] = -m.data[m_pos+j];
    // next column
    m_pos += m.no_rows;
    r_pos += r.no_rows;
  }

  return r;
}

template<class Num_T>
Mat<Num_T>& Mat<Num_T>::operator*=(const Mat<Num_T> &m)
{
  it_assert_debug(no_cols == m.no_rows, "Mat<>::operator*=(): Wrong sizes");
  Mat<Num_T> r(no_rows, m.no_cols);

  Num_T tmp;

  int i, j, k, r_pos = 0, pos = 0, m_pos = 0;

  for (i = 0; i < r.no_cols; i++) {
    for (j = 0; j < r.no_rows; j++) {
      tmp = Num_T(0);
      pos = 0;
      for (k = 0; k < no_cols; k++) {
        tmp += data[pos+j] * m.data[m_pos+k];
        pos += no_rows;
      }
      r.data[r_pos+j] = tmp;
    }
    r_pos += r.no_rows;
    m_pos += m.no_rows;
  }
  operator=(r); // time consuming
  return *this;
}

//! \cond
template<> ITPP_EXPORT mat& mat::operator*=(const mat &m);
template<> ITPP_EXPORT cmat& cmat::operator*=(const cmat &m);
//! \endcond

template<class Num_T> inline
Mat<Num_T>& Mat<Num_T>::operator*=(Num_T t)
{
  scal_vector(datasize, t, data);
  return *this;
}

//! Multiplication of two matrices
template<class Num_T>
Mat<Num_T> operator*(const Mat<Num_T> &m1, const Mat<Num_T> &m2)
{
  it_assert_debug(m1.cols() == m2.rows(),
                  "Mat<>::operator*(): Wrong sizes");
  Mat<Num_T> r(m1.rows(), m2.cols());

  Num_T tmp;
  int i, j, k;
  Num_T *tr = r._data();
  const Num_T *t1; const Num_T *t2 = m2._data();

  for (i = 0; i < r.cols(); i++) {
    for (j = 0; j < r.rows(); j++) {
      tmp = Num_T(0);
      t1 = m1._data() + j;
      for (k = m1.cols(); k > 0; k--) {
        tmp += *(t1) * *(t2++);
	    t1 += m1.rows();
      }
      *(tr++) = tmp;
	  t2 -= m2.rows();
    }
    t2 += m2.rows();
  }

  return r;
}

//! \cond
template<> ITPP_EXPORT mat operator*(const mat &m1, const mat &m2);
template<> ITPP_EXPORT cmat operator*(const cmat &m1, const cmat &m2);
//! \endcond

//! Multiplication of matrix \c m and vector \c v (column vector)
template<class Num_T>
Vec<Num_T> operator*(const Mat<Num_T> &m, const Vec<Num_T> &v)
{
  it_assert_debug(m.cols() == v.size(),
                  "Mat<>::operator*(): Wrong sizes");
  Vec<Num_T> r(m.rows());
  int i, k, m_pos;

  for (i = 0; i < m.rows(); i++) {
    r(i) = Num_T(0);
    m_pos = 0;
    for (k = 0; k < m.cols(); k++) {
      r(i) += m._data()[m_pos+i] * v(k);
      m_pos += m.rows();
    }
  }

  return r;
}

//! \cond
template<> ITPP_EXPORT vec operator*(const mat &m, const vec &v);
template<> ITPP_EXPORT cvec operator*(const cmat &m, const cvec &v);
//! \endcond

//! Multiplication of matrix and scalar
template<class Num_T>
Mat<Num_T> operator*(const Mat<Num_T> &m, Num_T t)
{
  Mat<Num_T> r(m.rows(), m.cols());

  const Num_T* m_data = m._data();
  Num_T* r_data = r._data();
  for (int i = 0; i < r._datasize(); i++)
    r_data[i] = m_data[i] * t;

  return r;
}

//! Multiplication of scalar and matrix
template<class Num_T> inline
Mat<Num_T> operator*(Num_T t, const Mat<Num_T> &m)
{
  return operator*(m, t);
}

template<class Num_T> inline
Mat<Num_T> elem_mult(const Mat<Num_T> &m1, const Mat<Num_T> &m2)
{
  Mat<Num_T> out;
  elem_mult_out(m1, m2, out);
  return out;
}

template<class Num_T>
void elem_mult_out(const Mat<Num_T> &m1, const Mat<Num_T> &m2,
                   Mat<Num_T> &out)
{
  it_assert_debug((m1.no_rows == m2.no_rows) && (m1.no_cols == m2.no_cols),
                  "Mat<>::elem_mult_out(): Wrong sizes");
  out.set_size(m1.no_rows, m1.no_cols);
  for (int i = 0; i < out.datasize; i++)
    out.data[i] = m1.data[i] * m2.data[i];
}

template<class Num_T>
void elem_mult_out(const Mat<Num_T> &m1, const Mat<Num_T> &m2,
                   const Mat<Num_T> &m3, Mat<Num_T> &out)
{
  it_assert_debug((m1.no_rows == m2.no_rows) && (m1.no_rows == m3.no_rows)
                  && (m1.no_cols == m2.no_cols) && (m1.no_cols == m3.no_cols),
                  "Mat<>::elem_mult_out(): Wrong sizes");
  out.set_size(m1.no_rows, m1.no_cols);
  for (int i = 0; i < out.datasize; i++)
    out.data[i] = m1.data[i] * m2.data[i] * m3.data[i];
}

template<class Num_T>
void elem_mult_out(const Mat<Num_T> &m1, const Mat<Num_T> &m2,
                   const Mat<Num_T> &m3, const Mat<Num_T> &m4,
                   Mat<Num_T> &out)
{
  it_assert_debug((m1.no_rows == m2.no_rows) && (m1.no_rows == m3.no_rows)
                  && (m1.no_rows == m4.no_rows) && (m1.no_cols == m2.no_cols)
                  && (m1.no_cols == m3.no_cols) && (m1.no_cols == m4.no_cols),
                  "Mat<>::elem_mult_out(): Wrong sizes");
  out.set_size(m1.no_rows, m1.no_cols);
  for (int i = 0; i < out.datasize; i++)
    out.data[i] = m1.data[i] * m2.data[i] * m3.data[i] * m4.data[i];
}

template<class Num_T>
#ifndef _MSC_VER
inline
#endif
void elem_mult_inplace(const Mat<Num_T> &m1, Mat<Num_T> &m2)
{
  it_assert_debug((m1.no_rows == m2.no_rows) && (m1.no_cols == m2.no_cols),
                  "Mat<>::elem_mult_inplace(): Wrong sizes");
  for (int i = 0; i < m2.datasize; i++)
    m2.data[i] *= m1.data[i];
}

template<class Num_T> inline
Num_T elem_mult_sum(const Mat<Num_T> &m1, const Mat<Num_T> &m2)
{
  it_assert_debug((m1.no_rows == m2.no_rows) && (m1.no_cols == m2.no_cols),
                  "Mat<>::elem_mult_sum(): Wrong sizes");
  Num_T acc = 0;

  for (int i = 0; i < m1.datasize; i++)
    acc += m1.data[i] * m2.data[i];

  return acc;
}

template<class Num_T> inline
Mat<Num_T>& Mat<Num_T>::operator/=(Num_T t)
{
  for (int i = 0; i < datasize; i++)
    data[i] /= t;
  return *this;
}

template<class Num_T> inline
Mat<Num_T>& Mat<Num_T>::operator/=(const Mat<Num_T> &m)
{
  it_assert_debug((m.no_rows == no_rows) && (m.no_cols == no_cols),
                  "Mat<>::operator/=(): Wrong sizes");
  for (int i = 0; i < datasize; i++)
    data[i] /= m.data[i];
  return *this;
}

template<class Num_T>
Mat<Num_T> operator/(const Mat<Num_T> &m, Num_T t)
{
  Mat<Num_T> r(m.no_rows, m.no_cols);
  for (int i = 0; i < r.datasize; ++i)
    r.data[i] = m.data[i] / t;
  return r;
}

template<class Num_T>
Mat<Num_T> operator/(Num_T t, const Mat<Num_T> &m)
{
  Mat<Num_T> r(m.no_rows, m.no_cols);
  for (int i = 0; i < r.datasize; ++i)
    r.data[i] = t / m.data[i];
  return r;
}

template<class Num_T> inline
Mat<Num_T> elem_div(const Mat<Num_T> &m1, const Mat<Num_T> &m2)
{
  Mat<Num_T> out;
  elem_div_out(m1, m2, out);
  return out;
}

template<class Num_T>
void elem_div_out(const Mat<Num_T> &m1, const Mat<Num_T> &m2,
                  Mat<Num_T> &out)
{
  it_assert_debug((m1.no_rows == m2.no_rows) && (m1.no_cols == m2.no_cols),
                  "Mat<>::elem_div_out(): Wrong sizes");

  if ((out.no_rows != m1.no_rows) || (out.no_cols != m1.no_cols))
    out.set_size(m1.no_rows, m1.no_cols);

  for (int i = 0; i < out.datasize; i++)
    out.data[i] = m1.data[i] / m2.data[i];
}

template<class Num_T> inline
Num_T elem_div_sum(const Mat<Num_T> &m1, const Mat<Num_T> &m2)
{
  it_assert_debug((m1.no_rows == m2.no_rows) && (m1.no_cols == m2.no_cols),
                  "Mat<>::elem_div_sum(): Wrong sizes");
  Num_T acc = 0;

  for (int i = 0; i < m1.datasize; i++)
    acc += m1.data[i] / m2.data[i];

  return acc;
}

template<class Num_T>
bool Mat<Num_T>::operator==(const Mat<Num_T> &m) const
{
  if (no_rows != m.no_rows || no_cols != m.no_cols) return false;
  for (int i = 0;i < datasize;i++) {
    if (data[i] != m.data[i]) return false;
  }
  return true;
}

template<class Num_T>
bool Mat<Num_T>::operator!=(const Mat<Num_T> &m) const
{
  if (no_rows != m.no_rows || no_cols != m.no_cols) return true;
  for (int i = 0;i < datasize;i++) {
    if (data[i] != m.data[i]) return true;
  }
  return false;
}

template <class Num_T>
std::ostream &operator<<(std::ostream &os, const Mat<Num_T> &m)
{
  int i;

  switch (m.rows()) {
  case 0 :
    os << "[]";
    break;
  case 1 :
    os << '[' << m.get_row(0) << ']';
    break;
  default:
    os << '[' << m.get_row(0) << std::endl;
    for (i = 1; i < m.rows() - 1; i++)
      os << ' ' << m.get_row(i) << std::endl;
    os << ' ' << m.get_row(m.rows() - 1) << ']';
  }

  return os;
}

template <class Num_T>
std::istream &operator>>(std::istream &is, Mat<Num_T> &m)
{
  std::ostringstream buffer;
  bool started = false;
  bool finished = false;
  bool brackets = false;
  bool within_double_brackets = false;
  char c;

  while (!finished) {
    if (is.eof()) {
      finished = true;
    }
    else {
      is.get(c);

      if (is.eof() || (c == '\n')) {
        if (brackets) {
          // Right bracket missing
          is.setstate(std::ios_base::failbit);
          finished = true;
        }
        else if (!((c == '\n') && !started)) {
          finished = true;
        }
      }
      else if ((c == ' ') || (c == '\t')) {
        if (started) {
          buffer << ' ';
        }
      }
      else if (c == '[') {
        if ((started && !brackets) || within_double_brackets) {
          // Unexpected left bracket
          is.setstate(std::ios_base::failbit);
          finished = true;
        }
        else if (!started) {
          started = true;
          brackets = true;
        }
        else {
          within_double_brackets = true;
        }
      }
      else if (c == ']') {
        if (!started || !brackets) {
          // Unexpected right bracket
          is.setstate(std::ios_base::failbit);
          finished = true;
        }
        else if (within_double_brackets) {
          within_double_brackets = false;
          buffer << ';';
        }
        else {
          finished = true;
        }
        while (!is.eof() && (((c = static_cast<char>(is.peek())) == ' ')
                             || (c == '\t'))) {
          is.get();
        }
        if (!is.eof() && (c == '\n')) {
          is.get();
        }
      }
      else {
        started = true;
        buffer << c;
      }
    }
  }

  if (!started) {
    m.set_size(0, false);
  }
  else {
    m.set(buffer.str());
  }

  return is;
}

//! \cond

// ---------------------------------------------------------------------
// Instantiations
// ---------------------------------------------------------------------

// class instantiations

ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Mat<double>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Mat<std::complex<double> >;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Mat<int>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Mat<short int>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Mat<bin>;

// addition operators

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat operator+(const mat &m1, const mat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat operator+(const cmat &m1, const cmat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat operator+(const imat &m1, const imat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat operator+(const smat &m1, const smat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat operator+(const bmat &m1, const bmat &m2);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat operator+(const mat &m, double t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat operator+(const cmat &m, std::complex<double> t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat operator+(const imat &m, int t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat operator+(const smat &m, short t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat operator+(const bmat &m, bin t);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat operator+(double t, const mat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat operator+(std::complex<double> t, const cmat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat operator+(int t, const imat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat operator+(short t, const smat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat operator+(bin t, const bmat &m);

// subtraction operators

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat operator-(const mat &m1, const mat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat operator-(const cmat &m1, const cmat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat operator-(const imat &m1, const imat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat operator-(const smat &m1, const smat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat operator-(const bmat &m1, const bmat &m2);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat operator-(const mat &m, double t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat operator-(const cmat &m, std::complex<double> t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat operator-(const imat &m, int t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat operator-(const smat &m, short t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat operator-(const bmat &m, bin t);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat operator-(double t, const mat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat operator-(std::complex<double> t, const cmat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat operator-(int t, const imat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat operator-(short t, const smat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat operator-(bin t, const bmat &m);

// unary minus

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat operator-(const mat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat operator-(const cmat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat operator-(const imat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat operator-(const smat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat operator-(const bmat &m);

// multiplication operators

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat operator*(const imat &m1, const imat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat operator*(const smat &m1, const smat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat operator*(const bmat &m1, const bmat &m2);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  ivec operator*(const imat &m, const ivec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  svec operator*(const smat &m, const svec &v);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bvec operator*(const bmat &m, const bvec &v);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat operator*(const mat &m, double t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat operator*(const cmat &m, std::complex<double> t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat operator*(const imat &m, int t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat operator*(const smat &m, short t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat operator*(const bmat &m, bin t);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat operator*(double t, const mat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat operator*(std::complex<double> t, const cmat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat operator*(int t, const imat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat operator*(short t, const smat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat operator*(bin t, const bmat &m);

// element-wise multiplication

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat elem_mult(const mat &m1, const mat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat elem_mult(const cmat &m1, const cmat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat elem_mult(const imat &m1, const imat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat elem_mult(const smat &m1, const smat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat elem_mult(const bmat &m1, const bmat &m2);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const mat &m1, const mat &m2, mat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const cmat &m1, const cmat &m2,
                                     cmat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const imat &m1, const imat &m2,
                                     imat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const smat &m1, const smat &m2,
                                     smat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const bmat &m1, const bmat &m2,
                                     bmat &out);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const mat &m1, const mat &m2,
                                     const mat &m3, mat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const cmat &m1, const cmat &m2,
                                     const cmat &m3, cmat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const imat &m1, const imat &m2,
                                     const imat &m3, imat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const smat &m1, const smat &m2,
                                     const smat &m3, smat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const bmat &m1, const bmat &m2,
                                     const bmat &m3, bmat &out);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const mat &m1, const mat &m2,
                                     const mat &m3, const mat &m4, mat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const cmat &m1, const cmat &m2,
                                     const cmat &m3, const cmat &m4,
                                     cmat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const imat &m1, const imat &m2,
                                     const imat &m3, const imat &m4,
                                     imat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const smat &m1, const smat &m2,
                                     const smat &m3, const smat &m4,
                                     smat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_out(const bmat &m1, const bmat &m2,
                                     const bmat &m3, const bmat &m4,
                                     bmat &out);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_inplace(const mat &m1, mat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_inplace(const cmat &m1, cmat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_inplace(const imat &m1, imat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_inplace(const smat &m1, smat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_mult_inplace(const bmat &m1, bmat &m2);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  double elem_mult_sum(const mat &m1, const mat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  std::complex<double> elem_mult_sum(const cmat &m1,
    const cmat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  int elem_mult_sum(const imat &m1, const imat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  short elem_mult_sum(const smat &m1, const smat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bin elem_mult_sum(const bmat &m1, const bmat &m2);

// division operator

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat operator/(double t, const mat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat operator/(std::complex<double> t, const cmat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat operator/(int t, const imat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat operator/(short t, const smat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat operator/(bin t, const bmat &m);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat operator/(const mat &m, double t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat operator/(const cmat &m, std::complex<double> t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat operator/(const imat &m, int t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat operator/(const smat &m, short t);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat operator/(const bmat &m, bin t);

// element-wise division

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat elem_div(const mat &m1, const mat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat elem_div(const cmat &m1, const cmat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat elem_div(const imat &m1, const imat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat elem_div(const smat &m1, const smat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat elem_div(const bmat &m1, const bmat &m2);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_div_out(const mat &m1, const mat &m2, mat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_div_out(const cmat &m1, const cmat &m2, cmat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_div_out(const imat &m1, const imat &m2, imat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_div_out(const smat &m1, const smat &m2, smat &out);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  void elem_div_out(const bmat &m1, const bmat &m2, bmat &out);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  double elem_div_sum(const mat &m1, const mat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  std::complex<double> elem_div_sum(const cmat &m1,
    const cmat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  int elem_div_sum(const imat &m1, const imat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  short elem_div_sum(const smat &m1, const smat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bin elem_div_sum(const bmat &m1, const bmat &m2);

// concatenation

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat concat_horizontal(const mat &m1, const mat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat concat_horizontal(const cmat &m1, const cmat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat concat_horizontal(const imat &m1, const imat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat concat_horizontal(const smat &m1, const smat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat concat_horizontal(const bmat &m1, const bmat &m2);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  mat concat_vertical(const mat &m1, const mat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  cmat concat_vertical(const cmat &m1, const cmat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  imat concat_vertical(const imat &m1, const imat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  smat concat_vertical(const smat &m1, const smat &m2);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  bmat concat_vertical(const bmat &m1, const bmat &m2);

// I/O streams

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  std::ostream &operator<<(std::ostream &os, const mat  &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  std::ostream &operator<<(std::ostream &os, const cmat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  std::ostream &operator<<(std::ostream &os, const imat  &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  std::ostream &operator<<(std::ostream &os, const smat  &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  std::ostream &operator<<(std::ostream &os, const bmat  &m);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  std::istream &operator>>(std::istream &is, mat  &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  std::istream &operator>>(std::istream &is, cmat &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  std::istream &operator>>(std::istream &is, imat  &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  std::istream &operator>>(std::istream &is, smat  &m);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT  std::istream &operator>>(std::istream &is, bmat  &m);

//! \endcond

} // namespace itpp

#endif // #ifndef MAT_H
