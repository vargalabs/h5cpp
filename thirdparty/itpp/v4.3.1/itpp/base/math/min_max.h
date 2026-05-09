/*!
 * \file
 * \brief Minimum and maximum functions on vectors and matrices
 * \author Tony Ottosson, Johan Bergman and Adam Piatyszek
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

#ifndef MIN_MAX_H
#define MIN_MAX_H

#include <itpp/base/mat.h>


namespace itpp
{

/*!
 * \addtogroup miscfunc
 * @{
 */

//! Maximum value of vector
template<class T>
T max(const Vec<T> &v)
{
  T maxdata = v(0);
  for (int i = 1; i < v.length(); i++)
    if (v(i) > maxdata)
      maxdata = v(i);
  return maxdata;
}

//! Maximum value of vector, also returns the index position of max value
template<class T>
T max(const Vec<T> &v, int& index)
{
  T maxdata = v(0);
  index = 0;
  for (int i = 1; i < v.length(); i++)
    if (v(i) > maxdata) {
      maxdata = v(i);
      index = i;
    }
  return maxdata;
}

/*!
 * Maximum values over each row/column in the matrix \c m
 *
 * <tt>max(m) = max(m, 1)</tt> returns a vector where the elements are
 * maximum over each column, whereas <tt>max(m, 2)</tt> returns a vector
 * where the elements are maximum over each row.
 */
template<class T>
Vec<T> max(const Mat<T> &m, int dim = 1)
{
  it_assert((dim == 1) || (dim == 2), "max(): dimension need to be 1 or 2");
  Vec<T> out;
  if (dim == 1) {
    out.set_size(m.cols(), false);
    for (int i = 0; i < m.cols(); i++)
      out(i) = max(m.get_col(i));
  }
  else {
    out.set_size(m.rows(), false);
    for (int i = 0; i < m.rows(); i++)
      out(i) = max(m.get_row(i));
  }
  return out;
}

/*!
 * Maximum values over each row/column in the matrix \c m
 *
 * <tt>max(m) = max(m, 1)</tt> returns a vector where the elements are
 * maximum over each column, whereas <tt>max(m, 2)</tt> returns a vector
 * where the elements are maximum over each row.
 *
 * Also returns a vector of indices with positions of maximum value within
 * a column/row.
 */
template<class T>
Vec<T> max(const Mat<T> &m, ivec &index, int dim = 1)
{
  it_assert((dim == 1) || (dim == 2), "max(): dimension need to be 1 or 2");
  Vec<T> out;
  if (dim == 1) {
    out.set_size(m.cols(), false);
    index.set_size(m.cols(), false);
    for (int i = 0; i < m.cols(); i++)
      out(i) = max(m.get_col(i), index(i));
  }
  else {
    out.set_size(m.rows(), false);
    index.set_size(m.rows(), false);
    for (int i = 0; i < m.rows(); i++)
      out(i) = max(m.get_row(i), index(i));
  }
  return out;
}

//! Minimum value of vector
template<class T>
T min(const Vec<T> &in)
{
  T mindata = in[0];
  for (int i = 1; i < in.length(); i++)
    if (in[i] < mindata)
      mindata = in[i];
  return mindata;
}

//! Minimum value of vector, also returns the index position of min value
template<class T>
T min(const Vec<T> &in, int& index)
{
  T mindata = in[0];
  index = 0;
  for (int i = 1; i < in.length(); i++)
    if (in[i] < mindata) {
      mindata = in[i];
      index = i;
    }
  return mindata;
}


/*!
 * Minimum values over each row/column in the matrix \c m
 *
 * <tt>min(m) = min(m, 1)</tt> returns a vector where the elements are
 * minimum over each column, whereas <tt>min(m, 2)</tt> returns a vector
 * where the elements are minimum over each row.
 */
template<class T>
Vec<T> min(const Mat<T> &m, int dim = 1)
{
  it_assert((dim == 1) || (dim == 2), "min(): dimension need to be 1 or 2");
  Vec<T> out;
  if (dim == 1) {
    out.set_size(m.cols(), false);
    for (int i = 0; i < m.cols(); i++)
      out(i) = min(m.get_col(i));
  }
  else {
    out.set_size(m.rows(), false);
    for (int i = 0; i < m.rows(); i++)
      out(i) = min(m.get_row(i));
  }
  return out;
}


/*!
 * Minimum values over each row/column in the matrix \c m
 *
 * <tt>min(m) = min(m, 1)</tt> returns a vector where the elements are
 * minimum over each column, whereas <tt>min(m, 2)</tt> returns a vector
 * where the elements are minimum over each row.
 *
 * Also returns a vector of indices with positions of minimum value within
 * a column/row.
 */
template<class T>
Vec<T> min(const Mat<T> &m,  ivec &index, int dim = 1)
{
  it_assert((dim == 1) || (dim == 2), "min(): dimension need to be 1 or 2");
  Vec<T> out;
  if (dim == 1) {
    out.set_size(m.cols(), false);
    index.set_size(m.cols(), false);
    for (int i = 0; i < m.cols(); i++)
      out(i) = min(m.get_col(i), index(i));
  }
  else {
    out.set_size(m.rows(), false);
    index.set_size(m.rows(), false);
    for (int i = 0; i < m.rows(); i++)
      out(i) = min(m.get_row(i), index(i));
  }
  return out;
}


//! Return the postion of the maximum element in the vector
template<class T>
int max_index(const Vec<T> &in)
{
  int maxindex = 0;
  for (int i = 1; i < in.length(); i++)
    if (in[i] > in[maxindex])
      maxindex = i;
  return maxindex;
}

//! Return the postion of the maximum element in the matrix
template<class T>
void max_index(const Mat<T> &m, int &row, int &col)
{
  T maxdata = m(0, 0);
  row = col = 0;
  for (int i = 0; i < m.rows(); i++)
    for (int j = 0; j < m.cols(); j++)
      if (m(i, j) > maxdata) {
        row = i;
        col = j;
        maxdata = m(i, j);
      }
}

//! Return the postion of the minimum element in the vector
template<class T>
int min_index(const Vec<T> &in)
{
  int minindex = 0;
  for (int i = 1; i < in.length(); i++)
    if (in[i] < in[minindex])
      minindex = i;
  return minindex;
}

//! Return the postion of the minimum element in the matrix
template<class T>
void min_index(const Mat<T> &m, int &row, int &col)
{
  T mindata = m(0, 0);
  row = col = 0;
  for (int i = 0; i < m.rows(); i++)
    for (int j = 0; j < m.cols(); j++)
      if (m(i, j) < mindata) {
        row = i;
        col = j;
        mindata = m(i, j);
      }
}

/*!
 * @}
 */

} //namespace itpp


#endif /* MIN_MAX_H */
