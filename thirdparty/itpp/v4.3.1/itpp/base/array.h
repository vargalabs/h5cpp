/*!
 * \file
 * \brief Definition of Array class (container)
 * \author Tobias Ringstrom and Adam Piatyszek
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

#ifndef ARRAY_H
#define ARRAY_H

#include <itpp/base/itassert.h>
#include <itpp/base/math/misc.h>
#include <itpp/base/factory.h>
#include <itpp/base/copy_vector.h>


namespace itpp
{

// Forward declarations
template<class T> class Array;
//! Append element \c e to the end of the Array \c a
template<class T> const Array<T> concat(const Array<T> &a, const T &e);
//! Append element \c e to the beginning of the Array \c a
template<class T> const Array<T> concat(const T &e, const Array<T> &a);
//! Concat Arrays \c a1 and \c a2
template<class T> const Array<T> concat(const Array<T> &a1,
                                        const Array<T> &a2);
//! Concat Arrays \c a1, \c a2 and \c a3
template<class T> const Array<T> concat(const Array<T> &a1,
                                        const Array<T> &a2,
                                        const Array<T> &a3);

/*!
  \ingroup arr_vec_mat
  \brief General array class
  \author Tobias Ringstrom and Adam Piatyszek

  This class is a general linear array class for arbitrary types. The
  operations and functions are the same as for the vector \c Vec class
  (except for the arithmetics).

  For rarely used types you will need to instantiate the class by
  \code
  template class Array<type>;
  \endcode

  The following example shows how to define an Array of vectors:
  \code
  vec a = randn(10);
  vec b = randn(20);
  vec c = randn(30);
  Array<vec> my_array(3);
  my_array(0) = a;
  my_array(1) = b;
  my_array(2) = c;
  \endcode

  For types T with istream \c operator>> defined special constructor or
  \c operator= or \c set_array functions (see Related Functions) can be
  used to assign a string literal to an Array. The string literal has the
  same format that is used by the istream/ostream operators:

  \code
  // Initialise an array with three bit vectors
  Array<bvec> B = "{[1 0 1] [0 0 1] [1 0 0 0 1]}";

  // Declare an Array of Arrays of vectors
  Array<Array<ivec> > an_array;

  // Assign with an Array containing 2 Arrays,
  // the first Array containing [1 2] and
  // the second Array containing [3 4 5] and [6 7]
  set_array(an_array, "{{[1 2]} {[3 4 5] [6 7]}}");
  \endcode

  By default, Array elements are created using the default constructor
  for the element type. This can be changed by specifying a suitable
  Factory in the Array constructor call (see Detailed Description for
  Factory).
*/
template<class T>
class Array
{
public:
  //! Default constructor. An element factory \c f can be specified.
  explicit Array(const Factory &f = DEFAULT_FACTORY);
  //! Create an Array of size \c n. An element factory \c f can be specified.
  Array(int n, const Factory &f = DEFAULT_FACTORY);
  //! Copy constructor. An element factory \c f can be specified.
  Array(const Array<T> &a, const Factory &f = DEFAULT_FACTORY);
  //! Create an Array from string. An element factory \c f can be specified.
  Array(const std::string& values, const Factory &f = DEFAULT_FACTORY);
  //! Create an Array from char*. An element factory \c f can be specified.
  Array(const char* values, const Factory &f = DEFAULT_FACTORY);

  //! Destructor
  virtual ~Array();

  //! Get the \c i element
  T &operator()(int i);
  //! Get the \c i element
  const T &operator()(int i) const;
  //! Sub-array from element \c i1 to element \c i2
  const Array<T> operator()(int i1, int i2) const;
  //! Sub-array with the elements given by the integer Array
  const Array<T> operator()(const Array<int> &indices) const;

  //! Get \c n left elements of the array
  Array<T> left(int n) const;
  //! Get \c n right elements of the array
  Array<T> right(int n) const;
  //! Get \c n elements of the array starting from \c pos
  Array<T> mid(int pos, int n) const;

  //! Assignment operator
  Array<T>& operator=(const T &e);
  //! Assignment operator
  Array<T>& operator=(const Array<T> &a);
  //! Assignment operator
  Array<T>& operator=(const char* values);

  //! Append element \c e to the end of the Array \c a
  friend const Array<T> concat <>(const Array<T> &a1, const T &e);
  //! Concat element \c e to the beginning of the Array \c a
  friend const Array<T> concat <>(const T &e, const Array<T> &a);
  //! Concat Arrays \c a1 and \c a2
  friend const Array<T> concat <>(const Array<T> &a1, const Array<T> &a2);
  //! Concat Arrays \c a1, \c a2 and \c a3
  friend const Array<T> concat <>(const Array<T> &a1, const Array<T> &a2,
                                  const Array<T> &a3);

  //! Returns the number of data elements in the array object
  int size() const { return ndata; }
  //! Returns the number of data elements in the array object
  int length() const { return ndata; }
  //! Resizing an Array<T>.
  void set_size(int n, bool copy = false);
  //! Resizing an Array<T>.
  void set_length(int n, bool copy = false) { set_size(n, copy); }

  //! Shift in data at position 0. Return data from the last position.
  T shift_right(const T& e);
  //! Shift in array at position 0. Return data from the last position.
  const Array<T> shift_right(const Array<T> &a);
  //! Shift in data at the last position. Return data from position 0.
  T shift_left(const T& e);
  //! Shift in array at the last position. Return data from position 0.
  const Array<T> shift_left(const Array<T> &a);
  //! Swap elements i and j.
  void swap(int i, int j);

  //! Set the subarray defined by indicies i1 to i2 to Array<T> a.
  void set_subarray(int i1, int i2, const Array<T> &a);
  //! Set the subarray defined by indicies i1 to i2 the element value t.
  void set_subarray(int i1, int i2, const T &t);

protected:
  //! Allocate storage for an array of length \c n
  void alloc(int n);
  //! Free the storage space allocated by the array
  void free();
  //! Check whether index \c i is in the allowed range
  bool in_range(int i) const { return ((i < ndata) && (i >= 0)); }
  //! The current number of elements in the Array
  int ndata;
  //! A pointer to the data area
  T *data;
  //! Element factory (by default set to DEFAULT_FACTORY)
  const Factory &factory;
};

// -------------------- Implementation starts here --------------------

template<class T> inline
void Array<T>::alloc(int n)
{
  if (n > 0) {
    create_elements(data, n, factory);
    ndata = n;
  }
  else {
    data = 0;
    ndata = 0;
  }
}

template<class T> inline
void Array<T>::free()
{
  destroy_elements(data, ndata);
  ndata = 0;
}

template<class T> inline
Array<T>::Array(const Factory &f) : ndata(0), data(0), factory(f) {}

template<class T> inline
Array<T>::Array(const int n, const Factory &f) : ndata(0), data(0), factory(f)
{
  alloc(n);
}

template<class T> inline
Array<T>::Array(const Array<T> &a, const Factory &f)
    : ndata(0), data(0), factory(f)
{
  alloc(a.ndata);
  for (int i = 0; i < a.ndata; i++)
    data[i] = a.data[i];
}

template<class T> inline
Array<T>::Array(const std::string& values, const Factory &f)
    : ndata(0), data(0), factory(f)
{
  std::istringstream buffer(values);
  buffer >> *this;
}

template<class T> inline
Array<T>::Array(const char* values, const Factory &f)
    : ndata(0), data(0), factory(f)
{
  std::istringstream buffer(values);
  buffer >> *this;
}

template<class T>
Array<T>::~Array()
{
  free();
}

template<class T>
void Array<T>::set_size(int size, bool copy)
{
  it_assert_debug(size >= 0, "Array::set_size(): New size must not be negative");
  if (ndata == size)
    return;
  if (copy) {
    // create a temporary pointer to the allocated data
    T* tmp = data;
    // store the current number of elements
    int old_ndata = ndata;
    // check how many elements we need to copy
    int min = (ndata < size) ? ndata : size;
    // allocate new memory
    alloc(size);
    // copy old elements into a new memory region
    for (int i = 0; i < min; ++i) {
      data[i] = tmp[i];
    }
    // initialize the rest of resized array
    for (int i = min; i < size; ++i) {
      data[i] = T();
    }
    // delete old elements
    destroy_elements(tmp, old_ndata);
  }
  else {
    free();
    alloc(size);
  }
}


template<class T> inline
T& Array<T>::operator()(int i)
{
  it_assert_debug(in_range(i), "Array::operator(): Improper index");
  return data[i];
}

template<class T> inline
const T& Array<T>::operator()(int i) const
{
  it_assert_debug(in_range(i), "Array::operator(): Improper index");
  return data[i];
}

template<class T> inline
const Array<T> Array<T>::operator()(int i1, int i2) const
{
  it_assert_debug(in_range(i1) && in_range(i2) && (i2 >= i1),
                  "Array::operator()(i1, i2): Improper indexes.");
  Array<T> s(i2 - i1 + 1);
  for (int i = 0; i < s.ndata; i++)
    s.data[i] = data[i1+i];
  return s;
}

template<class T> inline
const Array<T> Array<T>::operator()(const Array<int> &indices) const
{
  Array<T> a(indices.size());
  for (int i = 0; i < a.size(); i++) {
    it_assert_debug(in_range(indices(i)),
                    "Array::operator()(indices): Improper indices.");
    a(i) = data[indices(i)];
  }
  return a;
}

template<class T> inline
Array<T>& Array<T>::operator=(const Array<T> &a)
{
  if (this != &a) {
    set_size(a.ndata);
    for (int i = 0; i < ndata; i++)
      data[i] = a.data[i];
  }
  return *this;
}

template<class T> inline
Array<T>& Array<T>::operator=(const T &e)
{
  if (ndata == 0)
    set_size(1);
  for (int i = 0; i < ndata; i++)
    data[i] = e;
  return *this;
}

template<class T>
Array<T>& Array<T>::operator=(const char* values)
{
  std::istringstream buffer(values);
  buffer >> *this;
  return *this;
}


template<class T>
Array<T> Array<T>::left(int n) const
{
  it_assert_debug(in_range(n), "Array::left(): Index out of range");
  Array<T> tmp(n);
  for (int i = 0; i < n; ++i)
    tmp.data[i] = data[i];
  return tmp;
}

template<class T>
Array<T> Array<T>::right(int n) const
{
  it_assert_debug(in_range(n), "Array::right(): Index out of range");
  Array<T> tmp(n);
  for (int i = 0; i < n; ++i)
    tmp.data[i] = data[ndata-n+i];
  return tmp;
}

template<class T>
Array<T> Array<T>::mid(int pos, int n) const
{
  it_assert_debug((pos >= 0) && (n > 0) && (pos + n <= ndata), "Array::mid(): Indexing out of range");
  Array<T> tmp(n);
  for (int i = 0; i < n; ++i)
    tmp.data[i] = data[pos+i];
  return tmp;
}


template<class T>
T Array<T>::shift_right(const T& x)
{
  it_assert_debug(ndata > 0, "Array::shift_right(x): Array empty!");
  T ret;

  ret = data[ndata-1];
  for (int i = ndata - 1; i > 0; i--)
    data[i] = data[i-1];
  data[0] = x;

  return ret;
}


template<class T>
const Array<T> Array<T>::shift_right(const Array<T> &a)
{
  it_assert_debug(a.ndata <= ndata, "Array::shift_right(): Shift Array too large");
  Array<T> out(a.ndata);

  for (int i = 0; i < a.ndata; i++)
    out.data[i] = data[ndata-a.ndata+i];
  for (int i = ndata - 1; i >= a.ndata; i--)
    data[i] = data[i-a.ndata];
  for (int i = 0; i < a.ndata; i++)
    data[i] = a.data[i];

  return out;
}

template<class T>
T Array<T>::shift_left(const T& x)
{
  T temp = data[0];

  for (int i = 0; i < ndata - 1; i++)
    data[i] = data[i+1];
  data[ndata-1] = x;

  return temp;
}

template<class T>
const Array<T> Array<T>::shift_left(const Array<T> &a)
{
  it_assert_debug(a.ndata <= ndata, "Array::shift_left(): Shift Array too large");
  Array<T> out(a.ndata);

  for (int i = 0; i < a.ndata; i++)
    out.data[i] = data[i];
  for (int i = 0; i < ndata - a.ndata; i++)
    data[i] = data[i+a.ndata];
  for (int i = ndata - a.ndata; i < ndata; i++)
    data[i] = a.data[i-ndata+a.ndata];

  return out;
}

template<class T>
void Array<T>::swap(int i, int j)
{
  it_assert_debug(in_range(i) && in_range(j),
                  "Array::swap(): Indices out of range.");

  T temp = data[i];
  data[i] = data[j];
  data[j] = temp;
}

template<class T>
void Array<T>::set_subarray(int i1, int i2, const Array<T> &a)
{
  if (i1 == -1) i1 = ndata - 1;
  if (i2 == -1) i2 = ndata - 1;

  it_assert_debug(in_range(i1) && in_range(i2),
                  "Array<T>::set_subarray(): Indices out of range.");
  it_assert_debug(i2 >= i1, "Array<T>::set_subarray(): i2 >= i1 necessary.");
  it_assert_debug(i2 - i1 + 1 == a.ndata, "Array<T>::set_subarray(): Wrong sizes.");

  copy_vector(a.ndata, a.data, data + i1);
}

template<class T>
void Array<T>::set_subarray(int i1, int i2, const T &t)
{
  if (i1 == -1) i1 = ndata - 1;
  if (i2 == -1) i2 = ndata - 1;

  it_assert_debug(in_range(i1) && in_range(i2),
                  "Array<T>::set_subarray(): Indices out of range");
  it_assert_debug(i2 >= i1, "Array<T>::set_subarray(): i2 >= i1 necessary");

  for (int i = i1; i <= i2; i++)
    data[i] = t;
}

template<class T>
const Array<T> concat(const Array<T> &a, const T &e)
{
  Array<T> temp(a.size() + 1);

  for (int i = 0; i < a.size(); i++)
    temp(i) = a(i);
  temp(a.size()) = e;

  return temp;
}

template<class T>
const Array<T> concat(const T &e, const Array<T> &a)
{
  Array<T> temp(a.size() + 1);

  temp(0) = e;

  for (int i = 0; i < a.size(); i++)
    temp(i + 1) = a(i);

  return temp;
}

template<class T>
const Array<T> concat(const Array<T> &a1, const Array<T> &a2)
{
  Array<T> temp(a1.size() + a2.size());

  for (int i = 0; i < a1.size(); i++)
    temp(i) = a1(i);
  for (int i = 0; i < a2.size(); i++)
    temp(a1.size() + i) = a2(i);

  return temp;
}

template<class T>
const Array<T> concat(const Array<T> &a1, const Array<T> &a2,
                      const Array<T> &a3)
{
  // There should be some error control?
  Array<T> temp(a1.size() + a2.size() + a3.size());

  for (int i = 0; i < a1.size(); i++)
    temp(i) = a1(i);
  for (int i = 0; i < a2.size(); i++)
    temp(a1.size() + i) = a2(i);
  for (int i = 0; i < a3.size(); i++)
    temp(a1.size() + a2.size() + i) = a3(i);

  return temp;
}

/*!
  \relates Array
  \brief Output stream for Array<T>. T must have ostream operator<< defined.
*/
template<class T>
std::ostream &operator<<(std::ostream &os, const Array<T> &a)
{
  os << "{";
  for (int i = 0; i < a.size() - 1; i++)
    os << a(i) << " ";
  if (a.size() > 0)
    os << a(a.size() - 1);
  os << "}";

  return os;
}

/*!
  \relates Array
  \brief Input stream for Array<T>. T must have istream operator>> defined.
*/
template<class T>
std::istream &operator>>(std::istream &is, Array<T> &a)
{
  int nrof_elements = 0;
  char c;
  is >> c;
  if (c == '{') {
    is >> c;
    while (c != '}') {
      if (is.eof()) {
        is.setstate(std::ios_base::failbit);
        break;
      }
      if (c != ',') {  // Discard comma signs between elements
        is.putback(c);
      }
      if (++nrof_elements > a.size()) {
        a.set_size(nrof_elements, true);  // Too slow?
      }
      is >> a(nrof_elements - 1);
      is >> c;
    }
    if (a.size() > nrof_elements) {
      a.set_size(nrof_elements, true);
    }
  }
  else {
    is.setstate(std::ios_base::failbit);
  }

  return is;
}

/*!
  \relates Array
  \brief Assign a C-style string to an Array<T>. T must have istream
  operator>> defined.
*/
template<class T>
void set_array(Array<T> &a, const char *values)
{
  std::istringstream buffer(values);
  buffer >> a;
}

/*!
  \relates Array
  \brief Assign a string to an Array<T>. T must have istream operator>>
  defined.
*/
template<class T>
void set_array(Array<T> &a, const std::string &str)
{
  set_array(a, str.c_str());
}

} // namespace itpp

#endif // #ifndef ARRAY_H
