/*!
 * \file
 * \brief Sparse Vector Class definitions
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

#ifndef SVEC_H
#define SVEC_H

#include <itpp/base/vec.h>
#include <itpp/base/math/min_max.h>
#include <cstdlib>
#include <itpp/itexports.h>

namespace itpp
{

// Declaration of class Vec
template <class T> class Vec;
// Declaration of class Sparse_Vec
template <class T> class Sparse_Vec;

// ----------------------- Sparse_Vec Friends -------------------------------

//! v1+v2 where v1 and v2 are sparse vector
template <class T>
Sparse_Vec<T> operator+(const Sparse_Vec<T> &v1, const Sparse_Vec<T> &v2);

//! v1*v2 where v1 and v2 are sparse vectors
template <class T>
T operator*(const Sparse_Vec<T> &v1, const Sparse_Vec<T> &v2);

//! v1*v2 where v1 is a sparse vector and v2 is a dense vector
template <class T>
T operator*(const Sparse_Vec<T> &v1, const Vec<T> &v2);

//! v1*v2 where v1 is a dense vector and v2 is a sparse vector
template <class T>
T operator*(const Vec<T> &v1, const Sparse_Vec<T> &v2);

//! Elementwise multiplication of two sparse vectors returning a sparse vector
template <class T>
Sparse_Vec<T> elem_mult(const Sparse_Vec<T> &v1, const Sparse_Vec<T> &v2);

//! Elementwise multiplication of a sparse vector and a dense vector returning a dense vector
template <class T>
Vec<T> elem_mult(const Sparse_Vec<T> &v1, const Vec<T> &v2);

//! Elementwise multiplication of a sparse vector and a dense vector returning a sparse vector
template <class T>
Sparse_Vec<T> elem_mult_s(const Sparse_Vec<T> &v1, const Vec<T> &v2);

//! Elementwise multiplication of a dense vector and a sparse vector returning a dense vector
template <class T>
Vec<T> elem_mult(const Vec<T> &v1, const Sparse_Vec<T> &v2);

//! Elementwise multiplication of a dense vector and a sparse vector returning a sparse vector
template <class T>
Sparse_Vec<T> elem_mult_s(const Vec<T> &v1, const Sparse_Vec<T> &v2);

namespace details
{
	//this template selects appropriate type for Eps value used to remove small elements from
	//sparse containers
	template <typename NumT> struct Sparse_Eps_Type_Selector;
	template <> struct Sparse_Eps_Type_Selector<double> {typedef double eps_type;};
	template <> struct Sparse_Eps_Type_Selector<std::complex<double> > {typedef double eps_type;};
	template <> struct Sparse_Eps_Type_Selector<int> {typedef int eps_type;};
	template <> struct Sparse_Eps_Type_Selector<short> {typedef short eps_type;};
	template <> struct Sparse_Eps_Type_Selector<itpp::bin> {typedef int eps_type;};
}


/*!
  \brief Templated sparse vector class
  \author Tony Ottosson and Tobias Ringstrom

  A sparse vector is a vector where most elements are zero. The
  maximum number of none-zero elements is a parameter to the
  constructor. The elements are stored in random order, i.e. they
  are not sorted.

*/
template <class T>
class Sparse_Vec
{
public:

  //! Default constructor
  Sparse_Vec();

  /*!
    \brief Initiate an empty sparse vector

    \param sz Size of the sparse vector (i.e. maximum index is (\c sz - 1))
    \param data_init Maximum number of non-zero elements in the sparse vector (default value 200)
  */
  Sparse_Vec(int sz, int data_init = 200);

  /*!
    \brief Initiate a new sparse vector.

    \param v The elements of \c v are copied into the new sparse vector
  */
  Sparse_Vec(const Sparse_Vec<T> &v);

  /*!
    \brief Initiate a new sparse vector from a dense vector.

    \param v The elements of \c v are copied into the new sparse vector
  */
  Sparse_Vec(const Vec<T> &v);

  /*!
    \brief Initiate a new sparse vector from a dense vector. Elements of \c v larger than \c epsilon are copied into the new sparse vector.

    \note If the type T is complex<double>, then the elements of \c v larger than \c abs(epsilon) are copied into the new sparse vector.
  */
  Sparse_Vec(const Vec<T> &v, T epsilon);

  //! Destructor
  ~Sparse_Vec();

  /*!
    \brief Set the size \c sz of the sparse vector. Default value \c data_init=-1 \c => allocated size for the data is not changed.

    \param sz Size of the sparse vector (i.e. maximum index is (\c sz - 1))
    \param data_init Maximum number of non-zero elements in the sparse vector (default value -1 \c => allocated size for the data is not changed)
  */
  void set_size(int sz, int data_init = -1);

  //! Returns the size of the sparse vector
  int size() const { return v_size; }

  //! Number of non-zero elements in the sparse vector
  inline int nnz() {
    if (check_small_elems_flag) {
      remove_small_elements();
    }
    return used_size;
  }

  //! Returns the density of the sparse vector: (number of non-zero elements)/(size of the vector)
  double density();

  //! Set that all elements smaller than \a epsilon should be set to zero.
  void set_small_element(const T& epsilon);

  /*!
    Removes all elements that are smaller than \a epsilon from the non-zero elements.

    \note The small element \a epsilon can be set by the member function set_small_element. If no small value is set, the default value is always \c epsilon=0.
  */
  void remove_small_elements();

  /*!
    \brief Set the maximum number of non-zero elements to \c new_size

    \param new_size The new maximum number of non-zero elements.
  */
  void resize_data(int new_size);

  //! Set the maximum number of non-zero elements equal to the actual number of non-zero elements
  void compact();

  //! Returns a full, dense vector in \c v
  void full(Vec<T> &v) const;

  //! Returns a full, dense vector
  Vec<T> full() const;

  //! Returns the element with index \c i
  T operator()(int i) const;

  //! Set element \c i equal to \c v
  void set(int i, T v);

  //! Set the elements of the sparse vector with indices \c index_vec to the values in \c v
  void set(const ivec &index_vec, const Vec<T> &v);

  //! Set a new element with index \c i equal to \c v
  void set_new(int i, T v);

  //! Set new elements with indices \c index_vec equal to the values in \c v (no check whether the same index is used several times)
  void set_new(const ivec &index_vec, const Vec<T> &v);

  //! Add element \c i with \c v
  void add_elem(const int i, const T v);

  //! Add \c v to the elements specified by \c index_vec with \c v
  void add(const ivec &index_vec, const Vec<T> &v);

  //! Set the sparse vector to the all zero vector (removes all non-zero elements)
  void zeros();

  //! Set the i-th element to zero (i.e. clear that element if it contains a non-zero value)
  void zero_elem(const int i);

  //! Clear all non-zero elements of the sparse vector
  void clear();

  //! Clear the i-th element (if it contains a non-zero value)
  void clear_elem(const int i);

  /*!
    \brief Extract the reference to the p-th non-zero data element
  */
  inline void get_nz_data(int p, T& data_out) {
    if (check_small_elems_flag) {
      remove_small_elements();
    }
    data_out = data[p];
  }

  //! Returns the p-th non-zero data element
  inline T get_nz_data(int p) {
    if (check_small_elems_flag) {
      remove_small_elements();
    }
    return data[p];
  }

  //! Returns the vector index of the p-th non-zero element
  inline int get_nz_index(int p) {
    if (check_small_elems_flag) {
      remove_small_elements();
    }
    return index[p];
  }

  //! Returns the p-th non-zero value in \c dat and the corresponding vector index in \c idx
  inline void get_nz(int p, int &idx, T &dat) {
    if (check_small_elems_flag) {
      remove_small_elements();
    }
    idx = index[p];
    dat = data[p];
  }

  //! Return the indices of non-zero values
  ivec get_nz_indices();

  //! Return sparse subvector from index \c i1 to index \c i2
  Sparse_Vec<T> get_subvector(int i1, int i2) const;

  //! Returns the sum of all values squared
  T sqr() const;

  //! Assign sparse vector the value and length of the sparse vector \c v
  void operator=(const Sparse_Vec<T> &v);
  //! Assign sparse vector the value and length of the dense vector \c v
  void operator=(const Vec<T> &v);

  //! Returns the sign inverse of all elements in the sparse vector
  Sparse_Vec<T> operator-() const;

  //! Compare two sparse vectors. False if wrong sizes or different values
  bool operator==(const Sparse_Vec<T> &v);

  //! Add sparse vector \c v to all non-zero elements of the sparse vector
  void operator+=(const Sparse_Vec<T> &v);

  //! Add vector \c v to all non-zero elements of the sparse vector
  void operator+=(const Vec<T> &v);

  //! Subtract sparse vector \c v from all non-zero elements of the sparse vector
  void operator-=(const Sparse_Vec<T> &v);

  //! Subtract vector \c v from all non-zero elements of the sparse vector
  void operator-=(const Vec<T> &v);

  //! Multiply the scalar \c v to all non-zero elements of the sparse vector
  void operator*=(const T &v);

  //! Divide all non-zero elements of the sparse vector with the scalar \c v
  void operator/=(const T &v);

  //! Addition v1+v2 where v1 and v2 are sparse vector
  friend Sparse_Vec<T> operator+<>(const Sparse_Vec<T> &v1, const Sparse_Vec<T> &v2);
  //! Scalar product v1*v2 where v1 and v2 are sparse vectors
  friend T operator*<>(const Sparse_Vec<T> &v1, const Sparse_Vec<T> &v2);
  //! Scalar product v1*v2 where v1 is a sparse vector and v2 is a dense vector
  friend T operator*<>(const Sparse_Vec<T> &v1, const Vec<T> &v2);
  //! Scalar product v1*v2 where v1 is a dense vector and v2 is a sparse vector
  friend T operator*<>(const Vec<T> &v1, const Sparse_Vec<T> &v2);

  //! Element wise multiplication of two sparse vectors
  friend Sparse_Vec<T> elem_mult <>(const Sparse_Vec<T> &v1, const Sparse_Vec<T> &v2);

  //! Element wise multiplication of a sparse vector and a dense vector
  friend Vec<T> elem_mult <>(const Sparse_Vec<T> &v1, const Vec<T> &v2);

  //! Element wise multiplication of a sparse vector and a dense vector returning a sparse vector
  friend Sparse_Vec<T> elem_mult_s <>(const Sparse_Vec<T> &v1, const Vec<T> &v2);

  //! Element wise multiplication of a a dense vector and a sparse vector
  friend Vec<T> elem_mult <>(const Vec<T> &v1, const Sparse_Vec<T> &v2);

  //! Element wise multiplication of a a dense vector and a sparse vector returning a sparse vector
  friend Sparse_Vec<T> elem_mult_s <>(const Vec<T> &v1, const Sparse_Vec<T> &v2);

private:
  void init();
  void alloc();
  void free();

  int v_size, used_size, data_size;
  T *data;
  int *index;
  typename details::Sparse_Eps_Type_Selector<T>::eps_type eps;
  bool check_small_elems_flag;
};


/*!
  \relates Sparse_Vec
  \brief Type definition of an integer sparse vector
*/
typedef Sparse_Vec<int> sparse_ivec;

/*!
  \relates Sparse_Vec
  \brief Type definition of a double sparse vector
*/
typedef Sparse_Vec<double> sparse_vec;

/*!
  \relates Sparse_Vec
  \brief Type definition of a complex<double> sparse vector
*/
typedef Sparse_Vec<std::complex<double> > sparse_cvec;

// ----------------------- Implementation starts here --------------------------------

template <class T>
void Sparse_Vec<T>::init()
{
  v_size = 0;
  used_size = 0;
  data_size = 0;
  data = 0;
  index = 0;
  eps = 0;
  check_small_elems_flag = true;
}

template <class T>
void Sparse_Vec<T>::alloc()
{
  if (data_size != 0) {
    data = new T[data_size];
    index = new int[data_size];
  }
}

template <class T>
void Sparse_Vec<T>::free()
{
  delete [] data;
  data = 0;
  delete [] index;
  index = 0;
}

template <class T>
Sparse_Vec<T>::Sparse_Vec()
{
  init();
}

template <class T>
Sparse_Vec<T>::Sparse_Vec(int sz, int data_init)
{
  init();
  v_size = sz;
  used_size = 0;
  data_size = data_init;
  alloc();
}

template <class T>
Sparse_Vec<T>::Sparse_Vec(const Sparse_Vec<T> &v)
{
  init();
  v_size = v.v_size;
  used_size = v.used_size;
  data_size = v.data_size;
  eps = v.eps;
  check_small_elems_flag = v.check_small_elems_flag;
  alloc();

  for (int i = 0; i < used_size; i++) {
    data[i] = v.data[i];
    index[i] = v.index[i];
  }
}

template <class T>
Sparse_Vec<T>::Sparse_Vec(const Vec<T> &v)
{
  init();
  v_size = v.size();
  used_size = 0;
  data_size = std::min(v.size(), 10000);
  alloc();

  for (int i = 0; i < v_size; i++) {
    if (v(i) != T(0)) {
      if (used_size == data_size)
        resize_data(data_size*2);
      data[used_size] = v(i);
      index[used_size] = i;
      used_size++;
    }
  }
  compact();
}

template <class T>
Sparse_Vec<T>::Sparse_Vec(const Vec<T> &v, T epsilon)
{
  init();
  v_size = v.size();
  used_size = 0;
  data_size = std::min(v.size(), 10000);
  eps = std::abs(epsilon);
  alloc();

  for (int i = 0; i < v_size; i++) {
    if (std::abs(v(i)) > eps) {
      if (used_size == data_size)
        resize_data(data_size*2);
      data[used_size] = v(i);
      index[used_size] = i;
      used_size++;
    }
  }
  compact();
}

template <class T>
Sparse_Vec<T>::~Sparse_Vec()
{
  free();
}

template <class T>
void Sparse_Vec<T>::set_size(int new_size, int data_init)
{
  v_size = new_size;
  used_size = 0;
  if (data_init != -1) {
    free();
    data_size = data_init;
    alloc();
  }
}

template <class T>
double Sparse_Vec<T>::density()
{
  if (check_small_elems_flag) {
    remove_small_elements();
  }
  //return static_cast<double>(used_size) / v_size;
  return double(used_size) / v_size;
}

template <class T>
void Sparse_Vec<T>::set_small_element(const T& epsilon)
{
  eps = std::abs(epsilon);
  remove_small_elements();
}

template <class T>
void Sparse_Vec<T>::remove_small_elements()
{
  int i;
  int nrof_removed_elements = 0;

  //Remove small elements
  for (i = 0;i < used_size;i++) {
    if (std::abs(data[i]) <= eps) {
      nrof_removed_elements++;
    }
    else if (nrof_removed_elements > 0) {
      data[i-nrof_removed_elements] = data[i];
      index[i-nrof_removed_elements] = index[i];
    }
  }

  //Set new size after small elements have been removed
  used_size -= nrof_removed_elements;

  //Set the flag to indicate that all small elements have been removed
  check_small_elems_flag = false;
}


template <class T>
void Sparse_Vec<T>::resize_data(int new_size)
{
  it_assert(new_size >= used_size, "Sparse_Vec<T>::resize_data(int new_size): New size is to small");

  if (new_size != data_size) {
    if (new_size == 0)
      free();
    else {
      T *tmp_data = data;
      int *tmp_pos = index;
      data_size = new_size;
      alloc();
      for (int p = 0; p < used_size; p++) {
        data[p] = tmp_data[p];
        index[p] = tmp_pos[p];
      }
      delete [] tmp_data;
      delete [] tmp_pos;
    }
  }
}

template <class T>
void Sparse_Vec<T>::compact()
{
  if (check_small_elems_flag) {
    remove_small_elements();
  }
  resize_data(used_size);
}

template <class T>
void Sparse_Vec<T>::full(Vec<T> &v) const
{
  v.set_size(v_size);

  v = T(0);
  for (int p = 0; p < used_size; p++)
    v(index[p]) = data[p];
}

template <class T>
Vec<T> Sparse_Vec<T>::full() const
{
  Vec<T> r(v_size);
  full(r);
  return r;
}

// This is slow. Implement a better search
template <class T>
T Sparse_Vec<T>::operator()(int i) const
{
  it_assert_debug(i >= 0 && i < v_size, "The index of the element is out of range");

  bool found = false;
  int p;
  for (p = 0; p < used_size; p++) {
    if (index[p] == i) {
      found = true;
      break;
    }
  }
  return found ? data[p] : T(0);
}

template <class T>
void Sparse_Vec<T>::set(int i, T v)
{
  it_assert_debug(i >= 0 && i < v_size, "The index of the element is out of range");

  bool found = false;
  bool larger_than_eps;
  int p;

  for (p = 0; p < used_size; p++) {
    if (index[p] == i) {
      found = true;
      break;
    }
  }

  larger_than_eps = (std::abs(v) > eps);

  if (found && larger_than_eps)
    data[p] = v;
  else if (larger_than_eps) {
    if (used_size == data_size)
      resize_data(data_size*2 + 100);
    data[used_size] = v;
    index[used_size] = i;
    used_size++;
  }

  //Check if the stored element is smaller than eps. In that case it should be removed.
  if (std::abs(v) <= eps) {
    remove_small_elements();
  }

}

template <class T>
void Sparse_Vec<T>::set_new(int i, T v)
{
  it_assert_debug(v_size > i, "The index of the element exceeds the size of the sparse vector");

  //Check that the new element is larger than eps!
  if (std::abs(v) > eps) {
    if (used_size == data_size)
      resize_data(data_size*2 + 100);
    data[used_size] = v;
    index[used_size] = i;
    used_size++;
  }
}

template <class T>
void Sparse_Vec<T>::add_elem(const int i, const T v)
{
  bool found = false;
  int p;

  it_assert_debug(v_size > i, "The index of the element exceeds the size of the sparse vector");

  for (p = 0; p < used_size; p++) {
    if (index[p] == i) {
      found = true;
      break;
    }
  }
  if (found)
    data[p] += v;
  else {
    if (used_size == data_size)
      resize_data(data_size*2 + 100);
    data[used_size] = v;
    index[used_size] = i;
    used_size++;
  }

  check_small_elems_flag = true;

}

template <class T>
void Sparse_Vec<T>::add(const ivec& index_vec, const Vec<T>& v)
{
  bool found = false;
  int i, p, q;
  int nrof_nz = v.size();

  it_assert_debug(v_size > max(index_vec), "The indices exceeds the size of the sparse vector");

  //Elements are added if they have identical indices
  for (q = 0; q < nrof_nz; q++) {
    i = index_vec(q);
    for (p = 0; p < used_size; p++) {
      if (index[p] == i) {
        found = true;
        break;
      }
    }
    if (found)
      data[p] += v(q);
    else {
      if (used_size == data_size)
        resize_data(data_size*2 + 100);
      data[used_size] = v(q);
      index[used_size] = i;
      used_size++;
    }
    found = false;
  }

  check_small_elems_flag = true;

}

template <class T>
void Sparse_Vec<T>::zeros()
{
  used_size = 0;
  check_small_elems_flag = false;
}

template <class T>
void Sparse_Vec<T>::zero_elem(const int i)
{
  bool found = false;
  int p;

  it_assert_debug(v_size > i, "The index of the element exceeds the size of the sparse vector");

  for (p = 0; p < used_size; p++) {
    if (index[p] == i) {
      found = true;
      break;
    }
  }
  if (found) {
    data[p] = data[used_size-1];
    index[p] = index[used_size-1];
    used_size--;
  }
}

template <class T>
void Sparse_Vec<T>::clear()
{
  used_size = 0;
  check_small_elems_flag = false;
}

template <class T>
void Sparse_Vec<T>::clear_elem(const int i)
{
  bool found = false;
  int p;

  it_assert_debug(v_size > i, "The index of the element exceeds the size of the sparse vector");

  for (p = 0; p < used_size; p++) {
    if (index[p] == i) {
      found = true;
      break;
    }
  }
  if (found) {
    data[p] = data[used_size-1];
    index[p] = index[used_size-1];
    used_size--;
  }
}

template <class T>
void Sparse_Vec<T>::set(const ivec& index_vec, const Vec<T>& v)
{
  it_assert_debug(v_size > max(index_vec), "The indices exceeds the size of the sparse vector");

  //Clear all old non-zero elements
  clear();

  //Add the new non-zero elements
  add(index_vec, v);
}

template <class T>
void Sparse_Vec<T>::set_new(const ivec& index_vec, const Vec<T>& v)
{
  int q;
  int nrof_nz = v.size();

  it_assert_debug(v_size > max(index_vec), "The indices exceeds the size of the sparse vector");

  //Clear all old non-zero elements
  clear();

  for (q = 0; q < nrof_nz; q++) {
    if (std::abs(v[q]) > eps) {
      if (used_size == data_size)
        resize_data(data_size*2 + 100);
      data[used_size] = v(q);
      index[used_size] = index_vec(q);
      used_size++;
    }
  }
}

template <class T>
ivec Sparse_Vec<T>::get_nz_indices()
{
  int n = nnz();
  ivec r(n);
  for (int i = 0; i < n; i++) {
    r(i) = get_nz_index(i);
  }
  return r;
}

template <class T>
Sparse_Vec<T> Sparse_Vec<T>::get_subvector(int i1, int i2) const
{
  it_assert_debug(v_size > i1 && v_size > i2 && i1 <= i2 && i1 >= 0, "The index of the element exceeds the size of the sparse vector");

  Sparse_Vec<T> r(i2 - i1 + 1);

  for (int p = 0; p < used_size; p++) {
    if (index[p] >= i1 && index[p] <= i2) {
      if (r.used_size == r.data_size)
        r.resize_data(r.data_size*2 + 100);
      r.data[r.used_size] = data[p];
      r.index[r.used_size] = index[p] - i1;
      r.used_size++;
    }
  }
  r.eps = eps;
  r.check_small_elems_flag = check_small_elems_flag;
  r.compact();

  return r;
}

template <class T>
T Sparse_Vec<T>::sqr() const
{
  T sum(0);
  for (int p = 0; p < used_size; p++)
    sum += data[p] * data[p];

  return sum;
}

template <class T>
void Sparse_Vec<T>::operator=(const Sparse_Vec<T> &v)
{
  free();
  v_size = v.v_size;
  used_size = v.used_size;
  data_size = v.data_size;
  eps = v.eps;
  check_small_elems_flag = v.check_small_elems_flag;
  alloc();

  for (int i = 0; i < used_size; i++) {
    data[i] = v.data[i];
    index[i] = v.index[i];
  }
}

template <class T>
void Sparse_Vec<T>::operator=(const Vec<T> &v)
{
  free();
  v_size = v.size();
  used_size = 0;
  data_size = std::min(v.size(), 10000);
  eps = std::abs(T(0));
  check_small_elems_flag = false;
  alloc();

  for (int i = 0; i < v_size; i++) {
    if (v(i) != T(0)) {
      if (used_size == data_size)
        resize_data(data_size*2);
      data[used_size] = v(i);
      index[used_size] = i;
      used_size++;
    }
  }
  compact();
}

template <class T>
Sparse_Vec<T> Sparse_Vec<T>::operator-() const
{
  Sparse_Vec r(v_size, used_size);

  for (int p = 0; p < used_size; p++) {
    r.data[p] = -data[p];
    r.index[p] = index[p];
  }
  r.used_size = used_size;

  return r;
}

template <class T>
bool Sparse_Vec<T>::operator==(const Sparse_Vec<T> &v)
{
  int p, q;
  bool found = false;

  //Remove small elements before comparing the two sparse_vectors
  if (check_small_elems_flag)
    remove_small_elements();

  if (v_size != v.v_size) {
    //Return false if vector sizes are unequal
    return false;
  }
  else {
    for (p = 0;p < used_size;p++) {
      for (q = 0;q < v.used_size;q++) {
        if (index[p] == v.index[q]) {
          found = true;
          break;
        }
      }
      if (found == false)
        //Return false if non-zero element not found, or if elements are unequal
        return false;
      else if (data[p] != v.data[q])
        //Return false if non-zero element not found, or if elements are unequal
        return false;
      else
        //Check next non-zero element
        found = false;
    }
  }

  /*Special handling if sizes do not match.
  Required since v may need to do remove_small_elements() for true comparison*/
  if (used_size != v.used_size) {
    if (used_size > v.used_size) {
      //Return false if number of non-zero elements is less in v
      return false;
    }
    else {
      //Ensure that the remaining non-zero elements in v are smaller than v.eps
      int nrof_small_elems = 0;
      for (q = 0;q < v.used_size;q++) {
        if (std::abs(v.data[q]) <= v.eps)
          nrof_small_elems++;
      }
      if (v.used_size - nrof_small_elems != used_size)
        //Return false if the number of "true" non-zero elements are unequal
        return false;
    }
  }

  //All elements checks => return true
  return true;
}

template <class T>
void Sparse_Vec<T>::operator+=(const Sparse_Vec<T> &v)
{
  int i, p;
  T tmp_data;
  int nrof_nz_v = v.used_size;

  it_assert_debug(v_size == v.size(), "Attempted addition of unequal sized sparse vectors");

  for (p = 0; p < nrof_nz_v; p++) {
    i = v.index[p];
    tmp_data = v.data[p];
    //get_nz(p,i,tmp_data);
    add_elem(i, tmp_data);
  }

  check_small_elems_flag = true;
}

template <class T>
void Sparse_Vec<T>::operator+=(const Vec<T> &v)
{
  int i;

  it_assert_debug(v_size == v.size(), "Attempted addition of unequal sized sparse vectors");

  for (i = 0; i < v.size(); i++)
    if (v(i) != T(0))
      add_elem(i, v(i));

  check_small_elems_flag = true;
}


template <class T>
void Sparse_Vec<T>::operator-=(const Sparse_Vec<T> &v)
{
  int i, p;
  T tmp_data;
  int nrof_nz_v = v.used_size;

  it_assert_debug(v_size == v.size(), "Attempted subtraction of unequal sized sparse vectors");

  for (p = 0; p < nrof_nz_v; p++) {
    i = v.index[p];
    tmp_data = v.data[p];
    //v.get_nz(p,i,tmp_data);
    add_elem(i, -tmp_data);
  }

  check_small_elems_flag = true;
}

template <class T>
void Sparse_Vec<T>::operator-=(const Vec<T> &v)
{
  int i;

  it_assert_debug(v_size == v.size(), "Attempted subtraction of unequal sized sparse vectors");

  for (i = 0; i < v.size(); i++)
    if (v(i) != T(0))
      add_elem(i, -v(i));

  check_small_elems_flag = true;
}

template <class T>
void Sparse_Vec<T>::operator*=(const T &v)
{
  int p;

  for (p = 0; p < used_size; p++) {
    data[p] *= v;
  }

  check_small_elems_flag = true;
}

template <class T>
void Sparse_Vec<T>::operator/=(const T &v)
{
  int p;
  for (p = 0; p < used_size; p++) {
    data[p] /= v;
  }

  if (eps > 0) {
    check_small_elems_flag = true;
  }
}

template <class T>
T operator*(const Sparse_Vec<T> &v1, const Sparse_Vec<T> &v2)
{
  it_assert_debug(v1.v_size == v2.v_size, "Sparse_Vec<T> * Sparse_Vec<T>");

  T sum(0);
  Vec<T> v1f(v1.v_size);
  v1.full(v1f);
  for (int p = 0; p < v2.used_size; p++) {
    if (v1f[v2.index[p]] != T(0))
      sum += v1f[v2.index[p]] * v2.data[p];
  }

  return sum;
}

template <class T>
T operator*(const Sparse_Vec<T> &v1, const Vec<T> &v2)
{
  it_assert_debug(v1.size() == v2.size(), "Multiplication of unequal sized vectors attempted");

  T sum(0);
  for (int p1 = 0; p1 < v1.used_size; p1++)
    sum += v1.data[p1] * v2[v1.index[p1]];

  return sum;
}

template <class T>
T operator*(const Vec<T> &v1, const Sparse_Vec<T> &v2)
{
  it_assert_debug(v1.size() == v2.size(), "Multiplication of unequal sized vectors attempted");

  T sum(0);
  for (int p2 = 0; p2 < v2.used_size; p2++)
    sum += v1[v2.index[p2]] * v2.data[p2];

  return sum;
}

template <class T>
Sparse_Vec<T> elem_mult(const Sparse_Vec<T> &v1, const Sparse_Vec<T> &v2)
{
  it_assert_debug(v1.v_size == v2.v_size, "elem_mult(Sparse_Vec<T>, Sparse_Vec<T>)");

  Sparse_Vec<T> r(v1.v_size);
  ivec pos(v1.v_size);
  pos = -1;
  for (int p1 = 0; p1 < v1.used_size; p1++)
    pos[v1.index[p1]] = p1;
  for (int p2 = 0; p2 < v2.used_size; p2++) {
    if (pos[v2.index[p2]] != -1) {
      if (r.used_size == r.data_size)
        r.resize_data(r.used_size*2 + 100);
      r.data[r.used_size] = v1.data[pos[v2.index[p2]]] * v2.data[p2];
      r.index[r.used_size] = v2.index[p2];
      r.used_size++;
    }
  }
  r.compact();

  return r;
}

template <class T>
Vec<T> elem_mult(const Sparse_Vec<T> &v1, const Vec<T> &v2)
{
  it_assert_debug(v1.v_size == v2.size(), "elem_mult(Sparse_Vec<T>, Vec<T>)");

  Vec<T> r(v1.v_size);
  r = T(0);
  for (int p1 = 0; p1 < v1.used_size; p1++)
    r[v1.index[p1]] = v1.data[p1] * v2[v1.index[p1]];

  return r;
}

template <class T>
Sparse_Vec<T> elem_mult_s(const Sparse_Vec<T> &v1, const Vec<T> &v2)
{
  it_assert_debug(v1.v_size == v2.size(), "elem_mult(Sparse_Vec<T>, Vec<T>)");

  Sparse_Vec<T> r(v1.v_size);
  for (int p1 = 0; p1 < v1.used_size; p1++) {
    if (v2[v1.index[p1]] != T(0)) {
      if (r.used_size == r.data_size)
        r.resize_data(r.used_size*2 + 100);
      r.data[r.used_size] = v1.data[p1] * v2[v1.index[p1]];
      r.index[r.used_size] = v1.index[p1];
      r.used_size++;
    }
  }
  r.compact();

  return r;
}

template <class T>
Vec<T> elem_mult(const Vec<T> &v1, const Sparse_Vec<T> &v2)
{
  it_assert_debug(v1.size() == v2.v_size, "elem_mult(Vec<T>, Sparse_Vec<T>)");

  Vec<T> r(v2.v_size);
  r = T(0);
  for (int p2 = 0; p2 < v2.used_size; p2++)
    r[v2.index[p2]] = v1[v2.index[p2]] * v2.data[p2];

  return r;
}

template <class T>
Sparse_Vec<T> elem_mult_s(const Vec<T> &v1, const Sparse_Vec<T> &v2)
{
  it_assert_debug(v1.size() == v2.v_size, "elem_mult(Vec<T>, Sparse_Vec<T>)");

  Sparse_Vec<T> r(v2.v_size);
  for (int p2 = 0; p2 < v2.used_size; p2++) {
    if (v1[v2.index[p2]] != T(0)) {
      if (r.used_size == r.data_size)
        r.resize_data(r.used_size*2 + 100);
      r.data[r.used_size] = v1[v2.index[p2]] * v2.data[p2];
      r.index[r.used_size] = v2.index[p2];
      r.used_size++;
    }
  }
  r.compact();

  return r;
}

template <class T>
Sparse_Vec<T> operator+(const Sparse_Vec<T> &v1, const Sparse_Vec<T> &v2)
{
  it_assert_debug(v1.v_size == v2.v_size, "Sparse_Vec<T> + Sparse_Vec<T>");

  Sparse_Vec<T> r(v1);
  ivec pos(v1.v_size);
  pos = -1;
  for (int p1 = 0; p1 < v1.used_size; p1++)
    pos[v1.index[p1]] = p1;
  for (int p2 = 0; p2 < v2.used_size; p2++) {
    if (pos[v2.index[p2]] == -1) {// A new entry
      if (r.used_size == r.data_size)
        r.resize_data(r.used_size*2 + 100);
      r.data[r.used_size] = v2.data[p2];
      r.index[r.used_size] = v2.index[p2];
      r.used_size++;
    }
    else
      r.data[pos[v2.index[p2]]] += v2.data[p2];
  }
  r.check_small_elems_flag = true;  // added dec 7, 2006
  r.compact();

  return r;
}

//! Convert a dense vector \c v into its sparse representation
template <class T>
inline Sparse_Vec<T> sparse(const Vec<T> &v)
{
  Sparse_Vec<T> s(v);
  return s;
}

//! Convert a dense vector \c v into its sparse representation
template <class T>
inline Sparse_Vec<T> sparse(const Vec<T> &v, T epsilon)
{
  Sparse_Vec<T> s(v, epsilon);
  return s;
}

//! Convert a sparse vector \c s into its dense representation
template <class T>
inline Vec<T> full(const Sparse_Vec<T> &s)
{
  Vec<T> v;
  s.full(v);
  return v;
}

//! \cond

// ---------------------------------------------------------------------
// Instantiations
// ---------------------------------------------------------------------

ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Sparse_Vec<int>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Sparse_Vec<double>;
ITPP_EXPORT_TEMPLATE template class ITPP_EXPORT Sparse_Vec<std::complex<double> >;

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_ivec operator+(const sparse_ivec &,
                                        const sparse_ivec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_vec operator+(const sparse_vec &,
                                       const sparse_vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_cvec operator+(const sparse_cvec &,
                                        const sparse_cvec &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT int operator*(const sparse_ivec &, const sparse_ivec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT double operator*(const sparse_vec &, const sparse_vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT std::complex<double> operator*(const sparse_cvec &,
    const sparse_cvec &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT int operator*(const sparse_ivec &, const ivec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT double operator*(const sparse_vec &, const vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT std::complex<double> operator*(const sparse_cvec &,
    const cvec &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT int operator*(const ivec &, const sparse_ivec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT double operator*(const vec &, const sparse_vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT std::complex<double> operator*(const cvec &,
    const sparse_cvec &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_ivec elem_mult(const sparse_ivec &,
                                        const sparse_ivec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_vec elem_mult(const sparse_vec &, const sparse_vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_cvec elem_mult(const sparse_cvec &,
                                        const sparse_cvec &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec elem_mult(const sparse_ivec &, const ivec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec elem_mult(const sparse_vec &, const vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec elem_mult(const sparse_cvec &, const cvec &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_ivec elem_mult_s(const sparse_ivec &, const ivec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_vec elem_mult_s(const sparse_vec &, const vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_cvec elem_mult_s(const sparse_cvec &, const cvec &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT ivec elem_mult(const ivec &, const sparse_ivec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT vec elem_mult(const vec &, const sparse_vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT cvec elem_mult(const cvec &, const sparse_cvec &);

ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_ivec elem_mult_s(const ivec &, const sparse_ivec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_vec elem_mult_s(const vec &, const sparse_vec &);
ITPP_EXPORT_TEMPLATE template ITPP_EXPORT sparse_cvec elem_mult_s(const cvec &, const sparse_cvec &);

//! \endcond

} // namespace itpp

#endif // #ifndef SVEC_H

