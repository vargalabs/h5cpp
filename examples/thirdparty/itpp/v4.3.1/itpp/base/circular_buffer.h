/*!
 * \file
 * \brief Circular_Buffer class (container)
 * \author Tobias Tynderfeldt
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
 *
 * This file is not separated into .h and .cpp files. The reason is
 * to avoid problems with template initializations of this class.
 * A \c Circular_Buffer<type> can contain any type and it is not
 * possible to initialize and pre-compile all types that might be put
 * into a \c Circular_Buffer.
 */

#ifndef CIRCULAR_BUFFER_H
#define CIRCULAR_BUFFER_H

#include <itpp/base/vec.h>
#include <itpp/base/array.h>


namespace itpp
{

/*!
  \brief General circular buffer class

  This class is a general circular buffer class for arbitrary types.

  For rarely used types you will need to instantiate the class by
  \code
  template class Circular_Buffer<type>;
  \endcode

  The following example shows how to define a Circular_Buffer of doubles:
  \code
  vec a = randn(3);
  vec b = randn(8);
  vec out_vec;
  Array <double> out_array;

  Circular_Buffer<double> cb1(10);

  //Put the elements of a to the buffer
  cb1.put(a);

  //Vector output: Peek at the two oldest elements of the buffer (the two first extracted)
  cb1.peek(out_vec,2);
  cout << "peek(out_vec,2) = " << out_vec << ": display 2 first elements of the buffer, without affecting the content"  << endl;

  //Vector output: Peek at all elements of the buffer in reverse order
  cb1.peek_reverse(out_vec);
  cout << "peek_reverse(out_vec,-1) = " << out_vec << ": display buffer, without affecting the content"  << endl;

  //Array output: Peek at all elements of the buffer in reverse order
  cb1.peek_reverse(out_array);
  cout << "peek_reverse(out_array,-1) = " << out_array << ": display buffer, without affecting the content"  << endl;

  //Put the elements of \c b to the buffer
  cb1.put(b);

  //Extract the oldest element of the buffer
  cb1.get(out_vec,1);
  cout << "get(out_vec,1) = " << out_vec << endl ;

  //Extract all element of the buffer
  cb1.get(out_vec);
  cout << "get(out_vec) = " << out_vec << endl ;
  \endcode
*/
template<class T>
class Circular_Buffer
{
public:
  //! Default constructor
  Circular_Buffer();

  //! Create a Circular_Buffer of size \c n
  Circular_Buffer(int n);

  //! Create a copy of \c s
  Circular_Buffer(const Circular_Buffer<T> &s);

  //! Default destructor
  virtual ~Circular_Buffer();

  //! Write the element \a in to the buffer
  void put(const T& in);

  //! Write the vector of elements \a in to the circular buffer
  void put(const Vec<T>& in);

  //! Write the vector of elements \a in to the circular buffer
  void put(const Array<T>& in);

  //! Get the oldest element in the circular buffer
  void get(T& out);

  //! Get the oldest element in the circular buffer
  T get();

  //! Get the N oldest element in the circular buffer. \c N=-1 returns all elements in the buffer
  void get(Vec<T>& out, const int N = -1);

  //! Get the N oldest element in the circular buffer. \c N=-1 returns all elements in the buffer
  void get(Array<T>& out, const int N = -1);

  //! Peek at the oldest element in the circular buffer, without removing it.
  void peek(T& out) const;

  //! Peek at the oldest element in the circular buffer, without removing it.
  T peek() const;

  //! Peek at the element with index \a index in the circular buffer, without removing it.
  void peek(const int index, T& out) const;

  //! Peek at the N first elements of the circular buffer, without removing them. \c N=-1 peeks all elements in the buffer
  void peek(Vec<T>& out, const int N = -1) const;

  //! Peek at the elements with index \a index in the circular buffer, without removing them.
  void peek(const ivec& index, Vec<T>& out) const;

  //! Peek at the N first elements of the circular buffer, without removing them. \c N=-1 peeks all elements in the buffer
  void peek(Array<T>& out, const int N = -1) const;

  //! Peek at the elements with index \a index in the circular buffer, without removing them.
  void peek(const ivec& index, Array<T>& out) const;

  //! Peek at the latest element in the circular buffer, without removing it.
  void peek_reverse(T& out) const;

  //! Peek at the latest element in the circular buffer, without removing it.
  T peek_reverse() const;

  //! Peek at the N latest elements of the circular buffer in reverse order, without removing them. \c N=-1 returns all elements in the buffer
  void peek_reverse(Vec<T>& out, const int N = -1) const;

  //! Peek at the N latest elements of the circular buffer in reverse order, without removing them. \c N=-1 returns all elements in the buffer
  void peek_reverse(Array<T>& out, const int N = -1) const;

  //! Empty the circular buffer
  void clear();

  //! Assignment operator
  void operator=(const Circular_Buffer<T> &s);

  //! Returns the maximum number of data elements the circular buffer can store
  int size() const { return _ndata; }

  //! Returns the number of data elements currently stored in the circular buffer
  int nrof_elements() const { return _rw_dist; }

  //! Resizing a Circular_Buffer<T>.
  void set_size(int n, bool copy = false);

private:

  int _write;
  int _read;
  int _ndata;
  int _rw_dist;
  T *_data;

  void alloc(int n);
  void free();

};

// --------------------------- Implementation starts here ----------------------------------

template<class T>
Circular_Buffer<T>::Circular_Buffer()
{
  _data    = 0;
  _ndata   = 0;
  _rw_dist = 0;
  _read    = 0;
  _write   = 0;
}

template<class T>
Circular_Buffer<T>::Circular_Buffer(int n)
{
  alloc(n);
  _read    = 0;
  _write   = 0;
  _rw_dist = 0;
}

template<class T>
Circular_Buffer<T>::Circular_Buffer(const Circular_Buffer<T> &cb)
{
  _data    = NULL;
  _ndata   = 0;
  _read    = cb._read;
  _write   = cb._write;
  _rw_dist = cb._rw_dist;

  alloc(cb._ndata);
  for (int i = 0; i < cb._ndata; i++) { _data[i] = cb._data[i]; }
}

template<class T>
Circular_Buffer<T>::~Circular_Buffer()
{
  free();
}

template <class T>
void Circular_Buffer<T>::get(T& out)
{
  it_assert_debug(_rw_dist > 0, "Buffer empty. No data left to read from the buffer.");
  out = _data[_read];
  _read++;
  _rw_dist--;

  if (_read == _ndata) { _read = 0; }
}

template <class T>
T Circular_Buffer<T>::get()
{
  T out;

  get(out);
  return out;
}

template <class T>
void Circular_Buffer<T>::get(Vec<T>& out, const int N)
{
  int N_out;

  if (N == -1)
    N_out = _rw_dist;
  else
    N_out = N;

  out.set_size(N_out);

  for (int i = 0;i < N_out;i++) {
    it_assert_debug(_rw_dist > 0, "Buffer empty. No data left to read from the buffer.");
    out(i) = _data[_read];
    _read++;
    _rw_dist--;

    if (_read == _ndata)
      _read = 0;
  }
}

template <class T>
void Circular_Buffer<T>::get(Array<T>& out, const int N)
{
  int N_out;

  if (N == -1)
    N_out = _rw_dist;
  else
    N_out = N;

  out.set_size(N_out);

  for (int i = 0;i < N_out;i++) {
    it_assert_debug(_rw_dist > 0, "Buffer empty. No data left to read from the buffer.");
    out(i) = _data[_read];
    _read++;
    _rw_dist--;

    if (_read == _ndata)
      _read = 0;
  }
}

template <class T>
void Circular_Buffer<T>::peek(T& out) const
{
  it_assert_debug(_rw_dist > 0, "Attempted to peek at an empty buffer.");
  out = _data[_read];
}

template <class T>
T Circular_Buffer<T>::peek() const
{
  T out;

  peek(out);
  return out;
}

template <class T>
void Circular_Buffer<T>::peek(const int index, T& out) const
{
  it_assert_debug(_rw_dist > index && index >= 0, "The index exceeds the number of elements stored in the buffer.");
  out = _data[(_read+index)%_ndata];
}

template <class T>
void Circular_Buffer<T>::peek(Vec<T>& out, const int N) const
{
  int N_out;
  int read_tmp = _read;

  if (N == -1)
    N_out = _rw_dist;
  else
    N_out = N;

  it_assert_debug(_rw_dist >= N_out, "Attempted to peek at more elements than there are stored in the buffer.");
  out.set_size(N_out);

  for (int i = 0;i < N_out;i++) {
    out(i) = _data[read_tmp];
    read_tmp++;
    if (read_tmp == _ndata)
      read_tmp = 0;
  }
}

template <class T>
void Circular_Buffer<T>::peek(const ivec& index, Vec<T>& out) const
{
  out.set_size(index.size());

  for (int i = 0;i < index.size();i++) {
    it_assert_debug(_rw_dist >= index(i) && index(i) >= 0, "Attempted to peek at an element, whose index exceeds the number of buffered elements.");
    out(i) = _data[(_read+index(i))%_ndata];
  }
}

template <class T>
void Circular_Buffer<T>::peek(Array<T>& out, const int N) const
{
  int N_out;
  int read_tmp = _read;

  if (N == -1)
    N_out = _rw_dist;
  else
    N_out = N;

  it_assert_debug(_rw_dist >= N_out, "Attempted to peek at more elements than there are stored in the buffer.");
  out.set_size(N_out);

  for (int i = 0;i < N_out;i++) {
    out(i) = _data[read_tmp];
    read_tmp++;
    if (read_tmp == _ndata)
      read_tmp = 0;
  }
}

template <class T>
void Circular_Buffer<T>::peek(const ivec& index, Array<T>& out) const
{
  out.set_size(index.size());

  for (int i = 0;i < index.size();i++) {
    it_assert_debug(_rw_dist >= index(i) && index(i) >= 0, "Attempted to peek at an element, whose index exceeds the number of buffered elements.");
    out(i) = _data[(_read+index(i))%_ndata];
  }
}

template <class T>
void Circular_Buffer<T>::peek_reverse(T& out) const
{
  int read_tmp;

  it_assert_debug(_rw_dist > 0, "Attempted to peek at an empty buffer.");

  if (_write > 0)
    read_tmp = _write - 1;
  else
    read_tmp = _ndata - 1;

  out = _data[read_tmp];
}

template <class T>
T Circular_Buffer<T>::peek_reverse() const
{
  T out;

  peek_reverse(out);
  return out;
}

template <class T>
void Circular_Buffer<T>::peek_reverse(Vec<T>& out, const int N) const
{
  int N_out;
  int read_tmp;

  if (N == -1)
    N_out = _rw_dist;
  else
    N_out = N;

  it_assert_debug(_rw_dist >= N_out, "Attempted to peek at more elements than there are stored in the buffer.");
  out.set_size(N_out);

  if (_write > 0)
    read_tmp = _write - 1;
  else
    read_tmp = _ndata - 1;

  for (int i = 0;i < N_out;i++) {
    out(i) = _data[read_tmp];
    read_tmp--;
    if (read_tmp < 0)
      read_tmp = _ndata - 1;
  }
}

template <class T>
void Circular_Buffer<T>::peek_reverse(Array<T>& out, const int N) const
{
  int N_out;
  int read_tmp;

  if (N == -1)
    N_out = _rw_dist;
  else
    N_out = N;

  it_assert_debug(_rw_dist >= N_out, "Attempted to peek at more elements than there are stored in the buffer.");
  out.set_size(N_out);

  if (_write > 0)
    read_tmp = _write - 1;
  else
    read_tmp = _ndata - 1;

  for (int i = 0;i < N_out;i++) {
    out(i) = _data[read_tmp];
    read_tmp--;
    if (read_tmp < 0)
      read_tmp = _ndata - 1;
  }
}

template <class T>
void Circular_Buffer<T>::put(const T& in)
{
  //Remove the oldest element of the buffer if the buffer is full
  if (_rw_dist >= _ndata) {
    T dummy;
    get(dummy);
  }

  //Write data to the buffer and move the pointer to the next buffer slot
  _data[_write] = in;
  _write++;
  _rw_dist++;

  //Check if the pointer in the circular buffer should go back to zero
  if (_write >= _ndata)
    _write = 0;

}

template <class T>
void Circular_Buffer<T>::put(const Vec<T>& in)
{
  for (int i = 0;i < in.size();i++) {
    //Remove the oldest element of the buffer if the buffer is full
    if (_rw_dist >= _ndata) {
      T dummy;
      get(dummy);
    }

    //Write data to the buffer and move the pointer to the next buffer slot
    _data[_write] = in(i);
    _write++;
    _rw_dist++;

    //Check if the pointer in the circular buffer should go back to zero
    if (_write >= _ndata)
      _write = 0;
  }

}

template <class T>
void Circular_Buffer<T>::put(const Array<T>& in)
{
  for (int i = 0;i < in.size();i++) {
    //Remove the oldest element of the buffer if the buffer is full
    if (_rw_dist >= _ndata) {
      T dummy;
      get(dummy);
    }

    //Write data to the buffer and move the pointer to the next buffer slot
    _data[_write] = in(i);
    _write++;
    _rw_dist++;

    //Check if the pointer in the circular buffer should go back to zero
    if (_write >= _ndata)
      _write = 0;
  }
}

template <class T>
void Circular_Buffer<T>::clear()
{
  _write   = 0;
  _read    = 0;
  _rw_dist = 0;
}

template<class T>
void Circular_Buffer<T>::alloc(int n)
{
  if (n == 0) {
    _ndata = 0;
    _data  = NULL;
  }
  else if (n > 0) {
    _ndata = n;
    _data = new T[_ndata];
    it_assert(_data != 0, "Out of memory in Circular_Buffer::alloc");
  }
  else {
    it_error("Circular_Buffer<T>::alloc(int n): n must be positive");
  }
}

template<class T>
void Circular_Buffer<T>::free()
{
  delete [] _data;

  _data    = NULL;
  _ndata   = 0;
  _write   = 0;
  _read    = 0;
  _rw_dist = 0;
}

template<class T>
void Circular_Buffer<T>::operator=(const Circular_Buffer<T> &s)
{
  set_size(s._ndata);
  for (int i = 0; i < _ndata; i++)
    _data[i] = s._data[i];
  _read = s._read;
  _write = s._write;
  _rw_dist = _write - _read;
}

template<class T>
void Circular_Buffer<T>::set_size(int sz, bool copy)
{
  int i, min_nrof_elem;
  //T *tmp;
  Vec<T> tmp;

  if (_ndata == sz)
    return;

  if (copy) {
    peek_reverse(tmp, -1);
    min_nrof_elem = _rw_dist < sz ? _rw_dist : sz;
    alloc(sz);
    clear();
    for (i = 0; i < min_nrof_elem; i++)
      put(tmp(min_nrof_elem - 1 - i));
  }
  else {
    free();
    alloc(sz);
  }

  _ndata = sz;
}

} // namespace itpp

#endif // #ifndef CIRCULAR_BUFFER_H
