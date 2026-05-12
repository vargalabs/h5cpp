/*!
 * \file
 * \brief Stack class (container)
 * \author Thomas Eriksson
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
 * This file is not separated into a .h and a .cpp file. The reason is
 * to avoid problems with template initializations of this class.
 * An \c Stack<type> can contain any type and it is not possible to
 * initialize and pre-compile all types that might be put into an
 * \c Stack.
 */

#ifndef STACK_H
#define STACK_H

#include <itpp/base/itassert.h>


namespace itpp
{

/*!
  \brief General stack class

  This class is a general stack class for arbitrary types.

  For rarely used types you will need to instantiate the class by
  \code
  template class Stack<type>;
  \endcode

  The following example shows how to define a Stack of vectors:
  \code
  vec a = randn(10);
  vec b = randn(20);
  Stack<vec> my_stack(10);
  my_stack.push(a);
  my_stack.push(b);

  cout << my_stack.pop() << " " << my_stack.pop() << endl ;
  \endcode
*/
template<class T>
class Stack
{
public:
  //! Default constructor
  Stack();
  //! Create a Stack of size \c n
  Stack(int n);
  //! Create a copy of \c s
  Stack(const Stack<T> &s);
  //! Default destructor
  virtual ~Stack();

  //! Pop the topmost element of the stack
  T pop();
  //! Peek at the topmost element of the stack, without removing it
  T peek() const;
  //! Push an element at top of stack
  void push(T v);
  //! Empty the stack
  void clear();

  //! Assignment operator
  void operator=(const Stack<T> &s);

  //! Returns the maximum number of data elements the stack can store
  int size() const { return ndata; }
  //! Returns the number of data elements currently in the stack
  int no_elements() const { return valptr; }
  //! Resizing a Stack<T>.
  void set_size(int n, bool copy = false);

private:
  int valptr;
  int ndata;
  T *data;

private:
  void alloc(int n);
  void free();
};

// --------------------------- Implementation starts here ----------------------------------

template<class T>
Stack<T>::Stack()
{
  data = 0;
  ndata = 0;
  valptr = 0;
}

template<class T>
Stack<T>::Stack(int n)
{
  alloc(n);
  valptr = 0;
}

template<class T>
Stack<T>::Stack(const Stack<T> &s)
{
  data = NULL;
  ndata = 0;
  valptr = s.valptr;
  alloc(s.ndata);
  for (int i = 0; i < s.ndata; i++)
    data[i] = s.data[i];
}

template<class T>
Stack<T>::~Stack()
{
  free();
}

template <class T>
T Stack<T>::pop()
{
  it_error_if(valptr == 0, "Stack<T>::pop: Empty stack");
  valptr--;
  return data[valptr];
}

template <class T>
T Stack<T>::peek() const
{
  it_error_if(valptr == 0, "Stack<T>::peek: Empty stack");
  return data[valptr-1];
}

template <class T>
void Stack<T>::push(T v)
{
  it_error_if(valptr >= ndata, "Stack<T>::push: Full stack");
  data[valptr] = v;
  valptr++;
}

template <class T>
void Stack<T>::clear()
{
  valptr = 0;
}

template<class T>
void Stack<T>::alloc(int n)
{
  if (n == 0) {
    data = NULL;
    ndata = 0;
  }
  else {
    data = new T[n];
    it_assert_debug(data != 0, "Out of memory in Stack::alloc");
  }
  ndata = n;
}

template<class T>
void Stack<T>::free()
{

  delete [] data;

  data = 0;
  ndata = 0;
}

template<class T>
void Stack<T>::operator=(const Stack<T> &s)
{
  set_size(s.ndata);
  for (int i = 0; i < ndata; i++)
    data[i] = s.data[i];
  valptr = 0;
}

template<class T>
void Stack<T>::set_size(int sz, bool copy)
{
  int i, min;
  T *tmp;

  if (ndata == sz)
    return;

  if (copy) {
    tmp = data;
    min = ndata < sz ? ndata : sz;
    alloc(sz);
    for (i = 0; i < min; i++)
      data[i] = tmp[i];
    delete [] tmp;
  }
  else {
    free();
    alloc(sz);
  }
  ndata = sz;
}

} // namespace itpp

#endif // #ifndef STACK_H
