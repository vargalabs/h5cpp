/*!
 * \file
 * \brief Base class for class factories and memory allocation functions
 * \author Johan Bergman and Adam Piatyszek
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

#ifndef FACTORY_H
#define FACTORY_H

#include <complex>
#include <itpp/base/binary.h>
#include <itpp/itexports.h>

namespace itpp
{

// Forward declarations
template<class T> class Array;
template<class Num_T> class Mat;
template<class Num_T> class Vec;

/*!
  \brief Base class for class factories

  A class factory (or virtual constructor) is a class that can create instances
  of another class. Factory is a base class for such factories. When declaring
  an Array, Vec or Mat, a factory can be passed as an (optional) constructor
  argument:
  \code
  // Declare a Vec<type> with size=10 and factory=DEFAULT_FACTORY
  Vec<type> a(10);

  // Declare a Vec<type> with size=10 and factory=f
  Factory f;
  Vec<type> b(10, f);
  \endcode

  By default, the factory (\c DEFAULT_FACTORY and \c f in the above examples)
  is not used at all! However, by overloading a help function called
  \e create_elements we can force Array/Vec/Mat to use the factory for element
  creation (instead of using the default constructor for the element type).

  \note It is the \e numeric elements that will be created by the factory,
  i.e. for an Array<Mat<T> >, the factory will be used for creating the Mat
  elements rather than the Array elements.

  Here is an example that (partly) defines a user-defined numeric type My_Type,
  a corresponding factory My_Factory and a corresponding help function
  create_elements<My_Type> that will be used by Array, Vec and Mat for element
  creation.
  \code
  class My_Type {
  public:
  // Default constructor
  My_Type() : data(0) {}
  // Constructor
  My_Type(int d) : data(d) {}
  .
  .
  .
  protected:
  int data;
  };

  class My_Factory : public Factory {
  public:
  // Constructor
  explicit My_Factory(int d) : init_data(d) {}
  // Destructor
  virtual ~My_Factory() {}
  // Create an n-length array of My_Type
  virtual void create(My_Type* &ptr, int n) const {ptr = new My_Type[n](init_data);}
  protected:
  int init_data;
  };

  // Create an n-length array of My_Type using My_Factory f
  template<>
  void create_elements<My_Type>(My_Type* &ptr, int n, const Factory &f)
  {
  if (const My_Factory *my_factory_ptr = dynamic_cast<const My_Factory*>(&f)) {
  // Yes, f seems to be a My_Factory. Now call the My_Factory::create method
  my_factory_ptr->create(ptr, n);
  }
  else {
  // No, f does not seem to be a My_Factory. As a fallback solution,
  // assume that f is DEFAULT_FACTORY and use the default constructor
  ptr = new My_Type[n];
  }
  }
  \endcode

  Now,
  \code
  // Declare a My_Factory for init_data = 123
  My_Factory my123_factory(123);

  // Declare a Vec<My_Type> with size 10 that uses My_Type() for element creation
  Vec<My_Type> v1(10);

  // Declare a Vec<My_Type> with size 10 that uses My_Type(123) for element creation
  Vec<My_Type> v1(10, my123_factory);
  \endcode

  For a more interesting example, see Fix_Factory.
*/
class ITPP_EXPORT Factory
{
public:
  //! Default constructor
  Factory() {}
  //! Destructor
  virtual ~Factory() {}
};

//! Default (dummy) factory
const Factory DEFAULT_FACTORY;


//! Create an n-length array of T to be used as Array, Vec or Mat elements
template<class T> inline
void create_elements(T* &ptr, int n, const Factory &)
{
  void *p = operator new(sizeof(T) * n);
  ptr = reinterpret_cast<T*>(p);
  for (int i = 0; i < n; i++) {
    new(ptr + i) T();
  }
}


//! Specialization for unsigned char data arrays (used in GF2Mat)
template<> inline
void create_elements<unsigned char>(unsigned char* &ptr, int n,
                                    const Factory &)
{
  void *p = operator new(sizeof(unsigned char) * n);
  ptr = reinterpret_cast<unsigned char*>(p);
}

//! Specialization for binary data arrays
template<> inline
void create_elements<bin>(bin* &ptr, int n, const Factory &)
{
  void *p = operator new(sizeof(bin) * n);
  ptr = reinterpret_cast<bin*>(p);
}

//! Specialization for short integer data arrays
template<> inline
void create_elements<short int>(short int* &ptr, int n, const Factory &)
{
  void *p = operator new(sizeof(short int) * n);
  ptr = reinterpret_cast<short int*>(p);
}

//! Specialization for integer data arrays
template<> inline
void create_elements<int>(int* &ptr, int n, const Factory &)
{
  void *p = operator new(sizeof(int) * n);
  ptr = reinterpret_cast<int*>(p);
}

//! Specialization for 16-byte aligned double data arrays
template<> inline
void create_elements<double>(double* &ptr, int n, const Factory &)
{
  void *p0 = operator new(sizeof(double) * n + 16);
  void *p1 = reinterpret_cast<void*>((reinterpret_cast<std::size_t>(p0) + 16)
                                     & (~(std::size_t(15))));
  *(reinterpret_cast<void**>(p1) - 1) = p0;
  ptr = reinterpret_cast<double*>(p1);
}

//! Specialization for 16-byte aligned complex double data arrays
template<> inline
void create_elements<std::complex<double> >(std::complex<double>* &ptr,
    int n, const Factory &)
{
  void *p0 = operator new(sizeof(std::complex<double>) * n + 16);
  void *p1 = reinterpret_cast<void*>((reinterpret_cast<std::size_t>(p0) + 16)
                                     & (~(std::size_t(15))));
  *(reinterpret_cast<void**>(p1) - 1) = p0;
  ptr = reinterpret_cast<std::complex<double>*>(p1);
}



//! Destroy an array of Array, Vec or Mat elements
template<class T> inline
void destroy_elements(T* &ptr, int n)
{
  if (ptr) {
    for (int i = 0; i < n; ++i) {
      ptr[i].~T();
    }
    void *p = reinterpret_cast<void*>(ptr);
    operator delete(p);
    ptr = 0;
  }
}

//! Specialization for unsigned char data arrays (used in GF2Mat)
template<> inline
void destroy_elements<unsigned char>(unsigned char* &ptr, int)
{
  if (ptr) {
    void *p = reinterpret_cast<void*>(ptr);
    operator delete(p);
	ptr = 0;
  }
}

//! Specialization for binary data arrays
template<> inline
void destroy_elements<bin>(bin* &ptr, int)
{
  if (ptr) {
    void *p = reinterpret_cast<void*>(ptr);
    operator delete(p);
    ptr = 0;
  }
}
//! Specialization for short integer data arrays
template<> inline
void destroy_elements<short int>(short int* &ptr, int)
{
  if (ptr) {
    void *p = reinterpret_cast<void*>(ptr);
    operator delete(p);
    ptr = 0;
  }
}

//! Specialization for integer data arrays
template<> inline
void destroy_elements<int>(int* &ptr, int)
{
  if (ptr) {
    void *p = reinterpret_cast<void*>(ptr);
    operator delete(p);
    ptr = 0;
  }
}

//! Specialisation for 16-byte aligned double data arrays
template<> inline
void destroy_elements<double>(double* &ptr, int)
{
  if (ptr) {
    void *p = *(reinterpret_cast<void**>(ptr) - 1);
    operator delete(p);
    ptr = 0;
  }
}

//! Specialisation for 16-byte aligned complex double data arrays
template<> inline
void destroy_elements<std::complex<double> >(std::complex<double>* &ptr, int)
{
  if (ptr) {
    void *p = *(reinterpret_cast<void**>(ptr) - 1);
    operator delete(p);
    ptr = 0;
  }
}


//! Create an n-length array of Array<T> to be used as Array elements
template<class T>
void create_elements(Array<T>* &ptr, int n, const Factory &f)
{
  void *p = operator new(sizeof(Array<T>) * n);
  ptr = reinterpret_cast<Array<T>*>(p);
  for (int i = 0; i < n; ++i) {
    new(ptr + i) Array<T>(f);
  }
}

//! Create an n-length array of Mat<T> to be used as Array elements
template<class T>
void create_elements(Mat<T>* &ptr, int n, const Factory &f)
{
  void *p = operator new(sizeof(Mat<T>) * n);
  ptr = reinterpret_cast<Mat<T>*>(p);
  for (int i = 0; i < n; ++i) {
    new(ptr + i) Mat<T>(f);
  }
}

//! Create an n-length array of Vec<T> to be used as Array elements
template<class T>
void create_elements(Vec<T>* &ptr, int n, const Factory &f)
{
  void *p = operator new(sizeof(Vec<T>) * n);
  ptr = reinterpret_cast<Vec<T>*>(p);
  for (int i = 0; i < n; ++i) {
    new(ptr + i) Vec<T>(f);
  }
}

} // namespace itpp

#endif // #ifndef FACTORY_H
