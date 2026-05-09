//=================================================================================================
/*!
//  \file blazemark/classic/Vector.h
//  \brief Implementation of classic arbitrary sized vector
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZEMARK_CLASSIC_VECTOR_H_
#define _BLAZEMARK_CLASSIC_VECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <ostream>
#include <stdexcept>
#include <blaze/system/Restrict.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace classic {

//*************************************************************************************************
/*!\brief Classic implementation of an arbitrary sized vector.
//
// The Vector class is the representation of a vector with an arbitrary number N of dynamically
// allocated elements. The elements can be accessed directly with the subscript operator. The
// order of the elements is as following:

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right)\f]

// Vector can be used with any non-cv-qualified fundamental element type. The arithmetic operators
// for vector/vector and vector/element operations with the same element type work for any element
// type as long as the element type supports the arithmetic operation. Arithmetic operations
// between vectors and elements of different element types are not supported.
*/
template< typename Type >  // Data type of the vector
class Vector
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline Vector();
   inline Vector( size_t n );
   inline Vector( size_t n, const Type& init );
   inline Vector( const Vector& v );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~Vector();
   //@}
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   inline Vector&     operator= ( const Vector& rhs );
   inline Type&       operator[]( size_t index );
   inline const Type& operator[]( size_t index ) const;
   inline Vector&     operator+=( const Vector& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Operators */
   //@{
   inline size_t size() const;
   inline void   reset();
   inline void   resize( size_t n, bool preserve );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t n_;                //!< The current size/dimension of the vector.
   size_t capacity_;         //!< The maximum capacity of the vector.
   Type* BLAZE_RESTRICT v_;  //!< The dynamically allocated vector elements.
                             /*!< Access to the vector elements is gained via the subscript operator.
                                  The order of the elements is
                                  \f[\left(\begin{array}{*{5}{c}}
                                  0 & 1 & 2 & \cdots & N-1 \\
                                  \end{array}\right)\f] */
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST   ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE( Type );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor for Vector.
*/
template< typename Type >  // Data type of the vector
inline Vector<Type>::Vector()
   : n_       ( 0UL )      // The current size/dimension of the vector
   , capacity_( 0UL )      // The maximum capacity of the vector
   , v_       ( nullptr )  // The vector elements
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a vector of size \a n. No element initialization is performed!
//
// \param n The size of the vector.
//
// \note: This constructor is only responsible to allocate the required dynamic memory. No
// element initialization is performed!
*/
template< typename Type >  // Data type of the vector
inline Vector<Type>::Vector( size_t n )
   : n_       ( n )                    // The current size/dimension of the vector
   , capacity_( n_ )                   // The maximum capacity of the vector
   , v_       ( new Type[capacity_] )  // The vector elements
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a homogeneous initialization of all \a n vector elements.
//
// \param n The size of the vector.
// \param init The initial value of the vector elements.
//
// All vector elements are initialized with the specified value.
*/
template< typename Type >  // Data type of the vector
inline Vector<Type>::Vector( size_t n, const Type& init )
   : n_       ( n )                    // The current size/dimension of the vector
   , capacity_( n_ )                   // The maximum capacity of the vector
   , v_       ( new Type[capacity_] )  // The vector elements
{
   for( size_t i=0UL; i<n_; ++i )
      v_[i] = init;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for Vector.
//
// \param v Vector to be copied.
//
// The copy constructor is explicitly defined due to the required dynamic memory management
// and in order to enable/facilitate NRV optimization.
*/
template< typename Type >  // Data type of the vector
inline Vector<Type>::Vector( const Vector& v )
   : n_       ( v.n_ )          // The current size/dimension of the vector
   , capacity_( v.n_ )          // The maximum capacity of the vector
   , v_       ( new Type[n_] )  // The vector elements
{
   for( size_t i=0UL; i<n_; ++i )
      v_[i] = v.v_[i];
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for Vector.
*/
template< typename Type >  // Data type of the vector
inline Vector<Type>::~Vector()
{
   delete[] v_;
}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for Vector.
//
// \param rhs Vector to be copied.
// \return Reference to the assigned vector.
//
// The vector is resized according to the given N-dimensional vector and initialized as a
// copy of this vector.
*/
template< typename Type >  // Data type of the vector
inline Vector<Type>& Vector<Type>::operator=( const Vector& rhs )
{
   if( &rhs == this ) return *this;

   resize( rhs.n_, false );

   for( size_t i=0UL; i<n_; ++i )
      v_[i] = rhs.v_[i];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename Type >  // Data type of the vector
inline Type& Vector<Type>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < n_, "Invalid vector access index" );
   return v_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename Type >  // Data type of the vector
inline const Type& Vector<Type>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < n_, "Invalid vector access index" );
   return v_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a dense vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be added to the vector.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type >  // Data type of the vector
inline Vector<Type>& Vector<Type>::operator+=( const Vector& rhs )
{
   if( n_ != rhs.size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   for( size_t i=0UL; i<n_; ++i )
      v_[i] += rhs.v_[i];

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current size/dimension of the vector.
//
// \return The size of the vector.
*/
template< typename Type >  // Data type of the vector
inline size_t Vector<Type>::size() const
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename Type >  // Data type of the vector
inline void Vector<Type>::reset()
{
   for( size_t i=0UL; i<n_; ++i )
      v_[i] = Type();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the vector.
//
// \param n The new size of the vector.
// \param preserve \a true if the old values of the vector should be preserved, \a false if not.
// \return void
//
// This function resizes the vector using the given size to \a n. During this operation, new
// dynamic memory may be allocated in case the capacity of the vector is too small. Therefore
// this function potentially changes all vector elements. In order to preserve the old vector
// values, the \a preserve flag can be set to \a true. However, new vector elements are not
// initialized!\n
// The following example illustrates the resize operation of a vector of size 2 to a vector of
// size 4. The new, uninitialized elements are marked with \a x:

                              \f[
                              \left(\begin{array}{*{2}{c}}
                              1 & 2 \\
                              \end{array}\right)

                              \Longrightarrow

                              \left(\begin{array}{*{4}{c}}
                              1 & 2 & x & x \\
                              \end{array}\right)
                              \f]
*/
template< typename Type >  // Data type of the vector
inline void Vector<Type>::resize( size_t n, bool preserve )
{
   if( n == n_ ) return;

   if( preserve )
   {
      Type* BLAZE_RESTRICT v = new Type[n];
      const size_t minsize( ( n < n_ )?( n ):( n_ ) );

      for( size_t i=0UL; i<minsize; ++i )
         v[i] = v_[i];

      std::swap( v_, v );
      delete [] v;
      capacity_ = n;
   }
   else if( n > capacity_ ) {
      Type* BLAZE_RESTRICT v = new Type[n];
      std::swap( v_, v );
      delete [] v;
      capacity_ = n;
   }

   n_ = n;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Addition operator for the addition of two dense vectors (\f$ \vec{a}=\vec{b}+\vec{c} \f$).
//
// \param lhs The left-hand side dense vector for the vector addition.
// \param rhs The right-hand side dense vector for the vector addition.
// \return The sum of the two vectors.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the addition of two dense vectors:

   \code
   blazemark::classic::Vector<double> a, b, c;
   // ... Resizing and initialization
   c = a + b;
   \endcode

// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename Type >
inline const Vector<Type> operator+( const Vector<Type>& lhs, const Vector<Type>& rhs )
{
   if( lhs.size() != rhs.size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   Vector<Type> res( lhs.size() );
   const size_t end( lhs.size() & size_t(-2) );

   for( size_t i=0UL; i<end; i+=2UL ) {
      res[i    ] = lhs[i    ] + rhs[i    ];
      res[i+1UL] = lhs[i+1UL] + rhs[i+1UL];
   }
   if( end < lhs.size() ) {
      res[end] = lhs[end] + rhs[end];
   }

   //const size_t s( lhs.size() );
   //for( size_t i=0UL; i<s; ++i )
   //   res[i] = lhs[i] + rhs[i];

   return res;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction operator for the subtraction of two dense vectors (\f$ \vec{a}=\vec{b}-\vec{c} \f$).
//
// \param lhs The left-hand side dense vector for the vector subtraction.
// \param rhs The right-hand side dense vector to be subtracted from the vector.
// \return The difference of the two vectors.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the subtraction of two dense vectors:

   \code
   blazemark::classic::Vector<double> a, b, c;
   // ... Resizing and initialization
   c = a - b;
   \endcode

// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename Type >
inline const Vector<Type> operator-( const Vector<Type>& lhs, const Vector<Type>& rhs )
{
   if( lhs.size() != rhs.size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   Vector<Type> res( lhs.size() );
   const size_t end( lhs.size() & size_t(-2) );

   for( size_t i=0UL; i<end; i+=2UL ) {
      res[i    ] = lhs[i    ] - rhs[i    ];
      res[i+1UL] = lhs[i+1UL] - rhs[i+1UL];
   }
   if( end < lhs.size() ) {
      res[end] = lhs[end] - rhs[end];
   }

   return res;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the componentwise product of two dense vectors
//        (\f$ \vec{a}=\vec{b}*\vec{c} \f$).
//
// \param lhs The left-hand side dense vector for the component product.
// \param rhs The right-hand side dense vector for the component product.
// \return The product of the two vectors.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the component product of two dense vectors:

   \code
   blazemark::classic::Vector<double> a, b, c;
   // ... Resizing and initialization
   c = a * b;
   \endcode

// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename Type >
inline const Vector<Type> operator*( const Vector<Type>& lhs, const Vector<Type>& rhs )
{
   if( lhs.size() != rhs.size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   Vector<Type> res( lhs.size() );
   const size_t end( lhs.size() & size_t(-2) );

   for( size_t i=0UL; i<end; i+=2UL ) {
      res[i    ] = lhs[i    ] * rhs[i    ];
      res[i+1UL] = lhs[i+1UL] * rhs[i+1UL];
   }
   if( end < lhs.size() ) {
      res[end] = lhs[end] * rhs[end];
   }

   return res;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a dense vector and a scalar value
//        (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param vec The left-hand side dense vector for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return The scaled result vector.
//
// This operator represents the multiplication between a dense vector and a scalar value:

   \code
   blazemark::classic::Vector<double> a, b;
   // ... Resizing and initialization
   b = a * 1.25;
   \endcode

// Note that this operator only works for scalar values of built-in data type.
*/
template< typename Type >
inline const ::blaze::EnableIf_t< ::blaze::IsNumeric_v<Type>, Vector<Type> >
   operator*( const Vector<Type>& vec, Type scalar )
{
   Vector<Type> res( vec.size() );
   const size_t end( vec.size() & size_t(-2) );

   for( size_t i=0UL; i<end; i+=2UL ) {
      res[i    ] = vec[i    ] * scalar;
      res[i+1UL] = vec[i+1UL] * scalar;
   }
   if( end < vec.size() ) {
      res[end] = vec[end] * scalar;
   }

   return res;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a scalar value and a dense vector
//        (\f$ \vec{a}=s*\vec{b} \f$).
//
// \param scalar The left-hand side scalar value for the multiplication.
// \param vec The right-hand side vector for the multiplication.
// \return The scaled result vector.
//
// This operator represents the multiplication between a a scalar value and dense vector:

   \code
   blazemark::classic::Vector<double> a, b;
   // ... Resizing and initialization
   b = 1.25 * a;
   \endcode

// Note that this operator only works for scalar values of built-in data type.
*/
template< typename Type >
inline const ::blaze::EnableIf_t< ::blaze::IsNumeric_v<Type>, Vector<Type> >
   operator*( Type scalar, const Vector<Type>& vec )
{
   Vector<Type> res( vec.size() );
   const size_t end( vec.size() & size_t(-2) );

   for( size_t i=0UL; i<end; i+=2UL ) {
      res[i    ] = vec[i    ] * scalar;
      res[i+1UL] = vec[i+1UL] * scalar;
   }
   if( end < vec.size() ) {
      res[end] = vec[end] * scalar;
   }

   return res;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global output operator for arbitrary sized dense vectors.
//
// \param os Reference to the output stream.
// \param v Reference to a constant vector object.
// \return Reference to the output stream.
*/
template< typename Type >
inline std::ostream& operator<<( std::ostream& os, const Vector<Type>& v )
{
   for( size_t i=0UL; i<v.size(); ++i )
      os << v[i] << "\n";
   return os;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scalar product (inner product) of two dense vectors
//        (\f$ s=\vec{a}*\vec{b} \f$).
//
// \param lhs The left-hand side dense vector for the inner product.
// \param rhs The right-hand side dense vector for the inner product.
// \return The scalar product.
// \exception std::invalid_argument Vector sizes do not match.
//
// This function represents the scalar product (inner product) of two dense vectors:

   \code
   blazemark::classic::Vector<double> a, b;
   double res;
   // ... Resizing and initialization
   res = inner( a, b );
   \endcode

// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename Type >
inline double inner( const Vector<Type>& lhs, const Vector<Type>& rhs )
{
   if( lhs.size() != rhs.size() )
      throw std::invalid_argument( "Vector sizes do not match" );

   double s( 0.0 );

   for( size_t i=0UL; i<lhs.size(); ++i )
      s += lhs[i] * rhs[i];

   return s;
}
//*************************************************************************************************

} // namespace classic

} // namespace blazemark

#endif
