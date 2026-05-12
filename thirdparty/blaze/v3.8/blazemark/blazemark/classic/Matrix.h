//=================================================================================================
/*!
//  \file blazemark/classic/Matrix.h
//  \brief Implementation of classic arbitrary sized matrix
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

#ifndef _BLAZEMARK_CLASSIC_MATRIX_H_
#define _BLAZEMARK_CLASSIC_MATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <blaze/system/Restrict.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blazemark/classic/Vector.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace classic {

//*************************************************************************************************
/*!\brief Efficient implementation of a \f$ M \times N \f$ matrix.
//
// The Matrix class is the representation of a dynamic \f$ M \times N \f$ matrix with a total
// of \f$ M \cdot N \f$ dynamically allocated elements. These elements can be directly accessed
// with the 1D subscript operator or with the 2D function operator. The matrix is stored in a
// row-wise fashion:

               \f[\left(\begin{array}{*{5}{c}}
               0           & 1             & 2             & \cdots & N-1         \\
               N           & N+1           & N+2           & \cdots & 2 \cdot N-1 \\
               \vdots      & \vdots        & \vdots        & \ddots & \vdots      \\
               M \cdot N-N & M \cdot N-N+1 & M \cdot N-N+2 & \cdots & M \cdot N-1 \\
               \end{array}\right)\f]

// Matrix can be used with any non-cv-qualified fundamental element type. The arithmetic operators
// for matrix/matrix, matrix/vector and matrix/element operations with the same element type work
// for any element type as long as the element type supports the arithmetic operation. Arithmetic
// operations between matrices, vectors and elements of different element types are not supported.
*/
template< typename Type      // Data type of the matrix
        , bool SO = false >  // Storage order
class Matrix
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline Matrix();
   inline Matrix( size_t m, size_t n );
   inline Matrix( size_t m, size_t n, const Type& init );
   inline Matrix( const Matrix& m );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~Matrix();
   //@}
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   inline Matrix&     operator=( const Matrix& rhs );
   inline Type&       operator[]( size_t index );
   inline const Type& operator[]( size_t index ) const;
   inline Type&       operator()( size_t row, size_t col );
   inline const Type& operator()( size_t row, size_t col ) const;
   inline Matrix&     operator+=( const Matrix& rhs );
   //@}
   //**********************************************************************************************

    //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t rows   () const;
   inline size_t columns() const;
   inline void   reset  ();
   inline void   resize ( size_t m, size_t n, bool preserve );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;             //!< The current number of rows of the matrix.
   size_t n_;             //!< The current number of columns of the matrix.
   size_t capacity_;      //!< The maximum capacity of the matrix.
   Type* BLAZE_RESTRICT v_;  //!< The dynamically allocated matrix elements.
                          /*!< Access to the matrix elements is gained via the subscript or
                               function call operator. The order of the elements is
                               \f[\left(\begin{array}{*{5}{c}}
                               0            & 1             & 2             & \cdots & N-1         \\
                               N            & N+1           & N+2           & \cdots & 2 \cdot N-1 \\
                               \vdots       & \vdots        & \vdots        & \ddots & \vdots      \\
                               M \cdot N-N  & M \cdot N-N+1 & M \cdot N-N+2 & \cdots & M \cdot N-1 \\
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
/*!\brief The default constructor for Matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Matrix<Type,SO>::Matrix()
   : m_       ( 0UL )      // The current number of rows of the matrix
   , n_       ( 0UL )      // The current number of columns of the matrix
   , capacity_( 0UL )      // The maximum capacity of the matrix
   , v_       ( nullptr )  // The matrix elements
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a matrix of size \f$ m \times n \f$. No element initialization is performed!
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
//
// \note: This constructor is only responsible to allocate the required dynamic memory. No
// element initialization is performed!
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Matrix<Type,SO>::Matrix( size_t m, size_t n )
   : m_       ( m )                    // The current number of rows of the matrix
   , n_       ( n )                    // The current number of columns of the matrix
   , capacity_( m*n )                  // The maximum capacity of the matrix
   , v_       ( new Type[capacity_] )  // The matrix elements
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a homogenous initialization of all \f$ m \times n \f$ matrix elements.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param init The initial value of the matrix elements.
//
// All matrix elements are initialized with the specified value.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Matrix<Type,SO>::Matrix( size_t m, size_t n, const Type& init )
   : m_       ( m )                    // The current number of rows of the matrix
   , n_       ( n )                    // The current number of columns of the matrix
   , capacity_( m*n )                  // The maximum capacity of the matrix
   , v_       ( new Type[capacity_] )  // The matrix elements
{
   for( size_t i=0UL; i<capacity_; ++i )
      v_[i] = init;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for Matrix.
//
// \param m Matrix to be copied.
//
// The copy constructor is explicitly defined due to the required dynamic memory management
// and in order to enable/facilitate NRV optimization.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Matrix<Type,SO>::Matrix( const Matrix& m )
   : m_       ( m.m_  )                // The current number of rows of the matrix
   , n_       ( m.n_  )                // The current number of columns of the matrix
   , capacity_( m_*n_ )                // The maximum capacity of the matrix
   , v_       ( new Type[capacity_] )  // The matrix elements
{
   for( size_t i=0UL; i<capacity_; ++i )
      v_[i] = m.v_[i];
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for Matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Matrix<Type,SO>::~Matrix()
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
/*!\brief Copy assignment operator for Matrix.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
//
// The matrix is resized according to the given \f$ M \times N \f$ matrix and initialized as a
// copy of this matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Matrix<Type,SO>& Matrix<Type,SO>::operator=( const Matrix& rhs )
{
   if( &rhs == this ) return *this;

   resize( rhs.m_, rhs.n_, false );

   const size_t sqrsize( m_*n_ );
   for( size_t i=0UL; i<sqrsize; ++i )
      v_[i] = rhs.v_[i];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 1D-access to the matrix elements.
//
// \param index Access index. The index has to be in the range \f$[0..M \cdot N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Type& Matrix<Type,SO>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < m_*n_, "Invalid matrix access index" );
   return v_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 1D-access to the matrix elements.
//
// \param index Access index. The index has to be in the range \f$[0..M \cdot N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline const Type& Matrix<Type,SO>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < m_*n_, "Invalid matrix access index" );
   return v_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 2D-access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Type& Matrix<Type,SO>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i<m_ && j<n_, "Invalid matrix access index" );

   if( SO ) return v_[i+j*m_];
   else     return v_[i*n_+j];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 2D-access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline const Type& Matrix<Type,SO>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i<m_ && j<n_, "Invalid matrix access index" );

   if( SO ) return v_[i+j*m_];
   else     return v_[i*n_+j];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a dense matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side dense matrix to be added to the matrix.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline Matrix<Type,SO>& Matrix<Type,SO>::operator+=( const Matrix& rhs )
{
   if( rhs.rows() != m_ || rhs.columns() != n_ )
      throw std::invalid_argument( "Matrix sizes do not match" );

   const size_t sqrsize( m_ * n_ );
   for( size_t i=0UL; i<sqrsize; ++i )
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
/*!\brief Returns the current number of rows of the matrix.
//
// \return The number of rows of the matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline size_t Matrix<Type,SO>::rows() const
{
   return m_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of columns of the matrix.
//
// \return The number of columns of the matrix.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline size_t Matrix<Type,SO>::columns() const
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void Matrix<Type,SO>::reset()
{
   const size_t sqrsize( m_*n_ );
   for( size_t i=0UL; i<sqrsize; ++i )
      v_[i] = Type();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the matrix.
//
// \param m The new number of rows of the matrix.
// \param n The new number of columns of the matrix.
// \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
// \return void
//
// This function resizes the matrix using the given size to \f$ m \times n \f$. During this
// operation, new dynamic memory may be allocated in case the capacity of the matrix is too
// small. Therefore this function potentially changes all matrix elements. In order to preserve
// the old matrix values, the \a preserve flag can be set to \a true. However, new matrix
// elements are not initialized!\n
// The following example illustrates the resize operation of a \f$ 2 \times 4 \f$ matrix to a
// \f$ 4 \times 2 \f$ matrix. The new, uninitialized elements are marked with \a x:

                              \f[
                              \left(\begin{array}{*{4}{c}}
                              1 & 2 & 3 & 4 \\
                              5 & 6 & 7 & 8 \\
                              \end{array}\right)

                              \Longrightarrow

                              \left(\begin{array}{*{2}{c}}
                              1 & 2 \\
                              5 & 6 \\
                              x & x \\
                              x & x \\
                              \end{array}\right)
                              \f]
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline void Matrix<Type,SO>::resize( size_t m, size_t n, bool preserve )
{
   if( m == m_ && n == n_ ) return;

   if( preserve )
   {
      Type* BLAZE_RESTRICT v = new Type[m*n];
      const size_t min_m( ( m < m_ )?( m ):( m_ ) );
      const size_t min_n( ( n < n_ )?( n ):( n_ ) );

      for( size_t i=0UL; i<min_m; ++i )
         for( size_t j=0UL; j<min_n; ++j )
            v[i*n+j] = v_[i*n_+j];

      std::swap( v_, v );
      delete [] v;
      capacity_ = m*n;
   }
   else if( m*n > capacity_ ) {
      Type* BLAZE_RESTRICT v = new Type[m*n];
      std::swap( v_, v );
      delete [] v;
      capacity_ = m*n;
   }

   m_ = m;
   n_ = n;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Addition operator for the addition of two dense matrices (\f$ A=B+C \f$).
//
// \param lhs The left-hand side dense matrix for the matrix addition.
// \param rhs The right-hand side dense matrix to be added to the left-hand side matrix.
// \return The sum of the two matrices.
// \exception std::invalid_argument Matrix sizes do not match
//
// This operator represents the addition of two dense matrices:

   \code
   blazemark::classic::Matrix<double> A, B, C;
   // ... Resizing and initialization
   C = A + B;
   \endcode

// In case the current number of rows and columns of the two given  matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename Type  // Data type of the matrix
        , bool SO1       // Storage order of the left-hand side matrix
        , bool SO2 >     // Storage order of the right-hand side matrix
inline const Matrix<Type,SO1> operator+( const Matrix<Type,SO1>& A, const Matrix<Type,SO2>& B )
{
   if( A.rows() != B.rows() || A.columns() != B.columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   Matrix<Type,SO1> C( A.rows(), A.columns() );

   for( size_t i=0UL; i<A.rows(); ++i ) {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         C(i,j) = A(i,j) + B(i,j);
      }
   }

   return C;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction operator for the subtraction of two demse matrices (\f$ A=B-C \f$).
//
// \param lhs The left-hand side dense matrix for the matrix subtraction.
// \param rhs The right-hand side dense matrix to be subtracted from the left-hand side matrix.
// \return The difference of the two matrices.
// \exception std::invalid_argument Matrix sizes do not match
//
// This operator represents the subtraction of two dense matrices:

   \code
   blazemark::classic::Matrix<double> A, B, C;
   // ... Resizing and initialization
   C = A - B;
   \endcode

// In case the current number of rows and columns of the two given  matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename Type  // Data type of the matrix
        , bool SO1       // Storage order of the left-hand side matrix
        , bool SO2 >     // Storage order of the right-hand side matrix
inline const Matrix<Type,SO1> operator-( const Matrix<Type,SO1>& A, const Matrix<Type,SO2>& B )
{
   if( A.rows() != B.rows() || A.columns() != B.columns() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   Matrix<Type,SO1> C( A.rows(), A.columns() );

   for( size_t i=0UL; i<A.rows(); ++i ) {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         C(i,j) = A(i,j) - B(i,j);
      }
   }

   return C;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of two row-major matrices (\f$ A=B*C \f$).
//
// \param lhs The left-hand side matrix for the multiplication.
// \param rhs The right-hand side matrix for the multiplication.
// \return The resulting matrix.
//
// This operator represents the multiplication of two row-major dense matrices:

   \code
   blazemark::classic::Matrix<double,false> A, B, C;
   // ... Resizing and initialization
   C = A * B;
   \endcode

// In case the current number of columns of \a lhs and the current number of rows of \a rhs
// don't match, a \a std::invalid_argument is thrown.
*/
template< typename Type >  // Data type of the matrix
inline const Matrix<Type,false> operator*( const Matrix<Type,false>& A, const Matrix<Type,false>& B )
{
   if( A.columns() != B.rows() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   Matrix<Type,false> C( A.rows(), B.columns() );

   // Good loop ordering
   for( size_t i=0UL; i<A.rows(); ++i ) {
      for( size_t k=0UL; k<B.columns(); ++k ) {
         C(i,k) = A(i,0UL) * B(0UL,k);
      }
      for( size_t j=1UL; j<A.columns(); ++j ) {
         for( size_t k=0UL; k<B.columns(); ++k ) {
            C(i,k) += A(i,j) * B(j,k);
         }
      }
   }

   // Bad loop ordering
//    for( size_t i=0UL; i<A.rows(); ++i ) {
//       for( size_t j=0UL; j<B.columns(); ++j ) {
//          C(i,j) = A(i,0UL) * B(0UL,j);
//          for( size_t k=1UL; k<A.columns(); ++k ) {
//             C(i,j) += A(i,k) * B(k,j);
//          }
//       }
//    }

   return C;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of two column-major matrices (\f$ A=B*C \f$).
//
// \param lhs The left-hand side matrix for the multiplication.
// \param rhs The right-hand side matrix for the multiplication.
// \return The resulting matrix.
//
// This operator represents the multiplication of two column-major dense matrices:

   \code
   blazemark::classic::Matrix<double,true> A, B, C;
   // ... Resizing and initialization
   C = A * B;
   \endcode

// In case the current number of columns of \a lhs and the current number of rows of \a rhs
// don't match, a \a std::invalid_argument is thrown.
*/
template< typename Type >  // Data type of the matrix
inline const Matrix<Type,true> operator*( const Matrix<Type,true>& A, const Matrix<Type,true>& B )
{
   if( A.columns() != B.rows() )
      throw std::invalid_argument( "Matrix sizes do not match" );

   Matrix<Type,true> C( A.rows(), B.columns() );

   // Good loop ordering
   for( size_t i=0UL; i<B.columns(); ++i ) {
      for( size_t k=0UL; k<A.rows(); ++k ) {
         C(k,i) = A(k,0UL) * B(0UL,i);
      }
      for( size_t j=1UL; j<A.rows(); ++j ) {
         for( size_t k=0UL; k<A.rows(); ++k ) {
            C(k,i) += A(k,j) * B(j,i);
         }
      }
   }

   // Bad loop ordering
//    for( size_t i=0UL; i<A.rows(); ++i ) {
//       for( size_t j=0UL; j<B.columns(); ++j ) {
//          C(i,j) = A(i,0UL) * B(0UL,j);
//          for( size_t k=1UL; k<A.columns(); ++k ) {
//             C(i,j) += A(i,k) * B(k,j);
//          }
//       }
//    }

   return C;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a dense matrix and a dense vector
//        (\f$ \vec{a}=B*\vec{c} \f$).
//
// \param mat The left-hand side dense matrix for the multiplication.
// \param vec The right-hand side dense vector for the multiplication.
// \return The resulting vector.
// \exception std::invalid_argument Matrix and vector sizes do not match.
//
// This operator represents the multiplication between a row-major dense matrix and a dense vector:

   \code
   blazemark::classic::Matrix<double> A;
   blazemark::classic::Vector<double> x, y;
   // ... Resizing and initialization
   y = A * x;
   \endcode

// In case the current size of the vector \a vec doesn't match the current number of columns
// of the matrix \a mat, a \a std::invalid_argument is thrown.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline const Vector<Type> operator*( const Matrix<Type,SO>& A, const Vector<Type>& x )
{
   if( A.columns() != x.size() )
      throw std::invalid_argument( "Matrix and vector sizes do not match" );

   Vector<Type> y( A.rows() );

   for( size_t i=0UL; i<A.rows(); ++i ) {
      y[i] = Type();
      for( size_t j=0UL; j<A.columns(); ++j )
         y[i] += A(i,j) * x[j];
   }

   return y;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a dense vector and a row-major
//        dense matrix (\f$ \vec{a}^T=\vec{b}^T*C \f$).
//
// \param vec The left-hand side transpose dense vector for the multiplication.
// \param mat The right-hand side dense matrix for the multiplication.
// \return The resulting transpose vector.
// \exception std::invalid_argument Vector and matrix sizes do not match.
//
// This operator represents the multiplication between a transpose dense vector and a dense
// matrix:

   \code
   blazemark::classic::Matrix<double,false> A;
   blazemark::classic::Vector<double> x, y;
   // ... Resizing and initialization
   y = x * A;
   \endcode

// In case the current size of the vector \a vec doesn't match the current number of rows of
// the matrix \a mat, a \a std::invalid_argument is thrown.
*/
template< typename Type >
inline const Vector<Type> operator*( const Vector<Type>& x, const Matrix<Type,false>& A )
{
   if( x.size() != A.rows() )
      throw std::invalid_argument( "Vector and matrix sizes do not match" );

   Vector<Type> y( A.columns() );
   const size_t m  ( A.rows()    );
   const size_t n  ( A.columns() );
   const size_t end( n & size_t(-2) );

   for( size_t j=0UL; j<n; ++j ) {
      y[j] = x[0UL] * A(0UL,j);
   }
   for( size_t i=1UL; i<m; ++i ) {
      for( size_t j=0UL; j<end; j+=2UL ) {
         y[j    ] += x[i] * A(i,j    );
         y[j+1UL] += x[i] * A(i,j+1UL);
      }
      if( end < n ) {
         y[end] += x[i] * A(i,end);
      }
   }

   return y;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a dense matrix and a scalar value
//        (\f$ A=B*s \f$).
//
// \param mat The left-hand side dense matrix for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return The scaled result matrix.
//
// This operator represents the multiplication between a dense matrix and a scalar value:

   \code
   blazemark::classic::Matrix<double> A, B;
   // ... Resizing and initialization
   B = A * 1.25;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the involved data types \a T1::ElementType and \a T2. Note that this operator only
// works for scalar values of built-in data type.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline const ::blaze::EnableIf_t< ::blaze::IsNumeric_v<Type>, Matrix<Type,SO> >
   operator*( const Matrix<Type,SO>& mat, Type scalar )
{
   Matrix<Type,SO> res( mat.rows(), mat.columns() );
   const size_t sqrsize( mat.rows() * mat.columns() );

   for( size_t i=0UL; i<sqrsize; ++i ) {
      res[i] = mat[i] * scalar;
   }

   return res;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a scalar value and a dense matrix
//        (\f$ A=s*B \f$).
//
// \param scalar The left-hand side scalar value for the multiplication.
// \param mat The right-hand side dense matrix for the multiplication.
// \return The scaled result matrix.
//
// This operator represents the multiplication between a a scalar value and dense matrix:

   \code
   blazemark::classic::Matrix<double> A, B;
   // ... Resizing and initialization
   B = 1.25 * A;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the involved data types \a T1 and \a T2::ElementType. Note that this operator only
// works for scalar values of built-in data type.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline const ::blaze::EnableIf_t< ::blaze::IsNumeric_v<Type>, Matrix<Type,SO> >
   operator*( Type scalar, const Matrix<Type,SO>& mat )
{
   Matrix<Type,SO> res( mat.rows(), mat.columns() );
   const size_t sqrsize( mat.rows() * mat.columns() );

   for( size_t i=0UL; i<sqrsize; ++i ) {
      res[i] = mat[i] * scalar;
   }

   return res;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the outer product of two dense vectors
//        (\f$ A=\vec{b}*\vec{c}^T \f$).
//
// \param lhs The left-hand side dense vector for the outer product.
// \param rhs The right-hand side transpose dense vector for the outer product.
// \return The resulting matrix.
//
// This operator represents the outer product between a dense vector and a transpose dense
// vector:

   \code
   blazemark::classic::Vector<double> a, b;
   blazemark::classic::Matrix<double> A;
   // ... Resizing and initialization
   A = a * trans(b);
   \endcode

*/
template< typename Type >
inline const Matrix<Type,false> outer( const Vector<Type>& lhs, const Vector<Type>& rhs )
{
   Matrix<Type,false> A( lhs.size(), rhs.size() );

   for( size_t i=0UL; i<lhs.size(); ++i ) {
      for( size_t j=0UL; j<rhs.size(); ++j ) {
         A(i,j) = lhs[i] * rhs[j];
      }
   }

   return A;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global output operator for dense MxN matrices.
//
// \param os Reference to the output stream.
// \param dm Reference to a constant dense matrix object.
// \return Reference to the output stream.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
inline std::ostream& operator<<( std::ostream& os, const Matrix<Type,SO>& m )
{
   for( size_t i=0UL; i<m.rows(); ++i ) {
      for( size_t j=0UL; j<m.columns(); ++j ) {
         os << std::setw(14) << m(i,j);
      }
      os << "\n";
   }

   return os;
}
//*************************************************************************************************

} // namespace classic

} // namespace blazemark

#endif
