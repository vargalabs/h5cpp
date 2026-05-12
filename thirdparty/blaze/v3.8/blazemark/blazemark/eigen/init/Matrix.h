//=================================================================================================
/*!
//  \file blazemark/eigen/init/Matrix.h
//  \brief Header file for the Eigen dense matrix initialization functions
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

#ifndef _BLAZEMARK_EIGEN_INIT_MATRIX_H_
#define _BLAZEMARK_EIGEN_INIT_MATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <Eigen/Dense>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace eigen {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Eigen initialization functions */
//@{
template< typename Type >
void init( ::Eigen::Matrix<Type,3,1>& v );

template< typename Type >
void init( ::Eigen::Matrix<Type,6,1>& v );

template< typename Type >
void init( ::Eigen::Matrix<Type,::Eigen::Dynamic,1>& v );

template< typename Type >
void init( ::Eigen::Matrix<Type,3,3,::Eigen::RowMajor>& m );

template< typename Type >
void init( ::Eigen::Matrix<Type,3,3,::Eigen::ColMajor>& m );

template< typename Type >
void init( ::Eigen::Matrix<Type,6,6,::Eigen::RowMajor>& m );

template< typename Type >
void init( ::Eigen::Matrix<Type,6,6,::Eigen::ColMajor>& m );

template< typename Type >
void init( ::Eigen::Matrix<Type,::Eigen::Dynamic,::Eigen::Dynamic,::Eigen::RowMajor>& m );

template< typename Type >
void init( ::Eigen::Matrix<Type,::Eigen::Dynamic,::Eigen::Dynamic,::Eigen::ColMajor>& m );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given 3D vector.
//
// \param v The 3D vector to be initialized.
// \return void
//
// This function initializes the given 3D vector with random values.
*/
template< typename Type >  // Data type of the vector
void init( ::Eigen::Matrix<Type,3,1>& v )
{
   for( int i=0; i<3; ++i ) {
      v[i] = ::blaze::rand<Type>( 0, 10 );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given 6D vector.
//
// \param v The 6D vector to be initialized.
// \return void
//
// This function initializes the given 6D vector with random values.
*/
template< typename Type >  // Data type of the vector
void init( ::Eigen::Matrix<Type,6,1>& v )
{
   for( int i=0; i<6; ++i ) {
      v[i] = ::blaze::rand<Type>( 0, 10 );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given dense vector.
//
// \param v The dense vector to be initialized.
// \return void
//
// This function initializes the given dense vector with random values.
*/
template< typename Type >  // Data type of the vector
void init( ::Eigen::Matrix<Type,::Eigen::Dynamic,1>& v )
{
   const int N( v.size() );

   for( int i=0; i<N; ++i ) {
      v[i] = ::blaze::rand<Type>( 0, 10 );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given row-major 3x3 matrix.
//
// \param m The row-major 3x3 matrix to be initialized.
// \return void
//
// This function initializes the given row-major 3x3 matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( ::Eigen::Matrix<Type,3,3,::Eigen::RowMajor>& m )
{
   for( int i=0; i<3; ++i ) {
      for( int j=0; j<3; ++j ) {
         m(i,j) = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major 3x3 matrix.
//
// \param m The column-major 3x3 matrix to be initialized.
// \return void
//
// This function initializes the given column-major 3x3 matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( ::Eigen::Matrix<Type,3,3,::Eigen::ColMajor>& m )
{
   for( int i=0; i<3; ++i ) {
      for( int j=0; j<3; ++j ) {
         m(i,j) = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given row-major 6x6 matrix.
//
// \param m The row-major 6x6 matrix to be initialized.
// \return void
//
// This function initializes the given row-major 6x6 matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( ::Eigen::Matrix<Type,6,6,::Eigen::RowMajor>& m )
{
   for( int i=0; i<6; ++i ) {
      for( int j=0; j<6; ++j ) {
         m(i,j) = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major 6x6 matrix.
//
// \param m The column-major 6x6 matrix to be initialized.
// \return void
//
// This function initializes the given column-major 6x6 matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( ::Eigen::Matrix<Type,6,6,::Eigen::ColMajor>& m )
{
   for( int i=0; i<6; ++i ) {
      for( int j=0; j<6; ++j ) {
         m(i,j) = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given row-major dense matrix.
//
// \param m The row-major dense matrix to be initialized.
// \return void
//
// This function initializes the given row-major dense matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( ::Eigen::Matrix<Type,::Eigen::Dynamic,::Eigen::Dynamic,::Eigen::RowMajor>& m )
{
   const int M( m.rows() );
   const int N( m.cols() );

   for( int i=0; i<M; ++i ) {
      for( int j=0; j<N; ++j ) {
         m(i,j) = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major dense matrix.
//
// \param m The column-major dense matrix to be initialized.
// \return void
//
// This function initializes the given column-major dense matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( ::Eigen::Matrix<Type,::Eigen::Dynamic,::Eigen::Dynamic,::Eigen::ColMajor>& m )
{
   const int M( m.rows() );
   const int N( m.cols() );

   for( int j=0; j<N; ++j ) {
      for( int i=0; i<M; ++i ) {
         m(i,j) = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
//*************************************************************************************************

} // namespace eigen

} // namespace blazemark

#endif
