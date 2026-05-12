//=================================================================================================
/*!
//  \file blazemark/blaze/init/StaticUniLower.h
//  \brief Header file for the Blaze unilower static matrix initialization functions
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

#ifndef _BLAZEMARK_BLAZE_INIT_STATICUNILOWER_H_
#define _BLAZEMARK_BLAZE_INIT_STATICUNILOWER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <vector>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace blaze {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Blaze initialization functions */
//@{
template< typename Type, size_t N >
void init( ::blaze::UniLowerMatrix< ::blaze::StaticMatrix<Type,N,N,::blaze::rowMajor> >& m );

template< typename Type, size_t N >
void init( ::blaze::UniLowerMatrix< ::blaze::StaticMatrix<Type,N,N,::blaze::columnMajor> >& m );

template< typename Type, size_t N, bool SO, typename A >
void init( ::std::vector< ::blaze::UniLowerMatrix< ::blaze::StaticMatrix<Type,N,N,SO> >, A >& v );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given row-major unilower static matrix.
//
// \param m The row-major unilower static matrix to be initialized.
// \return void
//
// This function initializes the given row-major unilower static matrix with random values.
*/
template< typename Type  // Data type of the matrix
        , size_t N >     // Number of rows and columns
void init( ::blaze::UniLowerMatrix< ::blaze::StaticMatrix<Type,N,N,::blaze::rowMajor> >& m )
{
   for( size_t i=1UL; i<N; ++i ) {
      for( size_t j=0UL; j<i; ++j ) {
         m(i,j) = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major unilower static matrix.
//
// \param m The column-major unilower static matrix to be initialized.
// \return void
//
// This function initializes the given column-major unilower static matrix with random values.
*/
template< typename Type  // Data type of the matrix
        , size_t N >     // Number of rows and columns
void init( ::blaze::UniLowerMatrix< ::blaze::StaticMatrix<Type,N,N,::blaze::columnMajor> >& m )
{
   for( size_t j=0UL; j<N; ++j ) {
      for( size_t i=j+1UL; i<N; ++i ) {
         m(i,j) = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given vector of unilower static matrices.
//
// \param v The vector of unilower static matrices to be initialized.
// \return void
//
// This function initializes all unilower static matrices in the given vector with random values.
*/
template< typename Type  // Data type of the matrix
        , size_t N       // Number of rows and columns
        , bool SO        // Storage order
        , typename A >   // Type of the allocator
void init( ::std::vector< ::blaze::UniLowerMatrix< ::blaze::StaticMatrix<Type,N,N,SO> >, A >& v )
{
   const size_t size( v.size() );

   for( size_t i=0UL; i<size; ++i ) {
      init( v[i] );
   }
}
//*************************************************************************************************

} // namespace blaze

} // namespace blazemark

#endif
