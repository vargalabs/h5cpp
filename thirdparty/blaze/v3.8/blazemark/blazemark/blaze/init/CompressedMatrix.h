//=================================================================================================
/*!
//  \file blazemark/blaze/init/CompressedMatrix.h
//  \brief Header file for the Blaze compressed matrix initialization functions
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

#ifndef _BLAZEMARK_BLAZE_INIT_COMPRESSEDMATRIX_H_
#define _BLAZEMARK_BLAZE_INIT_COMPRESSEDMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <vector>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/util/Indices.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Config.h>
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
template< typename Type >
void init( ::blaze::CompressedMatrix<Type,::blaze::rowMajor>& m, size_t nonzeros );

template< typename Type >
void init( ::blaze::CompressedMatrix<Type,::blaze::columnMajor>& m, size_t nonzeros );

template< typename Type, bool SO >
void init( ::std::vector< ::blaze::CompressedMatrix<Type,SO> >& v, size_t nonzeros );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given row-major compressed matrix.
//
// \param m The row-major compressed matrix to be initialized.
// \param nonzeros The number of non-zero elements per row.
// \return void
//
// This function initializes the given row-major compressed matrix with random values.
// Each row will be filled with \a nonzeros non-zero elements, whose indices are randomly
// determined.
*/
template< typename Type >  // Data type of the matrix
void init( ::blaze::CompressedMatrix<Type,::blaze::rowMajor>& m, size_t nonzeros )
{
   const size_t M( m.rows()    );
   const size_t N( m.columns() );

   m.reserve( M * nonzeros );

   if( structure == band )
   {
      const size_t rrange( nonzeros / 2UL );
      const size_t lrange( ( nonzeros % 2UL )?( rrange ):( rrange-1UL ) );

      for( size_t i=0UL; i<M; ++i )
      {
         const size_t jbegin( ( i >= lrange )?( i-lrange ):( 0UL ) );
         const size_t jend  ( ( i+rrange+1UL < N )?( i+rrange+1UL ):( N ) );

         for( size_t j=jbegin; j<jend; ++j ) {
            m.append( i, j, ::blaze::rand<Type>( 0, 10 ) );
         }
         m.finalize( i );
      }
   }
   else
   {
      for( size_t i=0UL; i<M; ++i ) {
         ::blaze::Indices indices( 0UL, N-1UL, nonzeros );
         for( ::blaze::Indices::ConstIterator it=indices.begin(); it!=indices.end(); ++it ) {
            m.append( i, *it, ::blaze::rand<Type>( 0, 10 ) );
         }
         m.finalize( i );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major compressed matrix.
//
// \param m The column-major compressed matrix to be initialized.
// \param nonzeros The number of non-zero elements per column.
// \return void
//
// This function initializes the given column-major compressed matrix with random values.
// Each column will be filled with \a nonzeros non-zero elements, whose indices are randomly
// determined.
*/
template< typename Type >  // Data type of the matrix
void init( ::blaze::CompressedMatrix<Type,::blaze::columnMajor>& m, size_t nonzeros )
{
   const size_t M( m.rows()    );
   const size_t N( m.columns() );

   m.reserve( M * nonzeros );

   if( structure == band )
   {
      const size_t drange( nonzeros / 2UL );
      const size_t urange( ( nonzeros % 2UL )?( drange ):( drange-1UL ) );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( j >= urange )?( j-urange ):( 0UL ) );
         const size_t iend  ( ( j+drange+1UL < M )?( j+drange+1UL ):( M ) );

         for( size_t i=ibegin; i<iend; ++i ) {
            m.append( i, j, ::blaze::rand<Type>( 0, 10 ) );
         }
         m.finalize( j );
      }
   }
   else
   {
      for( size_t j=0UL; j<N; ++j ) {
         ::blaze::Indices indices( 0UL, M-1UL, nonzeros );
         for( ::blaze::Indices::ConstIterator it=indices.begin(); it!=indices.end(); ++it ) {
            m.append( *it, j, ::blaze::rand<Type>( 0, 10 ) );
         }
         m.finalize( j );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given vector of compressed matrices.
//
// \param v The vector of compressed matrices to be initialized.
// \param nonzeros The number of non-zero elements per row.
// \return void
//
// This function initializes the all compressed matrices in the given vector with random
// values. Each row will be filled with \a nonzeros non-zero elements, whose indices are
// randomly determined.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
void init( ::std::vector< ::blaze::CompressedMatrix<Type,SO> >& v, size_t nonzeros )
{
   const size_t size( v.size() );

   for( size_t i=0UL; i<size; ++i ) {
      init( v[i], nonzeros );
   }
}
//*************************************************************************************************

} // namespace blaze

} // namespace blazemark

#endif
