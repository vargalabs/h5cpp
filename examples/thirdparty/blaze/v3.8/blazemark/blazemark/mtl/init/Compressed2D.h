//=================================================================================================
/*!
//  \file blazemark/mtl/init/Compressed2D.h
//  \brief Header file for the MTL compressed matrix initialization functions
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

#ifndef _BLAZEMARK_MTL_INIT_COMPRESSED2D_H_
#define _BLAZEMARK_MTL_INIT_COMPRESSED2D_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <blaze/util/Indices.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Config.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace mtl {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name MTL initialization functions */
//@{
template< typename Type >
void init( ::mtl::compressed2D< Type, ::mtl::mat::parameters< ::mtl::tag::row_major > >& m
         , size_t nonzeros );

template< typename Type >
void init( ::mtl::compressed2D< Type, ::mtl::mat::parameters< ::mtl::tag::col_major > >& m
         , size_t nonzeros );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given row-major sparse matrix.
//
// \param m The row-major sparse matrix to be initialized.
// \param nonzeros The number of non-zero elements per row.
// \return void
//
// This function initializes the given row-major sparse matrix with random values. Each row
// will be filled with \a nonzeros non-zero elements, whose indices are randomly determined.
*/
template< typename Type >  // Data type of the matrix
void init( ::mtl::compressed2D< Type, ::mtl::mat::parameters< ::mtl::tag::row_major > >& m
         , size_t nonzeros )
{
   using row_major        = ::mtl::tag::row_major;
   using row_parameters   = ::mtl::mat::parameters<row_major>;
   using row_compressed2D = ::mtl::compressed2D<Type,row_parameters>;
   using row_inserter     = ::mtl::mat::inserter<row_compressed2D>;

   const size_t M( num_rows( m ) );
   const size_t N( num_cols( m ) );

   row_inserter ins( m );

   if( structure == band )
   {
      const size_t rrange( nonzeros / 2UL );
      const size_t lrange( ( nonzeros % 2UL )?( rrange ):( rrange-1UL ) );

      for( size_t i=0UL; i<M; ++i )
      {
         const size_t jbegin( ( i >= lrange )?( i-lrange ):( 0UL ) );
         const size_t jend  ( ( i+rrange+1UL < N )?( i+rrange+1UL ):( N ) );

         for( size_t j=jbegin; j<jend; ++j ) {
            ins[i][j] = ::blaze::rand<Type>( 0, 10 );
         }
      }
   }
   else
   {
      for( size_t i=0UL; i<M; ++i ) {
         ::blaze::Indices indices( 0UL, N-1UL, nonzeros );
         for( ::blaze::Indices::ConstIterator it=indices.begin(); it!=indices.end(); ++it ) {
            ins[i][*it] = ::blaze::rand<Type>( 0, 10 );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major sparse matrix.
//
// \param m The column-major sparse matrix to be initialized.
// \param nonzeros The number of non-zero elements per column.
// \return void
//
// This function initializes the given column-major sparse matrix with random values.
// Each column will be filled with \a nonzeros non-zero elements, whose indices are randomly
// determined.
*/
template< typename Type >  // Data type of the matrix
void init( ::mtl::compressed2D< Type, ::mtl::mat::parameters< ::mtl::tag::col_major > >& m
         , size_t nonzeros )
{
   using col_major        = ::mtl::tag::col_major;
   using col_parameters   = ::mtl::mat::parameters<col_major>;
   using col_compressed2D = ::mtl::compressed2D<Type,col_parameters>;
   using col_inserter     = ::mtl::mat::inserter<col_compressed2D>;

   const size_t M( num_rows( m ) );
   const size_t N( num_cols( m ) );

   col_inserter ins( m );

   if( structure == band )
   {
      const size_t drange( nonzeros / 2UL );
      const size_t urange( ( nonzeros % 2UL )?( drange ):( drange-1UL ) );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( j >= urange )?( j-urange ):( 0UL ) );
         const size_t iend  ( ( j+drange+1UL < M )?( j+drange+1UL ):( M ) );

         for( size_t i=ibegin; i<iend; ++i ) {
            ins[i][j] = ::blaze::rand<Type>( 0, 10 );
         }
      }
   }
   else
   {
      for( size_t j=0UL; j<N; ++j ) {
         ::blaze::Indices indices( 0UL, M-1UL, nonzeros );
         for( ::blaze::Indices::ConstIterator it=indices.begin(); it!=indices.end(); ++it ) {
            ins[*it][j] = ::blaze::rand<Type>( 0, 10 );
         }
      }
   }
}
//*************************************************************************************************

} // namespace mtl

} // namespace blazemark

#endif
