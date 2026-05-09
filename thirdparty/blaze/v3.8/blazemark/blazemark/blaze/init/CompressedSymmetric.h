//=================================================================================================
/*!
//  \file blazemark/blaze/init/CompressedSymmetric.h
//  \brief Header file for the Blaze symmetric compressed matrix initialization functions
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

#ifndef _BLAZEMARK_BLAZE_INIT_COMPRESSEDSYMMETRIC_H_
#define _BLAZEMARK_BLAZE_INIT_COMPRESSEDSYMMETRIC_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <vector>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/SymmetricMatrix.h>
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
void init( ::blaze::SymmetricMatrix< ::blaze::CompressedMatrix<Type,::blaze::rowMajor> >& m, size_t nonzeros );

template< typename Type >
void init( ::blaze::SymmetricMatrix< ::blaze::CompressedMatrix<Type,::blaze::columnMajor> >& m, size_t nonzeros );

template< typename Type, bool SO >
void init( ::std::vector< ::blaze::SymmetricMatrix< ::blaze::CompressedMatrix<Type,SO> > >& v, size_t nonzeros );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given row-major symmetric compressed matrix.
//
// \param m The row-major symmetric compressed matrix to be initialized.
// \param nonzeros The number of non-zero elements per row.
// \return void
//
// This function initializes the given row-major symmetric compressed matrix with random values.
// Each row will be filled with approximately \a nonzeros non-zero elements, whose indices are
// randomly determined.
*/
template< typename Type >  // Data type of the matrix
void init( ::blaze::SymmetricMatrix< ::blaze::CompressedMatrix<Type,::blaze::rowMajor> >& m, size_t nonzeros )
{
   const size_t N( m.rows() );

   m.reserve( N * nonzeros );
   for( size_t i=0UL; i<N; ++i ) {
      m.reserve( i, nonzeros );
   }

   if( structure == band )
   {
      const size_t half( ( nonzeros + 1UL ) / 2UL );

      for( size_t i=0UL; i<N; ++i )
      {
         const size_t jbegin( ( i+1UL > half )?( i+1UL-half ):( 0UL ) );

         for( size_t j=jbegin; j<=i; ++j ) {
            m.append( i, j, ::blaze::rand<Type>( 0, 10 ) );
         }
      }
   }
   else
   {
      for( size_t i=0UL; i<N; ++i ) {
         ::blaze::Indices indices( 0UL, N-1UL, nonzeros );
         for( ::blaze::Indices::ConstIterator it=indices.begin(); it!=indices.end(); ++it ) {
            if( *it > i ) break;
            m.append( i, *it, ::blaze::rand<Type>( 0, 10 ) );
         }
      }
   }

   m.trim();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major symmetric compressed matrix.
//
// \param m The column-major symmetric compressed matrix to be initialized.
// \param nonzeros The maximum number of non-zero elements per column.
// \return void
//
// This function initializes the given column-major symmetric compressed matrix with random values.
// Each column will be filled with approximately \a nonzeros non-zero elements, whose indices are
// randomly determined.
*/
template< typename Type >  // Data type of the matrix
void init( ::blaze::SymmetricMatrix< ::blaze::CompressedMatrix<Type,::blaze::columnMajor> >& m, size_t nonzeros )
{
   const size_t N( m.rows() );

   m.reserve( N * nonzeros );
   for( size_t j=0UL; j<N; ++j ) {
      m.reserve( j, nonzeros );
   }

   if( structure == band )
   {
      const size_t half( ( nonzeros + 1UL ) / 2UL );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( j+1UL > half )?( j+1UL-half ):( 0UL ) );

         for( size_t i=ibegin; i<=j; ++i ) {
            m.append( i, j, ::blaze::rand<Type>( 0, 10 ) );
         }
      }
   }
   else
   {
      for( size_t j=0UL; j<N; ++j ) {
         ::blaze::Indices indices( 0UL, N-1UL, nonzeros );
         for( ::blaze::Indices::ConstIterator it=indices.begin(); it!=indices.end(); ++it ) {
            if( *it > j ) break;
            m.append( *it, j, ::blaze::rand<Type>( 0, 10 ) );
         }
      }
   }

   m.trim();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given vector of symmetric compressed matrices.
//
// \param v The vector of symmetric compressed matrices to be initialized.
// \param nonzeros The number of non-zero elements per row.
// \return void
//
// This function initializes the all symmetric compressed matrices in the given vector with random
// values. Each row will be filled with approximately \a nonzeros non-zero elements, whose indices
// are randomly determined.
*/
template< typename Type  // Data type of the matrix
        , bool SO >      // Storage order
void init( ::std::vector< ::blaze::SymmetricMatrix< ::blaze::CompressedMatrix<Type,SO> > >& v, size_t nonzeros )
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
