//=================================================================================================
/*!
//  \file blazemark/eigen/init/SparseMatrix.h
//  \brief Header file for the Eigen sparse matrix initialization functions
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

#ifndef _BLAZEMARK_EIGEN_INIT_SPARSEMATRIX_H_
#define _BLAZEMARK_EIGEN_INIT_SPARSEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <Eigen/Sparse>
#include <blaze/util/Indices.h>
#include <blazemark/system/Config.h>
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
void init( ::Eigen::SparseMatrix<Type,::Eigen::RowMajor,EigenSparseIndexType>& m, size_t nonzeros );

template< typename Type >
void init( ::Eigen::SparseMatrix<Type,::Eigen::ColMajor,EigenSparseIndexType>& m, size_t nonzeros );
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
void init( ::Eigen::SparseMatrix<Type,::Eigen::RowMajor,EigenSparseIndexType>& m, size_t nonzeros )
{
   const int M( m.rows() );
   const int N( m.cols() );

   m.reserve( M*nonzeros );

   if( structure == band )
   {
      const int rrange( nonzeros / 2 );
      const int lrange( ( nonzeros % 2 )?( rrange ):( rrange-1 ) );

      for( int i=0; i<M; ++i )
      {
         m.startVec( i );

         const int jbegin( ( i >= lrange )?( i-lrange ):( 0 ) );
         const int jend  ( ( i+rrange+1 < N )?( i+rrange+1 ):( N ) );

         for( int j=jbegin; j<jend; ++j ) {
            m.insertBack(i,j) = ::blaze::rand<Type>( 0, 10 );
         }
      }
   }
   else
   {
      for( int i=0UL; i<M; ++i ) {
         m.startVec( i );
         ::blaze::Indices indices( 0UL, N-1UL, nonzeros );
         for( ::blaze::Indices::ConstIterator it=indices.begin(); it!=indices.end(); ++it ) {
            m.insertBack(i,*it) = ::blaze::rand<Type>( 0, 10 );
         }
      }
   }

   m.finalize();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major sparse matrix.
//
// \param m The column-major sparse matrix to be initialized.
// \param nonzeros The number of non-zero elements per column.
// \return void
//
// This function initializes the given column-major sparse matrix with random values. Each column
// will be filled with \a nonzeros non-zero elements, whose indices are randomly determined.
*/
template< typename Type >  // Data type of the matrix
void init( ::Eigen::SparseMatrix<Type,::Eigen::ColMajor,EigenSparseIndexType>& m, size_t nonzeros )
{
   const int M( m.rows() );
   const int N( m.cols() );

   m.reserve( N*nonzeros );

   if( structure == band )
   {
      const int drange( nonzeros / 2 );
      const int urange( ( nonzeros % 2 )?( drange ):( drange-1 ) );

      for( int j=0; j<N; ++j )
      {
         m.startVec( j );

         const int ibegin( ( j >= urange )?( j-urange ):( 0 ) );
         const int iend  ( ( j+drange+1 < M )?( j+drange+1 ):( M ) );

         for( int i=ibegin; i<iend; ++i ) {
            m.insertBack(i,j) = ::blaze::rand<Type>( 0, 10 );
         }
      }
   }
   else
   {
      for( int j=0UL; j<N; ++j ) {
         m.startVec( j );
         ::blaze::Indices indices( 0UL, M-1UL, nonzeros );
         for( ::blaze::Indices::ConstIterator it=indices.begin(); it!=indices.end(); ++it ) {
            m.insertBack(*it,j) = ::blaze::rand<Type>( 0, 10 );
         }
      }
   }

   m.finalize();
}
//*************************************************************************************************

} // namespace eigen

} // namespace blazemark

#endif
