//=================================================================================================
/*!
//  \file blazemark/flens/init/GeCRSMatrix.h
//  \brief Header file for the FLENS CRS matrix initialization functions
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

#ifndef _BLAZEMARK_FLENS_INIT_GECRSMATRIX_H_
#define _BLAZEMARK_FLENS_INIT_GECRSMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <flens/matrixtypes/general/impl/gecrsmatrix.h>
#include <blaze/util/Indices.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Config.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace flens {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name FLENS initialization functions */
//@{
template< typename Type, typename IndexType >
void init( ::flens::GeCRSMatrix< ::flens::CRS<Type,::flens::IndexBaseZero<IndexType> > >& m
         , size_t rows, size_t columns, size_t nonzeros );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given CRS matrix.
//
// \param m The CRS matrix to be initialized.
// \param rows The number of rows of the CRS matrix.
// \param columns The number of columns of the CRS matrix.
// \param nonzeros The number of non-zero elements per row.
// \return void
//
// This function initializes the given CRS matrix with random values. Each row will be filled
// with \a nonzeros non-zero elements, whose indices are randomly determined.
*/
template< typename Type         // Data type of the matrix
        , typename IndexType >  // Index type of the matrix
void init( ::flens::GeCRSMatrix< ::flens::CRS<Type,::flens::IndexBaseZero<IndexType> > >& m
         , size_t rows, size_t columns, size_t nonzeros )
{
   using IndexBase = ::flens::IndexBaseZero<IndexType>;
   using Coord     = ::flens::CoordStorage<Type,::flens::CoordRowColCmp,IndexBase>;

   const IndexType M( rows    );
   const IndexType N( columns );

   ::flens::GeCoordMatrix<Coord> tmp( M, N );

   if( structure == band )
   {
      const IndexType rrange( nonzeros / 2 );
      const IndexType lrange( ( nonzeros % 2 )?( rrange ):( rrange-1 ) );

      for( IndexType i=tmp.firstRow(); i<=tmp.lastRow(); ++i )
      {
         const IndexType jbegin( ( i >= lrange )?( i-lrange ):( 0 ) );
         const IndexType jend  ( ( i+rrange+1 < N )?( i+rrange+1 ):( N ) );

         for( IndexType j=jbegin; j<jend; ++j ) {
            tmp(i,j) += ::blaze::rand<Type>( 0, 10 );
         }
      }
   }
   else
   {
      for( IndexType i=tmp.firstRow(); i<=tmp.lastRow(); ++i ) {
         ::blaze::Indices indices( 0UL, columns-1UL, nonzeros );
         for( ::blaze::Indices::ConstIterator it=indices.begin(); it!=indices.end(); ++it ) {
            tmp(i,*it) += ::blaze::rand<Type>( 0, 10 );
         }
      }
   }

   m = tmp;
}
//*************************************************************************************************

} // namespace flens

} // namespace blazemark

#endif
